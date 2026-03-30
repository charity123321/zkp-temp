#ifndef BATCHING_H
#define BATCHING_H

namespace PIP_SIMD {
    extern "C" {
        #include "pip_threads.h"
    }
}
#include <iostream>
#include <memory>
#include<bit>
#include <map>
#include <unordered_map>

#include "subroutines/poly_iop/structs.h"
#include "subroutines/pcs/srs.h"
#include "field_and_polynomial/temp.h"
#include "transcript/transcript.h"
#include "arithmetic/virtual_polynomial.h"
#include "subroutines/poly_iop/sum_check/sumcheck.h"

//#include <libff/algebra/scalar_multiplication/multiexp.hpp>



// Pairing的第一个元素
template <typename P>
struct Commitment
{
    G1<P> point;
};

template <typename P>
struct MultilinearKzgProof
{
    std::vector<G1<P>> proofs;
};

template <typename P, typename scalar_field, typename PCS>
struct BatchProof
{
    /// A sum check proof proving tilde g's sum
    IOPProof<scalar_field> sum_check_proof;
    /// f_i(point_i)
    std::vector<scalar_field> f_i_eval_at_point_i;
    /// proof for g'(a_2)
    MultilinearKzgProof<P> g_prime_proof;
};

/// Steps:
/// 1. get challenge point t from transcript
/// 2. build eq(t,i) for i in [0..k]
/// 3. build \tilde g_i(b) = eq(t, i) * f_i(b)
/// 4. compute \tilde eq_i(b) = eq(b, point_i)
/// 5. run sumcheck on \sum_i=1..k \tilde eq_i * \tilde g_i
/// 6. build g'(X) = \sum_i=1..k \tilde eq_i(a2) * \tilde g_i(X) where (a2) is
///    the sumcheck's point
/// 7. open g'(X) at point (a2)

template <typename P, typename scalar_field, typename PCS>
BatchProof<P, scalar_field, PCS> multi_open_internal(
    const MultilinearProverParam<P> &prover_param,
    const std::vector<shared_ptr<DenseMultilinearExtension<scalar_field>>> &polynomials,
    const std::vector<std::vector<scalar_field>> &points,
    const std::vector<scalar_field> &evals,
    SimpleTranscript &transcript)
{

    for (std::vector<scalar_field> eval_point : points)
    {
        transcript.append_serializable_element("eval_point", eval_point);
        
    }
    for (scalar_field eval : evals)
    {
        transcript.append_field_element<scalar_field>("eval", eval);
        
    }

    // TODO: sanity checks
    size_t num_var = (*polynomials[0]).num_vars();
    size_t k = polynomials.size();

    // 注意k必须为2的幂次
    size_t ell =std::__bit_width(k-1);

   
    // challenge point t
    std::vector<scalar_field> t = transcript.get_and_append_challenge_vector<scalar_field>("t", ell);

    
    // eq(t,i) for i in [0..k]
    std::vector<scalar_field> eq_t_i_list = build_eq_x_r_vec<scalar_field>(t);

    
    // \tilde g_i(b) = eq(t,i) * f_i(b)

    // combine the polynomials that have same opening point first to reduce the
    // cost of sum check later.

    // points=[P1,P2,P1,P3]
    // point_indices[P1]=0
    // point_indices[P2]=1
    // point_indices[P3]=2
    
    struct vectorFieldComparator {
        bool operator()(const std::vector<scalar_field>& a, const std::vector<scalar_field>& b) const {
            if (a.size() != b.size()) {
                return a.size() < b.size();
            }
            for (size_t i = 0; i < a.size(); ++i) {
                if (a[i] != b[i]) {
                    // 直接依赖我们在 RelicFr 中新写的 operator<
                    return a[i] < b[i]; 
                }
            }
            return false;
        }
    };

    std::map<std::vector<scalar_field>, size_t, vectorFieldComparator> point_indices;
    for (std::vector<scalar_field> point : points)
    {
        if (point_indices.find(point) == point_indices.end())
        {
            point_indices[point] = point_indices.size();
        }
    }
    
    // deduped_points=[P1,P2,P3]
    std::vector<std::vector<scalar_field>> deduped_points(point_indices.size());
    for (const auto &[point, idx] : point_indices)
    {
        deduped_points[idx] = point;
    }

    std::vector<shared_ptr<DenseMultilinearExtension<scalar_field>>> merged_tilde_gs;
    merged_tilde_gs.reserve(point_indices.size());

    for (size_t i = 0; i < point_indices.size(); ++i)
    {
        merged_tilde_gs.push_back(make_shared<DenseMultilinearExtension<scalar_field>>(num_var));
    }
    
    // polynomials = [g1, g2, g3, g4]
    // points = [P1, P2, P1, P3]
    // eq_t_i_list = [c1, c2, c3, c4]

    // merge_tilde_gs[0] = c1*g1 + c3*g3
    // merge_tilde_gs[1] = c2*g2
    // merge_tilde_gs[2] = c4*g4

    for (size_t i = 0; i < k; ++i)
    {

        const auto &poly = polynomials[i];
        const auto &point = points[i];
        const scalar_field &coeff = eq_t_i_list[i];

        size_t idx = point_indices[point];
        auto temp = (*poly) * coeff;
        *merged_tilde_gs[idx] = *merged_tilde_gs[idx] + temp;
    }
        

    // Compute tilde_eqs
    std::vector<shared_ptr<DenseMultilinearExtension<scalar_field>>> tilde_eqs;
    tilde_eqs.reserve(deduped_points.size());

    for (const auto &point : deduped_points)
    {
        std::vector<scalar_field> eq_b_zi = build_eq_x_r_vec(point);
        auto eq_poly = make_shared<DenseMultilinearExtension<scalar_field>>(num_var, eq_b_zi);
        tilde_eqs.push_back(eq_poly);
    }
    
    // Build virtual polynomial for Sumcheck
    VirtualPolynomial<scalar_field> sum_check_vp(num_var);

    for (size_t i = 0; i < merged_tilde_gs.size(); ++i)
    {
        std::vector<shared_ptr<DenseMultilinearExtension<scalar_field>>> mle_list = {
            merged_tilde_gs[i], tilde_eqs[i]};
        sum_check_vp.add_mle_list(mle_list, scalar_field::one());
    }

    // Run Sumcheck protocol
    SumCheck<scalar_field> sumcheck;
    IOPProof<scalar_field> sumcheck_proof = sumcheck.prove(sum_check_vp, transcript);

    // a2:=sumcheck's point
    std::vector<scalar_field> a2;
    a2.reserve(num_var);
    for (size_t i = 0; i < num_var; ++i)
    {
        a2.push_back(sumcheck_proof.point[i]);
    }

    // build g'(X) = \sum_i=1..k \tilde eq_i(a2) * \tilde g_i(X) where (a2) is the
    // sumcheck's point \tilde eq_i(a2) = eq(a2, point_i)

    auto g_prime = make_shared<DenseMultilinearExtension<scalar_field>>(num_var);

    for (size_t i = 0; i < merged_tilde_gs.size(); ++i)
    {
        scalar_field eq_i_a2 = eq_eval(a2, deduped_points[i]);
        auto temp = (*merged_tilde_gs[i]) * eq_i_a2;
        *g_prime = *g_prime + temp;
    }

    // Open g' at point a2
    // 需要完善
    PCS pcs;
    auto [g_prime_proof, g_prime_eval] = pcs.open(prover_param, g_prime, a2);

    BatchProof<P, scalar_field, PCS> batch_proof;
    batch_proof.sum_check_proof = sumcheck_proof;
    batch_proof.f_i_eval_at_point_i = evals;
    batch_proof.g_prime_proof = g_prime_proof;

    return batch_proof;
}

// Steps:
// 1. get challenge point t from transcript
// 2. build g' commitment
// 3. ensure \sum_i eq(a2, point_i) * eq(t, <i>) * f_i_evals matches the sum
//    via SumCheck verification
// 4. verify commitment

template <typename P, typename scalar_field, typename PCS>
bool batch_verify_internal(
    const MultilinearVerifierParam<P> &verifier_param,
    const std::vector<Commitment<P>> &f_i_commitments,
    const std::vector<std::vector<scalar_field>> &points,
    const BatchProof<P, scalar_field, PCS> &proof,
    SimpleTranscript &transcript)
{
    for (const auto &eval_point : points)
    {
        transcript.append_serializable_element("eval_point", eval_point);
    }

    for (const auto &eval : proof.f_i_eval_at_point_i)
    {
        transcript.append_field_element("eval", eval);
    }

    // TODO: sanity checks

    // 注意k必须为2的幂次
    size_t k = f_i_commitments.size();
    //size_t ell = log2(k);
    size_t ell=std::__bit_width(k);
    size_t num_var = proof.sum_check_proof.point.size();

    // challenge point t
    std::vector<scalar_field> t = transcript.get_and_append_challenge_vector<scalar_field>("t", ell);

    // sum check point (a2)
    std::vector<scalar_field> a2;
    a2.reserve(num_var);
    for (size_t i = 0; i < num_var; ++i)
    {
        a2.push_back(proof.sum_check_proof.point[i]);
    }

    // build g' commitment
    std::vector<scalar_field> eq_t_list = build_eq_x_r_vec(t);

    std::vector<scalar_field> scalars;
    std::vector<G1<P>> bases;
    scalars.reserve(k);
    bases.reserve(k);

    for (size_t i = 0; i < k; ++i)
    {
        scalar_field eq_i_a2 = eq_eval(a2, points[i]);
        scalar_field scalar = eq_i_a2 * eq_t_list[i];
        scalars.push_back(scalar);
        bases.push_back(f_i_commitments[i].point);
    }

    // 使用pippenge计算msm
    // ==============================================================
    // 替换 libff::multi_exp 为 SimdMSM AVX-512 Pippenger
    // ==============================================================
    typename P::G1_type g_prime_commit;
    ep_set_infty(g_prime_commit.value); // 初始化为无穷远点

    size_t num_elements = scalars.size();
    if (num_elements > 0) {
        // 动态设置 Pippenger 窗口大小
        int WBITS = (num_elements > 1024) ? 16 : 12; 
        
        // 零拷贝指针强转：利用 vector 底层连续内存的特性
        const bn_t* scalars_ptr = reinterpret_cast<const bn_t*>(scalars.data());
        const ep_t* points_ptr  = reinterpret_cast<const ep_t*>(bases.data());
        
        // 调用多线程加速，注意使用 PIP_SIMD 隔离舱避免 pair 冲突
        PIP_SIMD::pip_threads(scalars_ptr, points_ptr, num_elements, g_prime_commit.value, WBITS);
        
        // 将投影坐标正则化为仿射坐标，防止后续 Pair 验证失败
        ep_norm(g_prime_commit.value, g_prime_commit.value);
    }
    // ==============================================================

    // Ensure \sum_i eq(t, <i>) * f_i_evals matches the sum via SumCheck
    scalar_field sum = scalar_field::zero();
    for (size_t i = 0; i < k; ++i)
    {
        sum = sum + eq_t_list[i] * proof.f_i_eval_at_point_i[i];
    }

    VPAuxInfo aux_info;
    aux_info.max_degree = 2;
    aux_info.num_variables = num_var;

    SumCheck<scalar_field> sumcheck;
    auto subclaim = sumcheck.verify(sum, proof.sum_check_proof, aux_info, transcript);

    scalar_field tilde_g_eval = subclaim.expected_evaluation;

    // Verify commitment
    Commitment<P> g_prime_commit_obj{g_prime_commit};
    PCS pcs;
    bool res = pcs.verify(
        verifier_param,
        g_prime_commit_obj,
        a2, tilde_g_eval,
        proof.g_prime_proof);

    return res;
}


#endif