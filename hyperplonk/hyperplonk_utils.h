#ifndef HYPERPLONK_UTILS_H
#define HYPERPLONK_UTILS_H

#include<iostream>
#include<memory>
#include<optional>
#include"field_and_polynomial/temp.h"
#include"subroutines/pcs/batching.h"
#include"arithmetic/multilinear_polynomial.h"
#include"hyperplonk/custom_gate.h"


using namespace std;

template<typename P,typename scalar_field,typename PCS>
class PcsAccumulator{
private:
    // sequence:
    // - prod(x) at 5 points
    // - w_merged at perm check point
    // - w_merged at zero check points (#witness points)
    // - selector_merged at zero check points (#selector points)
    // - w[0] at r_pi
    size_t num_var;
    vector<shared_ptr<DenseMultilinearExtension<scalar_field>>> polynomials;
    vector<Commitment<P>> commitments;
    vector<vector<scalar_field>> points;
    vector<scalar_field> evals;

public:
    PcsAccumulator(size_t nv):num_var(nv) {}

    void insert_poly_and_points(
        const shared_ptr<DenseMultilinearExtension<scalar_field>>& poly,
        Commitment<P>& commit,
        vector<scalar_field>& point
    ){
        assert((*poly).num_vars()==point.size());
        assert((*poly).num_vars()==num_var);

        scalar_field eval=evaluation_no_par<scalar_field>(*poly,point);

        evals.push_back(eval);
        polynomials.push_back(poly);
        points.push_back(point);
        commitments.push_back(commit);
    }


    BatchProof<P,scalar_field,PCS> multi_open(
        MultilinearProverParam<P>& prover_param,
        SimpleTranscript& transcript
    ){
        PCS pcs;
        return pcs.multi_open(
            prover_param,
            polynomials,
            points,
            evals,
            transcript
        );
    }

};

/// Build MLE from matrix of witnesses.
///
/// Given a matrix := [row1, row2, ...] where
/// row1:= (a1, a2, ...)
/// row2:= (b1, b2, ...)
/// row3:= (c1, c2, ...)
///
/// output mle(a1,b1,c1, ...), mle(a2,b2,c2, ...), ...
template<typename scalar_field,typename row>
vector<shared_ptr<DenseMultilinearExtension<scalar_field>>> build_mle(
    vector<row> rows
){
    size_t num_rows=rows.size();
    size_t num_vars=log2(num_rows);
    size_t num_mles=rows[0].data.size();

    vector<shared_ptr<DenseMultilinearExtension<scalar_field>>> res;
    res.reserve(num_mles);

    for(size_t i=0;i<num_mles;++i){
        vector<scalar_field> cur_coeffs;
        cur_coeffs.reserve(num_rows);

        for(const auto& r:rows){
            cur_coeffs.push_back(r.data[i]);
        }
        res.push_back(
            make_shared<DenseMultilinearExtension<scalar_field>>(num_vars,move(cur_coeffs))
        );
    }
    return res;
}

/// build `f(w_0(x),...w_d(x))` where `f` is the constraint polynomial
/// i.e., `f(a, b, c) = q_l a(x) + q_r b(x) + q_m a(x)b(x) - q_o c(x)` in
/// vanilla plonk
template<typename scalar_field>
VirtualPolynomial<scalar_field> build_f(
    CustomizedGates& gates,
    size_t num_vars,
    vector<shared_ptr<DenseMultilinearExtension<scalar_field>>>& selector_mles,
    vector<shared_ptr<DenseMultilinearExtension<scalar_field>>>& witness_mles
){
    // 确保所有selectors witness和num_vars一致

    VirtualPolynomial<scalar_field> res(num_vars);

    // 所有门
    for(GateTerm gate:gates.gates){
        // 计算系数
        scalar_field coeff_fr;
        if(gate.coefficient<0){
            coeff_fr=-scalar_field(static_cast<uint64_t>(-gate.coefficient));
        }else{
            coeff_fr=scalar_field(static_cast<uint64_t>(gate.coefficient));
        }
        // 构建MLE列表
        vector<shared_ptr<DenseMultilinearExtension<scalar_field>>> mle_list;

        if(gate.selector_index.has_value()){
            size_t s_idx=gate.selector_index.value();
            mle_list.push_back(selector_mles[s_idx]);
        }

        // 添加witness
        for(size_t witness_idx:gate.wire_indices){
            mle_list.push_back(witness_mles[witness_idx]);
        }

        res.add_mle_list(mle_list,coeff_fr);
    }


    return res;
}


template<typename scalar_field>
scalar_field eval_f(
    CustomizedGates& gates,
    vector<scalar_field>& selector_evals,
    vector<scalar_field>& witness_evals
){
    scalar_field res=scalar_field::zero();
    for(GateTerm gate:gates.gates){
        scalar_field cur_value;
        if(gate.coefficient<0){
            cur_value=-scalar_field(static_cast<uint64_t>(-gate.coefficient));
        }else{
            cur_value=scalar_field(static_cast<uint64_t>(gate.coefficient));
        }

        if(gate.selector_index.has_value()){
            size_t s_idx=gate.selector_index.value();
            cur_value*=selector_evals[s_idx];
        }else{
            cur_value*=scalar_field::one();
        }

        for(size_t witness_idx:gate.wire_indices){
            cur_value*=witness_evals[witness_idx];
        }
        res+=cur_value;
    }
    return res;
}

// check perm check subclaim:
// proof.witness_perm_check_eval ?= perm_check_sub_claim.expected_eval
// Q(x) := prod(x) - p1(x) * p2(x)
//     + alpha * frac(x) * g1(x) * ... * gk(x)
//     - alpha * f1(x) * ... * fk(x)
//
// where p1(x) = (1-x1) * frac(x2, ..., xn, 0)
//             + x1 * prod(x2, ..., xn, 0),
// and p2(x) = (1-x1) * frac(x2, ..., xn, 1)
//           + x1 * prod(x2, ..., xn, 1)
// and gi(x) = (wi(x) + beta * perms_i(x) + gamma)
// and fi(x) = (wi(x) + beta * s_id_i(x) + gamma)

template<typename scalar_field>
scalar_field eval_perm_gate(
    vector<scalar_field>& prod_evals,
    vector<scalar_field>& frac_evals,
    vector<scalar_field>& witness_perm_evals,
    vector<scalar_field>& id_evals,
    vector<scalar_field>& perm_evals,
    scalar_field& alpha,
    scalar_field& beta,
    scalar_field& gamma,
    scalar_field& x1
){
    scalar_field p1_eval=frac_evals[1]+x1*(prod_evals[1]-frac_evals[1]);
    scalar_field p2_eval=frac_evals[2]+x1*(prod_evals[2]-frac_evals[2]);

    scalar_field f_prod_eval=scalar_field::one();
    for(size_t i=0;i<witness_perm_evals.size();++i){
        f_prod_eval*=witness_perm_evals[i]+beta*id_evals[i]+gamma;
    }

    scalar_field g_prod_eval=scalar_field::one();
    for(size_t i=0;i<witness_perm_evals.size();++i){
        g_prod_eval*=witness_perm_evals[i]+beta*perm_evals[i]+gamma;
    }

    scalar_field res=prod_evals[0]-p1_eval*p2_eval+alpha*(frac_evals[0]*g_prod_eval-f_prod_eval);

    return res;
}
#endif