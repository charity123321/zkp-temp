#ifndef SNARK_H
#define SNARK_H

#include<utility>
#include"hyperplonk/hyperplonk_struct.h"
#include"subroutines/pcs/srs.h"
#include"subroutines/poly_iop/perm_check/permcheck.h"
#include"hyperplonk/witness.h"


template<typename P,typename scalar_field,typename PCS>
class HyperPlonkSNARK{
public: 
    std::pair<HyperPlonkProvingKey<P,scalar_field,PCS>,HyperPlonkVerifyingKey<P,scalar_field,PCS>> 
    preprocess(
        HyperPlonkIndex<scalar_field>& index,
        StructuredReferenceString<P,scalar_field>& pcs_srs
    );

    HyperPlonkProof<P,scalar_field,PCS>
    prove(
        HyperPlonkProvingKey<P,scalar_field,PCS>& pk,
        const std::vector<scalar_field>& pub_input,
        const std::vector<WitnessColumn<scalar_field>>& witnesses
    );

    bool verify(
        HyperPlonkVerifyingKey<P,scalar_field,PCS>& vk,
        std::vector<scalar_field>& pub_input,
        HyperPlonkProof<P,scalar_field,PCS>& proof
    );

};
template<typename P,typename scalar_field,typename PCS>
std::pair<
    HyperPlonkProvingKey<P,scalar_field,PCS>,
    HyperPlonkVerifyingKey<P,scalar_field,PCS>
> 
HyperPlonkSNARK<P,scalar_field,PCS>::preprocess(
        HyperPlonkIndex<scalar_field>& index,
        StructuredReferenceString<P,scalar_field>& pcs_srs
){
    size_t num_vars=index.num_variables();
    size_t supported_ml_degree=num_vars;

    PCS pcs;
    auto [pcs_prover_param,pcs_verifier_param]=pcs.trims(pcs_srs,supported_ml_degree);

    // build permutation oracles
    std::vector<shared_ptr<DenseMultilinearExtension<scalar_field>>> permutation_oracles;
    std::vector<Commitment<P>> perm_comms;
    size_t chunk_size=1ULL<<num_vars;

    std::vector<scalar_field> tmp;
    for(size_t i=0;i<index.num_witness_columns();++i){
        tmp.assign(index.permutation.begin()+i*chunk_size,index.permutation.begin()+(i+1)*chunk_size);
        auto perm_oracle=make_shared<DenseMultilinearExtension<scalar_field>>(num_vars,tmp);
    
        auto perm_comm=pcs.commit(pcs_prover_param,perm_oracle);
        permutation_oracles.push_back(perm_oracle);
        perm_comms.push_back(perm_comm);
    }

    // build selector oracles and commit to it
    std::vector<shared_ptr<DenseMultilinearExtension<scalar_field>>> selector_oracles;
    for(const auto& selector:index.selectors){
        SelectorColumn<scalar_field> sel;
        selector_oracles.push_back(sel.from(selector));
    }

    std::vector<Commitment<P>> selector_commitments;
    selector_commitments.reserve(selector_oracles.size());

    for(size_t i=0;i<selector_oracles.size();++i){
        auto commit=pcs.commit(pcs_prover_param,selector_oracles[i]);
        selector_commitments.push_back(commit);
    }

    HyperPlonkProvingKey<P,scalar_field,PCS> proving_key;
    proving_key.params=index.params;
    proving_key.permutation_oracle=move(permutation_oracles);
    proving_key.permutation_commitments=move(perm_comms);
    proving_key.selector_oracle=move(selector_oracles);
    proving_key.selector_commitments=move(selector_commitments);
    proving_key.pcs_param=move(pcs_prover_param);

    HyperPlonkVerifyingKey<P,scalar_field,PCS> verifying_key;
    verifying_key.params=index.params;
    verifying_key.pcs_param=move(pcs_verifier_param);
    verifying_key.selector_commitments=proving_key.selector_commitments;
    verifying_key.perm_commitments=proving_key.permutation_commitments;

    return std::make_pair(move(proving_key),move(verifying_key));
}


/// Generate HyperPlonk SNARK proof.
///
/// Inputs:
/// - `pk`: circuit proving key
/// - `pub_input`: online public input of length 2^\ell
/// - `witness`: witness assignment of length 2^n
///
/// Outputs:
/// - The HyperPlonk SNARK proof.
///
/// Steps:
///
/// 1. Commit Witness polynomials `w_i(x)` and append commitment to
///    transcript
///
/// 2. Run ZeroCheck on
///
///     `f(q_0(x),...q_l(x), w_0(x),...w_d(x))`  
///
/// where `f` is the constraint polynomial i.e.,
/// ```ignore
///     f(q_l, q_r, q_m, q_o, w_a, w_b, w_c)
///     = q_l w_a(x) + q_r w_b(x) + q_m w_a(x)w_b(x) - q_o w_c(x)
/// ```
/// in vanilla plonk, and obtain a ZeroCheckSubClaim
///
/// 3. Run permutation check on `\{w_i(x)\}` and `permutation_oracle`, and
///    obtain a PermCheckSubClaim.
///
/// 4. Generate evaluations and corresponding proofs
/// - 4.1. (deferred) batch opening prod(x) at
///   - [0, perm_check_point]
///   - [1, perm_check_point]
///   - [perm_check_point, 0]
///   - [perm_check_point, 1]
///   - [1,...1, 0]
///
/// - 4.2. permutation check evaluations and proofs
///   - 4.2.1. (deferred) wi_poly(perm_check_point)
///
/// - 4.3. zero check evaluations and proofs
///   - 4.3.1. (deferred) wi_poly(zero_check_point)
///   - 4.3.2. (deferred) selector_poly(zero_check_point)
///
/// - 4.4. public input consistency checks
///   - pi_poly(r_pi) where r_pi is sampled from transcript
///
/// - 5. deferred batch opening
template<typename P,typename scalar_field,typename PCS>
HyperPlonkProof<P,scalar_field,PCS>
HyperPlonkSNARK<P,scalar_field,PCS>::prove(
    HyperPlonkProvingKey<P,scalar_field,PCS>& pk,
    const std::vector<scalar_field>& pub_input,
    const std::vector<WitnessColumn<scalar_field>>& witnesses
){
    SimpleTranscript transcript("hyperplonk");

    size_t n=pub_input.size();
    if(!n){
        cout<<n<<endl;
        cout<<"这里为了使用pub_input"<<endl;
    }
    // prover_snaity_check  待实现

    size_t num_vars=pk.params.num_variables();

    size_t ell=log2(pk.params.num_pub_input);


    PcsAccumulator<P,scalar_field,PCS> pcs_acc(num_vars);
    // =======================================================================
    // 1. Commit Witness polynomials `w_i(x)` and append commitment to
    // transcript
    // ======================================================================= 
    std::vector<shared_ptr<DenseMultilinearExtension<scalar_field>>> witness_polys;
    witness_polys.reserve(witnesses.size());

    for(const auto& witness:witnesses){
        witness_polys.push_back(make_shared<DenseMultilinearExtension<scalar_field>>(witness.get_nv(),witness.data));
    }

    std::vector<Commitment<P>> witness_commits;
    witness_commits.reserve(witness_polys.size());


    PCS pcs;
    for(const auto& wpoly:witness_polys){
        auto commit=pcs.commit(pk.pcs_param,wpoly);
        witness_commits.push_back(commit);
    }

    for(const auto& w_com:witness_commits){
        transcript.append_serializable_element("w",w_com);
    }

    // =======================================================================
    // 2 Run ZeroCheck on
    //
    //     `f(q_0(x),...q_l(x), w_0(x),...w_d(x))`
    //
    // where `f` is the constraint polynomial i.e.,
    //
    //     f(q_l, q_r, q_m, q_o, w_a, w_b, w_c)
    //     = q_l w_a(x) + q_r w_b(x) + q_m w_a(x)w_b(x) - q_o w_c(x)
    //
    // in vanilla plonk, and obtain a ZeroCheckSubClaim
    // =======================================================================  

    VirtualPolynomial<scalar_field> fx=build_f<scalar_field>(
        pk.params.gate_func,
        pk.params.num_variables(),
        pk.selector_oracle,
        witness_polys
    );
           
    ZeroCheck<scalar_field> zero_check;
    IOPProof<scalar_field> zero_check_proof=zero_check.prove(fx,transcript);



    // =======================================================================
    // 3. Run permutation check on `\{w_i(x)\}` and `permutation_oracle`, and
    // obtain a PermCheckSubClaim.
    // =======================================================================
    PermutationCheck<P,scalar_field,PCS> perm_check;
    auto perm_check_result=perm_check.prove(
        pk.pcs_param,
        witness_polys,
        witness_polys,
        pk.permutation_oracle,
        transcript
    );

        
    
    auto [perm_check_proof,prod_x,frac_poly]=perm_check_result;
    std::vector<scalar_field> perm_check_point=perm_check_proof.zero_check_proof.point;
    // =======================================================================
    // 4. Generate evaluations and corresponding proofs
    // - permcheck
    //  1. (deferred) batch opening prod(x) at
    //   - [perm_check_point]
    //   - [perm_check_point[2..n], 0]
    //   - [perm_check_point[2..n], 1]
    //   - [1,...1, 0]
    //  2. (deferred) batch opening frac(x) at
    //   - [perm_check_point]
    //   - [perm_check_point[2..n], 0]
    //   - [perm_check_point[2..n], 1]
    //  3. (deferred) batch opening s_id(x) at
    //   - [perm_check_point]
    //  4. (deferred) batch opening perms(x) at
    //   - [perm_check_point]
    //  5. (deferred) batch opening witness_i(x) at
    //   - [perm_check_point]
    //
    // - zero check evaluations and proofs
    //   - 4.3.1. (deferred) wi_poly(zero_check_point)
    //   - 4.3.2. (deferred) selector_poly(zero_check_point)
    //
    // - 4.4. (deferred) public input consistency checks
    //   - pi_poly(r_pi) where r_pi is sampled from transcript
    // =======================================================================

    // (perm_check_point[2..n], 0)
    std::vector<scalar_field> perm_check_point_0;
    perm_check_point_0.reserve(num_vars);
    perm_check_point_0.push_back(scalar_field::zero());
    perm_check_point_0.insert(perm_check_point_0.end(),perm_check_point.begin(),perm_check_point.end()-1);

    // (perm_check_point[2..n], 1)
    std::vector<scalar_field> perm_check_point_1;
    perm_check_point_1.reserve(num_vars);
    perm_check_point_1.push_back(scalar_field::one());
    perm_check_point_1.insert(perm_check_point_1.end(),perm_check_point.begin(),perm_check_point.end()-1);

    // (1, ..., 1, 0)
    std::vector<scalar_field> prod_final_query_point(num_vars,scalar_field::one());
    prod_final_query_point[0]=scalar_field::zero();

    // prod(x)'s points
    pcs_acc.insert_poly_and_points(prod_x,perm_check_proof.prod_x_com,perm_check_point);
    pcs_acc.insert_poly_and_points(prod_x,perm_check_proof.prod_x_com,perm_check_point_0);
    pcs_acc.insert_poly_and_points(prod_x,perm_check_proof.prod_x_com,perm_check_point_1);
    pcs_acc.insert_poly_and_points(prod_x,perm_check_proof.prod_x_com,prod_final_query_point);

    
    // frac(x)'s points
    pcs_acc.insert_poly_and_points(frac_poly, perm_check_proof.frac_comm, perm_check_point);
    pcs_acc.insert_poly_and_points(frac_poly, perm_check_proof.frac_comm, perm_check_point_0);
    pcs_acc.insert_poly_and_points(frac_poly, perm_check_proof.frac_comm, perm_check_point_1);

    
    // perm(x)'s points
    for(size_t i=0;i<pk.permutation_oracle.size();++i){
        pcs_acc.insert_poly_and_points(pk.permutation_oracle[i],pk.permutation_commitments[i],perm_check_point);
    }

    // witnesses' points
    // TODO: refactor so it remains correct even if the order changed 
    for(size_t i=0;i<witness_polys.size();++i){
        pcs_acc.insert_poly_and_points(witness_polys[i],witness_commits[i],perm_check_point);
    }

    
    for(size_t i=0;i<witness_polys.size();++i){
        pcs_acc.insert_poly_and_points(witness_polys[i],witness_commits[i],zero_check_proof.point);
    }

    //   - 4.3.2. (deferred) selector_poly(zero_check_point)
    for(size_t i=0;i<pk.selector_oracle.size();++i){
        pcs_acc.insert_poly_and_points(pk.selector_oracle[i],pk.selector_commitments[i],zero_check_proof.point);
    }
    

    // - 4.4. public input consistency checks
    //   - pi_poly(r_pi) where r_pi is sampled from transcript
    std::vector<scalar_field> r_pi=transcript.get_and_append_challenge_vector<scalar_field>("r_pi",ell);
    std::vector<scalar_field> r_pi_padded=r_pi;
    std::vector<scalar_field> zeros(num_vars-ell,scalar_field::zero());
    r_pi_padded.insert(r_pi_padded.end(),zeros.begin(),zeros.end());
    // Evaluate witness_poly[0] at r_pi||0s which is equal to public_input evaluated
    // at r_pi. Assumes that public_input is a power of 2
    pcs_acc.insert_poly_and_points(witness_polys[0],witness_commits[0],r_pi_padded);

    
    // =======================================================================
    // 5. deferred batch opening
    // =======================================================================
    BatchProof<P,scalar_field,PCS> batch_openings=pcs_acc.multi_open(pk.pcs_param,transcript);
    
    HyperPlonkProof<P,scalar_field,PCS> proof;
    proof.witness_commits=move(witness_commits);
    proof.batch_openings=move(batch_openings);
    proof.zero_check_proof=move(zero_check_proof);
    proof.perm_check_proof=move(perm_check_proof);


    return proof;
}

/// Verify the HyperPlonk proof.
///
/// Inputs:
/// - `vk`: verification key
/// - `pub_input`: online public input
/// - `proof`: HyperPlonk SNARK proof
///
/// Outputs:
/// - Return a boolean on whether the verification is successful
///
/// 1. Verify zero_check_proof on
///
///     `f(q_0(x),...q_l(x), w_0(x),...w_d(x))`
///
/// where `f` is the constraint polynomial i.e.,
/// ```ignore
///     f(q_l, q_r, q_m, q_o, w_a, w_b, w_c)
///     = q_l w_a(x) + q_r w_b(x) + q_m w_a(x)w_b(x) - q_o w_c(x)
/// ```
/// in vanilla plonk, and obtain a ZeroCheckSubClaim
///
/// 2. Verify perm_check_proof on `\{w_i(x)\}` and `permutation_oracles`
///
/// 3. check subclaim validity
///
/// 4. Verify the opening against the commitment:
/// - check permutation check evaluations
/// - check zero check evaluations
/// - public input consistency checks
template<typename P,typename scalar_field,typename PCS>
bool 
HyperPlonkSNARK<P,scalar_field,PCS>:: verify(
    HyperPlonkVerifyingKey<P,scalar_field,PCS>& vk,
    std::vector<scalar_field>& pub_input,
    HyperPlonkProof<P,scalar_field,PCS>& proof
){
    SimpleTranscript transcript("hyperplonk");

    size_t num_selectors=vk.params.num_selector_columns();
    size_t num_witnesses=vk.params.num_witness_columns();
    size_t num_vars=vk.params.num_variables();

    size_t ell=log2(vk.params.num_pub_input);

    // =======================================================================
    // 0. sanity checks
    // =======================================================================
    // public input length    
    if(pub_input.size()!=vk.params.num_pub_input){
        cout<<"pub_input length mismatch"<<endl;
        throw;
    }

    // Extract evaluations from openings
    std::vector<scalar_field> prod_evals(
        proof.batch_openings.f_i_eval_at_point_i.begin(),
        proof.batch_openings.f_i_eval_at_point_i.begin()+4
    );

    std::vector<scalar_field> frac_evals(
        proof.batch_openings.f_i_eval_at_point_i.begin()+4,
        proof.batch_openings.f_i_eval_at_point_i.begin()+7
    );

    std::vector<scalar_field> perm_evals(
        proof.batch_openings.f_i_eval_at_point_i.begin()+7,
        proof.batch_openings.f_i_eval_at_point_i.begin()+7+num_witnesses
    );

    std::vector<scalar_field> witness_prem_evals(
        proof.batch_openings.f_i_eval_at_point_i.begin()+7+num_witnesses,
        proof.batch_openings.f_i_eval_at_point_i.begin()+7+2*num_witnesses
    );

    std::vector<scalar_field> witness_gate_evals(
        proof.batch_openings.f_i_eval_at_point_i.begin()+7+2*num_witnesses,
        proof.batch_openings.f_i_eval_at_point_i.begin()+7+3*num_witnesses
    );

    std::vector<scalar_field> selector_evals(
        proof.batch_openings.f_i_eval_at_point_i.begin()+7+3*num_witnesses,
        proof.batch_openings.f_i_eval_at_point_i.begin()+7+3*num_witnesses+num_selectors
    );

    scalar_field pi_eval=proof.batch_openings.f_i_eval_at_point_i.back();

    // =======================================================================
    // 1. Verify zero_check_proof on `f(q_0(x),...q_l(x), w_0(x),...w_d(x))`
    //
    // where `f` is the constraint polynomial i.e.,
    //
    //     f(q_l, q_r, q_m, q_o, w_a, w_b, w_c)
    //     = q_l w_a(x) + q_r w_b(x) + q_m w_a(x)w_b(x) - q_o w_c(x)
    //
    // =======================================================================

    // Zero check and perm check have different AuxInfo
    VPAuxInfo zero_check_aux_info;
    zero_check_aux_info.max_degree=vk.params.gate_func.degree();
    zero_check_aux_info.num_variables=num_vars;

    // push witness to transcript
    for(const auto& w_com:proof.witness_commits){
        transcript.append_serializable_element("w",w_com);
    }

    ZeroCheck<scalar_field> zero_check;
    ZeroCheckSubClaim<scalar_field> zero_check_sub_claim=zero_check.verify(
        proof.zero_check_proof,
        zero_check_aux_info,
        transcript
    );

    std::vector<scalar_field> zero_check_point=zero_check_sub_claim.point;

    // check zero check subclaim
    scalar_field f_eval=eval_f<scalar_field>(vk.params.gate_func,selector_evals,witness_gate_evals);
    if(f_eval!=zero_check_sub_claim.expected_evaluation){
        cout<<"zero check sub claim failed!"<<endl;
        throw;
    }

    // =======================================================================
    // 2. Verify perm_check_proof on `\{w_i(x)\}` and `permutation_oracle`
    // =======================================================================

    VPAuxInfo perm_check_aux_info;
    perm_check_aux_info.max_degree=proof.witness_commits.size()+1;
    perm_check_aux_info.num_variables=num_vars;

    PermutationCheck<P,scalar_field,PCS> perm_check;
    PermutationCheckSubclaim<P,scalar_field,PCS> perm_check_sub_claim=perm_check.verify(
        proof.perm_check_proof,
        perm_check_aux_info,
        transcript
    );

    std::vector<scalar_field> perm_check_point=perm_check_sub_claim.product_check_sub_claim.zero_check_sub_claim.point;
    scalar_field alpha=perm_check_sub_claim.product_check_sub_claim.alpha;
    auto [beta,gamma]=perm_check_sub_claim.challenge;

    std::vector<scalar_field> id_evals;
    for(size_t i=0;i<num_witnesses;++i){
        std::vector<scalar_field> ith_point=gen_eval_point<scalar_field>(i,log2(num_witnesses),perm_check_point);
        id_evals.push_back(vk.params.eval_id_oracle(ith_point));
    }

    // check evaluation subclaim
    scalar_field perm_gate_eval=eval_perm_gate<scalar_field>(
        prod_evals,
        frac_evals,
        witness_prem_evals,
        id_evals,
        perm_evals,
        alpha,
        beta,
        gamma,
        perm_check_point.back()
    );

    if(perm_gate_eval!=perm_check_sub_claim.product_check_sub_claim.zero_check_sub_claim.expected_evaluation){
        perm_gate_eval.print();
        cout<<endl;
        perm_check_sub_claim.product_check_sub_claim.zero_check_sub_claim.expected_evaluation.print();
        cout<<endl;
        cout<<"perm check failed!"<<endl;
        throw;
    }

    // =======================================================================
    // 3. Verify the opening against the commitment
    // ======================================================================= 

    std::vector<Commitment<P>> comms;
    std::vector<std::vector<scalar_field>> points;

    // (perm_check_point[2..n], 0)
    std::vector<scalar_field> perm_check_point_0;
    perm_check_point_0.reserve(num_vars);
    perm_check_point_0.push_back(scalar_field::zero());
    perm_check_point_0.insert(perm_check_point_0.end(),perm_check_point.begin(),perm_check_point.end()-1);

    // (perm_check_point[2..n], 1)
    std::vector<scalar_field> perm_check_point_1;
    perm_check_point_1.reserve(num_vars);
    perm_check_point_1.push_back(scalar_field::one());
    perm_check_point_1.insert(perm_check_point_1.end(),perm_check_point.begin(),perm_check_point.end()-1);

    // (1, ..., 1, 0)
    std::vector<scalar_field> prod_final_query_point(num_vars,scalar_field::one());
    prod_final_query_point[0]=scalar_field::zero();


    // prod(x)'s points
    for(int i=0;i<4;++i){
        comms.push_back(proof.perm_check_proof.prod_x_com);
    }
    points.push_back(perm_check_point);
    points.push_back(perm_check_point_0);
    points.push_back(perm_check_point_1);
    points.push_back(prod_final_query_point);

    // frac(x)'s points
    for(int i=0;i<3;++i){
        comms.push_back(proof.perm_check_proof.frac_comm);
    }
    points.push_back(perm_check_point);
    points.push_back(perm_check_point_0);
    points.push_back(perm_check_point_1);

    // perm's points
    for(const auto& pcom:vk.perm_commitments){
        comms.push_back(pcom);
        points.push_back(perm_check_point);
    }

    // witnesses's points
    // TODO: merge points
    for(const auto& wcom:proof.witness_commits){
        comms.push_back(wcom);
        points.push_back(perm_check_point);
    }
    for(const auto& wcom:proof.witness_commits){
        comms.push_back(wcom);
        points.push_back(zero_check_point);
    }

    // selector_poly(zero_check_point)
    for(const auto& com:vk.selector_commitments){
        comms.push_back(com);
        points.push_back(zero_check_point);
    }

    // - 4.4. public input consistency checks
    //   - pi_poly(r_pi) where r_pi is sampled from transcript

    std::vector<scalar_field> r_pi=transcript.get_and_append_challenge_vector<scalar_field>("r_pi",ell);

    // check public evaluation
    DenseMultilinearExtension<scalar_field> pi_poly(ell,pub_input);

    scalar_field expect_pi_eval=evaluation_no_par(pi_poly,r_pi);
    if(expect_pi_eval!=pi_eval){
        cout<<"Public input eval mismatch!"<<endl;
        throw;
    }

    std::vector<scalar_field> r_pi_padded=r_pi;
    std::vector<scalar_field> zeros(num_vars-ell,scalar_field::zero());
    r_pi_padded.insert(r_pi_padded.end(),zeros.begin(),zeros.end());

    comms.push_back(proof.witness_commits[0]);
    points.push_back(r_pi_padded);

    if(comms.size()!=proof.batch_openings.f_i_eval_at_point_i.size()){
        cout<<"Commitment and evaluations count mismatch!"<<endl;
        throw;
    }

    // check proof
    PCS pcs;
    bool res=pcs.batch_verify(
        vk.pcs_param,
        comms,
        points,
        proof.batch_openings,
        transcript
    );


    return res;
}






#endif