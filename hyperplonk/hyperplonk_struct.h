#ifndef HYPERPLONK_STRUCT_H
#define HYPERPLONK_STRUCT_H

#include"subroutines/pcs/batching.h"
#include"subroutines/poly_iop/structs.h"
#include"subroutines/poly_iop/prod_check/prodcheck.h"
#include"hyperplonk/custom_gate.h"
#include"hyperplonk/selectors.h"

/// The proof for the HyperPlonk PolyIOP, consists of the following:
///   - the commitments to all witness MLEs
///   - a batch opening to all the MLEs at certain index
///   - the zero-check proof for checking custom gate-satisfiability
///   - the permutation-check proof for checking the copy constraints
template<typename P,typename scalar_field,typename PCS>
struct HyperPlonkProof{
    // PCS commit for witnesses
    vector<Commitment<P>> witness_commits;
    BatchProof<P,scalar_field,PCS> batch_openings;
    // =======================================================================
    // IOP proofs
    // =======================================================================
    // the custom gate zerocheck proof
    IOPProof<scalar_field> zero_check_proof;
    ProductCheckProof<P,scalar_field> perm_check_proof;
};

/// The HyperPlonk instance parameters, consists of the following:
///   - the number of constraints
///   - number of public input columns
///   - the customized gate function
template<typename scalar_field>
class HyperPlonkParams{
public:
    /// the number of constraints
    size_t num_constraints;
    /// number of public input
    // public input is only 1 column and is implicitly the first witness column.
    // this size must not exceed number of constraints.
    size_t num_pub_input;
    // customized gate function
    CustomizedGates gate_func;

    HyperPlonkParams(): num_constraints(0),num_pub_input(0),gate_func() {}

    HyperPlonkParams(size_t num_constraints_,size_t num_pub_input_,CustomizedGates& gate_func_):
        num_constraints(num_constraints_),
        num_pub_input(num_pub_input_),
        gate_func(gate_func_)
    {}


    size_t num_variables(){
        return log2(num_constraints);
    }

    size_t num_selector_columns(){
        return gate_func.num_selector_columns();
    }

    size_t num_witness_columns(){
        return gate_func.num_witness_columns();
    }

    scalar_field eval_id_oracle(const vector<scalar_field>& point){
        size_t len=num_variables()+log2(num_witness_columns());
        if(point.size()!=len){
            cout<<"dismatch"<<endl;
        }

        scalar_field res=scalar_field::zero();
        scalar_field base=scalar_field::one();

        for(scalar_field v:point){
            res+=base*v;
            base+=base;
        }

        return res;

    }
    
};


/// The HyperPlonk index, consists of the following:
///   - HyperPlonk parameters
///   - the wire permutation
///   - the selector vectors
template<typename scalar_field>
class HyperPlonkIndex{
public:
    HyperPlonkParams<scalar_field> params;
    vector<scalar_field> permutation;
    vector<SelectorColumn<scalar_field>> selectors;

    HyperPlonkIndex()=default;

    HyperPlonkIndex(HyperPlonkParams<scalar_field>& params_,vector<scalar_field>& perm_,vector<SelectorColumn<scalar_field>>& selectors_):
        params(params_),
        permutation(perm_),
        selectors(selectors_)
    {}


    size_t num_variables(){
        return params.num_variables();
    }

    size_t num_selector_columns(){
        return params.num_selector_columns();
    }

    size_t num_witness_columns(){
        return params.num_witness_columns();
    }
};

/// The HyperPlonk proving key, consists of the following:
///   - the hyperplonk instance parameters
///   - the preprocessed polynomials output by the indexer
///   - the commitment to the selectors and permutations
///   - the parameters for polynomial commitment
template<typename P,typename scalar_field,typename PCS>
struct HyperPlonkProvingKey{
    HyperPlonkParams<scalar_field> params;
    vector<shared_ptr<DenseMultilinearExtension<scalar_field>>> permutation_oracle;
    vector<shared_ptr<DenseMultilinearExtension<scalar_field>>> selector_oracle;
    vector<Commitment<P>> selector_commitments;
    vector<Commitment<P>> permutation_commitments;
    MultilinearProverParam<P> pcs_param;
};

/// The HyperPlonk verifying key, consists of the following:
///   - the hyperplonk instance parameters
///   - the commitments to the preprocessed polynomials output by the indexer
///   - the parameters for polynomial commitment
template<typename P,typename scalar_field,typename PCS>
struct HyperPlonkVerifyingKey{
    HyperPlonkParams<scalar_field> params;
    MultilinearVerifierParam<P> pcs_param;
    vector<Commitment<P>> selector_commitments;
    vector<Commitment<P>> perm_commitments;
};


#endif