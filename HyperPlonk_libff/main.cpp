#include <iostream>
#include "test.h"

#include <libff/algebra/curves/bls12_381/bls12_381_pp.hpp> 
#include"subroutines/pcs/srs.h"
#include"arithmetic/util.h"

#include"hyperplonk/custom_gate.h"

/*
// 为VPAuxInfo添加序列化支持
std::vector<uint8_t> to_bytes(const VPAuxInfo& aux) {
    std::vector<uint8_t> bytes(sizeof(VPAuxInfo));
    memcpy(bytes.data(), &aux, sizeof(VPAuxInfo));
    return bytes;
}

// 为IOPProverMessage添加序列化支持
template<typename F>
std::vector<uint8_t> to_bytes(const IOPProverMessage<F>& msg) {
    std::vector<uint8_t> bytes;
    for (const auto& eval : msg.evaluations) {
        auto eval_bytes = to_bytes(eval);
        bytes.insert(bytes.end(), eval_bytes.begin(), eval_bytes.end());
    }
    return bytes;
}
*/


int main()
{
    //test_get_index();
    
    //test_sum_check_trivial();
    //test_zero_check_trivial();
    //test_sum_check_random();
    //test_zero_check_random();
    
    /*
    
    libff::bls12_381_pp::init_public_params();
    typedef libff::bls12_381_pp P;

    cout<<"test srs"<<endl;
    test_srs_for_check_h_mask<P,G1<P>::scalar_field>(3);
    */

    //test_multilinearKzgPCS();
    //debug_interpolate_uni_poly();
    //test_multilinearKzgPCS_for_batch();
   // test_productcheck_random();


    //test_permutationcheck_trivial();

    //test_vanilla_plonk_gates();
    //test_mock_gate();
    test_hyperplonk_snark();
   // test_hyperplonk_snark_for_mock_circuit();

    return 0;
}