#ifndef PERMCHECK_H
#define PERMCHECK_H

#include "subroutines/poly_iop/prod_check/prodcheck.h"
#include "transcript/transcript.h"
#include "subroutines/pcs/multilinear_KZG.h"
#include "subroutines/poly_iop/perm_check/perm_check_util.h"

/// A permutation subclaim consists of
/// - the SubClaim from the ProductCheck
/// - Challenges beta and gamma
template <typename P, typename scalar_field, typename PCS>
struct PermutationCheckSubclaim
{
    ProductCheckSubClaim<scalar_field> product_check_sub_claim;
    pair<scalar_field, scalar_field> challenge;
};

template <typename P, typename scalar_field, typename PCS>
class PermutationCheck
{
public:
    SimpleTranscript init_transcript()
    {
        return SimpleTranscript("Initializing PermutationCheck transcript");
    }

    tuple<ProductCheckProof<P, scalar_field>,
          shared_ptr<DenseMultilinearExtension<scalar_field>>,
          shared_ptr<DenseMultilinearExtension<scalar_field>>>
    prove(
        MultilinearProverParam<P> &pcs_param,
        const vector<shared_ptr<DenseMultilinearExtension<scalar_field>>> &fxs,
        const vector<shared_ptr<DenseMultilinearExtension<scalar_field>>> &gxs,
        const vector<shared_ptr<DenseMultilinearExtension<scalar_field>>> &perms,
        SimpleTranscript &transcript)
    {
        if (fxs.empty())
        {
            cout << " fxs is empty " << endl;
            throw;
        }

        if (fxs.size() != gxs.size() || fxs.size() != perms.size())
        {
            cout << "fxs.size, gxs.size and perms.size mismatch." << endl;
            throw;
        }

        size_t num_vars = (*fxs[0]).num_vars();

        for (size_t i = 0; i < fxs.size(); ++i)
        {
            if ((*fxs[i]).num_vars() != num_vars || (*gxs[i]).num_vars() != num_vars || (*perms[i]).num_vars() != num_vars)
            {
                throw;
            }
        }

        scalar_field beta = transcript.get_and_append_challenge<scalar_field>("beta");
        scalar_field gamma = transcript.get_and_append_challenge<scalar_field>("gamma");

        auto [numerators, denominators] = computer_nums_and_denoms<scalar_field>(beta, gamma, fxs, gxs, perms);

        ProductCheck<P, scalar_field, PCS> prod_check;
        auto [proof, prod_poly, frac_poly] = prod_check.prove(pcs_param, numerators, denominators, transcript);

        return make_tuple(
            proof,
            prod_poly,
            frac_poly);
    }

    PermutationCheckSubclaim<P, scalar_field, PCS> verify(
        ProductCheckProof<P, scalar_field> &proof,
        VPAuxInfo &aux_info,
        SimpleTranscript &transcript)
    {
        scalar_field beta = transcript.get_and_append_challenge<scalar_field>("beta");
        scalar_field gamma = transcript.get_and_append_challenge<scalar_field>("gamma");

        ProductCheck<P, scalar_field, PCS> prod_check;
        ProductCheckSubClaim<scalar_field> product_check_sub_claim = prod_check.verify(proof, aux_info, transcript);

        return PermutationCheckSubclaim<P, scalar_field, PCS>{
            product_check_sub_claim,
            {beta, gamma}};
    }
};

#endif