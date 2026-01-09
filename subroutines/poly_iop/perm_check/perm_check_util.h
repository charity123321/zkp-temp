#ifndef PERM_CHECK_UTIL_H
#define PERM_CHECK_UTIL_H

#include <memory>

#include "arithmetic/multilinear_polynomial.h"

/// Returns the evaluations of two list of MLEs:
/// - numerators = (a1, ..., ak)
/// - denominators = (b1, ..., bk)
///
///  where
///  - beta and gamma are challenges
///  - (f1, ..., fk), (g1, ..., gk),
///  - (s_id1, ..., s_idk), (perm1, ..., permk) are mle-s
///
/// - ai(x) is the MLE for `fi(x) + \beta s_id_i(x) + \gamma`
/// - bi(x) is the MLE for `gi(x) + \beta perm_i(x) + \gamma`
///
/// The caller is responsible for sanity-check
template <typename F>
pair<vector<shared_ptr<DenseMultilinearExtension<F>>>,
     vector<shared_ptr<DenseMultilinearExtension<F>>>>
computer_nums_and_denoms(
    const F &beta,
    const F &gamma,
    const vector<shared_ptr<DenseMultilinearExtension<F>>> &fxs,
    const vector<shared_ptr<DenseMultilinearExtension<F>>> &gxs,
    const vector<shared_ptr<DenseMultilinearExtension<F>>> &perms)
{
    size_t num_vars = (*fxs[0]).num_vars();
    vector<shared_ptr<DenseMultilinearExtension<F>>> numerators;
    vector<shared_ptr<DenseMultilinearExtension<F>>> denominators;

    vector<shared_ptr<DenseMultilinearExtension<F>>> s_id = identity_permutation_mles<F>(num_vars, fxs.size());

    for (size_t l = 0; l < fxs.size(); ++l)
    {
        vector<F> numerator_evals;
        vector<F> denominator_evals;

        size_t num_points = (*fxs[l]).get_evaluations().size();
        numerator_evals.reserve(num_points);
        denominator_evals.reserve(num_points);

        for (size_t i = 0; i < num_points; ++i)
        {
            F f_ev = (*fxs[l]).get_evaluations()[i];
            F g_ev = (*gxs[l]).get_evaluations()[i];
            F s_id_ev = (*s_id[l]).get_evaluations()[i];
            F perm_ev = (*perms[l]).get_evaluations()[i];

            F numerator = f_ev + beta * s_id_ev + gamma;
            F denominator = g_ev + beta * perm_ev + gamma;

            numerator_evals.push_back(numerator);
            denominator_evals.push_back(denominator);
        }

        numerators.push_back(make_shared<DenseMultilinearExtension<F>>(num_vars, numerator_evals));
        denominators.push_back(make_shared<DenseMultilinearExtension<F>>(num_vars, denominator_evals));
    }
    return {numerators, denominators};
}

#endif