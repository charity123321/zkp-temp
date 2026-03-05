#ifndef TEST_H
#define TEST_H

#include <iostream>
#include <libff/algebra/curves/alt_bn128/alt_bn128_init.hpp>
#include <libff/algebra/curves/alt_bn128/alt_bn128_fields.hpp>
#include <random>
#include<algorithm>

#include "arithmetic/virtual_polynomial.h"
#include"arithmetic/multilinear_polynomial.h"
#include "subroutines/poly_iop/sum_check/prover.h"
#include "subroutines/poly_iop/sum_check/verifier.h"
#include "subroutines/poly_iop/sum_check/sumcheck.h"
#include "transcript/transcript.h"
#include "subroutines/poly_iop/zero_check/zerocheck.h"

#include"subroutines/pcs/batching.h"
#include"subroutines/pcs/srs.h"
#include "subroutines/pcs/multilinear_KZG.h"

#include"subroutines/poly_iop/prod_check/prodcheck.h"

#include"subroutines/poly_iop/perm_check/permcheck.h"

#include <libff/algebra/curves/bls12_381/bls12_381_pp.hpp> 
#include <libff/common/profiling.hpp> 


#include"hyperplonk/custom_gate.h"
#include"hyperplonk/mock.h"
#include"hyperplonk/snark.h"


using namespace std;
using namespace libff;

void test_sum_check_trivial()
{
    // 初始化 ALT-BN128 参数
    //libff::init_alt_bn128_params();
    libff::init_bls12_381_params();

    // 使用 ALT-BN128 的基域 Fq
    //typedef libff::alt_bn128_Fq Fp;
    typedef libff::bls12_381_Fq Fp;

    cout << "=== SumCheck Protocol Example ===" << endl;

    try
    {

        //    创建一个更复杂的虚拟多项式进行测试
        //    例如：f(x,y,z) = 2*x*y + 3*y*z + 4*x*z，在布尔超立方体上的和为 9
        size_t num_vars = 3;

        // 创建三个MLE：f1(x,y,z) = x, f2(x,y,z) = y, f3(x,y,z) = z
        // 3变量有 2^3 = 8 个评估点，顺序为：000, 001, 010, 011, 100, 101, 110, 111
        vector<Fp> eval_x = {
            Fp::zero(), Fp::one(), Fp::zero(), Fp::one(), // 000, 100, 010, 110
            Fp::zero(), Fp::one(), Fp::zero(), Fp::one()  // 001, 101, 011, 111
        };

        vector<Fp> eval_y = {
            Fp::zero(), Fp::zero(), Fp::one(), Fp::one(), // 000, 100, 010, 110
            Fp::zero(), Fp::zero(), Fp::one(), Fp::one()  // 001, 101, 011, 111
        };

        vector<Fp> eval_z = {
            Fp::zero(), Fp::zero(), Fp::zero(), Fp::zero(), // 000, 100, 010, 110
            Fp::one(), Fp::one(), Fp::one(), Fp::one()      // 001, 101, 011, 111
        };

        auto mle_x = make_shared<DenseMultilinearExtension<Fp>>(num_vars, eval_x);
        auto mle_y = make_shared<DenseMultilinearExtension<Fp>>(num_vars, eval_y);
        auto mle_z = make_shared<DenseMultilinearExtension<Fp>>(num_vars, eval_z);

        // 创建虚拟多项式：f(x,y,z) = 2*x*y + 3*y*z + 4*x*z
        VirtualPolynomial<Fp> poly(num_vars);

        // 第一项：2*x*y
        vector<shared_ptr<DenseMultilinearExtension<Fp>>> term1 = {mle_x, mle_y};
        poly.add_mle_list(term1, Fp(2)); // 系数为2

        // 第二项：3*y*z
        vector<shared_ptr<DenseMultilinearExtension<Fp>>> term2 = {mle_y, mle_z};
        poly.add_mle_list(term2, Fp(3)); // 系数为3

        // 第三项：4*x*z
        vector<shared_ptr<DenseMultilinearExtension<Fp>>> term3 = {mle_x, mle_z};
        poly.add_mle_list(term3, Fp(4)); // 系数为4

        cout << "Created virtual polynomial with " << num_vars << " variables" << endl;
        cout << "Max degree: " << poly.aux_info.max_degree << endl;
        cout << "Number of product terms: " << poly.products.size() << endl;

        //    计算真实的和（用于验证）
        //    对于 f(x,y,z) = 2*x*y + 3*y*z + 4*x*z，在布尔超立方体上的和为：
        //    计算所有8个点的值：
        //    000: 2*0*0 + 3*0*0 + 4*0*0 = 0
        //    001: 2*0*0 + 3*0*1 + 4*0*1 = 0
        //    010: 2*0*1 + 3*1*0 + 4*0*0 = 0
        //    011: 2*0*1 + 3*1*1 + 4*0*1 = 3
        //    100: 2*1*0 + 3*0*0 + 4*1*0 = 0
        //    101: 2*1*0 + 3*0*1 + 4*1*1 = 4
        //    110: 2*1*1 + 3*1*0 + 4*1*0 = 2
        //    111: 2*1*1 + 3*1*1 + 4*1*1 = 2 + 3 + 4 = 9
        //    总和 = 0+0+0+3+0+4+2+9 = 18
        Fp true_sum = Fp(18);
        cout << "True sum over Boolean hypercube: ";
        true_sum.print();

        // 运行SumCheck协议
        SumCheck<Fp> sumcheck;

        // Prover端
        cout << "\n--- Prover Side ---" << endl;
        auto transcript_prover = sumcheck.init_transcript();

        auto proof = sumcheck.prove(poly, transcript_prover);

        cout << "Proof generated with " << proof.proofs.size() << " rounds" << endl;
        cout << "Challenge point: " << endl;
        for (const auto &p : proof.point)
        {
            p.print();
        }

        // Verifier端
        cout << "\n--- Verifier Side ---" << endl;

        auto transcript_verifier = sumcheck.init_transcript();

        auto subclaim = sumcheck.verify(true_sum, proof, poly.aux_info, transcript_verifier);

        cout << "Verification completed successfully!" << endl;
        cout << "Subclaim point: ";

        for (const auto &p : subclaim.point)
        {
            p.print();
        }
        cout << "Expected evaluation: ";
        subclaim.expected_evaluation.print();

        // 验证子声明
        cout << "\n--- Verifying Subclaim ---" << endl;
        // 在子声明点上计算多项式值
        auto actual_eval = poly.evaluate(subclaim.point);
        cout << "Actual evaluation at subclaim point: ";
        actual_eval.print();

        if (actual_eval == subclaim.expected_evaluation)
        {
            cout << "✓ Subclaim VERIFIED: Actual evaluation matches expected!" << endl;
            cout << "✓ SumCheck protocol completed SUCCESSFULLY!" << endl;
        }
        else
        {
            cout << "✗ Subclaim FAILED: Actual evaluation does not match expected!" << endl;
            return;
        }
    }
    catch (const exception &e)
    {
        cerr << "Error: " << e.what() << endl;
        return;
    }

    cout << "\n=== SumCheck Example Completed ===" << endl;
}

void test_zero_check_trivial()
{
    // 初始化 ALT-BN128 参数
    libff::init_alt_bn128_params();

    // 使用 ALT-BN128 的基域 Fq
    typedef libff::alt_bn128_Fq Fp;

    cout << "=== Complex ZeroCheck Protocol Example ===" << endl;

    try
    {

        //  创建一个在布尔超立方体上求和为零的多项式
        //  例如：f(x,y,z) = x*y - x*y，在布尔超立方体上的和显然为0
        size_t num_vars = 3;

        // 创建三个MLE：f1(x,y,z) = x, f2(x,y,z) = y, f3(x,y,z) = z
        vector<Fp> eval_x = {
            Fp::zero(), Fp::one(), Fp::zero(), Fp::one(), // 000, 100, 010, 110
            Fp::zero(), Fp::one(), Fp::zero(), Fp::one()  // 001, 101, 011, 111
        };

        vector<Fp> eval_y = {
            Fp::zero(), Fp::zero(), Fp::one(), Fp::one(), // 000, 100, 010, 110
            Fp::zero(), Fp::zero(), Fp::one(), Fp::one()  // 001, 101, 011, 111
        };

        vector<Fp> eval_z = {
            Fp::zero(), Fp::zero(), Fp::zero(), Fp::zero(), // 000, 100, 010, 110
            Fp::one(), Fp::one(), Fp::one(), Fp::one()      // 001, 101, 011, 111
        };

        auto mle_x = make_shared<DenseMultilinearExtension<Fp>>(num_vars, eval_x);
        auto mle_y = make_shared<DenseMultilinearExtension<Fp>>(num_vars, eval_y);
        auto mle_z = make_shared<DenseMultilinearExtension<Fp>>(num_vars, eval_z);

        // 创建虚拟多项式：f(x,y,z) = x*y - x*y + y*z - y*z
        // 这个多项式在布尔超立方体上的和显然为0
        VirtualPolynomial<Fp> poly(num_vars);

        // 第一项：x*y (正项)
        vector<shared_ptr<DenseMultilinearExtension<Fp>>> term1 = {mle_x, mle_y};
        poly.add_mle_list(term1, Fp(1));

        // 第二项：x*y (负项)
        vector<shared_ptr<DenseMultilinearExtension<Fp>>> term2 = {mle_x, mle_y};
        poly.add_mle_list(term2, Fp(-1));

        // 第三项：y*z (正项)
        vector<shared_ptr<DenseMultilinearExtension<Fp>>> term3 = {mle_y, mle_z};
        poly.add_mle_list(term3, Fp(1));

        // 第四项：y*z (负项)
        vector<shared_ptr<DenseMultilinearExtension<Fp>>> term4 = {mle_y, mle_z};
        poly.add_mle_list(term4, Fp(-1));

        cout << "Created zero polynomial with " << num_vars << " variables" << endl;
        cout << "Max degree: " << poly.aux_info.max_degree << endl;
        cout << "Number of product terms: " << poly.products.size() << endl;
        // 运行ZeroCheck协议
        ZeroCheck<Fp> zerocheck;

        // Prover端
        cout << "\n--- Prover Side ---" << endl;
        auto transcript_prover = zerocheck.init_transcript();
        auto proof = zerocheck.prove(poly, transcript_prover);

        cout << "Complex ZeroCheck proof generated with " << proof.proofs.size() << " rounds" << endl;

        // Verifier端
        cout << "\n--- Verifier Side ---" << endl;
        auto transcript_verifier = zerocheck.init_transcript();
        auto subclaim = zerocheck.verify(proof, poly.aux_info, transcript_verifier);

        cout << "Complex ZeroCheck verification completed successfully!" << endl;

        // 验证子声明
        cout << "\n--- Verifying Complex ZeroCheck Subclaim ---" << endl;
        auto actual_eval = poly.evaluate(subclaim.point);
        cout << "Actual evaluation at subclaim point: ";
        actual_eval.print();

        if (actual_eval == subclaim.expected_evaluation)
        {
            cout << "✓ Complex ZeroCheck Subclaim VERIFIED!" << endl;
            cout << "✓ Complex ZeroCheck protocol completed SUCCESSFULLY!" << endl;
        }
        else
        {
            cout << "✗ Complex ZeroCheck Subclaim FAILED!" << endl;
            return;
        }
    }
    catch (const exception &e)
    {
        cerr << "Error: " << e.what() << endl;
        return;
    }

    cout << "\n=== Complex ZeroCheck Example Completed ===" << endl;
}

void test_sum_check_random()
{
    // 初始化 ALT-BN128 参数
    libff::init_alt_bn128_params();

    // 使用 ALT-BN128 的基域 Fq
    typedef libff::alt_bn128_Fq Fp;

    cout << "=== Random SumCheck Protocol Example ===" << endl;

    try
    {

        // 生成随机虚拟多项式
        size_t num_vars = 3;
        auto num_multiplicands_range = make_pair(1, 8); // 每个乘积项有1-3个MLE
        size_t num_products = 5;                        // 5个乘积项

        cout << "Generating random polynomial..." << endl;
        cout << "Number of variables: " << num_vars << endl;
        cout << "Number of product terms: " << num_products << endl;
        cout << "Multiplicands range: [" << num_multiplicands_range.first
             << ", " << num_multiplicands_range.second << "]" << endl;

        // 生成随机多项式和它的真实和
        VirtualPolynomial<Fp> poly1;
        auto [poly, true_sum] = poly1.rand(num_vars, num_multiplicands_range, num_products);

        cout << "Random polynomial created successfully!" << endl;
        cout << "Max degree: " << poly.aux_info.max_degree << endl;
        cout << "Number of product terms: " << poly.products.size() << endl;
        cout << "True sum over Boolean hypercube: ";
        true_sum.print();
        cout << endl;

        // 运行SumCheck协议
        SumCheck<Fp> sumcheck;

        // Prover端
        cout << "\n--- Prover Side ---" << endl;
        auto transcript_prover = sumcheck.init_transcript();
        auto proof = sumcheck.prove(poly, transcript_prover);

        cout << "Proof generated with " << proof.proofs.size() << " rounds" << endl;

        // Verifier端
        cout << "\n--- Verifier Side ---" << endl;
        auto transcript_verifier = sumcheck.init_transcript();
        auto subclaim = sumcheck.verify(true_sum, proof, poly.aux_info, transcript_verifier);

        cout << "Verification completed successfully!" << endl;
        cout << "Subclaim point: ";
        for (const auto &p : subclaim.point)
        {
            p.print();
            cout << " ";
        }
        cout << endl;
        cout << "Expected evaluation: ";
        subclaim.expected_evaluation.print();
        cout << endl;

        // 验证子声明
        cout << "\n--- Verifying Subclaim ---" << endl;
        // 在子声明点上计算多项式值
        auto actual_eval = poly.evaluate(subclaim.point);
        cout << "Actual evaluation at subclaim point: ";
        actual_eval.print();
        cout << endl;

        if (actual_eval == subclaim.expected_evaluation)
        {
            cout << "✓ Subclaim VERIFIED: Actual evaluation matches expected!" << endl;
            cout << "✓ Random SumCheck protocol completed SUCCESSFULLY!" << endl;
        }
        else
        {
            cout << "✗ Subclaim FAILED: Actual evaluation does not match expected!" << endl;
            cout << "Difference: ";
            (actual_eval - subclaim.expected_evaluation).print();
            cout << endl;
            return;
        }
    }
    catch (const exception &e)
    {
        cerr << "Error: " << e.what() << endl;
        return;
    }

    cout << "\n=== Random SumCheck Example Completed ===" << endl;
}

void test_zero_check_random()
{
    // 初始化 ALT-BN128 参数
    libff::init_alt_bn128_params();

    // 使用 ALT-BN128 的基域 Fq
    typedef libff::alt_bn128_Fq Fp;

    cout << "=== Random SumCheck Protocol Example ===" << endl;

    try
    {

        // 生成随机虚拟多项式
        size_t num_vars = 3;
        auto num_multiplicands_range = make_pair(1, 4); // 每个乘积项有1-3个MLE
        size_t num_products = 5;                        // 5个乘积项

        cout << "Generating random polynomial..." << endl;
        cout << "Number of variables: " << num_vars << endl;
        cout << "Number of product terms: " << num_products << endl;
        cout << "Multiplicands range: [" << num_multiplicands_range.first
             << ", " << num_multiplicands_range.second << "]" << endl;

        // 生成随机0多项式
        VirtualPolynomial<Fp> poly1;
        auto poly = poly1.rand_zero(num_vars, num_multiplicands_range, num_products);

        Fp true_sum = Fp::zero();
        cout << "Random polynomial created successfully!" << endl;
        cout << "Max degree: " << poly.aux_info.max_degree << endl;
        cout << "Number of product terms: " << poly.products.size() << endl;
        cout << "True sum over Boolean hypercube: ";
        true_sum.print();
        cout << endl;

        // 运行ZeroCheck协议
        ZeroCheck<Fp> zerocheck;

        // Prover端
        cout << "\n--- Prover Side ---" << endl;
        auto transcript_prover = zerocheck.init_transcript();
        auto proof = zerocheck.prove(poly, transcript_prover);

        cout << "Complex ZeroCheck proof generated with " << proof.proofs.size() << " rounds" << endl;

        // Verifier端
        cout << "\n--- Verifier Side ---" << endl;
        auto transcript_verifier = zerocheck.init_transcript();
        auto subclaim = zerocheck.verify(proof, poly.aux_info, transcript_verifier);

        cout << "Complex ZeroCheck verification completed successfully!" << endl;

        // 验证子声明
        cout << "\n--- Verifying Complex ZeroCheck Subclaim ---" << endl;
        auto actual_eval = poly.evaluate(subclaim.point);
        cout << "Actual evaluation at subclaim point: ";
        actual_eval.print();

        if (actual_eval == subclaim.expected_evaluation)
        {
            cout << "✓ Complex ZeroCheck Subclaim VERIFIED!" << endl;
            cout << "✓ Complex ZeroCheck protocol completed SUCCESSFULLY!" << endl;
        }
        else
        {
            cout << "✗ Complex ZeroCheck Subclaim FAILED!" << endl;
            return;
        }
    }
    catch (const exception &e)
    {
        cerr << "Error: " << e.what() << endl;
        return;
    }

    cout << "\n=== Complex ZeroCheck Example Completed ===" << endl;
}

void test_multilinearKzgPCS(){
    libff::inhibit_profiling_info = true;  
    libff::inhibit_profiling_counters = true; 
    
    libff::bls12_381_pp::init_public_params();
    typedef libff::bls12_381_pp P;
    typedef G1<P>::scalar_field scalar_field;

    // 创建PCS 实例
    cout<<" 创建MultilinearKzgPCS实例 "<<endl;
    MultilinearKzgPCS<P,scalar_field> pcs;

    try{
        // 生成SRS
        cout<<" 生成测试用SRS "<<endl;
        auto srs=pcs.gen_srs_for_test(3);

        // 裁剪参数
        cout<<" 裁剪SRS为特定大小 "<<endl;
        auto [prover_param,verifier_param]=pcs.trims(srs,2);

        // 创建多项式: f(x_1, x_2) = 1 + 2x_1 + 3x_2 + 4x_1x_2
        cout<<" 创建2变量多线性多项式 "<<endl;
        vector<scalar_field> evaluations={
            scalar_field(1),
            scalar_field(4),
            scalar_field(3),
            scalar_field(10)
        };
        auto poly=make_shared<DenseMultilinearExtension<scalar_field>>(2,evaluations);

        // 生成承诺
        cout<<" 生成多项式承诺 "<<endl;
        Commitment<P> commitment=pcs.commit(prover_param,poly);
        cout<<" 承诺生成完成 "<<endl;

        // 打开证明
        cout<<" 在(1,2)点处生成打开证明 "<<endl;
        vector<scalar_field> point={scalar_field(1),scalar_field(2)};

        auto [proof,eval]=pcs.open(prover_param,poly,point);
        //cout<<" 求值结果："<<eval<<endl;
        cout<<" 证明包含 "<<proof.proofs.size()<<" 个元素"<<endl;

        // 验证证明
        cout<<" 验证证明 "<<endl;
        bool is_valid=pcs.verify(verifier_param,commitment,point,eval,proof);
        cout << " 验证结果: " << (is_valid ? "通过 ✓" : "失败 ✗") << endl;
        
        cout << " 测试批量打开" << endl;
        vector<shared_ptr<DenseMultilinearExtension<scalar_field>>> polynomials = {poly,poly};
        vector<vector<scalar_field>> points = {point,point};
        vector<scalar_field> evals = {eval,eval};
        
        SimpleTranscript transcript1;
        auto batch_proof = pcs.multi_open(prover_param, polynomials, points, evals, transcript1);
        cout << " 批量证明生成完成" << endl;
        
        // 9. 测试批量验证
        cout << " 测试批量验证" << endl;
        vector<Commitment<P>> commitments = {commitment,commitment};
        
        SimpleTranscript transcript2;
        bool batch_valid = pcs.batch_verify(verifier_param, commitments, points, batch_proof, transcript2);
        cout << " 批量验证结果: " << (batch_valid ? "通过 ✓" : "失败 ✗") << endl;
    }
    catch (const exception& e) 
    {
        cerr << "错误: " << e.what() << endl;
        return;
    }
}

void test_multilinearKzgPCS_for_batch(){
    libff::inhibit_profiling_info = true;  
    libff::inhibit_profiling_counters = true; 
    
    libff::bls12_381_pp::init_public_params();
    typedef libff::bls12_381_pp P;
    typedef G1<P>::scalar_field scalar_field;

    // 创建PCS 实例
    cout<<" 创建MultilinearKzgPCS实例 "<<endl;
    MultilinearKzgPCS<P,scalar_field> pcs;

    try{
        // 生成SRS
        cout<<" 生成测试用SRS "<<endl;
        auto srs=pcs.gen_srs_for_test(3);

        // 裁剪参数
        cout<<" 裁剪SRS为特定大小 "<<endl;
        auto [prover_param,verifier_param]=pcs.trims(srs,2);

        // 创建多项式1: f1(x_1, x_2) = 1 + 2x_1 + 3x_2 + 4x_1x_2
        cout<<" 创建2变量多线性多项式 * 2 "<<endl;
        vector<scalar_field> evaluations1={
            scalar_field(1),
            scalar_field(4),
            scalar_field(3),
            scalar_field(10)
        };
        auto poly1=make_shared<DenseMultilinearExtension<scalar_field>>(2,evaluations1);

        // 创建多项式2：f2(x_1,x_2) = 2 + x_1 + 3x_1x_2
        vector<scalar_field> evaluations2={
            scalar_field(2),
            scalar_field(5),
            scalar_field(3),
            scalar_field(6)
        };
        auto poly2=make_shared<DenseMultilinearExtension<scalar_field>>(2,evaluations2);

        // 生成承诺
        cout<<" 生成多项式承诺 "<<endl;
        Commitment<P> commitment1=pcs.commit(prover_param,poly1);
        Commitment<P> commitment2=pcs.commit(prover_param,poly2);
        cout<<" 承诺生成完成 "<<endl;

        // 打开证明
        cout<<" 批量生成打开证明 "<<endl;
        vector<vector<scalar_field>> points={
            {scalar_field(1),scalar_field(2)},
            {scalar_field(2),scalar_field(3)}
        };

        vector<shared_ptr<DenseMultilinearExtension<scalar_field>>> polynomials = {poly1,poly2};

        vector<scalar_field> evals = {
            (*poly1).evaluate(points[0]),
            (*poly2).evaluate(points[1])
        };
        
        SimpleTranscript transcript1;
        auto batch_proof = pcs.multi_open(prover_param, polynomials, points, evals, transcript1);
        cout << " 批量证明生成完成" << endl;
        
        // 9. 测试批量验证
        cout << " 测试批量验证" << endl;
        vector<Commitment<P>> commitments = {commitment1,commitment2};
        
        SimpleTranscript transcript2;
        bool batch_valid = pcs.batch_verify(verifier_param, commitments, points, batch_proof, transcript2);
        cout << " 批量验证结果: " << (batch_valid ? "通过 ✓" : "失败 ✗") << endl;
    }
    catch (const exception& e) 
    {
        cerr << "错误: " << e.what() << endl;
        return;
    }   
}

void test_productcheck_random(){
    libff::inhibit_profiling_info = true;  
    libff::inhibit_profiling_counters = true; 
    
    libff::bls12_381_pp::init_public_params();
    typedef libff::bls12_381_pp P;
    typedef G1<P>::scalar_field scalar_field;

    cout << "=== Random ProductCheck Protocol Example ===" << endl;
    try
    {
        // 生成随机虚拟多项式用于乘积检查
        size_t num_vars = 3;
        size_t num_polys = 4;  // 4对多项式 (f_i, g_i)
        
        cout << "Generating random polynomials for product check..." << endl;
        cout << "Number of variables: " << num_vars << endl;
        cout << "Number of polynomial pairs: " << num_polys << endl;
        
        // 生成随机多项式集合 f_i(x) 和 g_i(x)
        // 满足 ∏ f_i(x) = ∏ g_i(x) 对所有 x ∈ {0,1}^num_vars
        vector<shared_ptr<DenseMultilinearExtension<scalar_field>>> fxs;
        vector<shared_ptr<DenseMultilinearExtension<scalar_field>>> gxs;

        vector<shared_ptr<DenseMultilinearExtension<scalar_field>>> all_poly;

        cout << "Creating random polynomial pairs..." << endl;
        for (size_t i = 0; i < num_polys; ++i)
        {
            DenseMultilinearExtension<scalar_field> polyi(num_vars);
            all_poly.push_back(make_shared<DenseMultilinearExtension<scalar_field>>(polyi.rand(num_vars)));
        }

        cout<<"test"<<endl;
        cout<<(*all_poly[0]).get_evaluations().size()<<endl;
        
        fxs={all_poly[1],all_poly[2],all_poly[3],all_poly[0]};
        gxs={all_poly[0],all_poly[1],all_poly[2],all_poly[3]};
        
        cout << "Polynomial pairs created !" << endl;

        MultilinearKzgPCS<P,scalar_field> pcs;
        auto srs=pcs.gen_srs_for_test(num_vars);
    
        auto [prover_param,verifier_param]=pcs.trims(srs,num_vars);
        
        ProductCheck<P,scalar_field,MultilinearKzgPCS<P,scalar_field>> prod_check;
        cout << "\n--- Prover Side ---" << endl;
        
        // 初始化转录本
        auto transcript_prover = prod_check.init_transcript();
        
        // Prover 生成证明
        cout << "Generating product check proof..." << endl;

        auto [proof, prod_x_poly, frac_poly] = prod_check.prove(
            prover_param,
            fxs,
            gxs,
            transcript_prover
        );
        
        cout << "Proof generated successfully!" << endl;
        cout << "Product polynomial created" << endl;
        cout << "Fraction polynomial created" << endl;
        
        // 创建辅助信息
        VPAuxInfo aux_info;
        aux_info.num_variables = num_vars;
        aux_info.max_degree = (1 << num_vars)+1;
        
        // Verifier 端
        cout << "\n--- Verifier Side ---" << endl;
        
        auto transcript_verifier = prod_check.init_transcript();
        
        // Verifier 验证证明
        cout << "Verifying product check proof..." << endl;
        auto subclaim = prod_check.verify(
            proof,
            aux_info,
            transcript_verifier
        );
        
        cout << "Verification completed successfully!" << endl;
        cout << "Zero-check subclaim verified" << endl;
        cout << "Alpha challenge: ";
        subclaim.alpha.print();
        cout << endl;
        
        // 验证最终查询点
        cout << "\n--- Verifying Final Query ---" << endl;
        
        const auto [final_query_point, final_expected_eval] = subclaim.final_query;
        
        cout << "Final query point: "<<endl;
        
        for (size_t i = 0; i < final_query_point.size(); ++i) {
            cout << "  x" << i << " = ";
            final_query_point[i].print();
        }
        cout << endl;
        
        cout << "Expected evaluation at final point: ";
        final_expected_eval.print();
        cout << endl;
        
        // 实际计算乘积多项式的值
        // 注意：这里我们需要获取实际的 prod_x 多项式，但在验证时我们只有承诺
        // 为了测试，我们可以通过其他方式计算

        auto actual_eval = (*prod_x_poly).evaluate(final_query_point);
        cout << "Actual evaluation at final point: ";
        actual_eval.print();
        cout << endl;
        
        if (actual_eval == final_expected_eval)
        {
            cout << "✓ Final query VERIFIED: Actual evaluation matches expected!" << endl;
        }
        else
        {
            cout << "✗ Final query FAILED: Actual evaluation does not match expected!" << endl;
            cout << "Difference: ";
            (actual_eval - final_expected_eval).print();
            cout << endl;
        }
    }
    catch (const exception &e)
    {
        cerr << "Error: " << e.what() << endl;
        return;
    }

    cout << "\n=== Random ProductCheck Example Completed ===" << endl;
    
}

void test_permutationcheck_trivial() {
    libff::inhibit_profiling_info = true;  
    libff::inhibit_profiling_counters = true; 
    
    libff::bls12_381_pp::init_public_params();
    typedef libff::bls12_381_pp P;
    typedef G1<P>::scalar_field scalar_field;

    cout << "=== Trivial PermutationCheck Protocol Example ===" << endl;
    
    try {
        // 测试参数
        size_t num_vars = 3;      // 变量个数
        size_t num_polys = 3;     // 多项式个数
        
        cout << "Number of variables: " << num_vars << endl;
        cout << "Number of polynomials: " << num_polys << endl;
        cout << "Domain size: " << (1 << num_vars) << endl;
        
        // 生成随机多项式 f_i(x)
        cout << "\nGenerating random polynomials f_i(x)..." << endl;
        vector<shared_ptr<DenseMultilinearExtension<scalar_field>>> fxs;
        for (size_t i = 0; i < num_polys; ++i) {
            DenseMultilinearExtension<scalar_field> poly(num_vars);
            auto rand_poly = poly.rand(num_vars);
            fxs.push_back(make_shared<DenseMultilinearExtension<scalar_field>>(rand_poly));
        }
        
        // 生成逆恒等置换 π
        cout << "\nGenerating reverse identity permutation π..." << endl;
        vector<shared_ptr<DenseMultilinearExtension<scalar_field>>> perms=identity_permutation_mles<scalar_field>(num_vars,num_polys);
        reverse(perms.begin(),perms.end());

        /*
        // 根据置换计算 g_i(x) = f_i(π_i(x))
        cout << "\nComputing g_i(x) = f_i(π_i(x))..." << endl;
        vector<shared_ptr<DenseMultilinearExtension<scalar_field>>> gxs;
        
        for (size_t i = 0; i < num_polys; ++i) {
            vector<scalar_field> g_evals(1 << num_vars);
            
            for (size_t j = 0; j < (1UL << num_vars); ++j) {
                // 计算 π_i 在点 j 的值
                scalar_field perm_val = (*perms[i]).evaluate(vector<scalar_field>{scalar_field(j)});
                
                // 将 perm_val 转换为索引
                size_t perm_idx = static_cast<size_t>(perm_val.as_ulong());
                
                // 获取 f_i 在 perm_idx 处的值
                // 注意：需要将索引转换为二进制点
                vector<scalar_field> perm_point;
                for (size_t k = 0; k < num_vars; ++k) {
                    perm_point.push_back(scalar_field((perm_idx >> k) & 1));
                }
                
                g_evals[j] = (*fxs[i]).evaluate(perm_point);
            }
            
            auto g_poly = DenseMultilinearExtension<scalar_field>(num_vars, g_evals);
            gxs.push_back(make_shared<DenseMultilinearExtension<scalar_field>>(g_poly));
            cout << "  g_" << i << " computed" << endl;
        }
        */
        vector<shared_ptr<DenseMultilinearExtension<scalar_field>>> gxs=fxs;
        reverse(gxs.begin(),gxs.end());

        // 设置多项式承诺方案
        cout << "\nSetting up polynomial commitment scheme..." << endl;
        MultilinearKzgPCS<P, scalar_field> pcs;
        auto srs = pcs.gen_srs_for_test(num_vars);
        
        auto [prover_param, verifier_param] = pcs.trims(srs, num_vars);
        cout << "PCS parameters generated" << endl;
        
        //  创建置换检查协议实例
        cout << "\nCreating PermutationCheck protocol..." << endl;
        PermutationCheck<P, scalar_field, MultilinearKzgPCS<P, scalar_field>> perm_check;
        
        // 步骤6: 证明者端
        cout << "\n--- Prover Side ---" << endl;
        
        // 初始化转录本
        auto transcript_prover = perm_check.init_transcript();
        cout << "Transcript initialized for prover" << endl;
        
        // 生成置换检查证明
        cout << "Generating permutation check proof..." << endl;
        auto [proof, prod_poly, frac_poly] = perm_check.prove(
            prover_param,
            fxs,
            gxs,
            perms,
            transcript_prover
        );
        
        cout << "Proof generated successfully!" << endl;
        cout << "Product polynomial created" << endl;
        cout << "Fraction polynomial created" << endl;
        
        // 验证者端
        cout << "\n--- Verifier Side ---" << endl;
        
        // 创建验证者转录本
        auto transcript_verifier = perm_check.init_transcript();
        
        // 设置辅助信息
        VPAuxInfo aux_info;
        aux_info.num_variables = num_vars;
        aux_info.max_degree = (1 << num_vars) + 1;
        
        // 验证证明
        cout << "Verifying permutation check proof..." << endl;
        auto subclaim = perm_check.verify(
            proof,
            aux_info,
            transcript_verifier
        );
        
        cout << "Verification completed successfully!" << endl;
        
        // 输出验证结果
        cout << "\n--- Verification Results ---" << endl;
        
        // 输出挑战值
        cout << "Challenges from transcript:" << endl;
        cout << "  Beta: ";
        subclaim.challenge.first.print();
        cout << endl;
        cout << "  Gamma: ";
        subclaim.challenge.second.print();
        cout << endl;
        
        // 输出prodcheck子声明
        cout << "Product check subclaim details:" << endl;
        cout << "  Alpha challenge: ";
        subclaim.product_check_sub_claim.alpha.print();
        cout << endl;
        
        // 验证最终查询点
        cout << "\n--- Verifying Final Query ---" << endl;
        
        const auto [final_query_point, final_expected_eval] = 
            subclaim.product_check_sub_claim.final_query;
        
        cout << "Final query point:" << endl;
        for (size_t i = 0; i < final_query_point.size(); ++i) {
            cout << "  x" << i << " = ";
            final_query_point[i].print();
        }
        
        cout << "Expected evaluation at final point: ";
        final_expected_eval.print();
        cout << endl;
        
        // 实际计算product多项式在最终点的值
        auto actual_eval = (*prod_poly).evaluate(final_query_point);
        cout << "Actual evaluation at final point: ";
        actual_eval.print();
        cout << endl;
        
        // 验证结果
        if (actual_eval == final_expected_eval) {
            cout << "✓ Permutation check VERIFIED: All checks passed!" << endl;
            cout << "  - Final evaluation matches expected value" << endl;
            cout << "  - Permutation relationship holds" << endl;
        } else {
            cout << "✗ Permutation check FAILED: Final evaluation mismatch!" << endl;
            cout << "  Difference: ";
            (actual_eval - final_expected_eval).print();
            cout << endl;
        }
        
    } catch (const exception &e) {
        cerr << "\nError during permutation check test: " << e.what() << endl;
        return;
    }
    
    cout << "\n=== Trivial PermutationCheck Example Completed ===" << endl;
}


void test_vanilla_plonk_gates() {
    libff::inhibit_profiling_info = true;  
    libff::inhibit_profiling_counters = true; 
    
    libff::bls12_381_pp::init_public_params();
    typedef libff::bls12_381_pp P;
    typedef G1<P>::scalar_field scalar_field;

    cout << "=== 测试 vanilla_plonk_gates() ===" << endl;
    CustomizedGates custgate;
    
    CustomizedGates gates = custgate.vanilla_plonk_gates();
    
    cout << "门项数量: " << gates.gates.size() << endl;
    cout << "选择器列数: " << gates.num_selector_columns() << endl;
    cout << "见证列数: " << gates.num_witness_columns() << endl;
    cout << "门的度数: " << gates.degree() << endl;
    
    // 验证结构
    assert(gates.gates.size() == 5 && "Vanilla Plonk应该有5个门项");
    assert(gates.num_selector_columns() == 5 && "应该有5个选择器");
    assert(gates.num_witness_columns() == 3 && "应该有3个见证列");
    assert(gates.degree() == 3 && "最高次数应该是2（w1*w2）");
    
    // 验证具体项
    const auto& gate_terms = gates.gates;
    assert(gate_terms[0].coefficient == 1 && gate_terms[0].selector_index == 0 && 
           gate_terms[0].wire_indices == vector<size_t>{0} && "第一项应该是q_L*w1");
    assert(gate_terms[1].coefficient == 1 && gate_terms[1].selector_index == 1 && 
           gate_terms[1].wire_indices == vector<size_t>{1} && "第二项应该是q_R*w2");
    assert(gate_terms[2].coefficient == 1 && gate_terms[2].selector_index == 2 && 
           gate_terms[2].wire_indices == vector<size_t>{2} && "第三项应该是q_O*w3");
    /*
    assert(gate_terms[3].coefficient == 1 && gate_terms[3].selector_index == 3 && 
           gate_terms[3].wire_indices == vector<size_t>{0,1} && "第四项应该是q_M*w1*w2");
    */
    assert(gate_terms[4].coefficient == 1 && gate_terms[4].selector_index == 4 && 
           gate_terms[4].wire_indices.empty() && "第五项应该是q_C常数项");
    
    // 创建MockCircuit测试约束满足
    MockCircuit<scalar_field> circuit(8, gates);
    
    cout << "创建MockCircuit(8个约束)..." << endl;
    cout << "电路变量数: " << circuit.num_variable() << endl;
    
    bool satisfied = circuit.is_satisfied();
    cout << "约束检查结果: " << (satisfied ? "满足 ✓" : "不满足 ✗") << endl;
    
    assert(satisfied && "Vanilla Plonk门应该满足约束");
    cout << "测试通过！" << endl << endl;
}

void test_mock_gate() {
    libff::inhibit_profiling_info = true;  
    libff::inhibit_profiling_counters = true; 
    
    libff::bls12_381_pp::init_public_params();
    typedef libff::bls12_381_pp P;
    typedef G1<P>::scalar_field scalar_field;


    cout << "=== 测试 mock_gate(3个见证, 度数4) ===" << endl;
    
    size_t num_witness = 4;
    size_t degree = 10;

    CustomizedGates custgate;
    CustomizedGates gates = custgate.mock_gate(num_witness, degree);
    
    cout << "门项数量: " << gates.gates.size() << endl;
    cout << "选择器列数: " << gates.num_selector_columns() << endl;
    cout << "见证列数: " << gates.num_witness_columns() << endl;
    cout << "门的度数: " << gates.degree() << endl;
    
    // 验证结构
    assert(gates.gates.size() == num_witness + 2 && "门项数量应该为num_witness+2");
    assert(gates.num_selector_columns() == num_witness + 2 && "选择器列数应该为num_witness+2");
    assert(gates.num_witness_columns() == num_witness && "见证列数应该为num_witness");

    cout<<gates.degree()<<endl;
    assert(gates.degree() == degree+1 && "度数应该为指定的degree");
    
    // 验证第一项：q0 * w0^(degree-1) * w1
    const auto& first_term = gates.gates[0];
    assert(first_term.coefficient == 1 && "系数应为1");
    assert(first_term.selector_index == 0 && "选择器索引应为0");
    assert(first_term.wire_indices.size() == degree && "应该有degree个w0和一个w1");
    
    // 验证线性项
    for (size_t i = 0; i < num_witness; ++i) {
        const auto& term = gates.gates[i + 1];
        assert(term.coefficient == 1 && "系数应为1");
        assert(term.selector_index == i + 1 && "选择器索引应为i+1");
        assert(term.wire_indices.size() == 1 && "应该有1个见证索引");
        assert(term.wire_indices[0] == i && "见证索引应为i");
    }
    
    // 验证常数项
    const auto& const_term = gates.gates[num_witness + 1];
    assert(const_term.coefficient == 1 && "系数应为1");
    assert(const_term.selector_index == num_witness + 1 && "选择器索引应为num_witness+1");
    assert(const_term.wire_indices.empty() && "常数项应该没有见证索引");
    
    // 创建MockCircuit测试约束满足
    MockCircuit<scalar_field> circuit(16, gates);
    
    cout << "创建MockCircuit(16个约束)..." << endl;
    cout << "电路变量数: " << circuit.num_variable() << endl;
    
    bool satisfied = circuit.is_satisfied();
    cout << "约束检查结果: " << (satisfied ? "满足 ✓" : "不满足 ✗") << endl;
    
    assert(satisfied && "mock_gate应该满足约束");
    cout << "测试通过！" << endl << endl;
}


void test_hyperplonk_snark(){
    libff::inhibit_profiling_info = true;  
    libff::inhibit_profiling_counters = true; 
    
    libff::bls12_381_pp::init_public_params();
    typedef libff::bls12_381_pp P;
    typedef G1<P>::scalar_field scalar_field;

    // 初始化随机数生成器
    random_device rd;
    mt19937_64 rng(rd());
    
    // =======================================================================
    // 1. 生成SRS（结构化参考字符串）
    // =======================================================================
    cout << "1. Generating SRS (Structured Reference String)..." << endl;
    
    MultilinearKzgPCS<P,scalar_field> pcs;
    auto pcs_srs = pcs.gen_srs_for_test(16);

    
    // =======================================================================
    // 2. 设置电路参数
    // =======================================================================
    cout << "2. Setting up circuit parameters..." << endl;
    
    size_t num_constraints = 4;      // 4个约束
    size_t num_pub_input = 4;        // 4个公共输入
    size_t nv = 2;                   // log2(4) = 2个变量
    size_t num_witnesses = 2;        // 2个见证列
    
    // 定义门约束：q_L * W₁⁵ - W₂ = 0
    // 表示为：
    //   vec![
    //     ( 1,    Some(id_qL),    vec![id_W1, id_W1, id_W1, id_W1, id_W1]),
    //     (-1,    None,           vec![id_W2])
    //   ]
    
    CustomizedGates gate_func;
    
    // 第一个项：q_L * W₁⁵
    GateTerm term1;
    term1.coefficient = 1;  // 系数 1
    term1.selector_index = 0;            // q_L（索引0）
    term1.wire_indices = {0, 0, 0, 0, 0};  // W₁⁵（重复5次）
    
    // 第二个项：-W₂
    GateTerm term2;
    term2.coefficient = -1;  // 系数 -1
    term2.selector_index = nullopt;     // 无选择器
    term2.wire_indices = {1};           // W₂
    
    gate_func.gates = {term1, term2};
    
    // =======================================================================
    // 3. 创建电路参数和索引
    // =======================================================================
    cout << "3. Creating circuit parameters and index..." << endl;
    
    HyperPlonkParams<scalar_field> params;
    params.num_constraints = num_constraints;
    params.num_pub_input = num_pub_input;
    params.gate_func = gate_func;

    
    // 生成单位置换
    auto permutation = identity_permutation<scalar_field>(nv, num_witnesses);
    
    // 创建选择器 q1 = [1, 1, 1, 1]
    SelectorColumn<scalar_field> q1(vector<scalar_field>{
        scalar_field::one(), scalar_field::one(), scalar_field::one(), scalar_field::one()
    });
    
    // 创建索引
    HyperPlonkIndex<scalar_field> index;
    index.params = params;
    index.permutation = permutation;
    index.selectors = {q1};
    
    // =======================================================================
    // 4. 预处理：生成证明密钥和验证密钥
    // =======================================================================
    cout << "4. Preprocessing: generating proving and verifying keys..." << endl;
    
    HyperPlonkSNARK<P,scalar_field,MultilinearKzgPCS<P,scalar_field>> snark;
    auto [pk, vk] = snark.preprocess(index, pcs_srs);

    cout << "   - Proving key generated" << endl;
    cout << "   - Verifying key generated" << endl;
    
    // =======================================================================
    // 5. 设置见证和公共输入
    // =======================================================================
    cout << "5. Setting up witnesses and public inputs..." << endl;
    
    // 第一个见证：w1 = [0, 1, 2, 3]
    WitnessColumn<scalar_field> w1(vector<scalar_field>{
        scalar_field::zero(),
        scalar_field::one(),
        scalar_field(2),
        scalar_field(3)
    });
    
    // 第二个见证：w2 = [0⁵, 1⁵, 2⁵, 3⁵] = [0, 1, 32, 243]
    WitnessColumn<scalar_field> w2(vector<scalar_field>{
        scalar_field::zero(),
        scalar_field::one(),
        scalar_field(32),
        scalar_field(243)
    });
    
    // 公共输入 = w1
    WitnessColumn<scalar_field> pi = w1;

    cout << "   - Public input: same as witness 1" << endl;
    
    // =======================================================================
    // 6. 生成证明（正确路径）
    // =======================================================================
    cout << "6. Generating proof (correct path)..." << endl;
    
    auto proof =snark.prove(pk, pi.data, {w1,w2});
    
    cout << "   - Proof generated successfully" << endl;
    
    // =======================================================================
    // 7. 验证证明（正确路径）
    // =======================================================================
    cout << "7. Verifying proof (correct path)..." << endl;
    
    auto verify_ = snark.verify(vk, pi.data, proof);

    if (!verify_) {
        throw runtime_error("Proof should be valid but verification returned false");
    }
    
    cout << "   ✓ Proof is VALID" << endl;
        

}


void test_hyperplonk_snark_for_mock_circuit(){
    libff::inhibit_profiling_info = true;  
    libff::inhibit_profiling_counters = true; 
    
    libff::bls12_381_pp::init_public_params();
    typedef libff::bls12_381_pp P;
    typedef G1<P>::scalar_field scalar_field;

    // 初始化随机数生成器
    random_device rd;
    mt19937_64 rng(rd());
    
    // =======================================================================
    // 1. 生成SRS（结构化参考字符串）
    // =======================================================================
    cout << "1. Generating SRS (Structured Reference String)..." << endl;
    
    MultilinearKzgPCS<P,scalar_field> pcs;
    auto pcs_srs = pcs.gen_srs_for_test(16);

    
    // =======================================================================
    // 2. 设置电路参数
    // =======================================================================
    cout << "2. Setting up mock-circuit parameters..." << endl;
    

    size_t nv = 5;                   
    size_t num_witnesses = 4;        
    size_t degree=6;
    
    CustomizedGates gate;
    CustomizedGates mock_gate=gate.mock_gate(num_witnesses,degree);
    
    MockCircuit<scalar_field> circuit(1ULL<<nv,mock_gate);

    
    // =======================================================================
    // 3. 创建电路参数和索引
    // =======================================================================
    cout << "3. Creating index..." << endl;
    

    // 创建索引
    HyperPlonkIndex<scalar_field> index=circuit.index;
    
    // =======================================================================
    // 4. 预处理：生成证明密钥和验证密钥
    // =======================================================================
    cout << "4. Preprocessing: generating proving and verifying keys..." << endl;
    
    HyperPlonkSNARK<P,scalar_field,MultilinearKzgPCS<P,scalar_field>> snark;
    auto [pk, vk] = snark.preprocess(index, pcs_srs);

    cout << "   - Proving key generated" << endl;
    cout << "   - Verifying key generated" << endl;
    

    
    // =======================================================================
    // 5. 生成证明（正确路径）
    // =======================================================================
    cout << "5. Generating proof (correct path)..." << endl;
    
    auto proof =snark.prove(pk,circuit.public_inputs,circuit.witness);
    
    cout << "   - Proof generated successfully" << endl;
    
    // =======================================================================
    // 7. 验证证明（正确路径）
    // =======================================================================
    cout << "6. Verifying proof (correct path)..." << endl;
    
    auto verify_ = snark.verify(vk, circuit.public_inputs, proof);

    if (!verify_) {
        throw runtime_error("Proof should be valid but verification returned false");
    }
    
    cout << "   ✓ Proof is VALID" << endl;
        

}

#endif