#include <iostream>
#include <libff/algebra/curves/alt_bn128/alt_bn128_init.hpp>
#include <libff/algebra/curves/alt_bn128/alt_bn128_fields.hpp>
#include <random>
#include "arithmetic/virtual_polynomial.h"
#include "subroutines/poly_iop/sum_check/prover.h"
#include "subroutines/poly_iop/sum_check/verifier.h"
#include "subroutines/poly_iop/sum_check/sumcheck.h"
#include "transcript/transcript.h"

using namespace std;
using namespace libff;

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

void test_sum_check(){
        // 初始化 ALT-BN128 参数
    libff::init_alt_bn128_params();
    
    // 使用 ALT-BN128 的基域 Fq
    typedef libff::alt_bn128_Fq Fp;

    
    cout << "=== SumCheck Protocol Example ===" << endl;

    /*
    try {
        // 1. 创建随机数生成器
        std::random_device rd;
        std::mt19937_64 rng(rd());
        
        // 2. 创建一个简单的虚拟多项式进行测试
        //    例如：f(x,y) = x*y，在布尔超立方体上的和为 1
        size_t num_vars = 2;
        
        // 创建两个MLE：f1(x,y) = x, f2(x,y) = y
        vector<Fp> eval_f1 = {Fp::zero(), Fp::zero(), Fp::zero(), Fp::one()}; // x
        vector<Fp> eval_f2 = {Fp::zero(), Fp::zero(), Fp::one(), Fp::one()};  // y

        auto mle1 = std::make_shared<DenseMultilinearExtension<Fp>>(num_vars, eval_f1);
        auto mle2 = std::make_shared<DenseMultilinearExtension<Fp>>(num_vars, eval_f2);
        
        // 创建虚拟多项式：f(x,y) = x * y
        VirtualPolynomial<Fp> poly(num_vars);
        vector<shared_ptr<DenseMultilinearExtension<Fp>>> mle_list = {mle1, mle2};
        poly.add_mle_list(mle_list, Fp::one());
        
        cout << "Created virtual polynomial with " << num_vars << " variables" << endl;
        cout << "Max degree: " << poly.aux_info.max_degree << endl;
        
        // 3. 计算真实的和（用于验证）
        //    对于 f(x,y) = x*y，在布尔超立方体 {(0,0),(0,1),(1,0),(1,1)} 上的和为：
        //    f(0,0)=0, f(0,1)=0, f(1,0)=0, f(1,1)=1，总和 = 1
        Fp true_sum = Fp::one();
        cout << "True sum over Boolean hypercube: ";
        true_sum.print();
        
    */
    try {
        // 1. 创建随机数生成器
        std::random_device rd;
        std::mt19937_64 rng(rd());
        
        // 2. 创建一个更复杂的虚拟多项式进行测试
        //    例如：f(x,y,z) = 2*x*y + 3*y*z + 4*x*z，在布尔超立方体上的和为 9
        size_t num_vars = 3;
        
        // 创建三个MLE：f1(x,y,z) = x, f2(x,y,z) = y, f3(x,y,z) = z
        // 3变量有 2^3 = 8 个评估点，顺序为：000, 001, 010, 011, 100, 101, 110, 111
        vector <Fp> eval_x = {
            Fp::zero(), Fp::one(),  Fp::zero(), Fp::one(),   // 000, 100, 010, 110
            Fp::zero(), Fp::one(),  Fp::zero(), Fp::one()    // 001, 101, 011, 111
        };

        vector <Fp> eval_y = {
            Fp::zero(), Fp::zero(), Fp::one(),  Fp::one(),   // 000, 100, 010, 110
            Fp::zero(), Fp::zero(), Fp::one(),  Fp::one()    // 001, 101, 011, 111
        };

        vector <Fp> eval_z = {
            Fp::zero(), Fp::zero(), Fp::zero(), Fp::zero(),  // 000, 100, 010, 110
            Fp::one(),  Fp::one(),  Fp::one(),  Fp::one()    // 001, 101, 011, 111
        };

        auto mle_x = std::make_shared<DenseMultilinearExtension<Fp>>(num_vars, eval_x);
        auto mle_y = std::make_shared<DenseMultilinearExtension<Fp>>(num_vars, eval_y);
        auto mle_z = std::make_shared<DenseMultilinearExtension<Fp>>(num_vars, eval_z);
        
        // 创建虚拟多项式：f(x,y,z) = 2*x*y + 3*y*z + 4*x*z
        VirtualPolynomial<Fp> poly(num_vars);
        
        // 第一项：2*x*y
        vector<shared_ptr<DenseMultilinearExtension<Fp>>> term1 = {mle_x, mle_y};
        poly.add_mle_list(term1, Fp(2));  // 系数为2
        
        // 第二项：3*y*z  
        vector<shared_ptr<DenseMultilinearExtension<Fp>>> term2 = {mle_y, mle_z};
        poly.add_mle_list(term2, Fp(3));  // 系数为3
        
        // 第三项：4*x*z
        vector<shared_ptr<DenseMultilinearExtension<Fp>>> term3 = {mle_x, mle_z};
        poly.add_mle_list(term3, Fp(4));  // 系数为4
            
        cout << "Created virtual polynomial with " << num_vars << " variables" << endl;
        cout << "Max degree: " << poly.aux_info.max_degree << endl;
        cout << "Number of product terms: " << poly.products.size() << endl;
        
        // 3. 计算真实的和（用于验证）
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
        
        // 4. 运行SumCheck协议
        SumCheck<Fp> sumcheck;
        
        // Prover端
        cout << "\n--- Prover Side ---" << endl;
        auto transcript_prover = sumcheck.init_transcript();


        auto proof = sumcheck.prove(poly, transcript_prover);
        
        cout << "Proof generated with " << proof.proofs.size() << " rounds" << endl;
        cout << "Challenge point: "<<endl;
        for (const auto& p : proof.point) {
            p.print();
        }
        
        // 5. Verifier端
        cout << "\n--- Verifier Side ---" << endl;

        auto transcript_verifier = sumcheck.init_transcript();
        
        auto subclaim = sumcheck.verify(true_sum, proof, poly.aux_info, transcript_verifier);
        

        cout << "Verification completed successfully!" << endl;
        cout << "Subclaim point: ";
        
        for (const auto& p : subclaim.point) {
            p.print();
        }
        cout << "Expected evaluation: ";
        subclaim.expected_evaluation.print();
        
        // 6. 验证子声明
        cout << "\n--- Verifying Subclaim ---" << endl;
        // 在子声明点上计算多项式值
        auto actual_eval = poly.evaluate(subclaim.point);
        cout << "Actual evaluation at subclaim point: ";
        actual_eval.print();
        
        if (actual_eval == subclaim.expected_evaluation) {
            cout << "✓ Subclaim VERIFIED: Actual evaluation matches expected!" << endl;
            cout << "✓ SumCheck protocol completed SUCCESSFULLY!" << endl;
        } else {
            cout << "✗ Subclaim FAILED: Actual evaluation does not match expected!" << endl;
            return;

        }
        
        /*
        // 7. 测试错误情况（可选）
        cout << "\n--- Testing with Wrong Sum (should fail) ---" << endl;
        try {
            Fp wrong_sum = Fp::zero(); // 错误的和
            auto transcript_wrong = sumcheck.init_transcript();
            auto wrong_subclaim = sumcheck.verify(wrong_sum, proof, poly.aux_info, transcript_wrong);
            cout << "✗ ERROR: Verification should have failed with wrong sum!" << endl;
            return 1;
        } catch (const exception& e) {
            cout << "✓ Correctly rejected wrong sum: " << e.what() << endl;
        }
        */
        
    } catch (const exception& e) {
        cerr << "Error: " << e.what() << endl;
        return;
    }
    
    cout << "\n=== SumCheck Example Completed ===" << endl;

}

int main() {
    test_sum_check();
    return 0;
}