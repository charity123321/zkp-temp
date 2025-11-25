#include <iostream>
#include <unordered_map>
#include <vector>
#include <unordered_map>
#include <utility>
#include<memory>

#include"virtual_polynomial.h"
using namespace std;

template<typename F>
void build_eq_x_r_helper(
    const vector<F>& r,
    vector<F>& buf
){
    // r为空，则抛出错误
    if(r.empty()){
        throw;
    }else if(r.size()==1){
        buf.clear();
        buf.push_back(F::one()-r[0]);
        buf.push_back(r[0]);
    }else{
        vector<F> sub_r(r.begin()+1,r.end());
        build_eq_x_r_helper(sub_r,buf);

        // suppose at the previous step we received [b_1, ..., b_k]
        // for the current step we will need
        // if x_0 = 0:   (1-r0) * [b_1, ..., b_k]
        // if x_0 = 1:   r0 * [b_1, ..., b_k]
        // let mut res = vec![];
        // for &b_i in buf.iter() {
        //     let tmp = r[0] * b_i;
        //     res.push(b_i - tmp);
        //     res.push(tmp);
        // }
        // *buf = res;
        vector<F> res(buf.size()*2);

        for(size_t i=0;i<buf.size();++i){
            F tmp=r[0]*buf[i];
            res[2*i]=buf[i]-tmp;
            res[2*i+1]=tmp;
        }
        buf=move(res);
    }
}

template<typename F>
vector<F> build_eq_x_r_vec(vector<F>& r){
    // we build eq(x,r) from its evaluations
    // we want to evaluate eq(x,r) over x \in {0, 1}^num_vars
    // for example, with num_vars = 4, x is a binary vector of 4, then
    //  0 0 0 0 -> (1-r0)   * (1-r1)    * (1-r2)    * (1-r3)
    //  1 0 0 0 -> r0       * (1-r1)    * (1-r2)    * (1-r3)
    //  0 1 0 0 -> (1-r0)   * r1        * (1-r2)    * (1-r3)
    //  1 1 0 0 -> r0       * r1        * (1-r2)    * (1-r3)
    //  ....
    //  1 1 1 1 -> r0       * r1        * r2        * r3
    // we will need 2^num_var evaluations
    
    // r为空，则抛出错误
    if(r.empty()){
        throw;
    }
    vector<F> eval;
    build_eq_x_r_helper(r,eval);
    return eval;
}
/// This function build the eq(x, r) polynomial for any given r.
///
/// Evaluate
///      eq(x,y) = \prod_i=1^num_var (x_i * y_i + (1-x_i)*(1-y_i))
/// over r, which is
///      eq(x,y) = \prod_i=1^num_var (x_i * r_i + (1-x_i)*(1-r_i))
template<typename F>
shared_ptr<DenseMultilinearExtension<F>> build_eq_x_r(const vector<F>& r){
    auto evals=build_eq_x_r_vec<F>(r);
    return make_shared<DenseMultilinearExtension<F>>(
        DenseMultilinearExtension<F>::from_evaluations_vec(r.size(),move(evals))
    );
}

// 计算 \hat f(x) = \sum_{x_i \in eval_x} f(x_i) eq(x,r)
template<typename F>
VirtualPolynomial<F> build_f_hat(
    const VirtualPolynomial<F>& poly,
    const vector<F>& r
){
    // 变量不相等则抛出错误
    if (poly.aux_info.num_variables!=r.size()){
        throw; 
    }
    auto eq_x_r=build_eq_x_r<F>(r);
    VirtualPolynomial<F> res=poly;
    res.mul_by_mle(eq_x_r,F::one());
    
    return res;
}


// 计算 eq(x,y)=\prod (x_i * y_i + (1-x_i) * (1-y_i))
template<typename F>
F eq_eval(
    const vector<F>& x,
    const vector<F>& y
){
    // 两个向量长度不相等，抛出错误
    if(x.size()!=y.size()){
        throw;
    }
    F res =F::one();
    for(size_t i=0;i<x.size();++i){
        F xi_yi=x[i]*y[i];
        res=res*(F::one()+xi_yi+xi_yi-x[i]-y[i]);
    }
    return res;
}
