#ifndef RELICH
#define RELICH
extern "C"{
    #include<relic.h>
}

#include<iostream>

#include<vector>
class RelicFp{
public:
    fp_t value;

    RelicFp(){
        fp_new(value);
        fp_zero(value);
    }
    RelicFp(uint64_t val) {
        fp_new(value);
        if (val == 0) {
            fp_zero(value);
        } else {
            bn_t temp;
            bn_new(temp);
            bn_set_dig(temp, static_cast<dig_t>(val));
            fp_prime_conv(value, temp); 
            bn_free(temp);
        }
    }

    RelicFp(const RelicFp& other){
        fp_new(value);
        fp_copy(value,other.value);
    }

    ~RelicFp(){
        fp_free(value);
    }

    RelicFp& operator=(const RelicFp& other){
        if(this!=&other){
            fp_copy(value,other.value);
        }
        return *this;
    }

    RelicFp operator+(const RelicFp& other)const{
        RelicFp result;
        fp_add(result.value,this->value,other.value);
        return result;
    }

    RelicFp& operator+=(const RelicFp& other){
        fp_add(this->value, this->value, other.value);
        return *this; 
    }

    RelicFp operator-(const RelicFp& other)const{
        RelicFp result;
        fp_sub(result.value,this->value,other.value);
        return result;
    }
    RelicFp& operator-=(const RelicFp& other) {
        fp_sub(this->value, this->value, other.value);
        return *this;
    }
    RelicFp operator-()const{
        RelicFp result;
        fp_neg(result.value,this->value);
        return result;
    }

    RelicFp operator*(const RelicFp& other)const{
        RelicFp result;
        fp_mul(result.value,this->value,other.value);
        return result;
    }
    RelicFp& operator*=(const RelicFp& other) {
        fp_mul(this->value, this->value, other.value);
        return *this;
    }

    bool operator==(const RelicFp& other){
        return fp_cmp(this->value,other.value)==RLC_EQ;
    }
    bool operator!=(const RelicFp& other) const {
        return fp_cmp(this->value, other.value) != RLC_EQ;
    }
    RelicFp inverse(){
        RelicFp result;
        fp_inv(result.value,this->value);
        return result;
    }

    static RelicFp random_element(){
        RelicFp result;
        fp_rand(result.value);
        return result;
    }

    static RelicFp zero(){
        RelicFp result;
        fp_zero(result.value);
        return result;
    }

    static RelicFp one(){
        RelicFp result;
        bn_t temp;
        bn_new(temp);
        bn_set_dig(temp,1);
        fp_prime_conv(result.value,temp);
        bn_free(temp);
        return result;
    }
    bool is_zero() const {
        return fp_is_zero(this->value) == 1;
    }
    void from_words(const std::vector<uint64_t>& words) {
        bn_t temp;
        bn_new(temp);
        bn_zero(temp);

        for (int i = static_cast<int>(words.size()) - 1; i >= 0; --i) {

            bn_lsh(temp, temp, 64);
            
            bn_t w;
            bn_new(w);
            bn_set_dig(w, static_cast<dig_t>(words[i]));
            bn_add(temp, temp, w);
            bn_free(w);
        }

        fp_prime_conv(this->value, temp);
        bn_free(temp);
    }
    void print()const{
        fp_print(value);
    }

};

class RelicFr {
public:
    bn_t value;
    /*
    static constexpr size_t num_bits = 255;
    // 获取 BLS12-381 的群阶 r 
    static bn_st* get_order() {
        static bn_t r;
        static bool initialized = false;
        if (!initialized) {
            bn_new(r);
            ep_curve_get_ord(r);
            initialized = true;
        }
        return r;
    }
    */
    static bn_st* get_order() {
        struct OrderHolder {
            bn_t r;
            OrderHolder() { bn_new(r); ep_curve_get_ord(r); }
            ~OrderHolder() { bn_free(r); }
        };
        static OrderHolder holder; // C++11 保证这里是线程安全的
        return holder.r;
    }
    RelicFr() {
        bn_new(value);
        bn_zero(value);
    }

    RelicFr(uint64_t val) {
        bn_new(value);
        bn_set_dig(value, static_cast<dig_t>(val));
        bn_mod(value, value, get_order());
    }

    RelicFr(const RelicFr& other) {
        bn_new(this->value);
        bn_copy(this->value, other.value);
    }

    ~RelicFr() {
        bn_free(this->value);
    }

    RelicFr& operator=(const RelicFr& other) {
        if (this != &other) bn_copy(value, other.value);
        return *this;
    }

    // 核心：模 r 运算
    RelicFr operator+(const RelicFr& other) const {
        RelicFr res;
        bn_add(res.value, this->value, other.value);
        bn_mod(res.value, res.value, get_order());
        return res;
    }

    RelicFr& operator+=(const RelicFr& other) {
        bn_add(this->value, this->value, other.value);
        bn_mod(this->value, this->value, get_order());
        return *this;
    }
    /*
    RelicFr operator-(const RelicFr& other) const {
        RelicFr res;
        bn_sub(res.value, this->value, other.value);
        if (bn_sign(res.value) == RLC_NEG) {
            bn_add(res.value, res.value, get_order());
        }
        return res;
    }
    */
    RelicFr operator-(const RelicFr& b) const {
        RelicFr res;
        // 1. 执行大数减法
        bn_sub(res.value, this->value, b.value);
        
        // 2. 【核心修复】如果结果是负数，必须加回模数 r
        if (bn_sign(res.value) == RLC_NEG) {
            // 获取你的标量场阶 r (73eda7...)
            bn_add(res.value, res.value, RelicFr::get_order());
            // 确保符号位转正
            res.value->sign = RLC_POS; 
        }
        return res;
    }
    RelicFr operator*(const RelicFr& other) const {
        RelicFr res;
        bn_mul(res.value, this->value, other.value);
        bn_mod(res.value, res.value, get_order());
        return res;
    }

    RelicFr& operator*=(const RelicFr& other) {
        bn_mul(this->value, this->value, other.value);
        bn_mod(this->value, this->value, get_order());
        return *this;
    }

    RelicFr inverse() const {
        RelicFr res;
        bn_t r_order, u, v;
        bn_new(r_order); bn_new(u); bn_new(v);

        // 1. 获取 BLS12-381 的标量场阶 r
        g1_get_ord(r_order); 

        // 2. 确保输入在 [0, r-1] 范围内 (处理负号位)
        bn_t temp_val;
        bn_new(temp_val);
        bn_copy(temp_val, this->value);
        if (bn_sign(temp_val) == RLC_NEG) {
            bn_add(temp_val, temp_val, r_order);
            temp_val->sign = RLC_POS;
        }
        bn_mod(temp_val, temp_val, r_order);

        // 3. 计算逆元: u * temp_val + v * r_order = 1
        bn_gcd_ext(u, res.value, v, temp_val, r_order);

        // 4. 处理负数结果并确保符号正确
        if (bn_sign(res.value) == RLC_NEG) {
            bn_add(res.value, res.value, r_order);
        }
        res.value->sign = RLC_POS;

        // 释放临时变量
        bn_free(r_order); bn_free(u); bn_free(v); bn_free(temp_val);
        
        return res;
    }
    bool operator==(const RelicFr& other) const {
        return bn_cmp(this->value, other.value) == RLC_EQ;
    }

    // 允许比较两个标量的大小（用于 map 的自动排序）
    bool operator<(const RelicFr& other) const {
        // RLC_LT 代表 Less Than (小于)
        return bn_cmp(this->value, other.value) == RLC_LT;
    }
    // 1. 不等号
    bool operator!=(const RelicFr& other) const {
        return !(*this == other);
    }

    // 2. 一元负号 (取反)
    RelicFr operator-() const {
        return RelicFr::zero() - *this; 
    }

    static RelicFr random_element() {
        RelicFr res;
        bn_rand_mod(res.value, get_order());
        return res;
    }
    void from_words(const std::vector<uint64_t>& words) {
            bn_zero(value);
            // 小端序逻辑：从最高索引开始处理
            for (int i = static_cast<int>(words.size()) - 1; i >= 0; --i) {
                bn_lsh(value, value, 64);
                bn_add_dig(value, value, static_cast<dig_t>(words[i]));
            }
            // 注意：如果你发现两端挑战还是对不上，请尝试下面这个真正的“小端序”写法：
            /*
            bn_zero(value);
            bn_t base, temp;
            bn_new(base); bn_new(temp);
            bn_set_dig(base, 1);
            for (uint64_t w : words) {
                bn_set_dig(temp, static_cast<dig_t>(w));
                bn_mul(temp, temp, base); // 这里逻辑较复杂，建议直接对齐 libff 的字节序
                bn_add(value, value, temp);
                bn_lsh(base, base, 64);
            }
            */
            bn_mod(value, value, get_order());
    }

    static RelicFr zero() { RelicFr res; bn_zero(res.value); return res; }
    static RelicFr one() { RelicFr res; bn_set_dig(res.value, 1); return res; }


    bool is_zero() const {
        bn_t temp;
        bn_new(temp);
        bn_mod(temp, this->value, get_order());
        bool res = bn_is_zero(temp);
        bn_free(temp);
        return res;
    }
    void print() const { bn_print(value); }
};


class RelicG1{
public:
    using scalar_field = RelicFr;
    ep_t value;

    RelicG1(){
        ep_new(value);
        ep_set_infty(value);
    }

    RelicG1(const RelicG1& other){
        ep_new(value);
        ep_copy(value,other.value);
    }

    ~RelicG1() {ep_free(value);}

    RelicG1& operator=(const RelicG1& other) {
        if (this != &other) ep_copy(value, other.value);
        return *this;
    }

    RelicG1 operator+(const RelicG1& other) const {
        RelicG1 res;
        ep_add(res.value, this->value, other.value);
        
        return res;
    }
    RelicG1 operator-(const RelicG1& other) const {
        RelicG1 res;
        ep_sub(res.value, this->value, other.value);
        return res;
    }
    RelicG1 operator*(const RelicFr& s) const {
        RelicG1 res;
        // RelicFr 内部是 bn_t，直接调用 ep_mul
        ep_mul(res.value, this->value, s.value);
        return res;
    }
    friend RelicG1 operator*(const RelicFr& s, const RelicG1& p) {
        return p * s;
    }
    static RelicG1 one() { 
        RelicG1 res; 
        ep_curve_get_gen(res.value); // 获取生成元
        return res; 
    }

    static RelicG1 zero() { 
        RelicG1 res; 
        ep_set_infty(res.value); 
        return res; 
    }
    
    bool operator==(const RelicG1& other) const {
        return ep_cmp(this->value, other.value) == RLC_EQ;
    }

    void print() const { ep_print(value); }
};

class RelicG2 {
public:
    using scalar_field = RelicFr;
    ep2_t value;

    RelicG2() { 
        ep2_new(value); 
        ep2_set_infty(value); // 修正：G2 对应的无穷远点
    }
    
    RelicG2(const RelicG2& other) { 
        ep2_new(value); 
        ep2_copy(value, other.value); 
    }
    
    ~RelicG2() { ep2_free(value); }

    RelicG2& operator=(const RelicG2& other) {
        if (this != &other) ep2_copy(value, other.value);
        return *this;
    }

    RelicG2 operator+(const RelicG2& other) const {
        RelicG2 res;
        ep2_add(res.value, this->value, other.value);
        return res;
    }
    RelicG2 operator-(const RelicG2& other) const {
        RelicG2 res;
        ep2_sub(res.value, this->value, other.value);
        return res;
    }
    RelicG2 operator*(const RelicFr& s) const {
        RelicG2 res;
        ep2_mul(res.value, this->value, s.value);
        return res;
    }
    friend RelicG2 operator*(const RelicFr& s, const RelicG2& p) {
        return p * s;
    }
    static RelicG2 one() { 
        RelicG2 res; 
        ep2_curve_get_gen(res.value); 
        return res; 
    }

    static RelicG2 zero() { 
        RelicG2 res; 
        ep2_set_infty(res.value); 
        return res; 
    }
};

// =================================================================
// 目标群包装类 (GT - Target Group on Fp12)
// =================================================================
class RelicGT {
public:
    fp12_t value;

    RelicGT() {
        fp12_new(value);
        // 初始化为乘法单位元 1
        fp12_zero(value);
        fp_set_dig(value[0][0][0], 1); 
    }

    RelicGT(const RelicGT& other) {
        fp12_new(value);
        fp12_copy(value, other.value);
    }

    ~RelicGT() {
        fp12_free(value);
    }

    RelicGT& operator=(const RelicGT& other) {
        if (this != &other) fp12_copy(value, other.value);
        return *this;
    }

    RelicGT operator*(const RelicGT& other) const {
        RelicGT res;
        fp12_mul(res.value, this->value, other.value);
        return res;
    }

    bool operator==(const RelicGT& other) const {
        return fp12_cmp(this->value, other.value) == RLC_EQ;
    }

    static RelicGT one() {
        return RelicGT(); // 构造函数默认就是单位元
    }
};

struct relic_bls12_381_pp {
    using scalar_field = RelicFr;
    using base_field   = RelicFp;
    using G1_type      = RelicG1;
    using G2_type      = RelicG2;
    using GT_type      = RelicGT;

    static void init_public_params() {
        if (core_init() != RLC_OK) {
            std::cerr << "RELIC core_init failed!" << std::endl;
            exit(1);
        }
        // 显式加载 BLS12-381 曲线
        ep_param_set(B12_P381);

        // 预热群阶缓存，确保 RelicFr 逻辑正确
        RelicFr::get_order();
    }
    static GT_type multi_pairing(const std::vector<G1_type>& g1_vec, const std::vector<G2_type>& g2_vec) {
        GT_type res;
        if (g1_vec.empty() || g1_vec.size() != g2_vec.size()) return res;

        // 零拷贝指针转换
        const ep_t* p_ptr = reinterpret_cast<const ep_t*>(g1_vec.data());
        const ep2_t* q_ptr = reinterpret_cast<const ep2_t*>(g2_vec.data());

        // RELIC 底层的同时配对函数：合并 Miller Loop，只需 1 次 Final Exp
        pc_map_sim(res.value, p_ptr, q_ptr, g1_vec.size());
        
        return res;
    }
};

/*
template<typename pp_t>
struct G1 {
    using scalar_field = typename pp_t::scalar_field;
    using group_type   = typename pp_t::G1_type;
};

template<typename pp_t>
struct G2 {
    using scalar_field = typename pp_t::scalar_field;
    using group_type   = typename pp_t::G2_type;
};
*/
template<typename P> using G1 = typename P::G1_type;
template<typename P> using G2 = typename P::G2_type;
template<typename P> using GT = typename P::GT_type;
#endif