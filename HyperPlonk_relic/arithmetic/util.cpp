#include <iostream>
#include <vector>
#include <cmath>
#include <cstdint>
#include <tuple>



// 将整数转化为小端存储的二进制向量
// 5 --> [1,0,1,0]
std::vector<bool> bit_decompose(uint64_t input, size_t num_var)
{
    std::vector<bool> res;
    res.reserve(num_var);
    uint64_t i = input;
    for (size_t j = 0; j < num_var; j++)
    {
        res.push_back((i & 1) == 1);
        i >>= 1;
    }
    return res;
}

// 返回合并MLE后的变量数
size_t get_batched_nv(size_t num_var, size_t polynomials_len)
{
    return num_var + static_cast<size_t>(log2(polynomials_len));
}

// 将小段存储的二进制向量转为整数
uint64_t project(const std::vector<bool> &input)
{
    uint64_t res = 0;
    for (auto it = input.rbegin(); it != input.rend(); ++it)
    {
        res <<= 1;
        res += *it ? 1 : 0;
    }
    return res;
}

/// Input index
/// - `i := (i_0, ...i_{n-1})`,
/// - `num_vars := n`
/// return three elements:
/// - `x0 := (i_1, ..., i_{n-1}, 0)`
/// - `x1 := (i_1, ..., i_{n-1}, 1)`
/// - `sign := i_0`
std::tuple<size_t, size_t, bool> get_index(size_t i, size_t num_vars)
{
    std::vector<bool> bit_sequence = bit_decompose(static_cast<uint64_t>(i), num_vars);

    // 构建 base_bits: bit_sequence[1..num_vars-1]
    std::vector<bool> base_bits(bit_sequence.begin(), bit_sequence.end() - 1);

    // 构建 x0_bits
    std::vector<bool> x0_bits = base_bits;
    x0_bits.insert(x0_bits.begin(), false);

    // 构建 x1_bits
    std::vector<bool> x1_bits = base_bits;
    x1_bits.insert(x1_bits.begin(), true);

    size_t x0 = static_cast<size_t>(project(x0_bits));
    size_t x1 = static_cast<size_t>(project(x1_bits));
    bool sign = bit_sequence[num_vars - 1];

    return std::make_tuple(x0, x1, sign);
}
