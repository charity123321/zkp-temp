#ifndef CUSTOM_GATE_H
#define CUSTOM_GATE_H

#include<iostream>
#include<vector>
#include<optional>
using namespace std;

/// Customized gate is a list of tuples of
///     (coefficient, selector_index, wire_indices)
///
/// Example:
///     q_L(X) * W_1(X)^5 - W_2(X) = 0
/// is represented as
/// vec![
///     ( 1,    Some(id_qL),    vec![id_W1, id_W1, id_W1, id_W1, id_W1]),
///     (-1,    None,           vec![id_W2])
/// ]
///
/// CustomizedGates {
///     gates: vec![
///         (1, Some(0), vec![0, 0, 0, 0, 0]),
///         (-1, None, vec![1])
///     ],
/// };
/// where id_qL = 0 // first selector
/// id_W1 = 0 // first witness
/// id_w2 = 1 // second witness
///
/// here coeff is a signed integer, instead of a field element



struct GateTerm{
    int64_t coefficient;
    optional<size_t> selector_index;
    vector<size_t> wire_indices;
};
class CustomizedGates{
public:
    vector<GateTerm> gates;

    CustomizedGates() =default;
    CustomizedGates(const vector<GateTerm>& gate):gates(gate) {}

    size_t degree(){
        size_t res=0;
        for(GateTerm term:gates){
            size_t term_degree=term.wire_indices.size();
            if(term.selector_index.has_value()){
                term_degree+=1;
            }
            res=max(res,term_degree);
        }
        return res;
    }

    size_t num_selector_columns(){
        size_t res=0;
        for(GateTerm term:gates){
            if(term.selector_index.has_value()){
                res+=1;
            }
        }
        return res;
    }

    size_t num_witness_columns(){
        size_t res=0;
        for(GateTerm term:gates){
            if(!term.wire_indices.empty()){
                res=max(res,term.wire_indices.back());
            }
        }
        return res+1;
    }
    /// Return a vanilla plonk gate:
    /// ``` ignore
    ///   q_L w_1 + q_R w_2 + q_O w_3 + q_M w1w2 + q_C = 0
    /// ```
    /// which is
    /// ``` ignore
    ///     (1,    Some(id_qL),     vec![id_W1]),
    ///     (1,    Some(id_qR),     vec![id_W2]),
    ///     (1,    Some(id_qO),     vec![id_W3]),
    ///     (1,    Some(id_qM),     vec![id_W1, id_w2]),
    ///     (1,    Some(id_qC),     vec![]),
    /// ```
    CustomizedGates vanilla_plonk_gates(){
        return CustomizedGates({
            {1,0,{0}},
            {1,1,{1}},
            {1,2,{2}},
            {1,3,{0,1}},
            {1,4,{}}
        });
    }

    /// Generate a random gate for `num_witness` with a highest degree =
    /// `degree`
    CustomizedGates mock_gate(size_t num_witness,size_t degree){
        vector<GateTerm> gates;

        vector<size_t> high_degree_term(degree-1,0);
        // q_0 w_0^{degree-1} * w_1
        high_degree_term.push_back(1);
        gates.push_back({1,0,move(high_degree_term)});

        // q_{i+1} w_i
        for(size_t i=0;i<num_witness;++i){
            gates.push_back({1,i+1,{i}});
        }
        // 常数项
        gates.push_back({1,num_witness+1,{}});

        return CustomizedGates(move(gates));
    }
    /// Return a plonk gate where #selector > #witness * 2
    /// ``` ignore
    ///   q_1 w_1   + q_2 w_2   + q_3 w_3   +
    ///   q_4 w1w2  + q_5 w1w3  + q_6 w2w3  +
    ///   q_7 = 0
    /// ```
    /// which is
    /// ``` ignore
    ///     (1,    Some(id_qL),     vec![id_W1]),
    ///     (1,    Some(id_qR),     vec![id_W2]),
    ///     (1,    Some(id_qO),     vec![id_W3]),
    ///     (1,    Some(id_qM),     vec![id_W1, id_w2]),
    ///     (1,    Some(id_qC),     vec![]),
    /// ```
    CustomizedGates super_long_selector_gate() {
        return CustomizedGates({
            {1, 0, {0}},
            {1, 1, {1}},
            {1, 2, {2}},
            {1, 3, {0, 1}},
            {1, 4, {0, 2}},
            {1, 5, {1, 2}},
            {1, 6, {}}
        });
    }
};
#endif