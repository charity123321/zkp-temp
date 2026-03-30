#ifndef MOCK_H
#define MOCK_H

#include"hyperplonk/witness.h"
#include"hyperplonk/hyperplonk_struct.h"
#include"hyperplonk/custom_gate.h"
#include"arithmetic/multilinear_polynomial.h"

template<typename scalar_field>
class MockCircuit{
public: 
    vector<scalar_field> public_inputs;
    vector<WitnessColumn<scalar_field>> witness;
    HyperPlonkIndex<scalar_field> index;

    size_t num_variable(){
        return index.num_variables();
    }

    size_t num_selector_columns(){
        return index.num_selector_columns();
    }

    size_t num_witness_columns(){
        return index.num_witness_columns();
    }

    MockCircuit(size_t num_contraints,CustomizedGates& gate){
        /*
        random_device rd;
        mt19937_64 rng(rd());
        */
        size_t nv=log2(num_contraints);
        size_t num_selectors=gate.num_selector_columns();
        size_t num_witness=gate.num_witness_columns();
        size_t log_n_wires=log2(num_witness);
        size_t merge_nv=nv+log_n_wires;

        vector<SelectorColumn<scalar_field>> selectors(num_selectors);
        vector<WitnessColumn<scalar_field>> witnesses(num_witness);

        // 为每个约束生成数据
        for(size_t cs_counter=0;cs_counter<num_contraints;++cs_counter){
            vector<scalar_field> cur_selectors;
            vector<scalar_field> cur_witness;

            // 生成随机选择器
            for(size_t i=0;i<num_selectors-1;++i){
                cur_selectors.push_back(scalar_field::random_element());
            }
            // 生成随机witness
            for(size_t i=0;i<num_witness;++i){
                cur_witness.push_back(scalar_field::random_element());
            }

            // 计算最后一个选择器
            scalar_field last_selector=scalar_field::zero();
            size_t index_counter=0;
            for(GateTerm gate_term:gate.gates){
                if(index_counter!=num_selectors-1){
                    scalar_field cur_monomial;
                    if(gate_term.coefficient<0){
                        cur_monomial=-scalar_field(static_cast<uint64_t>(-gate_term.coefficient));
                    }else{
                        cur_monomial=scalar_field(static_cast<uint64_t>(gate_term.coefficient));
                    }
                    if(gate_term.selector_index.has_value()){
                        size_t p=gate_term.selector_index.value();
                        cur_monomial*=cur_selectors[p];
                    }
                    for(size_t wit_index:gate_term.wire_indices){
                        cur_monomial*=cur_witness[wit_index];
                    }
                    last_selector+=cur_monomial;
                }else{
                    scalar_field cur_monomial;
                    if(gate_term.coefficient<0){
                        cur_monomial=-scalar_field(static_cast<uint64_t>(-gate_term.coefficient));
                    }else{
                        cur_monomial=scalar_field(static_cast<uint64_t>(gate_term.coefficient));
                    }
                    for(size_t wit_index:gate_term.wire_indices){
                        cur_monomial*=cur_witness[wit_index];
                    }
                    last_selector=last_selector*(-cur_monomial.inverse());
                }
                index_counter++;
            }

            cur_selectors.push_back(last_selector);

            // 将当前约束数据添加到各列
            for(size_t i=0;i<num_selectors;++i){
                selectors[i].append(cur_selectors[i]);
            }

            for(size_t i=0;i<num_witness;++i){
                witnesses[i].append(cur_witness[i]);
            }
        }
        // 设置公共输入
        size_t pub_input_len=min<size_t>(4,num_contraints);

        vector<scalar_field> public_inputs_vec;
        for(size_t i=0;i<pub_input_len;++i){
            public_inputs_vec.push_back(witnesses[0].data[i]);
        }

        HyperPlonkParams<scalar_field> params(num_contraints,pub_input_len,gate);
        vector<scalar_field> permutation=identity_permutation<scalar_field>(merge_nv,1);

        this->index=HyperPlonkIndex<scalar_field>(params,permutation,selectors);
        this->witness=witnesses;
        this->public_inputs=public_inputs_vec;

    }

    bool is_satisfied(){
        size_t num_vars=num_variable();

        for(size_t current_row=0;current_row<num_vars;++current_row){
            scalar_field cur=scalar_field::zero();

            for(GateTerm gate_term:index.params.gate_func.gates){
                scalar_field cur_monomial;
                if(gate_term.coefficient<0){
                    cur_monomial=-scalar_field(static_cast<uint64_t>(-gate_term.coefficient));
                }else{
                    cur_monomial=scalar_field(static_cast<uint64_t>(gate_term.coefficient));
                }

                if(gate_term.selector_index.has_value()){
                    size_t p=gate_term.selector_index.value();
                    cur_monomial*=index.selectors[p].data[current_row];
                }

                for(size_t wit_index:gate_term.wire_indices){
                    cur_monomial*=witness[wit_index].data[current_row];
                }
                cur+=cur_monomial;
            }
            /*
            cout<<"check "<<current_row<<endl;
            cur.print();
            cout<<endl;
            */
            if(!cur.is_zero()){
                return false;
            }
        }

        return true;
    }
};
#endif