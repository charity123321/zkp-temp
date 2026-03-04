#ifndef WITNESS_H
#define WITNESS_H

#include"hyperplonk/hyperplonk_utils.h"
#include<iostream>

using namespace std;


template<typename scalar_field>
struct WitnessRow{
    vector<scalar_field> data;
};


template<typename scalar_field>
class WitnessColumn{
public:
    vector<scalar_field> data;

    WitnessColumn()=default;

    WitnessColumn(const vector<scalar_field>& data_):data(data_) {}

    size_t get_nv()const {
        return log2(data.size());
    }

    void append(const scalar_field& new_element){
        data.push_back(new_element);
    }

    vector<WitnessColumn<scalar_field>> form_witness_rows(
        const vector<WitnessRow<scalar_field>>& witness_rows
    ){
        if(witness_rows.empty()){
            cout<<"empty"<<endl;
            throw;
        }

        vector<WitnessColumn<scalar_field>> res;
        size_t num_columns=witness_rows[0].data.size();

        for(size_t i=0;i<num_columns;++i){
            vector<scalar_field> cur_column;
            for(WitnessRow<scalar_field> row:witness_rows){
                cur_column.push_back(row.data[i]);
            }
            res.push_back(WitnessColumn<scalar_field>(move(cur_column)));
        }
        return res;
    }
};



#endif