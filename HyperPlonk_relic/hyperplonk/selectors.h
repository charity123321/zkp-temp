#ifndef SELECTORS_H
#define SELECTORS_H

#include"hyperplonk/hyperplonk_utils.h"
#include<iostream>



template<typename scalar_field>
struct SelectorRow{
    std::vector<scalar_field> data;
};
/// Build MLE from matrix of selectors.
///
/// Given a matrix := [row1, row2, ...] where
/// row1:= (a1, a2, ...)
/// row2:= (b1, b2, ...)
/// row3:= (c1, c2, ...)
///
/// output mle(a1,b1,c1, ...), mle(a2,b2,c2, ...), ...
/*
template<typename scalar_field>
std::vector<shared_ptr<DenseMultilinearExtension<scalar_field>>> build_mles(
    std::vector<SelectorRow<scalar_field>>& matrix
){
    return build_mle<scalar_field,SelectorRow>(matrix);
}
*/
template<typename scalar_field>
class SelectorColumn{
public:
    std::vector<scalar_field> data;

    SelectorColumn()=default;

    SelectorColumn(const std::vector<scalar_field>& data_):data(data_) {}

    size_t get_nv() const{
        return log2(data.size());
    }

    void append(const scalar_field& new_element){
        data.push_back(new_element);
    }

    // build selector columns from rows
    std::vector<SelectorColumn<scalar_field>> from_selector_rows(
        std::vector<SelectorRow<scalar_field>>& selector_rows
    ){
        if(selector_rows.empty()){
            cout<<"empty"<<endl;
            throw;
        }

        std::vector<SelectorColumn<scalar_field>> res;
        size_t num_columns=selector_rows[0].data.size();

        for(size_t i=0;i<num_columns;++i){
            std::vector<scalar_field> cur_column;
            for(SelectorRow<scalar_field> row:selector_rows){
                cur_column.push_back(row.data[i]);
            }
            res.push_back(SelectorColumn<scalar_field>(move(cur_column)));
        }
        return res;
    }

    shared_ptr<DenseMultilinearExtension<scalar_field>> from(
       const SelectorColumn<scalar_field>& witness
    ){
        size_t nv=witness.get_nv();
        return make_shared<DenseMultilinearExtension<scalar_field>>(nv,witness.data);
    }
};



#endif