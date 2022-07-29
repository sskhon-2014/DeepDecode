/*******************************************************
 * Copyright (C) 2022 California Institute of Technology - All Rights Reserved
 *
 * Author/Inventor: Harshaan Sekhon (sskhon2014@berkeley.edu), July 2022
 *
 * Unauthorized distribution or modification of this file via any medium is strictly prohibited.
 * All commerical applications of this software - or its excerpts - are strictly prohibited.
 * Avalibility of this program on GitHub is not permission to use this program in the aforementioned ways.
 * All materials here held proprietary.
 *
 * See Caltech Copyright Policy for Further Details (https://innovation.caltech.edu/patents-licensing/policies/caltech-copyright-and-software-policy)
 */


#include "utils.h"
#include <iterator>
#include <bits/stdc++.h>
#include <cmath>
#include <limits>
#include <string.h>


class valuePair {
  public:
    valuePair(int idx, double value);
    int idx;
    double value;
};
valuePair::valuePair(int idx, double value){
    this->idx = idx;
    this->value = value;
};

vector<vector<double>> create_2d_vector(int from, int to, int n){
    vector<vector<double>> retval;
    int rowsize = (to-from)/n;
    for(int i=0;i<n;i++){
        
        retval.push_back(create_vector(from, from+rowsize));
        from += rowsize;
    }
    return retval;
}

vector<double> create_vector(int from, int to){
    vector<double> row;

    for(int i=from;i<to;i++){
        row.push_back(i);
    }
    return row;
}


vector<vector<double>> csv_extract(string fname){
    
    
    vector<vector<double>> content;
    vector<double> row;
    string line, val;
    
    fstream file (fname, ios::in);
    file.ignore(10000, '\n');
    if(file.is_open())
    {
        while(getline(file, line))
        {
            row.clear();
            
            stringstream str(line);
            
            while(getline(str, val, ','))
            {
                
                if(val == "inf"){
                    row.push_back(INFINITY);
                }
                else if(val == "-inf"){
                    row.push_back(-INFINITY);
                }
                else{
                    row.push_back(stod(val));
                }
            }
            content.push_back(row);
        }
    }
    else
        cout<<"Could not open the file\n";
    
    
    
    return content;
}

void print_2d_vector(vector<vector<double>> content){
    if(content.size() <= 10){
        for(int i=0;i<content.size();i++)
        {
            for(int j=0;j<content[i].size();j++)
            {
                cout<<content[i][j]<<" ";
            }
            cout<<"\n";
        }
    }
    else{
        for(int i=0;i<5;i++)
        {
            for(int j=0;j<content[i].size();j++)
            {
                cout<<content[i][j]<<" ";
            }
            cout<<"\n";
        }
        cout << "\t........" << endl;
        for(int i=content.size()-5;i<content.size();i++)
        {
            for(int j=0;j<content[i].size();j++)
            {
                cout<<content[i][j]<<" ";
            }
            cout<<"\n";
        }
    }
    cout << content.size() << " rows" << endl;

}

void print_vector(vector<double> content){
    for(int i=0;i<content.size();i++){
        cout << content[i] << " ";
    }
    printf("\n");
}

string vector_to_str(vector<double> content){
    string retval = "[";
    retval += to_string((int)content[0]);
    for(int i=1;i<content.size();i++){
        retval += ", " + to_string((int)content[i]);
    }
    retval += "]";
    return retval;
}


bool compareVector(vector<double> v1, vector<double> v2)
{
    return (v1[PEAK] > v2[PEAK]);
}
bool compareVector2(valuePair* v1, valuePair* v2)
{
    return (v1->value > v2->value);
}



vector<vector<double>> remove_ignorable_rows(vector<vector<double>> content){
    for(int i=0;i<content.size();i++){
        if(content[i][XCENTROID] <= 0.1 || content[i][YCENTROID] <= 0.1){
            content.erase(content.begin() + i);
            i--;
        }
        else if(content[i][XCENTROID] >= MAX_DIM - 1.1 || content[i][YCENTROID] >= MAX_DIM - 1.1){
            content.erase(content.begin() + i);
            i--;
        }
    }
    return content;
}

vector<double> get_column(vector<vector<double>> content, int idx, bool rounded){
    vector<double> retval;
    for(int i=0;i<content.size();i++){
        if(rounded)
            retval.push_back(round(content[i][idx]));
        else
            retval.push_back(content[i][idx]);
    }
    return retval;
}

double* get_column_arr(vector<vector<double>> content, int idx, bool rounded){
    double* retval =(double*) malloc(sizeof(double) * content.size());
    for(int i=0;i<content.size();i++){
        if(rounded)
            retval[i] = round(content[i][idx]);
        else
            retval[i] = content[i][idx];
    }
    return retval;
}

vector<double> get_address_values(vector<double>*** address_book, double left_edge_x, double right_edge_x, double left_edge_y, double right_edge_y){
    vector<double> retval;
    for(int i=left_edge_x; i<right_edge_x;i++){
        if(address_book[i] != NULL){
            for(int j=left_edge_y;j<right_edge_y;j++){
                if(address_book[i][j] != NULL){
                    for(int k=0;k<address_book[i][j]->size();k++){
                        retval.push_back((*address_book[i][j])[k]);
                    }
                }
            }
        }
    }
    return retval;
}

vector<double> get_vector_indices(vector<double> vctr, vector<double> indices){
    vector<double> retval;
    for(int i=0;i<indices.size();i++){
        retval.push_back(vctr[indices[i]]);
    }
    return retval;
}

vector<double> get_distances(double x_i, double y_i, vector<double> x_short, vector<double> y_short){
    vector<double> retval;
    for(int i=0;i<x_short.size();i++){
        retval.push_back(((x_short[i] - x_i) * (x_short[i] - x_i)) + ((y_short[i] - y_i) * (y_short[i] - y_i)));
    }
    return retval;
}

vector<double> get_vector_cleaned(vector<double> vctr, vector<bool> indices){
    vector<double> retval;
    for(int i =0;i<vctr.size();i++){
        if(indices[i]){
            retval.push_back(vctr[i]);
        }
    }
    return retval;
}

vector<double> get_peak_sorted_indices(vector<double> peak_cleaned, vector<double> distances_cleaned){
    vector<valuePair*> pairs;
    for(int i=peak_cleaned.size()-1;i>=0;i--){
        pairs.push_back(new valuePair(i, (peak_cleaned[i] / (2+distances_cleaned[i]))));
    }
    sort(pairs.begin(), pairs.end(), compareVector2);
    vector<double> retval;
    for(int i=0;i<pairs.size();i++){
        retval.push_back(pairs[i]->idx);
    }

    return retval;
}

vector<double> get_n_elm(vector<double> v, int n){
    vector<double> retval;
    if(n> v.size()) n = v.size();
    for(int i=0;i<n;i++){
        retval.push_back(v[i]);
    }
    return retval;
}

double average(std::vector<double> const& v){
    if(v.empty()){
        return 0;
    }
    double total = 0;
    for(int i=0;i<v.size();i++){
        total+=v[i];
    }
    double count = (double)v.size();
    return total / count;
}

double stddev(vector<double> v){
    double mean = average(v);
    double result = 0;
    for(int i=0;i<v.size();i++){
        result += pow(v[i]-mean, 2);
    }
    result = sqrt(result/v.size());
    return result;
}

vector<double> get_vector_indices_from_array(double* arr, vector<double> indices){
    vector<double> retval;
    for(int i=0;i<indices.size();i++){
        retval.push_back(arr[(int)indices[i]]);
    }
    return retval;
}
