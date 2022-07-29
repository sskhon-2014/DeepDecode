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

#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <list>
#include <cmath>
#include "utils.h"
#include <chrono>


using namespace std;

vector<vector<double>> barcoding_sets;
double convert_edge(double edge){
    if(edge <= 0) return 0;
    else if(edge >= MAX_DIM-1) return MAX_DIM -1;
    else return edge;
}

vector<bool> repeat_truth(vector<double> v){
    vector<bool> back_propogation;
    vector<double> mem;
    vector<double> mem_set;

    for(int i=0;i<v.size();i++){
        double hyb = v[i];
        if(find(mem.begin(), mem.end(), hyb) != mem.end()){
            back_propogation.push_back(false);
        }
        else{
            if(find(mem_set.begin(), mem_set.end(), hyb) != mem_set.end()){
                back_propogation.push_back(false);
            }
            else{
                back_propogation.push_back(true);
            }
        }
        mem.push_back(hyb);
        for(int i=0;i<barcoding_sets.size();i++){
            if(find(barcoding_sets[i].begin(), barcoding_sets[i].end(), hyb) != barcoding_sets[i].end()){
                for(int j=0;j<barcoding_sets[i].size();j++){
                    mem_set.push_back(barcoding_sets[i][j]);
                }
            }
        }

    }

    return back_propogation;
}


int main() {
    auto s = std::chrono::high_resolution_clock::now();
    barcoding_sets = create_2d_vector(45,81,4);
    print_2d_vector(barcoding_sets);

    vector<vector<double>>all_dots;
    
    all_dots = csv_extract("input.csv");

    sort(all_dots.begin(), all_dots.end(), compareVector);
    print_2d_vector(all_dots);


    all_dots = remove_ignorable_rows(all_dots);
    print_2d_vector(all_dots);
    vector<double> x_rounded = get_column(all_dots, XCENTROID, true);
    vector<double> y_rounded = get_column(all_dots, YCENTROID, true);
    vector<double> indices = create_vector(0, y_rounded.size());
    
    //print_vector(indices);

    vector<double>** address_book[MAX_DIM];
    for(int i=0;i<MAX_DIM;i++){
        address_book[i] = nullptr;
    }
    for(int i=0;i<indices.size()-1;i++){
        
        if(address_book[(int)x_rounded[i]] == nullptr){
            address_book[(int)x_rounded[i]] = (vector<double>**)malloc(sizeof(vector<double>*) * MAX_DIM);
            for(int j=0;j<MAX_DIM;j++){
                address_book[(int)x_rounded[i]][j] = nullptr;
            }
        }
        if(address_book[(int)x_rounded[i]] != nullptr && address_book[(int)x_rounded[i]][(int)y_rounded[i]] == nullptr){
            address_book[(int)x_rounded[i]][(int)y_rounded[i]] = new vector<double>();
        }
        if(address_book[(int)x_rounded[i]] != nullptr && address_book[(int)x_rounded[i]][(int)y_rounded[i]] != nullptr){
            address_book[(int)x_rounded[i]][(int)y_rounded[i]]->push_back(indices[i]);
        }
        
        
        //cout << "----------------------------------------------------------" << endl;
    }

    int size = all_dots.size();
    double *x = get_column_arr(all_dots, XCENTROID, false);
    double *y = get_column_arr(all_dots, YCENTROID, false);
    double *peak = get_column_arr(all_dots, PEAK, false);
    double *hyb = get_column_arr(all_dots, HYB, false);

    bool indices_run[size];
    for(int i=0;i<size;i++)indices_run[i] = false;

    vector<vector<double>> barcodes;
    vector<double> mean_peak;
    vector<double> mean_sd;
    vector<double> sd_xpos;
    vector<double> sd_ypos;
    vector<double> mean_xpos;
    vector<double> mean_ypos;

    double part1 = 0, part2 = 0, part3 = 0, part4 = 0, part5=0, part6;


    //print_vector(indices_to_span);
    auto start = std::chrono::high_resolution_clock::now();
    auto current = std::chrono::high_resolution_clock::now();

    for(int i=0;i<size;i++){
        if(i%1000 == 0){
            cout << i << endl;
            /*
            cout << "part1: " << (int)part1 << endl;
            cout << "part2: " << (int)part2 << endl;
            cout << "part3: " << (int)part3 << endl;
            cout << "part4: " << (int)part4 << endl;
            cout << "part5: " << (int)part5 << endl;
            */
        }
        //print_vector(indices_run);
        if(indices_run[i] == false){
            double x_i = x[i], y_i = y[i], peak_i = peak[i], hyb_i = hyb[i], index_i = i;

            double x_r = round(x_i);
            double y_r = round(y_i);

            double left_edge_x = convert_edge(x_r - 2);
            double right_edge_x = convert_edge(x_r +2);

            double left_edge_y = convert_edge(y_r - 2);
            double right_edge_y = convert_edge(y_r +2);
            
            vector<double> nearby_indices = get_address_values(address_book, left_edge_x, right_edge_x,left_edge_y, right_edge_y);
            //print_vector(nearby_indices);
            for(int i=0;i<nearby_indices.size();i++){
                if(indices_run[(int)nearby_indices[i]]){ 
                    nearby_indices.erase(nearby_indices.begin()+i);
                    i--;
                }
            }
            //print_vector(nearby_indices);
            
            current = std::chrono::high_resolution_clock::now();
            part1 += std::chrono::duration_cast<std::chrono::milliseconds>(current-start).count();
            start = current;

            

            vector<double> x_short = get_vector_indices_from_array(x, nearby_indices);
            vector<double> y_short = get_vector_indices_from_array(y, nearby_indices);
            vector<double> peak_short = get_vector_indices_from_array(peak, nearby_indices);
            vector<double> hyb_short = get_vector_indices_from_array(hyb, nearby_indices);
/*
            print_vector(x_short);
            print_vector(y_short);
*/
            auto current = std::chrono::high_resolution_clock::now();
            part2 += std::chrono::duration_cast<std::chrono::milliseconds>(current-start).count();
            start = current;

            vector<double> distances = get_distances(x_i, y_i, x_short, y_short);
            vector<bool> within_1px;
            for(int i=0;i<distances.size();i++){
                if(distances[i] < DIST_CUTOFF)within_1px.push_back(true);
                else within_1px.push_back(false);
            }

            
            vector<double> peak_cleaned = get_vector_cleaned(peak_short, within_1px);
            //print_vector(peak_cleaned);
            vector<double> hyb_cleaned = get_vector_cleaned(hyb_short, within_1px);
            vector<double> x_cleaned = get_vector_cleaned(x_short, within_1px);
            vector<double> y_cleaned = get_vector_cleaned(y_short, within_1px);
            vector<double> distances_cleaned = get_vector_cleaned(distances, within_1px);
            vector<double> indices_cleaned = get_vector_cleaned(nearby_indices, within_1px);

/*
            current = std::chrono::high_resolution_clock::now();
            part3 += std::chrono::duration_cast<std::chrono::milliseconds>(current-start).count();
            start = current;

            cout << "part3" << endl;
            print_vector(peak_cleaned);
            print_vector(hyb_cleaned);
            print_vector(x_cleaned);
            print_vector(y_cleaned);
            print_vector(indices_cleaned);
*/

            vector<double> peak_sorted_indices = get_peak_sorted_indices(peak_cleaned, distances_cleaned);
/*
            for(int i=0;i<peak_cleaned.size();i++){
                cout <<  peak_sorted_indices[i] << ":"<< peak_cleaned[i] << endl;
            }

*/
            vector<double> hyb_sorted = get_vector_indices(hyb_cleaned, peak_sorted_indices);
            vector<double> peak_sorted = get_vector_indices(peak_cleaned, peak_sorted_indices);
            vector<double> distance_sorted = get_vector_indices(distances_cleaned, peak_sorted_indices);
            for(int i=0;i<peak_sorted.size();i++){
                peak_sorted[i] = peak_sorted[i]/(1+distance_sorted[i]);
            }
            vector<double> x_sorted = get_vector_indices(x_cleaned, peak_sorted_indices);
            vector<double> y_sorted = get_vector_indices(y_cleaned, peak_sorted_indices);
            vector<double> indices_sorted = get_vector_indices(indices_cleaned, peak_sorted_indices);
            
/*
            current = std::chrono::high_resolution_clock::now();
            part4 += std::chrono::duration_cast<std::chrono::milliseconds>(current-start).count();
            start = current;

            cout << "part4" << endl;
            print_vector(hyb_short);
            print_vector(hyb_sorted);
            print_vector(peak_sorted);
            print_vector(x_sorted);
            print_vector(y_sorted);
            print_vector(indices_sorted);
          
    */        
            vector<bool> unrepeated_indices = repeat_truth(hyb_sorted);
  /*          print_vector(hyb_sorted);
            for(int i=0;i<unrepeated_indices.size();i++){
                cout << unrepeated_indices[i] << " ";
            }
            
            cout << endl;
    
            getchar();
*/
            vector<double> hybs_unrepeat = get_vector_cleaned(hyb_sorted, unrepeated_indices);
            vector<double> peak_unrepeat = get_vector_cleaned(peak_sorted, unrepeated_indices);
            vector<double> x_unrepeat = get_vector_cleaned(x_sorted, unrepeated_indices);
            vector<double> y_unrepeat = get_vector_cleaned(y_sorted, unrepeated_indices);
            vector<double> indices_unrepeat = get_vector_cleaned(indices_sorted, unrepeated_indices);

/*
            current = std::chrono::high_resolution_clock::now();
            part5 += std::chrono::duration_cast<std::chrono::milliseconds>(current-start).count();
            start = current;
                        cout << "part5" << endl;

            print_vector(hybs_unrepeat);
            print_vector(peak_unrepeat);
            print_vector(x_unrepeat);
            print_vector(y_unrepeat);
            print_vector(indices_unrepeat);
  */          
            
            if(hybs_unrepeat.size() > 1){
                barcodes.push_back(get_n_elm(hybs_unrepeat, 4));
                mean_peak.push_back(average(get_n_elm(peak_unrepeat, 4)));
                mean_sd.push_back(stddev(get_n_elm(peak_unrepeat, 4)));

                sd_xpos.push_back(stddev(get_n_elm(x_unrepeat, 4)));
                sd_ypos.push_back(stddev(get_n_elm(y_unrepeat, 4)));

                mean_xpos.push_back(average(get_n_elm(x_unrepeat, 4)));
                mean_ypos.push_back(average(get_n_elm(y_unrepeat, 4)));
   /*             
                cout << "--------//////////---------" << endl;
                
                print_2d_vector(barcodes);
                print_vector(mean_peak);
                print_vector(mean_sd);
                print_vector(sd_xpos);
                print_vector(sd_ypos);
                print_vector(mean_xpos);
                print_vector(mean_ypos);
                
*/          
                //print_vector(indices_unrepeat);
                for(int i=0;i<indices_unrepeat.size();i++){
                    indices_run[(int)indices_unrepeat[i]] = true;
                }
                //print_vector(mean_peak);


            }
/*
            current = std::chrono::high_resolution_clock::now();
            part6 += std::chrono::duration_cast<std::chrono::milliseconds>(current-start).count();
            start = current;
            getchar();

  */          
        }
    }

        
    ofstream oFile ("out.csv"); 
    oFile << ",code,Mean-Peak,Mean-SD,SD-X,SD-Y,Mean-X,Mean-y\n";
    for(int i=0;i<mean_peak.size();i++){
        oFile << i << ",\"" << vector_to_str(barcodes[i]) << "\"," << to_string(mean_peak[i]) << ","  << to_string(mean_sd[i]) << "," <<to_string(sd_xpos[i]) << "," << to_string(sd_ypos[i])<< "," << to_string(mean_xpos[i]) << "," << to_string(mean_ypos[i])<<endl; 
    }
    oFile.close();


    return 0;

}
