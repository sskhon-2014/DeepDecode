/*******************************************************
 * Copyright (C) 2022 Harshaan Sekhon, sskhon2014@berkeley.edu - All Rights Reserved
 * Unauthorized distribution or modification of this file via any medium is strictly prohibited.
 * All commerical applications of this software or its excerpts are strictly prohibited.
 * Avalibility of this program on GitHub is not permission to use this program in the aforementioned ways.
 * All materials here held proprietary.
 *
 * Author: Harshaan Sekhon (sskhon2014@berkeley.edu), July 2022
 */

#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <string>

#define UNNAMED 0
#define ID 1
#define XCENTROID 2
#define YCENTROID 3
#define SHARPNESS 4
#define ROUNDESS1 5
#define ROUNDNESS2 6
#define NPIX 7
#define SKY 8
#define PEAK 9
#define FLUX 10
#define MAG 11
#define HYB 12

#define MAX_DIM 2048
#define DIST_CUTOFF 1


using namespace std;


vector<vector<double>> create_2d_vector(int from, int to, int n);
vector<vector<double>> csv_extract(string filename);
void print_2d_vector(vector<vector<double>> content);
vector<double> create_vector(int from, int to);

void print_vector(vector<double> content);

bool compareVector(vector<double> v1, vector<double> v2);

vector<vector<double>> remove_ignorable_rows(vector<vector<double>> content);

double* get_column_arr(vector<vector<double>> content, int idx, bool rounded);
vector<double> get_column(vector<vector<double>> content, int idx, bool rounded);

vector<double> get_vector_indices(vector<double> vctr, vector<double> indices);
vector<double> get_vector_indices_from_array(double* arr, vector<double> indices);

vector<double> get_address_values(vector<double>*** address_book, double left_edge_x, double right_edge_x, double left_edge_y, double right_edge_y);
vector<double> get_distances(double x_i, double y_i, vector<double> x_short, vector<double> y_short);
vector<double> get_vector_cleaned(vector<double> vctr, vector<bool> indices);
vector<double> get_peak_sorted_indices(vector<double> peak_cleaned, vector<double> distances_cleaned);

vector<double> get_n_elm(vector<double> v, int n);
double average(std::vector<double> const& v);
double stddev(vector<double> v);
string vector_to_str(vector<double> content);
