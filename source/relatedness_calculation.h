#ifndef RELATEDNESS_CALCULATION_H
#define RELATEDNESS_CALCULATION_H

#include <deque>
#include <string>
#include <thread>
#include <set>
#include <cmath>
#include <iomanip>
#include <map>
#include <vector>
#include <thread>


#include "basic_functions.h"
#include "node.h"
#include "dyad.h"
#include "general_class_functions.h"

using namespace std; 

double calculate_f_xy(node * indiv_1,node * indiv_2,std::deque<std::deque<double>> *f_matrix,int ncol,int nrow,std::deque <node> * all_nodes_ptr,std::deque <dyad> * all_dyads_ptr,int dyad_index,std::deque <int> * path_indices_ptr,int path_index,map <string,int> * dyad_dict_ptr,int current_generation_1,int current_generation_2,int generation_limit, std::set<string>* node_space_ptr);
void fill_f_matrix(std::deque<std::deque<double>> *f_matrix,int ncol, int nrow,std::deque <node> *all_nodes_ptr,std::deque <dyad> *all_dyads_ptr,map<string,int>*dyad_dict_ptr,int generation_limit,bool multithreading = false, int dyads_start = 0,int dyads_end = 0,int thread = 0);
void add_missing_parent_in_dyad_list(std::deque<node>*all_nodes_ptr,std::string all_node_names_str,std::deque<string>*mom_names_ptr,std::deque<string>*sire_names_ptr);
string all_nodes_info_file(std::deque <node> *all_nodes_ptr,std::map<string,int> * dyad_dict,std::deque<dyad>* all_dyads_ptr, string filename);
void write_dyad_list(std::deque <dyad> *all_dyads_ptr,std::deque<std::deque<double>> &f_mat,int ncol,string filename, string write_dyadlist = "full");
void set_min_f(std::deque<std::deque<double>> &f_matrix,std::deque<dyad>*all_dyads_ptr,std::deque<node>*all_nodes_ptr);
void set_all_min_DGD(std::deque<node>*all_nodes_ptr, std::deque<dyad>*all_dyads_ptr);
void pip_forward(string file,string output_file,string input_dyadlist,int maturation_age_f,int maturation_age_m,int gestation_length,string write_dyadlist,int generation_limit = -1,int n_cores = 1);
void reduce_node_space(node* indiv_1,node* indiv_2,std::set<string>*node_space);
bool indivs_in_node_space(node* indiv_1, node* indiv_2,std::set<string>*node_space_ptr);
#endif