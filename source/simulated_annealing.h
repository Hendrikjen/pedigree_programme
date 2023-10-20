#ifndef SIMULATED_ANNEALING_H
#define SIMULATED_ANNEALING_H

#include <string>
#include <deque>
#include <map>
#include <set>
#include <cmath>
#include <thread>
#include <vector>
#include <iomanip>
#include <random>
#include <fstream>
#include <filesystem>

#include "node.h"
#include "dyad.h"
#include "gap.h"
#include "general_class_functions.h"
#include "relatedness_calculation.h"

using namespace std; 
namespace fs = std::filesystem;

double randomnumber();
double load_data_for_sim_annealing(string file_full,string file_gaps, std::deque<dyad>*all_dyads_ptr,std::deque<int>*subset_idx_ptr, std::deque<node>*all_nodes_ptr,map<string, int> *dyad_dict_ptr,int maturation_age_f,int maturation_age_m,int gestation_length,string file_dyads="NA");
std::string get_next_solution_all_gaps(int count,double temperature,std::deque<gap>*gaps_ptr,int default_year,int maturation_age_m,int maturation_age_f,int gestation_length,std::deque<node>*all_nodes_ptr,std::deque<node>*new_nodes_ptr);
double compare_f_matrix(std::deque<std::deque<double>> &f_mat_1,std::deque<std::deque<double>> & f_mat_2,std::deque<node>*all_nodes_ptr,std::deque<node>*nodes_ptr_2, std::deque<int>*subset_idx_ptr, std::deque<dyad>*all_dyads_ptr);
double get_init_temp_factor(std::deque<std::deque<double>> &f_mat,std::deque<node>*all_nodes_ptr);
void fill_pure_f_matrix(std::deque<dyad>* all_dyads_ptr,std::deque<int>* subset_idx_ptr,std::deque<node>* all_nodes_ptr,std::deque<node>* current_nodes_ptr,std::deque<std::deque<double>>* f_matrix,bool full_ped,int dyads_start = 0,int dyads_end = 0,bool multithreading = false,int thread_no = 0);
void simulated_annealing_complete_pedigree(double init_temp,double stop_temp,double temp_decay,std::deque<node>*all_nodes_ptr,std::deque<dyad>*all_dyads_ptr,std::deque<int>*subset_idx_ptr,int gestation_length,int maturation_age_f,int maturation_age_m,int default_year,string filepath,bool visualization,bool multithreading = false, int n_cores = 1);
#endif