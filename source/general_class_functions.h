#ifndef GENERAL_CLASS_FUNCTIONS_H
#define GENERAL_CLASS_FUNCTIONS_H

#include <string>
#include <deque>
#include <map>
#include <iostream>
#include <fstream>

#include "basic_functions.h"
#include "node.h"
#include "dyad.h"
#include "gap.h"

using namespace std; 

void set_parent_pool(std::deque <node> *all_nodes_ptr, int maturation_age_m,int maturation_age_f,int gestation_length,bool twins);
void get_ancestor_of_focal_indiv(std::deque<node*>*anc_nodes_ptr,node*indiv,std::map<string,std::deque<node*>>*offspring_dict);
bool contains_ancestor(node* pot_offspring,node* pot_ancestor,std::deque <node> *all_nodes_ptr);
double calculate_pure_f_xy(node*indiv_1,node*indiv_2,std::deque<dyad>*dyads_ptr,std::deque <node>*nodes_ptr,std::deque<std::deque<double>>* f_matrix,bool fill_in);
std::deque <node*> get_parent_pool(string sex,int mat_age_m,int mat_age_f,int gest_length, std::deque <node>*all_nodes_ptr, node* indiv, bool twins,string nonparents = "");
string all_nodes_to_population_file(std::deque <node> *all_nodes_ptr,string filename);
std::deque<double> get_all_gaps(std::deque<gap>*gaps_ptr,std::deque<node>*all_nodes_ptr);
void get_start_solution(double temperature,std::deque<gap>*gaps_ptr,std::deque<node>*all_nodes_ptr,std::deque<node>*current_nodes_ptr,int maturation_age_m,int maturation_age_f,int gestation_length,bool twins);
void get_random_full_ped(string file,std::deque<node>*all_nodes_ptr,int maturation_age_f,int maturation_age_m,int gestation_length,bool twins);
int get_full_generation_depth(node* indiv);
#endif