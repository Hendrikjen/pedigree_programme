#ifndef POPULATION_SIMULATION_H
#define POPULATION_SIMULATION_H

#include <string>
#include <deque>
#include <map>
#include <iostream>
#include <fstream>

#include "datefmt.h"
#include "node.h"
#include "general_class_functions.h"
#include "simpath.h"

using namespace std; 

string generate_name(int number);
std::deque <int> filter_for_pot_anc(string sex,int maturation_age,int current_year,std::deque<node>*all_nodes_ptr,int max_age);
void expand_paths(std::deque<simpath>parent_paths,std::string parent_sex,node* offspring,std::deque<simpath>*all_dyadic_paths_ptr,std::map<string,std::deque<simpath>>*path_dict_ptr);
void write_simulated_r_values_by_dyads(int simulation_duration,int start_individuals,std::deque<simpath>*all_dyadic_paths_ptr);
string population_simulation(int simulation_duration, int start_individuals,int gestation_length,int maturation_age_f,int maturation_age_m,int max_age, int default_year);
string add_parental_gaps(string file,double mat_gaps, double pat_gaps, int start_individuals);
#endif