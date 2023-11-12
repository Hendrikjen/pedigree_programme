#ifndef DYAD_H
#define DYAD_H

#include <deque>
#include <string>
#include <cmath>
#include <set>

#include "basic_functions.h"
#include "node.h"

using namespace std;

class dyad{
protected:
    double r_value,real_r_value,gap_r_value; 
    int indiv_1_idx,indiv_2_idx,min_dyadic_genealogical_depth;
    string indiv_1,indiv_2,dyad_name,path_str;
    std::deque <deque <node*>> paths;
    std::set <int> filtered_path_indices;
    std::deque <string> path_characteristics;
public:
    dyad(string indiv_1,string indiv_2,double r_value=-1);
    string get_indiv_1_name();
    string get_indiv_2_name();
    string get_dyad_name();
    double get_r_value();
    double get_real_r_value();
    double get_gap_r_value();
    std::set <int>* get_filtered_path_indices();
    double get_r_diff(string kind_of_r);
    int get_indiv_1_idx();
    int get_indiv_2_idx();
    int get_min_dyadic_genealogical_depth();
    std::deque <deque <node*>>* get_paths_ptr();
    std::deque <string>* get_path_characteristics_ptr();
    string get_depths(int generation_limit);
    string get_common_ancs(int generation_limit);
    string get_kinlines(int generation_limit);
    string get_pathlines(int generation_limit);
    void set_min_DGD(int min_DGD);
    void set_r_value(double r_val);
    void set_real_r_value(double r_val);
    void set_gap_r_value(double r_val);
    void set_path_str(string paths);
    void set_filtered_path_indices(std::set<int>indices);
    string get_path_str();
    string print_paths();
    string get_kinlabel(int generation_limit);
    string get_full_half(std::deque<std::deque<node*>> paths, string kinlabel_str,string lca_str);
    void compute_path_characteristics(int generation_limit);
    int add_new_path(deque<node*>* pathdeq_to_add = nullptr);
    void set_idx(std::deque <node> *all_nodes_ptr);
    std::tuple<int, int> get_dyad_idx_in_f_matrix();
    string get_dyad_infos(int generation_limit);
    string get_r_infos();
};
#endif