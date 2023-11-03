#ifndef NODE_H
#define NODE_H

#include <string>
#include <deque>  

#include "datefmt.h"
#include "basic_functions.h"

using namespace std;

class node {
protected:
    string name = "",sex = "",mom = "",sire = "",nonsires = "",nondams = "",real_mom = "",real_sire = "",pot_sires = "",pot_moms = "";
    int birthseason = 0,matidx = 0,full_generations = 0;
    double min_f = 0;
    datefmt DOB, DOD;
    std::deque<node*> affected_nodes = {};
    std::deque<int>mom_pool = {}, sire_pool = {};
public:
    node* mom_node=nullptr;
    node* sire_node=nullptr;
    node(string name="NA", string sex="NA",string mom="NA",string sire="NA",string nonsires="NA",string nondams="NA",string DOB = "NA",string DOD = "NA",int birthseason=0,int matidx = 0);
    string get_sex();
    string get_name();
    string get_mom(string unk = "unkn_f");
    string get_sire(string unk = "unkn_m");
    string get_parent(string type);
    string get_real_mom();
    string get_real_sire();
    string get_nonsires();
    string get_nondams();
    int get_birthseason();
    int get_matidx();
    string get_pot_sires();
    string get_pot_moms();
    int get_full_generations();
    double get_min_f();
    std::deque<int> get_mom_pool();
    std::deque<int> get_sire_pool();
    void set_min_f(double min_f);
    void set_sire(string new_sire);
    void set_mom(string new_mom);
    void set_real_mom(string real_mom);
    void set_real_sire(string real_sire);
    void set_sire_node(node* new_sire);
    void set_mom_node(node* new_mom);
    void set_pot_sires(string pot_sires);
    void set_pot_moms(string pot_moms);
    void push_back_pot_parent(string sex,int idx);
    void set_full_generations(int full_generations);
    node* get_sire_node();
    node* get_mom_node();
    datefmt get_DOB();
    datefmt get_DOD();
    void set_DOD(string DOD);
    void set_unkn_x();
    string compare_pedigree_infos();
    std::deque<node*>* get_affected_nodes_ptr();
    string get_affected_nodes_str();
    string get_infos(bool nonparent = false,bool potparent = false);
    void create_parent_ptr(std::deque <node> *all_nodes_ptr,bool reassigning = false);
};

#endif
