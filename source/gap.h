#ifndef GAP_H
#define GAP_H

#include <deque>
#include <string>

#include "node.h"

using namespace std;

class gap{
protected:
    node* indiv;
    node* pot_candidate;
    std::deque<string> sexseq;
    bool was_altered;
    int idx;
    std::deque<int> dependend_gaps_idx;
public:
    gap(node* indiv=nullptr,std::deque<string> sexseq = {},int idx = 0,node* pot_candidate = nullptr);
    node* get_indiv();
    node* get_pot_candidate();
    std::deque<string> get_sexseq();
    int get_idx();
    std::deque<int>* get_dependend_gaps_ptr();
    string get_dependend_gaps_str();
    bool get_was_altered();
    void set_indiv(node* new_indiv);
    void set_pot_candidate(node* new_candidate);
    void set_was_altered(bool altered);
    string print_sexseq();
    string get_gap_infos();
};
#endif