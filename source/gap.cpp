#include "gap.h"

gap::gap(node* indiv,std::deque<string> sexseq,int idx,node* pot_candidate){
    gap::indiv = indiv;
    gap::sexseq = sexseq;
    gap::pot_candidate = pot_candidate;
    gap::idx = idx;
    gap::was_altered = false;
}
node* gap::get_indiv(){
    return gap::indiv;
}
node* gap::get_pot_candidate(){
    return gap::pot_candidate;
}
std::deque<string> gap::get_sexseq(){
    return gap::sexseq;
}
int gap::get_idx(){
    return gap::idx;
}
std::deque<int>* gap::get_dependend_gaps_ptr(){
    return &(gap::dependend_gaps_idx);
}
string gap::get_dependend_gaps_str(){
    if(gap::dependend_gaps_idx.size()>0){
        string dependend_gaps = to_string(dependend_gaps_idx[0]);
        if(gap::dependend_gaps_idx.size()>1){
            for(int i = 1;i<gap::dependend_gaps_idx.size();i++){
                dependend_gaps += ("," + to_string(dependend_gaps_idx[i]));
            }
        }
        return dependend_gaps;
    }else{
        return "no dependend gaps";
    }
}
bool gap::get_was_altered(){
    return gap::was_altered;
}
void gap::set_indiv(node* new_indiv){
    gap::indiv = new_indiv;
}
void gap::set_pot_candidate(node* new_candidate){
    gap::pot_candidate = new_candidate;
}
void gap::set_was_altered(bool altered){
    gap::was_altered = altered;
}
string gap::print_sexseq(){
    string sexseq_str = "";
    for(int i = 0;i<gap::sexseq.size();i++){
        sexseq_str += gap::sexseq[i];
    }
    return sexseq_str;
}
string gap::get_gap_infos(){
    if(gap::get_pot_candidate()==nullptr){
        return "("+get_indiv()->get_name()+", "+print_sexseq()+", gap)";
    }else{
        return "("+get_indiv()->get_name()+", "+print_sexseq()+", "+get_pot_candidate()->get_name()+")";
    }
}