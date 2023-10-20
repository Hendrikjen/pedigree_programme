#include "simpath.h"

simpath::simpath(string path_name, double r_val){
    simpath::path_name=path_name;
    simpath::r_val=r_val;
    std::deque<string> path_indivs = str_split(simpath::path_name,"_");
    simpath::first = path_indivs[0];
    simpath::last = path_indivs[path_indivs.size()-1];
    if(simpath::first<simpath::last){
        simpath::dyad_name = simpath::first+"_"+simpath::last;
    }else{
        simpath::dyad_name = simpath::last+"_"+simpath::first;
    }
}
string simpath::get_path_name(){
    return simpath::path_name;
}
double simpath::get_r_val(){
    return simpath::r_val;
}
string simpath::get_first(){
    return simpath::first;
}
string simpath::get_last(){
    return simpath::last;
}
string simpath::get_dyad_name(){
    return simpath::dyad_name;
}
string simpath::get_info(){
    return "[path] "+simpath::path_name+"\t[r] "+to_string(simpath::r_val);
}