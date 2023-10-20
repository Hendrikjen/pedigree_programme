#include "dyad.h"

dyad::dyad(string indiv_1,string indiv_2,double r_value){ 
    if(indiv_1<indiv_2){
        dyad::indiv_1=indiv_1;
        dyad::indiv_2=indiv_2;
    }else{
        dyad::indiv_1=indiv_2;
        dyad::indiv_2=indiv_1;
    }
    dyad::dyad_name=dyad::indiv_1+"_"+dyad::indiv_2;
    dyad::r_value=r_value;
    dyad::real_r_value=-1;
    dyad::gap_r_value=-1;
    dyad::indiv_1_idx=-1;
    dyad::indiv_2_idx=-1;
    dyad::min_dyadic_genealogical_depth = -1;
}
string dyad::get_indiv_1_name(){
    return indiv_1;
}
string dyad::get_indiv_2_name(){
    return indiv_2;
}
string dyad::get_dyad_name(){
    return dyad_name;
}
double dyad::get_r_value(){
    return r_value;
}
double dyad::get_real_r_value(){
    return dyad::real_r_value;
}
double dyad::get_gap_r_value(){
    return dyad::gap_r_value;
}
double dyad::get_r_diff(string kind_of_r){ // calc or gap
    if(kind_of_r == "real_minus_r"){
        if(dyad::real_r_value!=-1 and dyad::r_value!=-1){
            return sqrt((dyad::real_r_value - dyad::r_value)*(dyad::real_r_value - dyad::r_value));
        }else{
            return -1;
        }
    }else if(kind_of_r == "real_minus_gap"){
        if(dyad::real_r_value!=-1 and dyad::gap_r_value!=-1){
            return sqrt((dyad::real_r_value - dyad::gap_r_value)*(dyad::real_r_value - dyad::gap_r_value));
        }else{
            return -1;
        }
    }else{
        cout << "kind_of_r error in get_r_diff(). please choose between 'real_minus_r' or 'real_minus_gap'"<<endl;
        return -1;
    }
}
int dyad::get_indiv_1_idx(){
    return indiv_1_idx;
}
int dyad::get_indiv_2_idx(){
    return indiv_2_idx;
}
int dyad::get_min_dyadic_genealogical_depth(){
    return dyad::min_dyadic_genealogical_depth;
}
std::deque <deque <node*>>* dyad::get_paths_ptr(){
    return &(dyad::paths);
}
std::deque <string>* dyad::get_path_characteristics_ptr(){
    return &(dyad::path_characteristics);
}
string dyad::get_depths(){
    if(dyad::path_characteristics.empty()){
        dyad::compute_path_characteristics();
    }
    return dyad::path_characteristics.at(3);
}
string dyad::get_common_ancs(){
    if(dyad::path_characteristics.empty()){
        dyad::compute_path_characteristics();
    }
    return dyad::path_characteristics.at(2);
}
string dyad::get_kinlines(){
    if(dyad::path_characteristics.empty()){
        dyad::compute_path_characteristics();
    }
    return dyad::path_characteristics.at(1);
}
string dyad::get_pathlines(){
    if(dyad::path_characteristics.empty()){
        dyad::compute_path_characteristics();
    }
    return dyad::path_characteristics.at(0);
}
void dyad::set_min_DGD(int min_DGD){
    dyad::min_dyadic_genealogical_depth = min_DGD;
}
void dyad::set_r_value(double r_val){
    dyad::r_value=r_val;
}
void dyad::set_real_r_value(double r_val){
    dyad::real_r_value = r_val;
}
void dyad::set_gap_r_value(double r_val){
    dyad::gap_r_value = r_val;
}
void dyad::set_path_str(string paths){
    dyad::path_str = paths;
}
string dyad::get_path_str(){
    return dyad::path_str;
}
string dyad::print_paths(){
    if(!dyad::paths.empty()){ 
        string output_paths = "";
        for (int i = 0; i<dyad::paths.size(); i++) { // interate over all paths
            string current_path = ""; // prepare output of current path
            if(dyad::paths[i].size()>0){ // if atleast one indiv is saved at path (should be)
                if(dyad::paths[i].size()>1){
                    for(int j = 0;j<((dyad::paths[i]).size()-1);j++){
                        current_path = current_path + dyad::paths[i].at(j)->get_name() + "@";
                    }
                }
                if(!dyad::paths[i].empty()){
                    output_paths = output_paths+ current_path + dyad::paths[i].back()->get_name();
                }
            }
            if(i<(dyad::paths.size()-1)){
                output_paths = output_paths + "/@/";
            }
        }
        return output_paths;
    }else{
        return "NA";
    }
}
string dyad::get_kinlabel(){
    if(dyad::path_characteristics.empty()){
        dyad::compute_path_characteristics();
    }
    if(!dyad::paths.empty() && dyad::paths[0].size()>0){ 
        string kinlabel = "";
        std::deque <string> depths = str_split(dyad::get_depths(),"/@/");
        std::deque <string> grand_counter = {"grand","great","parent","child","mother","daughter","father","son","parsib","nibling","aunt","niece","uncle","nephew","sister","brother"};
        std::deque <string> cousin_counter = {"1st","2nd","3rd","th"};
        std::deque <string> removed_counter = {"once","twice","thrice","times"};
        for (int i = 0; i<dyad::paths.size(); i++) { 
            string current_kinlabel = "";
            std::deque <string> current_depth = str_split(depths[i],"/");
            int depth_1 = 0, depth_2 = 0;
            string sex_1 = "", sex_2 = "";
            if(std::stoi(current_depth[0])<std::stoi(current_depth[1])){
                depth_1 = std::stoi(current_depth[1]);
                depth_2 = std::stoi(current_depth[0]);
                sex_1 = dyad::paths[i].at(dyad::paths[i].size()-1)->get_sex();
                sex_2 = dyad::paths[i].at(0)->get_sex();
            }else{
                depth_2 = std::stoi(current_depth[1]);
                depth_1 = std::stoi(current_depth[0]);
                sex_2 = dyad::paths[i].at(dyad::paths[i].size()-1)->get_sex();
                sex_1 = dyad::paths[i].at(0)->get_sex();
            }
            int difference = depth_1-depth_2;
            if(depth_2 == 0){
                int counter[2] = {3,2};
                if(sex_1=="m"){
                    counter[0] = counter[0]+4;
                }
                if(sex_2=="m"){
                    counter[1] = counter[1]+4;
                }
                if(sex_1=="f"){
                    counter[0] = counter[0]+2;
                }
                if(sex_2=="f"){
                    counter[1] = counter[1]+2;
                }
                while(difference>2){
                    current_kinlabel = current_kinlabel + grand_counter[1] + "-";
                    difference--;
                }
                if(difference==2){
                    current_kinlabel = current_kinlabel + grand_counter[0] + "-";
                }
                current_kinlabel = current_kinlabel + grand_counter[counter[0]] + "|" + current_kinlabel + grand_counter[counter[1]];
            }else if(depth_1 == 1 | depth_2== 1){
                if (depth_1==depth_2) {
                    if (sex_1==sex_2 & sex_2=="m") {
                        current_kinlabel = "brothers";
                    } else if (sex_1==sex_2 & sex_2=="f"){
                        current_kinlabel = "sisters";
                    } else {
                        current_kinlabel = "siblings";
                    }
                }else{
                    int counter[2] = {9,8};
                    if(sex_1=="m"){
                        counter[0] = counter[0]+4;
                    }
                    if(sex_2=="m"){
                        counter[1] = counter[1]+4;
                    }
                    if(sex_1=="f"){
                        counter[0] = counter[0]+2;
                    }
                    if(sex_2=="f"){
                        counter[1] = counter[1]+2;
                    }
                    while(difference>2){
                        current_kinlabel = current_kinlabel + grand_counter[1] + "-";
                        difference--;
                    }
                    if(difference==2){
                        current_kinlabel = current_kinlabel + grand_counter[0];
                    }
                    current_kinlabel = current_kinlabel + grand_counter[counter[0]] + "|" + current_kinlabel + grand_counter[counter[1]];
                }
            }else{
                if(depth_1==depth_2) {
                    if(depth_1>=5){
                        current_kinlabel = to_string(depth_1-1)+cousin_counter[3]+"-cousins";
                    }else{
                        current_kinlabel = cousin_counter[depth_1-2]+"-cousins";
                    }
                }else{
                    if(depth_2>=5){
                        current_kinlabel = to_string(depth_2-1) + cousin_counter[3] + "-cousins";
                    }else{
                        current_kinlabel = cousin_counter[depth_2-2]+"-cousins";
                    }
                    if(difference >=4){
                        current_kinlabel = current_kinlabel + "-" + to_string(difference) + removed_counter[3]+"-removed";
                    }else{
                        current_kinlabel = current_kinlabel + "-" + removed_counter[difference-1]+"-removed"; 
                    }
                }
            }
            if(i<(dyad::paths.size()-1)){
                kinlabel = kinlabel + current_kinlabel + "/@/";
            }else{
                kinlabel = kinlabel + current_kinlabel;
            }
        }
        return kinlabel;
    }else{
        return "nonkin";
    }
}
string dyad::get_full_half(std::deque<std::deque<node*>> paths, string kinlabel_str,string lca_str){
    string output_full_half = "";
    if(paths.size()<1){ // no path at all
        output_full_half = "NA";
    }else if(paths.size()==1){ // only one path
        output_full_half = "half";
    }else{ // at least two paths ( possibly full)
        std::deque<string> kinlabel = str_split(kinlabel_str,"/@/");
        std::deque<string> lca = str_split(lca_str,"/@/");
        std::deque<string> vec_full_half = {};
        for(int i = 0;i<kinlabel.size();i++){ // initialize vec_full_half with "half" for each path
            vec_full_half.push_back("half");
        }
        for(int i = 1;i<kinlabel.size();i++){ // iterate through paths (i == right one; j == left one, since j<i)
            if(vec_full_half[i]!="full"){
                for(int j = 0;j<i;j++){
                    if(kinlabel[i]==kinlabel[j] && lca[i]!=lca[j]){
                        std::deque<node*> pot_full_path_1 = paths[j];
                        std::deque<node*> pot_full_path_2 = paths[i];
                        if(pot_full_path_1.size()==pot_full_path_2.size()){
                            if(pot_full_path_2[0]->get_name() == pot_full_path_1[pot_full_path_1.size()-1]->get_name() && pot_full_path_2[pot_full_path_2.size()-1]->get_name() == pot_full_path_1[0]->get_name()){
                                std::deque<node*> pot_full_path_temp = pot_full_path_1;
                                for(int k = 0;k<pot_full_path_1.size();k++){
                                    pot_full_path_temp[pot_full_path_1.size()-1-k] = pot_full_path_1[k];
                                }
                                swap(pot_full_path_1,pot_full_path_temp);
                            }
                            bool full_path = true;
                            int pot_lca = -1;
                            for(int l = 0;l<pot_full_path_2.size();l++){
                                if(pot_full_path_2[l]->get_name()!=pot_full_path_1[l]->get_name()){// && pot_full_path_2[l]->get_name()!=lca[j] && pot_full_path_1[l]->get_name()!=lca[i]){
                                    if(pot_lca==-1){
                                        pot_lca = l;
                                    }else{
                                        full_path = false;
                                        break;
                                    }
                                }
                            }
                            if(pot_lca == -1 || lca[j]==lca[i] || pot_full_path_2[pot_lca]->get_name()!=lca[i] || pot_full_path_1[pot_lca]->get_name()!=lca[j]){
                                full_path = false;
                            }
                            if(full_path==true){
                                if(vec_full_half[j] == "full" && vec_full_half[i] != "full"){
                                    cout << "three full paths. error."<<endl;
                                }else{
                                    vec_full_half[i] = "full";
                                    vec_full_half[j] = "full";
                                }
                            }
                        }
                    }
                }
            }   
        }
        for (int i = 0; i<vec_full_half.size(); i++) {
            if(i<(vec_full_half.size()-1)){
                output_full_half = output_full_half + vec_full_half[i] + "/@/";
            }else{
                output_full_half = output_full_half + vec_full_half[i];
            }
        }
    }
    return output_full_half;
}
void dyad::compute_path_characteristics(){
    if(!dyad::paths.empty()){ 
        string output_pathlines = "";
        string output_kinlines = "";
        string output_commonancs = "";
        string output_depths = "";
        for (int i = 0; i<dyad::paths.size(); i++) { 
            string current_pathline = ""; 
            string current_kinline = "";
            string current_commonanc = "";
            string current_depth = "";
            std::deque <int> current_BS; 
            if(dyad::paths[i].size()>0){
                for(node* n: dyad::paths[i]){
                    current_pathline = current_pathline + n->get_sex();
                    current_BS.push_back(n->get_birthseason());
                }
                int oldest = 0;
                for(int k = 0;k<current_BS.size();k++){
                    if(current_BS[k]<current_BS[oldest]){
                        oldest = k;
                    }
                }
                current_commonanc = paths[i].at(oldest)->get_name();
                current_depth = to_string(oldest)+"/"+to_string((paths[i].size()-1)-oldest);
                bool mixed = false;
                string comp_c = "";
                if(oldest == 0||oldest == (paths[i].size()-1)){
                    if(oldest==0){
                        comp_c = paths[i].at(0)->get_sex();
                        for(int k = 0;k<(paths[i].size()-1);k++){
                            if(comp_c != paths[i].at(k)->get_sex()){
                                mixed = true;
                            }
                        }
                    }else if(oldest == (paths[i].size()-1)){
                        comp_c = paths[i].at(1)->get_sex();
                        for(int k = 1;k<paths[i].size();k++){
                            if(comp_c != paths[i].at(k)->get_sex()){
                                mixed = true;
                            }
                        }
                    }
                }else{
                    comp_c = paths[i].at(1)->get_sex();
                    for(int k = 1;k<(paths[i].size()-1);k++){
                        if(comp_c != paths[i].at(k)->get_sex()){
                            mixed = true;
                        }
                    }
                }
                if(mixed){
                    current_kinline = "mixed";
                }else if(comp_c=="f"){
                    current_kinline = "mat";
                }else{
                    current_kinline = "pat";
                }
            }
            if(i<(dyad::paths.size()-1)){
                output_pathlines = output_pathlines + current_pathline + "/@/";
                output_depths = output_depths + current_depth + "/@/";
                output_kinlines = output_kinlines + current_kinline + "/@/";
                output_commonancs= output_commonancs + current_commonanc + "/@/";
            }else{
                output_pathlines = output_pathlines + current_pathline;
                output_depths = output_depths + current_depth;
                output_kinlines = output_kinlines + current_kinline;
                output_commonancs = output_commonancs + current_commonanc;
            }
        }
        dyad::path_characteristics = {output_pathlines,output_kinlines,output_commonancs,output_depths};
    }else{
        dyad::path_characteristics = {"NA","NA","NA","NA"};
    }
}
int dyad::add_new_path(deque<node*>* pathdeq_to_add){
    deque <node*> path = {};
    if(pathdeq_to_add!=nullptr){
        path = *pathdeq_to_add; 
    }
    dyad::paths.push_back(path);
    return (dyad::paths.size())-1;
}
void dyad::set_idx(std::deque <node> *all_nodes_ptr){
    for(int i = 0;(dyad::indiv_1_idx==-1||dyad::indiv_2_idx==-1)&&i<all_nodes_ptr->size();i++){
        //cout << "-"<<all_nodes_ptr->at(i).get_name()<<endl;
        if(all_nodes_ptr->at(i).get_name()==indiv_1){
            //cout << ":"<<i<<endl;
            dyad::indiv_1_idx = all_nodes_ptr->at(i).get_matidx();                  
        }else if(all_nodes_ptr->at(i).get_name()==indiv_2){
            //cout << ";"<<i<<endl;
            dyad::indiv_2_idx = all_nodes_ptr->at(i).get_matidx();
        }
    }
    //cout << indiv_1<<", "<<indiv_2<<" | idx: "<<dyad::get_indiv_1_idx()<<","<<dyad::get_indiv_2_idx()<<endl;
}
std::tuple<int, int> dyad::get_dyad_idx_in_f_matrix(){
    if(dyad::indiv_1_idx < dyad::indiv_2_idx){
        return {dyad::indiv_2_idx,dyad::indiv_1_idx};
    }else{
        return {dyad::indiv_1_idx,dyad::indiv_2_idx};
    }
}
string dyad::get_dyad_infos(){
    dyad::compute_path_characteristics();
    std::deque <string> * path_characteristics = dyad::get_path_characteristics_ptr();
    if(r_value>0){
        return get_indiv_1_name()+"\t"+get_indiv_2_name()+"\t"+get_dyad_name()+"\t"+to_string_with_precision(get_r_value(),15)+"\t"+print_paths()+"\t"+get_pathlines()+"\t"+get_kinlines()+"\t"+get_common_ancs()+"\t"+get_depths()+"\t"+get_kinlabel()+"\t"+get_full_half(dyad::paths,get_kinlabel(),get_common_ancs())+"\t"+to_string(get_min_dyadic_genealogical_depth())+"\n";
    }else if (r_value==0){
        return get_indiv_1_name()+"\t"+get_indiv_2_name()+"\t"+get_dyad_name()+"\t0\tNA\tNA\tNA\tNA\tNA\t"+get_kinlabel()+"\tNA\t"+to_string(get_min_dyadic_genealogical_depth())+"\n";
    }else{
        return get_indiv_1_name()+"\t"+get_indiv_2_name()+"\t"+get_dyad_name()+"\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\t"+to_string(get_min_dyadic_genealogical_depth())+"\n";
    }
}
string dyad::get_r_infos(){
    return get_dyad_name()+"\tr   "+to_string_with_precision(get_r_value(),5)+"\treal   "+to_string_with_precision(get_real_r_value(),5)+"\tgap   "+to_string_with_precision(get_gap_r_value(),5)+"\treal-r   "+to_string_with_precision(get_r_diff("real_minus_r"),5)+"\treal-gap   "+to_string_with_precision(get_r_diff("real_minus_gap"),5)+"\n";
}