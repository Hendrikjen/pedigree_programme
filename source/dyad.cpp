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
double dyad::get_r_diff(string kind_of_r){  //calculates the difference between two relatedness coefficients based on the given 'kind_of_r'
// options kind_of_r: The type of comparison to perform -> 'real_minus_r' (|realized r - current r|) or 'real_minus_gap' (|realized r - pedigree r with gaps|)
    try{
        if(kind_of_r == "real_minus_r"){
            if(dyad::real_r_value!=-1 and dyad::r_value!=-1){
                return sqrt((dyad::real_r_value - dyad::r_value)*(dyad::real_r_value - dyad::r_value));
            }else{
                throw std::runtime_error("One or more relatedness values are undefined.");
            }
        }else if(kind_of_r == "real_minus_gap"){
            if(dyad::real_r_value!=-1 and dyad::gap_r_value!=-1){
                return sqrt((dyad::real_r_value - dyad::gap_r_value)*(dyad::real_r_value - dyad::gap_r_value));
            }else{
                throw std::runtime_error("One or more relatedness values are undefined.");
            }
        }else{
            throw std::invalid_argument("Invalid kind_of_r value. Please choose between 'real_minus_r' or 'real_minus_gap'");
        }
    }catch(const std::exception &ex) {
        std::cerr << "Error in dyad::get_r_diff(): " << ex.what() << std::endl;
        return 999.0;
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
void limit_path_characteristics_to_genealogical_depth(int generation){
    
}
string dyad::get_depths(int generation_limit){
    if(dyad::path_characteristics.empty()){
        dyad::compute_path_characteristics(generation_limit);
    }
    if((dyad::path_characteristics.at(3).length()>3)&&(dyad::path_characteristics.at(3).substr(dyad::path_characteristics.at(3).length() - 3) == "/@/")){
        return dyad::path_characteristics.at(3).substr(0, dyad::path_characteristics.at(3).length() - 3);
    }else{
        return dyad::path_characteristics.at(3);
    }
}
string dyad::get_common_ancs(int generation_limit){
    if(dyad::path_characteristics.empty()){
        dyad::compute_path_characteristics(generation_limit);
    }
    if((dyad::path_characteristics.at(2).length()>3)&&(dyad::path_characteristics.at(2).substr(dyad::path_characteristics.at(2).length() - 3) == "/@/")){
        return dyad::path_characteristics.at(2).substr(0, dyad::path_characteristics.at(2).length() - 3);
    }else{
        return dyad::path_characteristics.at(2);
    }
}
string dyad::get_kinlines(int generation_limit){
    if(dyad::path_characteristics.empty()){
        dyad::compute_path_characteristics(generation_limit);
    }
    if((dyad::path_characteristics.at(1).length()>3)&&(dyad::path_characteristics.at(1).substr(dyad::path_characteristics.at(1).length() - 3) == "/@/")){
        return dyad::path_characteristics.at(1).substr(0, dyad::path_characteristics.at(1).length() - 3);
    }else{
        return dyad::path_characteristics.at(1);
    }
}
string dyad::get_pathlines(int generation_limit){
    if(dyad::path_characteristics.empty()){
        dyad::compute_path_characteristics(generation_limit);
    }
    if((dyad::path_characteristics.at(0).length()>3)&&(dyad::path_characteristics.at(0).substr(dyad::path_characteristics.at(0).length() - 3) == "/@/")){
        return dyad::path_characteristics.at(0).substr(0, dyad::path_characteristics.at(0).length() - 3);
    }else{
        return dyad::path_characteristics.at(0);
    }
}
string dyad::get_path_str(){
    return dyad::path_str;
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
string dyad::print_paths(){ // generates a string representation of paths (A,B and A,D,B to A@B/@/A@D@B), if paths is empty, return NA
    try{
        if(!dyad::paths.empty()){ 
            string output_paths = "";
            for (int i = 0; i<dyad::paths.size(); i++) { // interate over all paths (A@B, A@D@B)
                string current_path = ""; 
                if(dyad::paths[i].size()>0){ // if at least one indiv is saved within current path (should be)
                    if(dyad::paths[i].size()>1){ // if multiple individuals are saved within current path
                        for(int j = 0;j<((dyad::paths[i]).size()-1);j++){
                            current_path = current_path + dyad::paths[i].at(j)->get_name() + "@"; // string  them together with @ as delimiter
                        }
                    }
                    if(!dyad::paths[i].empty()){
                        output_paths = output_paths+ current_path + dyad::paths[i].back()->get_name(); // string current path to prior processed paths 
                    }
                }
                if(i<(dyad::paths.size()-1)){
                    output_paths = output_paths + "/@/"; // use /@/ as delimiter between multiple paths
                }
            }
            return output_paths;
        }else{
            return "NA";
        }
    }catch(const std::exception &ex) {
        std::cerr << "Error in dyad::print_paths(): " << ex.what() << std::endl;
        return "unable to print path";
    }
}
string dyad::get_kinlabel(int generation_limit){ // generates kinlabel based on depth & sex
    try{
        if(dyad::path_characteristics.empty()){
            dyad::compute_path_characteristics(generation_limit); // ensure path characteristics were calculation previously
        }
        if(!dyad::paths.empty() && dyad::paths[0].size()>0){ // if paths to analyze exist
            string kinlabel = ""; // to collect all final kinlabels
            std::deque <string> depths = str_split(dyad::get_depths(generation_limit),"/@/");
            std::deque <string> grand_counter = {"grand","great","parent","child","mother","daughter","father","son","parsib","nibling","aunt","niece","uncle","nephew","sister","brother"};
            std::deque <string> cousin_counter = {"1st","2nd","3rd","th"};
            std::deque <string> removed_counter = {"once","twice","thrice","times"};
            for (int i = 0; i<dyad::paths.size(); i++) { // iterate through paths
                string current_kinlabel = ""; // to collect infos about the current kinlabel
                std::deque <string> current_depth = str_split(depths[i],"/");
                int depth_1 = 0, depth_2 = 0; // init depths
                string sex_1 = "", sex_2 = ""; // init sexes
                bool swapped = false;
                if(std::stoi(current_depth[0])<std::stoi(current_depth[1])){ //get depths & sexes in the correct order so that depth_2 <= depth_1
                    depth_1 = std::stoi(current_depth[1]);
                    depth_2 = std::stoi(current_depth[0]);
                    sex_1 = dyad::paths[i].at(dyad::paths[i].size()-1)->get_sex();
                    sex_2 = dyad::paths[i].at(0)->get_sex();
                }else{
                    depth_2 = std::stoi(current_depth[1]);
                    depth_1 = std::stoi(current_depth[0]);
                    sex_2 = dyad::paths[i].at(dyad::paths[i].size()-1)->get_sex();
                    sex_1 = dyad::paths[i].at(0)->get_sex();
                    swapped = true;
                }
                int difference = depth_1-depth_2; // get difference between both depths
                if(depth_2 == 0){ // if one focal is the ancestor of the other
                    int counter[2] = {3,2}; // set both focal labels on the sex independent version: 3 and 2 which later will codes for child and parent; if sex info is available update the kinlabel counter
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
                    while(difference>2){ // add "great" as often as the difference <= 2
                        current_kinlabel = current_kinlabel + grand_counter[1] + "-";
                        difference--;
                    }
                    if(difference==2){ // add "grand"
                        current_kinlabel = current_kinlabel + grand_counter[0] + "-";
                    }
                    if(swapped){ // string everything together in the correct order with | as delimiter between both focals (e.g. great-grand-parent|great-grand-child)
                        current_kinlabel = current_kinlabel + grand_counter[counter[1]] + "|" + current_kinlabel + grand_counter[counter[0]];
                    }else{
                        current_kinlabel = current_kinlabel + grand_counter[counter[0]] + "|" + current_kinlabel + grand_counter[counter[1]];
                    }
                }else if(depth_1 == 1 | depth_2 == 1){ 
                    if (depth_1==depth_2) { // if the focals are siblings
                        if (sex_1==sex_2 & sex_2=="m") {
                            current_kinlabel = "brothers";
                        } else if (sex_1==sex_2 & sex_2=="f"){
                            current_kinlabel = "sisters";
                        } else {
                            current_kinlabel = "siblings";
                        }
                    }else{ // if depth_1 > 1
                        int counter[2] = {9,8}; // set both focal labels on the sex independent version: 9 and 8 which later will codes for nibling (niece/nephew) and parsib (parent's sibling: aunt/uncle); if sex info is available update the kinlabel counter
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
                        while(difference>2){// add "great" as often as the difference <= 2
                            current_kinlabel = current_kinlabel + grand_counter[1] + "-";
                            difference--;
                        }
                        if(difference==2){ // add "grand"
                            current_kinlabel = current_kinlabel + grand_counter[0];
                        }
                        if(swapped){ // string everything together in the correct order with | as delimiter between both focals (e.g. great-grand-aunt|great-grand-niece)
                            current_kinlabel = current_kinlabel + grand_counter[counter[1]] + "|" + current_kinlabel + grand_counter[counter[0]];
                        }else{
                            current_kinlabel = current_kinlabel + grand_counter[counter[0]] + "|" + current_kinlabel + grand_counter[counter[1]];
                        }
                        
                    }
                }else{
                    if(depth_1==depth_2) { // if focals are cousins (unreduced)
                        if(depth_1>=5){ // from 5 on (because 2 = first cousins, 3 = second cousins, 4 = third cousins) use number + 'th'
                            current_kinlabel = to_string(depth_1-1)+cousin_counter[3]+"-cousins";
                        }else{ // if 2 = first cousins, 3 = second cousins, 4 = third cousins
                            current_kinlabel = cousin_counter[depth_1-2]+"-cousins";
                        }
                    }else{ // if focals are reduced cousins
                        if(depth_2>=5){ // as above, but the main cousin characterization is determined by the lower depth
                            current_kinlabel = to_string(depth_2-1) + cousin_counter[3] + "-cousins";
                        }else{ // if 2 = first cousins, 3 = second cousins, 4 = third cousins
                            current_kinlabel = cousin_counter[depth_2-2]+"-cousins";
                        }
                        if(difference >=4){ // add the difference as 'reduced'; from 4 on use number +'times'
                            current_kinlabel = current_kinlabel + "-" + to_string(difference) + removed_counter[3]+"-removed";
                        }else{ //add the difference as 'reduced'; 1 = once, 2 = twice, 3 = thrice
                            current_kinlabel = current_kinlabel + "-" + removed_counter[difference-1]+"-removed"; 
                        }
                    }
                }
                if(i<(dyad::paths.size()-1)){// add current kinlabel to the previous kinlabels of the dyad, delimited by /@/
                    kinlabel = kinlabel + current_kinlabel + "/@/"; 
                }else{
                    kinlabel = kinlabel + current_kinlabel;
                }
            }
            return kinlabel;
        }else{ // if no relatedness path exist (and r = 0)
            return "nonkin";
        }
    }catch(const std::exception &ex) {
        std::cerr << "Error in dyad::get_kinlabel(): " << ex.what() << std::endl;
        return "unable to get kinlabel";
    }
}
string dyad::get_full_half(std::deque<std::deque<node*>> paths, string kinlabel_str,string lca_str){ // detemine if identical paths exist, except for the common ancestor; e.g. to determine whether siblings are full- or half-siblings
    try{
        string output_full_half = "";
        if(paths.size()<1){ // no path at all
            output_full_half = "NA";
        }else if(paths.size()==1){ // only one path --> half
            output_full_half = "half";
        }else{ // at least two paths ( possibly full)
            std::deque<string> kinlabel = str_split(kinlabel_str,"/@/");
            std::deque<string> lca = str_split(lca_str,"/@/");
            std::deque<string> vec_full_half = {};
            for(int i = 0;i<kinlabel.size();i++){ // initialize vec_full_half with "half" for each path
                vec_full_half.push_back("half");
            }
            for(int i = 1;i<kinlabel.size();i++){ // iterate through paths (i == 'right' one; j == 'left' one, since j<i)
                if(vec_full_half[i]!="full"){ 
                    for(int j = 0;j<i;j++){ // iterate through all paths from the beginning to i-path
                        if(kinlabel[i]==kinlabel[j] && lca[i]!=lca[j]){ // if kinlabels match, but not the LCA 
                            std::deque<node*> pot_full_path_1 = paths[j]; // get boths paths as deque
                            std::deque<node*> pot_full_path_2 = paths[i];
                            if(pot_full_path_1.size()==pot_full_path_2.size()){ // if both paths are of the same length
                                if(pot_full_path_2[0]->get_name() == pot_full_path_1[pot_full_path_1.size()-1]->get_name() && pot_full_path_2[pot_full_path_2.size()-1]->get_name() == pot_full_path_1[0]->get_name()){ // if both paths are in reverse orders -> inverse one path (backwards)
                                    std::deque<node*> pot_full_path_temp = pot_full_path_1;
                                    for(int k = 0;k<pot_full_path_1.size();k++){
                                        pot_full_path_temp[pot_full_path_1.size()-1-k] = pot_full_path_1[k];
                                    }
                                    swap(pot_full_path_1,pot_full_path_temp);
                                }
                                bool full_path = true; // default: full
                                int pot_lca = -1;
                                for(int l = 0;l<pot_full_path_2.size();l++){ // iterate through the paths and get the first node which does not match (should be the lca)
                                    if(pot_full_path_2[l]->get_name()!=pot_full_path_1[l]->get_name()){
                                        if(pot_lca==-1){
                                            pot_lca = l;
                                        }else{ // if no mismatching node can be found
                                            full_path = false;
                                            break;
                                        }
                                    }
                                } 
                                if(pot_lca == -1 || lca[j]==lca[i] || pot_full_path_2[pot_lca]->get_name()!=lca[i] || pot_full_path_1[pot_lca]->get_name()!=lca[j]){// if LCAs does not differ, or if the potential LCA is not the actual LCA -> not full paths
                                    full_path = false;
                                }
                                if(full_path==true){
                                    if(vec_full_half[j] == "full" && vec_full_half[i] != "full"){ // if one of the paths is already claimed as "full" means that there are three "full" paths, which throws an error
                                        throw std::runtime_error("Three 'full' paths are not possible");
                                    }else{ // set both paths as "full"
                                        vec_full_half[i] = "full";
                                        vec_full_half[j] = "full";
                                    }
                                }
                            }
                        }
                    }
                }   
            }
            for (int i = 0; i<vec_full_half.size(); i++) { // string the full/half info alltogether, delimited by /@/
                if(i<(vec_full_half.size()-1)){
                    output_full_half = output_full_half + vec_full_half[i] + "/@/";
                }else{
                    output_full_half = output_full_half + vec_full_half[i];
                }
            }
        }
        return output_full_half;
    }catch(const std::exception &ex) {
        std::cerr << "Error in dyad::get_full_half(): " << ex.what() << std::endl;
        return "unable to get full/half";
    }
}
void dyad::compute_path_characteristics(int generation_limit){ // set pathline, kinline, lowest common ancestor and depth as dyad attributes from path
    try{
        if(!dyad::paths.empty()){ // if paths to analyze exist
            std::set<int> paths_to_remove;
            string output_pathlines = "";
            string output_kinlines = "";
            string output_commonancs = "";
            string output_depths = "";
            for (int i = 0; i<dyad::paths.size(); i++) {  // iterate through all paths
                string current_pathline = ""; 
                string current_kinline = "";
                string current_commonanc = "";
                string current_depth = "";
                std::deque <int> current_BS; 
                if(dyad::paths[i].size()>0){ // if path consists of at least one individual (should be)
                    for(node* n: dyad::paths[i]){ // iterate throuch all node along the path to get the sexes for pathline 
                        current_pathline = current_pathline + n->get_sex();
                        current_BS.push_back(n->get_birthseason()); // get also all birth years/birth season
                    }
                    int oldest = 0;
                    for(int k = 0;k<current_BS.size();k++){ // determine the index of the oldest individual of all nodes along the path based on birth season data
                        if(current_BS[k]<current_BS[oldest]){
                            oldest = k;
                        }
                    }
                    current_commonanc = paths[i].at(oldest)->get_name(); // the oldest individual == lowest common ancestor
                    current_depth = to_string(oldest)+"/"+to_string((paths[i].size()-1)-oldest); // depth == index of LCA/(pathlength-LCAindex)
                    bool mixed = false; 
                    string comp_c = "";
                    if(oldest == 0||oldest == (paths[i].size()-1)){ // if one focal is the ancestor of the other one
                        if(oldest==0){ // if the LCA is the first individual in the path
                            comp_c = paths[i].at(0)->get_sex(); // get sex of LCA
                            for(int k = 0;k<(paths[i].size()-1);k++){ // iterate through all individuals except the last one (focal) 
                                if(comp_c != paths[i].at(k)->get_sex()){ // if another sex than the sex of LCA occur, set mixed == true
                                    mixed = true;
                                }
                            }
                        }else if(oldest == (paths[i].size()-1)){ // the LCA is the last individual in the path
                            comp_c = paths[i].at(1)->get_sex(); // get sex of LCA
                            for(int k = 1;k<paths[i].size();k++){// iterate through all individuals except the first one (focal) 
                                if(comp_c != paths[i].at(k)->get_sex()){// if another sex than the sex of LCA occurs, set mixed == true
                                    mixed = true;
                                }
                            }
                        }
                    }else{ // if both focal's sexes have to be excluded (because noone is the ancestor of the other and therefore not part of the kinline)
                        comp_c = paths[i].at(1)->get_sex(); // get sex of the second individual along the path
                        for(int k = 1;k<(paths[i].size()-1);k++){ //iterate through all individuals except the first and the last one (focals) 
                            if(comp_c != paths[i].at(k)->get_sex()){// if another sex than the sex of the second individual occurs, set mixed == true
                                mixed = true;
                            }
                        }
                    }
                    if(mixed){ // ancestors are female(s) as well as male(s)
                        current_kinline = "mixed";
                    }else if(comp_c=="f"){ // ancestors are solely female(s)
                        current_kinline = "mat";
                    }else if(comp_c=="m"){ // ancestors are solely male(s)
                        current_kinline = "pat";
                    }else{ // an unknown sex occurs along the path; no certain declaration is possible
                        current_kinline = "unk";
                    }
                }
                if(generation_limit<0 || (stoi(str_split(current_depth,"/")[0]) < generation_limit && stoi(str_split(current_depth,"/")[1]) < generation_limit )){ // if required save only these path characteristics of paths within the given generation_limit (no generation limit == -1)
                    if(i<(dyad::paths.size()-1)){ //add current attribute to the previous attributes of the dyad, delimited by /@/
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
                }else{
                    paths_to_remove.insert(i);
                }
            }
            dyad::path_characteristics = {output_pathlines,output_kinlines,output_commonancs,output_depths};
            if (!paths_to_remove.empty()) {
                std::deque<std::deque<node*>> temp_paths;
                for (int i = 0; i < dyad::paths.size(); i++) {
                    if (paths_to_remove.find(i) == paths_to_remove.end()) {
                        temp_paths.push_back(dyad::paths[i]);
                    } else {
                        dyad::set_r_value(dyad::get_r_value() - std::pow(0.5, (dyad::paths[i].size() - 1)));
                    }
                }
                dyad::paths = temp_paths;
            }
        }else{ // no paths at all
            dyad::path_characteristics = {"NA","NA","NA","NA"};
        }
    }catch(const std::exception &ex) {
        std::cerr << "Error in dyad::compute_path_characteristics(): " << ex.what() << std::endl;
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
void dyad::set_idx(std::deque <node> *all_nodes_ptr){ // link node index from all_nodes to dyad individuals
    try{
        for(int i = 0;(dyad::indiv_1_idx==-1||dyad::indiv_2_idx==-1)&&i<all_nodes_ptr->size();i++){
            if(all_nodes_ptr->at(i).get_name()==indiv_1){
                dyad::indiv_1_idx = all_nodes_ptr->at(i).get_matidx();                  
            }else if(all_nodes_ptr->at(i).get_name()==indiv_2){
                dyad::indiv_2_idx = all_nodes_ptr->at(i).get_matidx();
            }
        }
    }catch(const std::exception &ex) {
        std::cerr << "Error in dyad::set_idx(): " << ex.what() << std::endl;
    }
}
std::tuple<int, int> dyad::get_dyad_idx_in_f_matrix(){
    if(dyad::indiv_1_idx < dyad::indiv_2_idx){
        return {dyad::indiv_2_idx,dyad::indiv_1_idx};
    }else{
        return {dyad::indiv_1_idx,dyad::indiv_2_idx};
    }
}
string dyad::get_dyad_infos(int generation_limit){ // generate string representation of all dyad attributes (e.g. to write dyadfile)
    try{
        dyad::compute_path_characteristics(generation_limit); // ensure all path characteristics are computed
        std::deque <string> * path_characteristics = dyad::get_path_characteristics_ptr();
        if(r_value>0){ // if dyad is related -> return general dyad information + all attributes
            return get_indiv_1_name()+"\t"+get_indiv_2_name()+"\t"+get_dyad_name()+"\t"+to_string_with_precision(get_r_value(),15)+"\t"+print_paths()+"\t"+get_pathlines(generation_limit)+"\t"+get_kinlines(generation_limit)+"\t"+get_common_ancs(generation_limit)+"\t"+get_depths(generation_limit)+"\t"+get_kinlabel(generation_limit)+"\t"+get_full_half(dyad::paths,get_kinlabel(generation_limit),get_common_ancs(generation_limit))+"\t"+to_string(get_min_dyadic_genealogical_depth())+"\n";
        }else if (r_value==0){ // if dyad is not related -> only dyad info + kinlabel & minDGD
            return get_indiv_1_name()+"\t"+get_indiv_2_name()+"\t"+get_dyad_name()+"\t0\tNA\tNA\tNA\tNA\tNA\t"+get_kinlabel(generation_limit)+"\tNA\t"+to_string(get_min_dyadic_genealogical_depth())+"\n";
        }else{ // if no information is available whether dyad is related or not -> only dyad info + minDGD
            return get_indiv_1_name()+"\t"+get_indiv_2_name()+"\t"+get_dyad_name()+"\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\t"+to_string(get_min_dyadic_genealogical_depth())+"\n";
        }
    }catch(const std::exception &ex) {
        std::cerr << "Error in dyad::get_dyad_infos(): " << ex.what() << std::endl;
        return "unable to get dyad info";
    }
}
string dyad::get_r_infos(){ // generate string representation of available relatedness coefficient information (r value from current pedigree, realized r value, original r value from pedigree with gaps) + their difference
    return get_dyad_name()+"\tr   "+to_string_with_precision(get_r_value(),5)+"\treal   "+to_string_with_precision(get_real_r_value(),5)+"\tgap   "+to_string_with_precision(get_gap_r_value(),5)+"\treal-r   "+to_string_with_precision(get_r_diff("real_minus_r"),5)+"\treal-gap   "+to_string_with_precision(get_r_diff("real_minus_gap"),5)+"\n";
}