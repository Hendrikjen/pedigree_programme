#include "node.h" 

node::node(string name, string sex,string mom,string sire,string nonsires,string nondams,string DOB,string DOD,int birthseason,int matidx){
    node::name = name;
    node::sex = sex;
    node::birthseason = birthseason;
    if(mom=="unknown"||mom=="NA"||mom=="UNK"||mom=="unkn_f"){
        node::mom="unkn_f";
    }else{
        node::mom = mom;
    }
    if(sire=="unknown"||sire=="NA"||sire=="UNK"||sire=="unkn_m"){
        node::sire="unkn_m";
    }else{
        node::sire = sire;
    }            
    node::matidx = matidx;
    node::nonsires = nonsires;
    node::nondams = nondams;
    if(DOB.find("-")!=-1){
        std::deque <string> DOB_vec = str_split(DOB,"-");
        node::DOB = datefmt(stoi(DOB_vec[0]),stoi(DOB_vec[1]),stoi(DOB_vec[2]));
    }
    if(DOD.find("-")!=-1){
        std::deque <string> DOD_vec = str_split(DOD,"-");
        node::DOD = datefmt(stoi(DOD_vec[0]),stoi(DOD_vec[1]),stoi(DOD_vec[2]));
    }
} 
string node::get_sex(){
    return sex;
}
string node::get_name(){
    return name;
}
string node::get_mom(){
    return mom;
}
string node::get_sire(){
    return sire;
}
string node::get_parent(string type){// "f" or "m" or "dyad"
    if(type == "f"){
        return node::mom;
    }else if(type == "m"){
        return node::sire;
    }else if(type == "dyad"){
        if(node::mom<node::sire){
            return node::mom+"_"+node::sire;
        }else{
            return node::sire+"_"+node::mom;
        }
    }else{
        return "";
    }
}
string node::get_real_mom(){
    return node::real_mom;
}
string node::get_real_sire(){
    return node::real_sire;
}
string node::get_nonsires(){
    return node::nonsires;
}
string node::get_nondams(){
    return node::nondams;
}
int node::get_birthseason(){
    return node::birthseason;
}
int node::get_matidx(){
    return node::matidx;
}
string node::get_pot_sires(){
    if(node::pot_sires==""){
        return "NA";
    }else{
        return node::pot_sires;
    }
}
string node::get_pot_moms(){
    if(node::pot_moms==""){
        return "NA";
    }else{
        return node::pot_moms;
    }
}
int node::get_full_generations(){
    return node::full_generations;
}
double node::get_min_f(){
    return node::min_f;
}
std::deque<int> node::get_mom_pool(){
    return node::mom_pool;
}
std::deque<int> node::get_sire_pool(){
    return node::sire_pool;
}
void node::set_min_f(double min_f){
    node::min_f = min_f;
}
void node::set_sire(string new_sire){
    node::sire = new_sire;
}
void node::set_mom(string new_mom){
    node::mom = new_mom;
}
void node::set_real_mom(string real_mom){
    node::real_mom = real_mom;
}
void node::set_real_sire(string real_sire){
    node::real_sire = real_sire;
}
void node::set_sire_node(node* new_sire){
    node::sire_node = new_sire;
}
void node::set_mom_node(node* new_mom){
    node::mom_node = new_mom;
}
void node::set_pot_sires(string pot_sires){
    node::pot_sires = pot_sires;
}
void node::set_pot_moms(string pot_moms){
    node::pot_moms = pot_moms;
}
void node::push_back_pot_parent(string sex,int idx){
    if(sex=="f"){
        node::mom_pool.push_back(idx);
    }else if(sex=="m"){
        node::sire_pool.push_back(idx);
    }
}
void node::set_full_generations(int full_generations){
    node::full_generations = full_generations;
}
node* node::get_sire_node(){
    return node::sire_node;
}
node* node::get_mom_node(){
    return node::mom_node;
}
datefmt node::get_DOB(){
    return node::DOB;
}
datefmt node::get_DOD(){
    return node::DOD;
}
void node::set_DOD(string DOD){
    if(DOD.find("-")!=-1){
        std::deque <string> DOD_vec = str_split(DOD,"-");
        node::DOD = datefmt(stoi(DOD_vec[0]),stoi(DOD_vec[1]),stoi(DOD_vec[2]));
    }
}
void node::set_unkn_x(){
    if(node::mom=="unknown"){
        node::mom = "unkn_f";
    }
    if(node::real_mom=="unknown"){
        node::real_mom = "unkn_f";
    }
    if(node::sire=="unknown"){
        node::sire = "unkn_m";
    }
    if(node::real_sire=="unknown"){
        node::real_sire = "unkn_m";
    }
}
string node::compare_pedigree_infos(){
    try{
        node::set_unkn_x();
        if(node::mom==node::real_mom && node::sire==node::real_sire){
            return "true";
        }else{
            return ("false\t[real_mom/mom] " + node::real_mom+"/"+node::mom + "\t[real_sire/sire] "+node::real_sire+"/"+node::sire);
        }
    }catch(const std::exception &ex) {
        std::cerr << "Error in node::compare_pedigree_infos(): " << ex.what() << std::endl;
        return "unable to compare pedigree infos";
    }
}
std::deque<node*>* node::get_affected_nodes_ptr(){ 
    return &(node::affected_nodes); 
}
string node::get_affected_nodes_str(){
    try{
        string affected_nodes_str = "";
        for(int i = 0;i<affected_nodes.size();i++){
            if(i==0){
                affected_nodes_str += affected_nodes[i]->get_name();
            }else{
                affected_nodes_str = affected_nodes_str + "@" + affected_nodes[i]->get_name();
            }
        }
        return affected_nodes_str;
    }catch(const std::exception &ex) {
        std::cerr << "Error in node::get_affected_nodes_str(): " << ex.what() << std::endl;
        return "unable to get affected nodes as string";
    }
}
string node::get_infos(bool nonparent,bool potparent){
    try{
        string info = name;
        if(sex=="f"||sex=="m"){
            info.append(", " + sex);
        }
        else{
            info.append(", unknown sex");
        }
        if(get_DOB().get_date()!="0-0-0"){//birthseason!=0){
            info.append(", born: " + get_DOB().get_date());//to_string(birthseason));
        }
        if(get_DOD().get_date()!="0-0-0"){
            info.append(", died: " + get_DOD().get_date());
        }
        info.append(", mom: " + mom+ ", sire: " + sire+ "\n");
        if(nonparent == true){
            info.append("nonsire: " + nonsires + "\nnondam: " + nondams + "\n");
        }
        if(potparent == true){
            info.append("potential sire: " + pot_sires + "\npotential dam: " + pot_moms + "\n");
        }
        return info;
    }catch(const std::exception &ex) {
        std::cerr << "Error in node::get_infos(): " << ex.what() << std::endl;
        return "unable to get node infos";
    }
}
void node::create_parent_ptr(std::deque <node> *all_nodes_ptr,bool reassigning){ 
    try{
        if(reassigning==false){
            if(node::name==(all_nodes_ptr->at(0)).get_name()||(node::name==(all_nodes_ptr->at(1)).get_name())){ // for imaginary nodes
                node::mom_node = &(all_nodes_ptr->at(0));
                node::sire_node = &(all_nodes_ptr->at(1));
            }
            else{
                for(int i = 0;node::mom_node==nullptr||node::sire_node==nullptr;i++){ // as long as one parent_node is a nullptr, try next i in all_nodes_ptr
                    if((node::mom=="unknown"||node::mom=="unkn_f")&&(all_nodes_ptr->at(i)).get_name()=="unkn_f"){
                        node::mom_node = &(all_nodes_ptr->at(i));
                        continue;
                    }
                    if((all_nodes_ptr->at(i)).get_name()==node::mom){
                        node::mom_node = &(all_nodes_ptr->at(i));
                        continue;
                    }
                    if((node::sire=="unknown"||node::sire=="unkn_m")&&(all_nodes_ptr->at(i)).get_name()=="unkn_m"){
                        node::sire_node = &(all_nodes_ptr->at(i));
                        continue;
                    }
                    if((all_nodes_ptr->at(i)).get_name()==node::sire){
                        node::sire_node = &(all_nodes_ptr->at(i));
                        continue;
                    }
                    if(i>=all_nodes_ptr->size()){
                        cout << "ERROR. Not able to assign a sire_node as well as a mom_node pointer. please check."<<endl;
                        break;
                    }
                }
            }
        }else{
            if(node::name==(all_nodes_ptr->at(0)).get_name()||(node::name==(all_nodes_ptr->at(1)).get_name())){ // for imaginary nodes
                node::mom_node = &(all_nodes_ptr->at(0));
                node::sire_node = &(all_nodes_ptr->at(1));
            }else{
                //cout << node::get_infos()<< "reassigning (mom = "<<node::mom_node->get_name()<<"): "<<node::mom_node << " --> ";
                if(node::mom==all_nodes_ptr->at(0).get_name()){
                    node::mom_node = &(all_nodes_ptr->at(0));
                }else if(node::mom_node != &(all_nodes_ptr->at(0)) && node::mom_node != nullptr){
                    node::mom_node = &(all_nodes_ptr->at(node::mom_node->get_matidx()));
                }else{
                    cout << "WARNING. No priorly assigned mom_node pointer. please check."<<endl;
                }
                //cout << node::mom_node<<" ("<<node::mom_node->get_name()<<")\nand (sire = "<<node::sire_node->get_name()<<"): "<<node::sire_node<<" --> ";
                if(node::sire==all_nodes_ptr->at(1).get_name()){
                    node::sire_node = &(all_nodes_ptr->at(1));
                }else if(node::sire_node != &(all_nodes_ptr->at(1)) && node::sire_node != nullptr){
                    node::sire_node = &(all_nodes_ptr->at(node::sire_node->get_matidx()));
                }else{
                    cout << "WARNING. No priorly assigned sire_node pointer. please check."<<endl;
                }
                //cout << node::sire_node<<" ("<<node::sire_node->get_name()<<indiv_1_idx")"<<endl;
            }  
        }
    }catch(const std::exception &ex) {
        std::cerr << "Error in node::create_parent_ptr(): " << ex.what() << std::endl;
    }
}