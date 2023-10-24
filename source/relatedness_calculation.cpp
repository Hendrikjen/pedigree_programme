#include "relatedness_calculation.h"

double calculate_f_xy(node * indiv_1,node * indiv_2,std::deque<std::deque<double>> *f_matrix,int ncol,int nrow,std::deque <node> * all_nodes_ptr,std::deque <dyad> * all_dyads_ptr,int dyad_index,std::deque <int> * path_indices_ptr,int path_index,map <string,int> * dyad_dict_ptr, std::set<string>* node_space_ptr){ // calculates dyadic relatedness coeffiecients recursively and collects path characteristics  
    try{
        int i = (*indiv_1).get_matidx();
        int j = (*indiv_2).get_matidx();
        if(i<j){
            swap(i,j);
        }
        double f_xy = f_matrix->at(i)[j];// get (default) r value from f_matrix
        int next_path_index = path_index;
        std::deque <deque<node*>>* paths_ptr = all_dyads_ptr->at(dyad_index).get_paths_ptr();
        if(f_xy==1){ // indiv_1 == indiv_2
            if(path_index>path_indices_ptr->back()){ // if new path is needed -> add new path & new path_index
                path_indices_ptr->push_back(all_dyads_ptr->at(dyad_index).add_new_path());
            }
            paths_ptr->at(path_index).push_back(indiv_1); // indiv_1/2 correlates with common ancestor -> push
        }else if(f_xy<1&f_xy>0){ // i.e. dyad of interest was already computed for prior dyad -> info & paths can be found in dyadlist
            string dyad_to_search = "";
            if((*indiv_1).get_name()<(*indiv_2).get_name()){ // determine dyad_name in the right order of individuals
                dyad_to_search = (*indiv_1).get_name()+"_"+(*indiv_2).get_name();
            }else{
                dyad_to_search = (*indiv_2).get_name()+"_"+(*indiv_1).get_name();
            }
            if(dyad_dict_ptr->count(dyad_to_search)>0){ // only if dyadlist contains dyad already
                int used_dyad_index = (*dyad_dict_ptr)[dyad_to_search]; // determine at which index in dyadlist dyad_to_search can be found              
                std::deque<std::deque<node*>> * used_paths_ptr = all_dyads_ptr->at(used_dyad_index).get_paths_ptr(); // get paths from dyad_to_search
                for(int k=0;k<used_paths_ptr->size();k++){ // iterate through paths to copy paths from used_dyad to current paths_ptr
                    if(k==0){
                        if(path_index>path_indices_ptr->back()){// if new path is needed -> add new path & new path_index
                            path_indices_ptr->push_back(all_dyads_ptr->at(dyad_index).add_new_path());
                        }
                        if((*indiv_1).get_name()==used_paths_ptr->at(k).at(0)->get_name()){ // decide whether the nodes has to be added ahead or behind
                            for(node* n: used_paths_ptr->at(k)){
                                paths_ptr->at(path_index).push_back(n);
                            }
                        }else{
                            for(node* n: used_paths_ptr->at(k)){
                                paths_ptr->at(path_index).push_front(n);
                            }
                        }
                    }else{
                        int new_index = all_dyads_ptr->at(dyad_index).add_new_path();
                        path_indices_ptr->push_back(new_index);
                        if((*indiv_1).get_name()==used_paths_ptr->at(k).at(0)->get_name()){// decide whether the nodes has to be added ahead or behind
                            for(node* n: used_paths_ptr->at(k)){
                                paths_ptr->at(new_index).push_back(n);
                            }
                        }else{
                            for(node* n: used_paths_ptr->at(k)){
                                paths_ptr->at(new_index).push_front(n);
                            }
                        }
                    }
                }
            }
        }else if(f_xy!=0){ // no entry in f_matrix ([i][j] == 999):
            //if(indivs_in_node_space(indiv_1,indiv_2,node_space_ptr)==false){
            //    f_xy = 0;
            //}else{
                if((*indiv_1).get_mom()==(*indiv_2).get_name()||(*indiv_1).get_sire()==(*indiv_2).get_name()||(*indiv_2).get_mom()==(*indiv_1).get_name()||(*indiv_2).get_sire()==(*indiv_1).get_name()){ // if parent/offspring 
                    if(path_index>path_indices_ptr->back()){// if new path is needed -> add new path & new path_index
                        path_indices_ptr->push_back(all_dyads_ptr->at(dyad_index).add_new_path());
                    }
                    int commit_index = path_indices_ptr->back()+1;
                    double xp = 0;
                    if((*indiv_1).get_mom()==(*indiv_2).get_name()){  // offspring|mom 
                        xp = calculate_f_xy((*indiv_1).sire_node,indiv_2,f_matrix,ncol,nrow,all_nodes_ptr,all_dyads_ptr,dyad_index,path_indices_ptr,commit_index,dyad_dict_ptr, node_space_ptr);
                        if(xp!=0){ 
                            for(int k = commit_index;k<=path_indices_ptr->back();k++){
                                paths_ptr->at(path_indices_ptr->at(k)).push_front(indiv_1);
                            }
                        }
                    }else if((*indiv_1).get_sire()==(*indiv_2).get_name()){  // offspring|sire
                        xp = calculate_f_xy((*indiv_1).mom_node,indiv_2,f_matrix,ncol,nrow,all_nodes_ptr,all_dyads_ptr,dyad_index,path_indices_ptr,commit_index,dyad_dict_ptr,node_space_ptr);
                        if(xp!=0){ 
                            for(int k = commit_index;k<=path_indices_ptr->back();k++){
                                paths_ptr->at(path_indices_ptr->at(k)).push_front(indiv_1);
                            }
                        }
                    }else if((*indiv_1).get_name()==(*indiv_2).get_mom()){  // mom|offspring
                        xp = calculate_f_xy(indiv_1,(*indiv_2).sire_node,f_matrix,ncol,nrow,all_nodes_ptr,all_dyads_ptr,dyad_index,path_indices_ptr,commit_index,dyad_dict_ptr,node_space_ptr);
                        if(xp!=0){
                            for(int k = commit_index;k<=path_indices_ptr->back();k++){
                                paths_ptr->at(path_indices_ptr->at(k)).push_back(indiv_2);
                            }
                        }
                    }else if((*indiv_1).get_name()==(*indiv_2).get_sire()){ // sire|offspring
                        xp = calculate_f_xy(indiv_1,(*indiv_2).mom_node,f_matrix,ncol,nrow,all_nodes_ptr,all_dyads_ptr,dyad_index,path_indices_ptr,commit_index,dyad_dict_ptr,node_space_ptr);
                        if(xp!=0){ 
                            for(int k = commit_index;k<=path_indices_ptr->back();k++){
                                paths_ptr->at(path_indices_ptr->at(k)).push_back(indiv_2);
                            }
                        }
                    }
                    f_xy = 0.5*(1+xp);
                    for(int k = path_index;k<commit_index;k++){
                        paths_ptr->at(path_indices_ptr->at(k)).push_back(indiv_1);
                        paths_ptr->at(path_indices_ptr->at(k)).push_back(indiv_2);
                    }
                }else if(contains_ancestor(indiv_2,indiv_1,all_nodes_ptr)){ // indiv 1 is ancestor of indiv 2
                    if(path_index>path_indices_ptr->back()){ // if new path is needed -> add new path & new path_index
                        path_indices_ptr->push_back(all_dyads_ptr->at(dyad_index).add_new_path());
                    }
                    double xs = calculate_f_xy(indiv_1,(*indiv_2).sire_node,f_matrix,ncol,nrow,all_nodes_ptr,all_dyads_ptr,dyad_index,path_indices_ptr,next_path_index,dyad_dict_ptr,node_space_ptr);
                    if(xs!=0){
                        next_path_index = path_indices_ptr->back()+1;
                    }
                    double xm = calculate_f_xy(indiv_1,(*indiv_2).mom_node,f_matrix,ncol,nrow,all_nodes_ptr,all_dyads_ptr,dyad_index,path_indices_ptr,next_path_index,dyad_dict_ptr,node_space_ptr);
                    f_xy = 0.5*(xs + xm);
                    for(int k = path_index;k<=path_indices_ptr->back();k++){
                        paths_ptr->at(path_indices_ptr->at(k)).push_back(indiv_2);
                    }
                }else if(contains_ancestor(indiv_1,indiv_2,all_nodes_ptr)){ // indiv 2 is ancestor of indiv 1
                    if(path_index>path_indices_ptr->back()){
                        path_indices_ptr->push_back(all_dyads_ptr->at(dyad_index).add_new_path());
                    }
                    double sx = calculate_f_xy((*indiv_1).sire_node,indiv_2,f_matrix,ncol,nrow,all_nodes_ptr,all_dyads_ptr,dyad_index,path_indices_ptr,next_path_index,dyad_dict_ptr,node_space_ptr);
                    if(sx!=0){
                        next_path_index = path_indices_ptr->back()+1;
                    }   
                    double mx = calculate_f_xy((*indiv_1).mom_node,indiv_2,f_matrix,ncol,nrow,all_nodes_ptr,all_dyads_ptr,dyad_index,path_indices_ptr,next_path_index,dyad_dict_ptr,node_space_ptr);
                    f_xy = 0.5*(sx+mx);
                    for(int k = path_index;k<=path_indices_ptr->back();k++){
                        paths_ptr->at(path_indices_ptr->at(k)).push_front(indiv_1);
                    }
                }else{ // indiv 1 and indiv 2 are not (in)direct ancestors of each other, but might be related otherwise -> check
                    double mm, ms, sm, ss;
                    mm = calculate_f_xy((*indiv_1).mom_node,(*indiv_2).mom_node,f_matrix,ncol,nrow,all_nodes_ptr,all_dyads_ptr,dyad_index,path_indices_ptr,next_path_index,dyad_dict_ptr,node_space_ptr);
                    if(mm!=0){
                        next_path_index = path_indices_ptr->back()+1; // needs to be raised in case another path parts from here
                    }
                    ms = calculate_f_xy((*indiv_1).mom_node,(*indiv_2).sire_node,f_matrix,ncol,nrow,all_nodes_ptr,all_dyads_ptr,dyad_index,path_indices_ptr,next_path_index,dyad_dict_ptr,node_space_ptr);
                    if(ms!=0){
                        next_path_index = path_indices_ptr->back()+1;
                    }
                    sm = calculate_f_xy((*indiv_1).sire_node,(*indiv_2).mom_node,f_matrix,ncol,nrow,all_nodes_ptr,all_dyads_ptr,dyad_index,path_indices_ptr,next_path_index,dyad_dict_ptr,node_space_ptr);
                    if(sm!=0){
                        next_path_index = path_indices_ptr->back()+1;
                    }
                    ss = calculate_f_xy((*indiv_1).sire_node,(*indiv_2).sire_node,f_matrix,ncol,nrow,all_nodes_ptr,all_dyads_ptr,dyad_index,path_indices_ptr,next_path_index,dyad_dict_ptr,node_space_ptr);
                    f_xy=0.25*(mm+ms+sm+ss);
                    if(f_xy!=0){  
                        for(int k=path_index;k<=path_indices_ptr->back();k++){ // push indiv_1 and indiv_2 to all previously created paths (embracing) 
                            paths_ptr->at(k).push_front(indiv_1);
                            paths_ptr->at(k).push_back(indiv_2);
                        }
                    }
                }
            //}
            if(f_xy==0 && f_matrix->at(i)[j]==-1){ // fill f_matrix
                f_matrix->at(i)[j]=f_xy;
            }
        }
        return f_xy;
    }catch(const std::exception &ex) {
        std::cerr << "Error in calculate_f_xy(): " << ex.what() << std::endl;
    }
}
void fill_f_matrix(std::deque<std::deque<double>> *f_matrix,int ncol, int nrow,std::deque <node> *all_nodes_ptr,std::deque <dyad> *all_dyads_ptr,map<string,int>*dyad_dict_ptr,bool reduce_node_space_tf,bool multithreading, int dyads_start,int dyads_end,int thread){ 
    // function to start the dyadic relatedness calculation & add r values to relatedness matrix f
    try{
        if(multithreading == false){ // determines the range of dyads to process if multi-threading is deactivated
            dyads_start = 0;
            dyads_end = all_dyads_ptr->size();
        }else if(multithreading == true && dyads_end == 0){
            throw runtime_error("Unable to fill the 'f_matrix' with multiple threads (dyads_end == 0)");
        }
        for(int i = dyads_start;i<dyads_end;i++){ // Iterate over dyads of interest
            int progress = static_cast<int>((static_cast<double>(i-dyads_start)/(dyads_end-dyads_start))*100);
            if (progress % 10 == 0) {
                string printer = "["+ to_string(thread)+"] "+to_string(progress)+"% \t("+ to_string(i) +" of "+to_string(dyads_start)+".."+to_string(dyads_end)+")";
                cout <<printer<<endl;
            }
            auto [x,y] = all_dyads_ptr->at(i).get_dyad_idx_in_f_matrix();// get the dyadic matrix indices [x][y] in the relatedness matrix f
            if(all_dyads_ptr->at(i).get_indiv_1_name()!=all_nodes_ptr->at(0).get_name() // check that the dyad consists of non-imaginary nodes and that x, y are valid
                    &&all_dyads_ptr->at(i).get_indiv_2_name()!=all_nodes_ptr->at(0).get_name()
                    &&all_dyads_ptr->at(i).get_indiv_1_name()!=all_nodes_ptr->at(1).get_name()
                    &&all_dyads_ptr->at(i).get_indiv_2_name()!=all_nodes_ptr->at(1).get_name()
                    && x >= 0 && y >= 0){ 
                // add a new empty path (index) and initialize path_indices
                int path_index = all_dyads_ptr->at(i).add_new_path();
                std::deque <int> path_indices = {path_index}; 
                std::set<string> node_space;
                if (reduce_node_space_tf){ // if indicated, reduce the node space
                    reduce_node_space(&all_nodes_ptr->at(x),&all_nodes_ptr->at(y),&node_space);
                }
                if(reduce_node_space_tf && node_space.size()==0){ // if node space is empty, no common relatives exist
                    f_matrix->at(x)[y]=0;
                }else{ 
                    node_space.insert(all_dyads_ptr->at(i).get_indiv_1_name()); // include focal IDs in the node space
                    node_space.insert(all_dyads_ptr->at(i).get_indiv_2_name());
                    
                    // calculate dyadic relatedness and assign it to the relatedness matrix f
                    double f_xy = calculate_f_xy(&all_nodes_ptr->at(all_dyads_ptr->at(i).get_indiv_1_idx()),&all_nodes_ptr->at(all_dyads_ptr->at(i).get_indiv_2_idx()),f_matrix,ncol,nrow,all_nodes_ptr,all_dyads_ptr,i,&path_indices,path_index,dyad_dict_ptr,&node_space);
                    f_matrix->at(x)[y]=f_xy;
                }
            }
        }
    }catch(const std::exception &ex) {
        std::cerr << "Error in fill_f_matrix(): " << ex.what() << std::endl;
    }
}
void add_missing_parent_in_dyad_list(std::deque<node>*all_nodes_ptr,std::string all_node_names_str,std::deque<string>*mom_names_ptr,std::deque<string>*sire_names_ptr){
    try{
        int last_matidx = all_nodes_ptr->at(all_nodes_ptr->size()-1).get_matidx() + 1;
        for(int i = 0;i<mom_names_ptr->size();i++){
            if(all_node_names_str.find(("|"+mom_names_ptr->at(i))+"|") == string::npos){
                node n(mom_names_ptr->at(i),"f","unkn_f","unkn_m","NA","NA","NA","NA",0,last_matidx);
                all_nodes_ptr->push_back(n);
                all_node_names_str.append(n.get_name()+"|");
                last_matidx += 1;
            }
        }
        //last_matidx = all_nodes_ptr->at(all_nodes_ptr->size()-1).get_matidx() +1;
        for(int i = 0;i<sire_names_ptr->size();i++){
            if(all_node_names_str.find(("|"+sire_names_ptr->at(i))+"|") == string::npos){
                node n(sire_names_ptr->at(i),"m","unkn_f","unkn_m","NA","NA","NA","NA",0,last_matidx);
                all_nodes_ptr->push_back(n);
                all_node_names_str.append(n.get_name()+"|");
                last_matidx+=1;
            }
        }
    }catch(const std::exception &ex) {
        std::cerr << "Error in add_missing_parent_in_dyad_list(): " << ex.what() << std::endl;
    }
}
string all_nodes_info_file(std::deque <node> *all_nodes_ptr,std::map<string,int> * dyad_dict,std::deque<dyad>* all_dyads_ptr, string filename){ // save all_nodes as file
    try{
        ofstream out; 
        out.open(filename+".txt"); 
        for(int i = 0;i<all_nodes_ptr->size();i++){
            if(i==0){
                out << "ID\tsex\tBS\tmom\tsire\tDOB\tDOD\tpot_sire\tpot_mom\tfull_generations\tmin_f"<<endl;
            }
            if(all_nodes_ptr->at(i).get_name()!="unkn_f"&&all_nodes_ptr->at(i).get_name()!="unkn_m"){ // no imaginary nodes within population file
                out << all_nodes_ptr->at(i).get_name()<<"\t"<<all_nodes_ptr->at(i).get_sex()<<"\t"<<all_nodes_ptr->at(i).get_birthseason()<<"\t"<<all_nodes_ptr->at(i).get_mom()<<"\t"<<all_nodes_ptr->at(i).get_sire()<<"\t"<<all_nodes_ptr->at(i).get_DOB().get_date()<<"\t"<<all_nodes_ptr->at(i).get_DOD().get_date()<<"\t";
                out << all_nodes_ptr->at(i).get_pot_sires()<<"\t"<<all_nodes_ptr->at(i).get_pot_moms()<<"\t"<<all_nodes_ptr->at(i).get_full_generations()<<"\t"<<to_string_with_precision(all_nodes_ptr->at(i).get_min_f(),15)<<endl;
            }
        }
        out.close();
        cout<<"write "<<filename+".txt"<<endl;
        return filename;
    }catch(const std::exception &ex) {
        std::cerr << "Error in all_nodes_info_file(): " << ex.what() << std::endl;
    }
}
void write_dyad_list(std::deque <dyad> *all_dyads_ptr,std::deque<std::deque<double>> &f_mat,int ncol,string filename, string write_dyadlist){// writes assigned dyadlist into file
    try{
        ofstream out; 
        out.open(filename+"_dyad_list.txt"); 
        for(int i = 0;i<all_dyads_ptr->size();i++){
            auto [x,y] = all_dyads_ptr->at(i).get_dyad_idx_in_f_matrix();
            if(x >= 0 && y >= 0){
                all_dyads_ptr->at(i).set_r_value(f_mat[x][y]);//*(f_mat_ptr + (all_dyads_ptr->at(i).get_indiv_1_idx()*ncol+all_dyads_ptr->at(i).get_indiv_2_idx())));
                if(write_dyadlist == "full"){
                    out <<setprecision(15)<< all_dyads_ptr->at(i).get_dyad_infos();
                }else if (write_dyadlist == "reduced"){
                    out << setprecision(15)<<all_dyads_ptr->at(i).get_dyad_name() << "\t"<<all_dyads_ptr->at(i).get_r_value()<<endl;
                }
            }else{
                cout << "set_idx() error: dyad index == -1"<<endl;
            }
        }
        out.close();
        cout<<"write "<<filename+"_dyad_list.txt"<<endl;
    }catch(const std::exception &ex) {
        std::cerr << "Error in write_dyad_list(): " << ex.what() << std::endl;
    }
}
void set_min_f(std::deque<std::deque<double>> &f_matrix,std::deque<dyad>*all_dyads_ptr,std::deque<node>*all_nodes_ptr){
    try{
        for(int i = 0;i<all_nodes_ptr->size();i++){
            if(all_nodes_ptr->at(i).get_sire() == "unkn_m"||all_nodes_ptr->at(i).get_mom() == "unkn_f"){
                all_nodes_ptr->at(i).set_min_f(0);
            }else{
                int idx_1 = all_nodes_ptr->at(i).get_sire_node()->get_matidx();
                int idx_2 = all_nodes_ptr->at(i).get_mom_node()->get_matidx();
                if(idx_1 < idx_2){
                    swap(idx_1,idx_2);
                }
                if(f_matrix[idx_1][idx_2] == -1){
                    all_nodes_ptr->at(i).set_min_f(calculate_pure_f_xy(all_nodes_ptr->at(i).get_mom_node(),all_nodes_ptr->at(i).get_sire_node(),all_dyads_ptr,all_nodes_ptr,&f_matrix,false) * 0.5);
                }else{
                    all_nodes_ptr->at(i).set_min_f(f_matrix[idx_1][idx_2] * 0.5);
                }
            }
        }
    }catch(const std::exception &ex) {
        std::cerr << "Error in set_min_f(): " << ex.what() << std::endl;
    }
}
void set_all_min_DGD(std::deque<node>*all_nodes_ptr, std::deque<dyad>*all_dyads_ptr){
    try{
        int dgd_1 = -1;
        int dgd_2 = -1;
        for (int i = 0; i < all_dyads_ptr->size(); i++){
            dgd_1 = all_nodes_ptr->at(all_dyads_ptr->at(i).get_indiv_1_idx()).get_full_generations();
            dgd_2 = all_nodes_ptr->at(all_dyads_ptr->at(i).get_indiv_2_idx()).get_full_generations();
            if(dgd_1 != -1 && dgd_2 != -1 && all_dyads_ptr->at(i).get_min_dyadic_genealogical_depth() == -1){
                if(dgd_1 > dgd_2){
                    all_dyads_ptr->at(i).set_min_DGD(dgd_2);
                }else{
                    all_dyads_ptr->at(i).set_min_DGD(dgd_1);
                }
            }
        }
    }catch(const std::exception &ex) {
        std::cerr << "Error in set_all_min_DGD(): " << ex.what() << std::endl;
    }
}
void pip_forward(string file,string output_file,string input_dyadlist,int maturation_age_f,int maturation_age_m,int gestation_length,string write_dyadlist,int generation_limit,int n_cores,bool reduce_node_space_tf){ // pipeline or forward simulation -> calculates r and path characteristics for given dyads in a given, imperfect pedigree (pipeline for masterthesis part one)
    try{
        // INIT BASIC OBJECTS 
        string sex,name,mom,sire,DOB,DOD,nonsires,nondams;
        int birthseason;
        node imaginary_female(name="unkn_f","f","NA","NA","NA","NA","NA","NA",0,0); // create male & female imaginary nodes as dummy in case of an unknown parent
        node imaginary_male(name="unkn_m","m","NA","NA","NA","NA","NA","NA",0,1);
        std::deque <node> all_nodes;
        all_nodes.push_back(imaginary_female);
        all_nodes.push_back(imaginary_male);
        std::deque <string> all_node_names;
        std::string all_node_names_str = "|";
        std::deque <string> mom_node_names;
        std::deque <string> sire_node_names;
        all_node_names.push_back(imaginary_female.get_name());
        all_node_names.push_back(imaginary_male.get_name());

        // load pedigree data and transfer information into node objects
        ifstream data(file);
        if (! data) {
            cout << "unable to open "<<file<<" for reading" << endl;
        }else{
            cout << "load pedigree data ("<<file<<")"<<endl;
        }
        int matidx = 2;
        while(data >> name >> sex >> birthseason >> mom >> sire >> DOB >> DOD >> nonsires >> nondams){
            all_node_names.push_back(name);
            if(mom=="unknown"||mom=="NA"||mom=="unkn_f"||mom=="UNK"){
                mom="unkn_f";
            }else{
                mom_node_names.push_back(mom);
            }
            if(sire=="unknown"||sire=="NA"||sire=="unkn_m"||sire=="UNK"){
                sire="unkn_m";
            }else{
                sire_node_names.push_back(sire);
            }
            node x(name,sex,mom,sire,nonsires,nondams,DOB,DOD,birthseason,matidx);
            all_nodes.push_back(x);
            all_node_names.push_back(x.get_name());
            all_node_names_str.append(x.get_name()+"|");
            matidx += 1;
        }
        data.close(); 
        
        add_missing_parent_in_dyad_list(&all_nodes,all_node_names_str,&mom_node_names,&sire_node_names);
        
        // load dyad list with dyads of interest or if file is not available -> create all possible dyads including all_nodes
        ifstream data_dyadlist(input_dyadlist); 
        std::deque <dyad> all_dyads;
        map<string, int> dyad_dict;
        if (! data_dyadlist or input_dyadlist == "") {
            cout << "Unable to find dyad selection file for reading. Dyad list will be generated for all dyad combinations." << endl;
            int index = 0;
            for(int i = 0;i<all_nodes.size();i++){
                for(int j = i+1;j<all_nodes.size();j++){
                    if (all_nodes[i].get_name()!="unkn_f" && all_nodes[i].get_name()!="unkn_m" && all_nodes[j].get_name()!="unkn_f" && all_nodes[j].get_name()!="unkn_m"){
                        dyad x(all_nodes[i].get_name(),all_nodes[j].get_name());
                        all_dyads.push_back(x);
                        dyad_dict[x.get_dyad_name()]=index;
                        index += 1;
                    }
                }
            }
        }else{
            cout << "load dyad data ("<<input_dyadlist<<")"<<endl;
            string indiv_1,indiv_2;
            int index = 0;
            while(data_dyadlist>>indiv_1>>indiv_2){
                dyad x(indiv_1,indiv_2);
                if(dyad_dict.count(x.get_dyad_name())==0){
                    dyad_dict[x.get_dyad_name()]=index;
                    all_dyads.push_back(x);
                    index += 1;
                }
            }
            data_dyadlist.close();
        }
        // START MAIN PIPELINE
        cout << "create parent pointer"<<endl;
        for(int i = 0;i<all_nodes.size();i++){
            all_nodes[i].create_parent_ptr(&all_nodes);
        }
        cout << "create relatedness matrix"<<endl;
        int nrow=all_nodes.size();
        int ncol=all_nodes.size();
        std::deque<std::deque<double>> f_mat = {};
        for(int i = 0;i<ncol;i++){
            f_mat.push_back(std::deque<double>(i+1, -1.0f));
            f_mat[i][i] = 1;
            f_mat[i][0] = 0;
            if(i>0){
                f_mat[i][1] = 0;
            }
        }
        for(int i = 0;i<all_dyads.size();i++){
            all_dyads[i].set_idx(&all_nodes);
        }
        
        if(n_cores > 1){
            std::vector<std::thread> threads = {};
            for(int i = 0;i<n_cores;i++){ 
                int dyads_start = (int) (i*floor(all_dyads.size()/n_cores));
                int dyads_end = (int) ((i+1)*floor(all_dyads.size()/n_cores));
                if(i==(n_cores-1)){
                    dyads_end = all_dyads.size();
                }
                cout << "calculate dyadic relatedness: thread "<<i<<" from dyad "<<dyads_start<<" to "<<dyads_end<<endl;
                thread th1(fill_f_matrix,&f_mat,ncol,nrow,&all_nodes,&all_dyads,&dyad_dict,reduce_node_space_tf,true,dyads_start,dyads_end,i+1);
                threads.push_back(std::move(th1));
            }
            for(int i = 0;i<threads.size();i++){
                threads[i].join();
                //cout << "join thread "<<i<<endl;
            }
        }else{
            cout << "calculate dyadic relatedness"<<endl;
            fill_f_matrix(&f_mat,ncol,nrow,&all_nodes,&all_dyads,&dyad_dict,reduce_node_space_tf,false);
        }
        cout<<"set info attributes (parent pool, min_DGD, min_f)"<<endl;
        set_parent_pool(&all_nodes,maturation_age_m,maturation_age_f,gestation_length);
        set_all_min_DGD(&all_nodes,&all_dyads);
        set_min_f(f_mat,&all_dyads,&all_nodes);
        write_dyad_list(&all_dyads,f_mat,ncol,output_file,write_dyadlist);
        all_nodes_info_file(&all_nodes,&dyad_dict,&all_dyads,output_file+"_info");// potential parent information as new column in masterfile   
    }catch(const std::exception &ex) {
        std::cerr << "Error in pip_forward(): " << ex.what() << std::endl;
    }
}
void reduce_node_space(node* indiv_1,node* indiv_2,std::set<string>*node_space){
    try{
        std::deque<node*> ancestor_1;
        std::deque<node*> ancestor_2;
        std::deque<node*> common_ancs;
        std::map<string,std::deque<node*>> offspring_dict;
        get_ancestor_of_focal_indiv(&ancestor_1,indiv_1,&offspring_dict); // ptr-return: all_ancestors of indiv_1 (ancestor_1) && updated offspring_dict with all parent-offspring(s) "pairs"
        get_ancestor_of_focal_indiv(&ancestor_2,indiv_2,&offspring_dict);// ptr-return: all_ancestors of indiv_2 (ancestor_2) && updated offspring_dict with all parent-offspring(s) "pairs"
        
        for(int i = 0;i<ancestor_1.size();i++){ // only ancestors indiv1/2 have in common 
            for(int j = 0;j<ancestor_2.size();j++){
                if(ancestor_1.at(i)->get_name()==ancestor_2.at(j)->get_name()){
                    common_ancs.push_back(ancestor_1[i]);
                }
            }
        }
        std::deque<node*>* common_ancs_ptr = &common_ancs;
        for(int i = 0;i<common_ancs_ptr->size();i++){ // add all descendents of the common ancestors to the set (by dictionary) to generate the search space for paths
            node_space->insert(common_ancs_ptr->at(i)->get_name());
            if(offspring_dict.count(common_ancs_ptr->at(i)->get_name())>0){// || common_ancs_ptr->at(i)->get_name()!=indiv_1->get_name() && common_ancs_ptr->at(i)->get_name()!=indiv_2->get_name()){
                std::deque<node*> offspr = offspring_dict.at(common_ancs_ptr->at(i)->get_name());
                common_ancs_ptr->insert(common_ancs_ptr->end(),offspr.begin(),offspr.end()); // expand common_ancs so that the offspring of the added offspring will be added too
            }
        }
    }catch(const std::exception &ex) {
        std::cerr << "Error in reduce_node_space(): " << ex.what() << std::endl;
    }
}
bool indivs_in_node_space(node* indiv_1, node* indiv_2,std::set<string>*node_space_ptr){
    try{
        bool i1_in_node_space = false;
        bool i2_in_node_space = false;
        for (auto& node_name : *node_space_ptr) {
            if(indiv_1->get_name() == node_name){
                i1_in_node_space = true;																	  
            }						   
            if(indiv_2->get_name() == node_name){
                i2_in_node_space = true;																 
            }
        }
        return (i1_in_node_space & i2_in_node_space);
    }catch(const std::exception &ex) {
        std::cerr << "Error in indivs_in_node_space(): " << ex.what() << std::endl;
    }
}