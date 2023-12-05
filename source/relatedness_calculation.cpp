#include "relatedness_calculation.h"

double calculate_f_xy(node * indiv_1,node * indiv_2,std::deque<std::deque<double>> *f_matrix,int ncol,int nrow,std::deque <node> * all_nodes_ptr,std::deque <dyad> * all_dyads_ptr,int dyad_index,std::deque <int> * path_indices_ptr,int path_index,map <string,int> * dyad_dict_ptr, std::set<string>* node_space_ptr,bool reduce_node_space){ // calculate dyadic relatedness coeffiecients recursively and collect path characteristics  
    try{
        int i = (*indiv_1).get_matidx(); // get matrix indices of the dyad
        int j = (*indiv_2).get_matidx();
        if(i<j){
            swap(i,j);
        }
        double f_xy = f_matrix->at(i)[j];// get (default) relatedness value from f_matrix
        int next_path_index = path_index;
        std::deque <deque<node*>>* paths_ptr = all_dyads_ptr->at(dyad_index).get_paths_ptr();
        if(f_xy==1){ // indiv_1 == indiv_2
            if(path_index>path_indices_ptr->back()){ // if new path is needed -> add new path & new path_index
                path_indices_ptr->push_back(all_dyads_ptr->at(dyad_index).add_new_path());
            }
            paths_ptr->at(path_index).push_back(indiv_1); // indiv_1/2 correlates with common ancestor (because it is the same individual) -> push
        }else if(f_xy<1&f_xy>0){ // for instance if dyad of interest was already computed for prior dyad -> info & paths can be found in dyadlist
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
        }else if(f_xy!=0){ // no entry in f_matrix ([i][j] == 999): calculate the dyadic relatedness coefficient
            if(reduce_node_space == true && indivs_in_node_space(indiv_1,indiv_2,node_space_ptr)==false){ // if node space is reduced, check first if a calculation is necessary
                f_xy = 0;
            }else{
                if((*indiv_1).get_mom()==(*indiv_2).get_name()||(*indiv_1).get_sire()==(*indiv_2).get_name()||(*indiv_2).get_mom()==(*indiv_1).get_name()||(*indiv_2).get_sire()==(*indiv_1).get_name()){ // if parent/offspring 
                    if(path_index>path_indices_ptr->back()){// if new path is needed -> add new path & new path_index
                        path_indices_ptr->push_back(all_dyads_ptr->at(dyad_index).add_new_path());
                    }
                    int commit_index = path_indices_ptr->back()+1;
                    double xp = 0;
                    if((*indiv_1).get_mom()==(*indiv_2).get_name()){  // offspring|mom 
                        xp = calculate_f_xy((*indiv_1).sire_node,indiv_2,f_matrix,ncol,nrow,all_nodes_ptr,all_dyads_ptr,dyad_index,path_indices_ptr,commit_index,dyad_dict_ptr, node_space_ptr,reduce_node_space);
                        if(xp!=0){ 
                            for(int k = commit_index;k<=path_indices_ptr->back();k++){
                                paths_ptr->at(path_indices_ptr->at(k)).push_front(indiv_1);
                            }
                        }
                    }else if((*indiv_1).get_sire()==(*indiv_2).get_name()){  // offspring|sire
                        xp = calculate_f_xy((*indiv_1).mom_node,indiv_2,f_matrix,ncol,nrow,all_nodes_ptr,all_dyads_ptr,dyad_index,path_indices_ptr,commit_index,dyad_dict_ptr,node_space_ptr,reduce_node_space);
                        if(xp!=0){ 
                            for(int k = commit_index;k<=path_indices_ptr->back();k++){
                                paths_ptr->at(path_indices_ptr->at(k)).push_front(indiv_1);
                            }
                        }
                    }else if((*indiv_1).get_name()==(*indiv_2).get_mom()){  // mom|offspring
                        xp = calculate_f_xy(indiv_1,(*indiv_2).sire_node,f_matrix,ncol,nrow,all_nodes_ptr,all_dyads_ptr,dyad_index,path_indices_ptr,commit_index,dyad_dict_ptr,node_space_ptr,reduce_node_space);
                        if(xp!=0){
                            for(int k = commit_index;k<=path_indices_ptr->back();k++){
                                paths_ptr->at(path_indices_ptr->at(k)).push_back(indiv_2);
                            }
                        }
                    }else if((*indiv_1).get_name()==(*indiv_2).get_sire()){ // sire|offspring
                        xp = calculate_f_xy(indiv_1,(*indiv_2).mom_node,f_matrix,ncol,nrow,all_nodes_ptr,all_dyads_ptr,dyad_index,path_indices_ptr,commit_index,dyad_dict_ptr,node_space_ptr,reduce_node_space);
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
                    double xs = calculate_f_xy(indiv_1,(*indiv_2).sire_node,f_matrix,ncol,nrow,all_nodes_ptr,all_dyads_ptr,dyad_index,path_indices_ptr,next_path_index,dyad_dict_ptr,node_space_ptr,reduce_node_space);
                    if(xs!=0){
                        next_path_index = path_indices_ptr->back()+1;
                    }
                    double xm = calculate_f_xy(indiv_1,(*indiv_2).mom_node,f_matrix,ncol,nrow,all_nodes_ptr,all_dyads_ptr,dyad_index,path_indices_ptr,next_path_index,dyad_dict_ptr,node_space_ptr,reduce_node_space);
                    f_xy = 0.5*(xs + xm); 
                    for(int k = path_index;k<=path_indices_ptr->back();k++){
                        paths_ptr->at(path_indices_ptr->at(k)).push_back(indiv_2);
                    }
                }else if(contains_ancestor(indiv_1,indiv_2,all_nodes_ptr)){ // indiv 2 is ancestor of indiv 1
                    if(path_index>path_indices_ptr->back()){
                        path_indices_ptr->push_back(all_dyads_ptr->at(dyad_index).add_new_path());
                    }
                    double sx = calculate_f_xy((*indiv_1).sire_node,indiv_2,f_matrix,ncol,nrow,all_nodes_ptr,all_dyads_ptr,dyad_index,path_indices_ptr,next_path_index,dyad_dict_ptr,node_space_ptr,reduce_node_space);
                    if(sx!=0){
                        next_path_index = path_indices_ptr->back()+1;
                    }   
                    double mx = calculate_f_xy((*indiv_1).mom_node,indiv_2,f_matrix,ncol,nrow,all_nodes_ptr,all_dyads_ptr,dyad_index,path_indices_ptr,next_path_index,dyad_dict_ptr,node_space_ptr,reduce_node_space);
                    f_xy = 0.5*(sx+mx);
                    for(int k = path_index;k<=path_indices_ptr->back();k++){
                        paths_ptr->at(path_indices_ptr->at(k)).push_front(indiv_1);
                    }
                }else{ // indiv 1 and indiv 2 are not (in)direct ancestors of each other, but might be related otherwise -> check
                    double mm, ms, sm, ss;
                    mm = calculate_f_xy((*indiv_1).mom_node,(*indiv_2).mom_node,f_matrix,ncol,nrow,all_nodes_ptr,all_dyads_ptr,dyad_index,path_indices_ptr,next_path_index,dyad_dict_ptr,node_space_ptr,reduce_node_space);
                    if(mm!=0){
                        next_path_index = path_indices_ptr->back()+1; // needs to be raised in case another path parts from here
                    }
                    ms = calculate_f_xy((*indiv_1).mom_node,(*indiv_2).sire_node,f_matrix,ncol,nrow,all_nodes_ptr,all_dyads_ptr,dyad_index,path_indices_ptr,next_path_index,dyad_dict_ptr,node_space_ptr,reduce_node_space);
                    if(ms!=0){
                        next_path_index = path_indices_ptr->back()+1;
                    }
                    sm = calculate_f_xy((*indiv_1).sire_node,(*indiv_2).mom_node,f_matrix,ncol,nrow,all_nodes_ptr,all_dyads_ptr,dyad_index,path_indices_ptr,next_path_index,dyad_dict_ptr,node_space_ptr,reduce_node_space);
                    if(sm!=0){
                        next_path_index = path_indices_ptr->back()+1;
                    }
                    ss = calculate_f_xy((*indiv_1).sire_node,(*indiv_2).sire_node,f_matrix,ncol,nrow,all_nodes_ptr,all_dyads_ptr,dyad_index,path_indices_ptr,next_path_index,dyad_dict_ptr,node_space_ptr,reduce_node_space);
                    f_xy=0.25*(mm+ms+sm+ss);
                    if(f_xy!=0){  
                        for(int k=path_index;k<=path_indices_ptr->back();k++){ // push indiv_1 and indiv_2 to all previously created paths (embracing) 
                            paths_ptr->at(k).push_front(indiv_1);
                            paths_ptr->at(k).push_back(indiv_2);
                        }
                    }
                }
            }
            if(f_xy==0 && f_matrix->at(i)[j]==-1){ // fill f_matrix
                f_matrix->at(i)[j]=f_xy;
            }
        }
        return f_xy;
    }catch(const std::exception &ex) {
        std::cerr << "Error in calculate_f_xy(): " << ex.what() << std::endl;
        return 999;
    }
}
void fill_f_matrix(std::deque<std::deque<double>> *f_matrix,int ncol, int nrow,std::deque <node> *all_nodes_ptr,std::deque <dyad> *all_dyads_ptr,map<string,int>*dyad_dict_ptr,bool reduce_node_space_tf,bool multithreading, int dyads_start,int dyads_end,int thread){ // function to start the dyadic relatedness calculation & add r values to relatedness matrix f (with path characteristics)
    try{
        if(multithreading == false){ // determines the range of dyads to process if multi-threading is not requested
            dyads_start = 0;
            dyads_end = all_dyads_ptr->size();
        }else if(multithreading == true && dyads_end == 0){
            throw runtime_error("Unable to fill the 'f_matrix' with multiple threads (dyads_end == 0)");
        }
        int last_printed_progress = 0;
        for(int i = dyads_start;i<dyads_end;i++){ // Iterate over dyads of interest
            int progress = static_cast<int>((static_cast<double>(i-dyads_start)/(dyads_end-dyads_start))*100);
            if (progress % 10 == 0 && progress >  last_printed_progress) {
                string printer = "["+ to_string(thread)+"] "+to_string(progress)+"% \t("+ to_string(i) +" of "+to_string(dyads_start)+".."+to_string(dyads_end)+")";
                cout <<printer<<endl;
                last_printed_progress = progress;
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
                    double f_xy = calculate_f_xy(&all_nodes_ptr->at(all_dyads_ptr->at(i).get_indiv_1_idx()),&all_nodes_ptr->at(all_dyads_ptr->at(i).get_indiv_2_idx()),f_matrix,ncol,nrow,all_nodes_ptr,all_dyads_ptr,i,&path_indices,path_index,dyad_dict_ptr,&node_space,reduce_node_space_tf); // calculate dyadic relatedness and assign it to the relatedness matrix f
                    f_matrix->at(x)[y]=f_xy;
                }
            }
        }
    }catch(const std::exception &ex) {
        std::cerr << "Error in fill_f_matrix(): " << ex.what() << std::endl;
    }
}
void add_missing_parent_nodes(std::deque<node>*all_nodes_ptr,std::string all_node_names_str,std::deque<string>*mom_names_ptr,std::deque<string>*sire_names_ptr){ // check if all parents are also included in all_nodes, else they will be added
    try{
        int last_matidx = all_nodes_ptr->at(all_nodes_ptr->size()-1).get_matidx() + 1;// get index of last individual in all nodes and add 1 (next potential index)
        for(int i = 0;i<mom_names_ptr->size();i++){ // iterate through all moms
            if(all_node_names_str.find(("|"+mom_names_ptr->at(i))+"|") == string::npos){ // if all_nodes does not include mom
                node n(mom_names_ptr->at(i),"f","unkn_f","unkn_m","NA","NA","NA","NA",0,last_matidx); // create mom node and add it to all_nodes
                all_nodes_ptr->push_back(n);
                all_node_names_str.append(n.get_name()+"|");
                last_matidx += 1;
            }
        }
        for(int i = 0;i<sire_names_ptr->size();i++){ // iterate through all sires
            if(all_node_names_str.find(("|"+sire_names_ptr->at(i))+"|") == string::npos){ //if all_nodes does not included sire
                node n(sire_names_ptr->at(i),"m","unkn_f","unkn_m","NA","NA","NA","NA",0,last_matidx);//create sire node and add it to all_nodes
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
        for(int i = 0;i<all_nodes_ptr->size();i++){ // iterate through all nodes
            if(i==0){
                out << "ID\tfull_generations\tmin_f"<<endl; // header
            }
            if(all_nodes_ptr->at(i).get_name()!="unkn_f"&&all_nodes_ptr->at(i).get_name()!="unkn_m"){ // no imaginary nodes within pedigree file
                out << all_nodes_ptr->at(i).get_name()<<"\t"<<all_nodes_ptr->at(i).get_full_generations()<<"\t"<<to_string_with_precision(all_nodes_ptr->at(i).get_min_f(),15)<<endl; // make file with additional (min_f & generation_depth) per individual
            }
        }
        out.close();
        cout<<"write "<<filename+".txt"<<endl;
        return filename;
    }catch(const std::exception &ex) {
        std::cerr << "Error in all_nodes_info_file(): " << ex.what() << std::endl;
        return "unable to write nodes info file";
    }
}
void write_dyad_list(std::deque <dyad> *all_dyads_ptr,std::deque<std::deque<double>> &f_mat,int ncol,string filename, string write_dyadlist,int generation_limit){// writes given dyads and their information into output file
    try{
        ofstream out; 
        out.open(filename+"_dyad_list.txt"); 
        for(int i = 0;i<all_dyads_ptr->size();i++){ // iterate through all dyads
            auto [x,y] = all_dyads_ptr->at(i).get_dyad_idx_in_f_matrix(); // get relatedness matrix indices to get relatedness value
            if(x >= 0 && y >= 0){
                all_dyads_ptr->at(i).set_r_value(f_mat[x][y]); // set r value as attribute
                if(write_dyadlist == "full"){ // write all information into file
                    out <<setprecision(15)<< all_dyads_ptr->at(i).get_dyad_infos(generation_limit);
                }else if (write_dyadlist == "reduced"){ // write only r value and dyadname into output file
                    out << setprecision(15)<<all_dyads_ptr->at(i).get_dyad_name() << "\t"<<all_dyads_ptr->at(i).get_r_value()<<endl;
                }
            }else{
                throw runtime_error("set_idx() error: dyad index == -1");
            }
        }
        out.close();
        cout<<"write "<<filename+"_dyad_list.txt"<<endl;
    }catch(const std::exception &ex) {
        std::cerr << "Error in write_dyad_list(): " << ex.what() << std::endl;
    }
}
void set_min_f(std::deque<std::deque<double>> &f_matrix,std::deque<dyad>*all_dyads_ptr,std::deque<node>*all_nodes_ptr){ // calculate and set minimal inbreeding value for each node/individual
    try{
        for(int i = 0;i<all_nodes_ptr->size();i++){ // iterate through all nodes
            if(all_nodes_ptr->at(i).get_sire() == "unkn_m"||all_nodes_ptr->at(i).get_mom() == "unkn_f"){ // if at least one parent is unknown, parents are considered as unrelated -> min_f = 0
                all_nodes_ptr->at(i).set_min_f(0);
            }else{
                int idx_1 = all_nodes_ptr->at(i).get_sire_node()->get_matidx(); // get index of sire
                int idx_2 = all_nodes_ptr->at(i).get_mom_node()->get_matidx(); // get index of mom
                if(idx_1 < idx_2){
                    swap(idx_1,idx_2);
                }
                if(f_matrix[idx_1][idx_2] == -1){ // mom_sire dyad relatedness coefficient is not yet calculated -> calculate & set
                    all_nodes_ptr->at(i).set_min_f(calculate_pure_f_xy(all_nodes_ptr->at(i).get_mom_node(),all_nodes_ptr->at(i).get_sire_node(),all_dyads_ptr,all_nodes_ptr,&f_matrix,false) * 0.5);
                }else{ // get dyadic relatedness coefficient of parents to calculate minimal inbreeding value
                    all_nodes_ptr->at(i).set_min_f(f_matrix[idx_1][idx_2] * 0.5);
                }
            }
        }
    }catch(const std::exception &ex) {
        std::cerr << "Error in set_min_f(): " << ex.what() << std::endl;
    }
}
void set_all_min_DGD(std::deque<node>*all_nodes_ptr, std::deque<dyad>*all_dyads_ptr){ // set minimal genetic depth (how many generations are completely known) for each dyad (get it from node attribute)
    try{
        int dgd_1 = -1; // default
        int dgd_2 = -1;
        for (int i = 0; i < all_dyads_ptr->size(); i++){ // iterate through each dyad and set min_DGD
            dgd_1 = all_nodes_ptr->at(all_dyads_ptr->at(i).get_indiv_1_idx()).get_full_generations();
            dgd_2 = all_nodes_ptr->at(all_dyads_ptr->at(i).get_indiv_2_idx()).get_full_generations();
            if(dgd_1 != -1 && dgd_2 != -1 && all_dyads_ptr->at(i).get_min_dyadic_genealogical_depth() == -1){
                if(dgd_1 > dgd_2){ // determine which one is the lower one
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
void relatedness_calculation(string file,string output_file,string input_dyadlist,int maturation_age_f,int maturation_age_m,int gestation_length,bool twins,string write_dyadlist,int generation_limit,int n_cores,bool reduce_node_space_tf){ // calculates r and path characteristics for given dyads in a given, imperfect pedigree (pipeline for masterthesis part one)
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
            throw runtime_error("Unable to open "+file+" for reading");
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
            node x(name,sex,mom,sire,nonsires,nondams,DOB,DOD,birthseason,matidx); // create node object for each row/individual in the pedigree file and add it to all nodes
            all_nodes.push_back(x);
            all_node_names.push_back(x.get_name());
            all_node_names_str.append(x.get_name()+"|");
            matidx += 1;
        }
        data.close(); 
        
        add_missing_parent_nodes(&all_nodes,all_node_names_str,&mom_node_names,&sire_node_names); // ensure that all parents are included in all_nodes too
        
        // load dyad list with dyads of interest or if file is not available -> create all possible dyads including all_nodes
        ifstream data_dyadlist(input_dyadlist); 
        std::deque <dyad> all_dyads;
        map<string, int> dyad_dict;
        if (! data_dyadlist or input_dyadlist == "NA") {
            cout << "No file of selected dyads. Dyad list will be generated for all dyad combinations." << endl;
            int index = 0;
            for(int i = 0;i<all_nodes.size();i++){ //iterate through all nodes systematically to find all combinations
                for(int j = i+1;j<all_nodes.size();j++){ 
                    if (all_nodes[i].get_name()!="unkn_f" && all_nodes[i].get_name()!="unkn_m" && all_nodes[j].get_name()!="unkn_f" && all_nodes[j].get_name()!="unkn_m"){ // no dyads created with imaginary nodes
                        dyad x(all_nodes[i].get_name(),all_nodes[j].get_name()); // create dyad object and add to all_dyads
                        all_dyads.push_back(x);
                        dyad_dict[x.get_dyad_name()]=index;
                        index += 1;
                    }
                }
            }
        }else{ // dyad selection exists
            cout << "load dyad data ("<<input_dyadlist<<")"<<endl;
            string indiv_1,indiv_2;
            int index = 0;
            while(data_dyadlist>>indiv_1>>indiv_2){
                dyad x(indiv_1,indiv_2); // create dyad object and add to all_dyads
                if(dyad_dict.count(x.get_dyad_name())==0){
                    dyad_dict[x.get_dyad_name()]=index;
                    all_dyads.push_back(x);
                    index += 1;
                }
            }
            data_dyadlist.close();
        }
        
        cout << "create parent pointer"<<endl;
        for(int i = 0;i<all_nodes.size();i++){
            all_nodes[i].create_parent_ptr(&all_nodes);
        }
        cout << "initialize relatedness matrix"<<endl;
        int nrow=all_nodes.size();
        int ncol=all_nodes.size();
        std::deque<std::deque<double>> f_mat = {};
        for(int i = 0;i<ncol;i++){
            f_mat.push_back(std::deque<double>(i+1, -1.0f)); // relatedness matrix (n*n == nodes.size()) with -1 as default value
            f_mat[i][i] = 1; // relatedness of an individual to itself
            f_mat[i][0] = 0; // relatedness to imaginary female
            if(i>0){
                f_mat[i][1] = 0; // relatedness to imaginary male
            }
        }
        for(int i = 0;i<all_dyads.size();i++){ // link index from all_nodes & relatedness matrix to dyad
            all_dyads[i].set_idx(&all_nodes);
        }
        if(n_cores > 1){ // prepare multiprocessing
            std::vector<std::thread> threads = {};
            for(int i = 0;i<n_cores;i++){ 
                int dyads_start = (int) (i*floor(all_dyads.size()/n_cores));// determine an almost uniform number of dyads to process for each core
                int dyads_end = (int) ((i+1)*floor(all_dyads.size()/n_cores));
                if(i==(n_cores-1)){// last core
                    dyads_end = all_dyads.size();
                }
                cout << "calculate dyadic relatedness: thread "+to_string(i)+" from dyad "+to_string(dyads_start)+" to "+to_string(dyads_end)<<endl;
                thread th1(fill_f_matrix,&f_mat,ncol,nrow,&all_nodes,&all_dyads,&dyad_dict,reduce_node_space_tf,true,dyads_start,dyads_end,i+1); // calculate dyadic relatedness coefficients and path characteristics
                threads.push_back(std::move(th1));
            }
            for(int i = 0;i<threads.size();i++){ // wait for all threads to finish
                threads[i].join();
            }
        }else{ // no multiprocessing
            cout << "calculate dyadic relatedness"<<endl;
            fill_f_matrix(&f_mat,ncol,nrow,&all_nodes,&all_dyads,&dyad_dict,reduce_node_space_tf,false);// calculate dyadic relatedness coefficients and path characteristics without multiprocessing
        }
        cout<<"set info attributes (parent pool, min_DGD, min_f)"<<endl;
        set_parent_pool(&all_nodes,maturation_age_m,maturation_age_f,gestation_length,twins); //parent pool as node attribute
        set_all_min_DGD(&all_nodes,&all_dyads); //minimal completely known generations as dyad attribute
        set_min_f(f_mat,&all_dyads,&all_nodes); //minimal inbreeding value as node attribute
        write_dyad_list(&all_dyads,f_mat,ncol,output_file,write_dyadlist,generation_limit); // save relatedness & path characteristics in output file
        all_nodes_info_file(&all_nodes,&dyad_dict,&all_dyads,output_file+"_info");// save pedigree with extra information in additional output file 
    }catch(const std::exception &ex) {
        std::cerr << "Error in relatedness_calculation(): " << ex.what() << std::endl;
    }
}
void reduce_node_space(node* indiv_1,node* indiv_2,std::set<string>*node_space){ // return set of node names that includes only common ancestor of two focal and their descandants to ultimately reduce the number of individuals to be part of the recursive relatedness calculation/path search
    try{
        std::deque<node*> ancestor_1; // lists all ancestors of individual 1
        std::deque<node*> ancestor_2; // lists all ancestors of individual 2
        std::deque<node*> common_ancs; // filtered common ancestors of both individuals
        std::map<string,std::deque<node*>> offspring_dict;
        get_ancestor_of_focal_indiv(&ancestor_1,indiv_1,&offspring_dict); // ptr-return: all_ancestors of indiv_1 (ancestor_1) && updated offspring_dict with all parent-offspring(s) "pairs"
        get_ancestor_of_focal_indiv(&ancestor_2,indiv_2,&offspring_dict);// ptr-return: all_ancestors of indiv_2 (ancestor_2) && updated offspring_dict with all parent-offspring(s) "pairs"
        
        for(int i = 0;i<ancestor_1.size();i++){ // iterate through ancestors of individual 1 and filter for common ancestors
            for(int j = 0;j<ancestor_2.size();j++){
                if(ancestor_1.at(i)->get_name()==ancestor_2.at(j)->get_name()){
                    common_ancs.push_back(ancestor_1[i]); // add as common ancestor
                }
            }
        }
        std::deque<node*>* common_ancs_ptr = &common_ancs; // convert to pointer to add elements during iteration
        for(int i = 0;i<common_ancs_ptr->size();i++){ // iterate through common ancestors to add all of their descendants  == node space
            node_space->insert(common_ancs_ptr->at(i)->get_name()); // add common ancestor to node space
            if(offspring_dict.count(common_ancs_ptr->at(i)->get_name())>0){ // other offspring from common ancestor are irrelevant because they are already excluded to be an ancestor of one of the individuals 
                std::deque<node*> offspr = offspring_dict.at(common_ancs_ptr->at(i)->get_name());//get all offspring rom offspring dict for the respective common ancestor
                common_ancs_ptr->insert(common_ancs_ptr->end(),offspr.begin(),offspr.end()); // expand common_ancs so that the offspring of the added offspring will be added too
            }
        }
    }catch(const std::exception &ex) {
        std::cerr << "Error in reduce_node_space(): " << ex.what() << std::endl;
    }
}
bool indivs_in_node_space(node* indiv_1, node* indiv_2,std::set<string>*node_space_ptr){ // check if both individuals are included in the node space
    try{
        bool i1_in_node_space = false; // default
        bool i2_in_node_space = false;
        for (auto& node_name : *node_space_ptr) { // iterate through node space
            if(indiv_1->get_name() == node_name){ // if indiv 1 is included in node space, first true
                i1_in_node_space = true;																	  
            }						   
            if(indiv_2->get_name() == node_name){ // if indiv 2 is also included in node space, second true
                i2_in_node_space = true;																 
            }
        }
        return (i1_in_node_space & i2_in_node_space); // true, only if boths individuals are in node space
    }catch(const std::exception &ex) {
        std::cerr << "Error in indivs_in_node_space(): " << ex.what() << std::endl;
        return false;
    }
}

// generation limitation implementieren
// why pot_mom/sire empty in example?