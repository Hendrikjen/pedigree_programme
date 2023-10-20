#include "simulated_annealing.h"

double randomnumber() {
    static std::default_random_engine randomnumber;
    std::uniform_real_distribution<double> dist(0.0, 1.0); 
    return dist(randomnumber); 
}
double load_data_for_sim_annealing(string file_full,string file_gaps, std::deque<dyad>*all_dyads_ptr,std::deque<int>*subset_idx_ptr, std::deque<node>*all_nodes_ptr,map<string, int> *dyad_dict_ptr,int maturation_age_f,int maturation_age_m,int gestation_length,string file_dyads){ // Loads and prepares all relevant data for simulated annealing (full pedigree, gap pedigree and their respectively dyadic information + calculates basic statistic (ratio, mean, max and min of gap-resulting discrepancy)
    // LOAD PEDIGREE INFORMATION
    cout << "load data for simulated annealing ..."<<endl;
    string sex,name,mom,sire,DOB,DOD,nonsires,nondams,indiv1, indiv2,dyad_name,paths,pathline,kinline,common_anc,depth,kinlabel,fullhalf;
    int birthseason;
    double r_value,real_r;
    node imaginary_female(name="unkn_f","f","NA","NA","NA","NA","NA","NA",0,0); // create male & female imaginary nodes as dummy in case of an unknown parent
    node imaginary_male(name="unkn_m","m","NA","NA","NA","NA","NA","NA",0,1);
    all_nodes_ptr->push_back(imaginary_female);
    all_nodes_ptr->push_back(imaginary_male);
    ifstream data(file_gaps+".txt"); // load pedigree with gaps
    if (! data) {
        cout << "unable to open file '"<<file_gaps <<".txt' for reading" << endl;
    }
    int matidx = 2;
    while(data >> name >> sex >> birthseason >> mom >> sire >> DOB >> DOD >> nonsires >> nondams){
        if(mom=="unknown"||mom=="unkn_f"||mom=="UNK"||mom=="NA"){
            mom=="unkn_f";
        }
        if(sire=="unknown"||sire=="unkn_m"||sire=="UNK"||sire=="NA"){
            sire=="unkn_m";
        }
        node x(name,sex,mom,sire,nonsires,nondams,DOB,DOD,birthseason,matidx);
        all_nodes_ptr->push_back(x);
        matidx += 1;
    }
    data.close(); 
    cout<<" ... gaps pedigree loaded ..."<<endl;
    if(fs::exists(file_full+".txt")){
        ifstream data_full(file_full+".txt"); // load completely known pedigree
        if (! data_full) {
            cout << "file '"<<file_full <<".txt' exists but program is unable to open file for reading" << endl;
        }
        matidx = 2;
        while(data_full >> name >> sex >> birthseason >> mom >> sire >> DOB >> DOD >> nonsires >> nondams){
            if(mom=="unknown"||mom=="unkn_f"||mom=="UNK"||mom=="NA"){
                mom=="unkn_f";
            }
            if(sire=="unknown"||sire=="unkn_m"||sire=="UNK"||sire=="NA"){
                sire=="unkn_m";
            }
            all_nodes_ptr->at(matidx).set_real_mom(mom);
            all_nodes_ptr->at(matidx).set_real_sire(sire);
            matidx += 1;
        }
        data_full.close(); 
    }else{
        cout << "file '"<<file_full <<".txt' (real & completely filled pedigree) does not exist. no pedigree comparison at the end possible."<<endl;
    }
    for(int i = 0;i<all_nodes_ptr->size();i++){
        all_nodes_ptr->at(i).create_parent_ptr(all_nodes_ptr);
    }
    cout<<" ... complete pedigree (add only real parent) loaded ..."<<endl;
    // LOAD DYAD INFOMATION
    double normalize_factor = 0;
    std::deque<int> dyads_with_erroneous_r_values = {};
    if (file_dyads == "NA"){
        int dl_count = 0;
        cout << "try to load both dyad_lists and concatenate..."<<endl;
        if(fs::exists(file_gaps+"_dyad_list.txt") == false){
            cout << "dyad_list file does not exist. calculate."<<endl;
            pip_forward(file_gaps,maturation_age_f,maturation_age_m,gestation_length);
        }
        ifstream data_dyads(file_gaps+ "_dyad_list.txt");// load dyad information which belongs to pedigree with gaps -> gap r values // please ensure a complete dyad_list with all dyads of interest (from WGS)
        if (! data_dyads) {
            cout << "unable to open file '"<<file_gaps <<"_dyad_list.txt' for reading" << endl;
        }else{
            while(data_dyads>>indiv1>>indiv2>>dyad_name>>r_value>>paths>>pathline>>kinline>>common_anc>>depth>>kinlabel>>fullhalf){
                dyad x(indiv1,indiv2);
                if(r_value>=0 && r_value<=1){
                    x.set_gap_r_value(r_value);
                    x.set_path_str(paths);
                }else{
                    cout << "Erroneous r_value ("<<r_value<<") for dyad "<<x.get_dyad_name()<<endl;
                    dyads_with_erroneous_r_values.push_back(dl_count);
                }
                all_dyads_ptr->push_back(x);
                dl_count += 1;
            }
        }
        data_dyads.close(); 
        cout<<" ... gaps dyad list loaded..."<<endl;

        ifstream data_full_dyads(file_full+ "_dyad_list.txt"); // load dyad information which belongs to completely known pedigree -> real r values
        if (! data_full_dyads) {
            cout << "unable to open file '"<<file_full <<"_dyad_list.txt' for reading" << endl;
        }
        dl_count = 0;
        while(data_full_dyads>>dyad_name>>r_value){
            if(dl_count%1000000==0){
                cout << dl_count<<endl;
            }
            if(r_value>=0 && r_value<=1){
                bool found = false;
                for(int i = 0;i<all_dyads_ptr->size();i++){
                    if(all_dyads_ptr->at(i).get_dyad_name()==dyad_name){
                        all_dyads_ptr->at(i).set_real_r_value(r_value);
                        dyad_dict_ptr->insert({dyad_name,matidx});
                        normalize_factor += r_value;
                        subset_idx_ptr->push_back(i);
                        found = true;
                        break;
                    }
                }
                if (found == false){
                    std::deque<string> indivs = str_split(dyad_name,"_");
                    dyad x(indivs[0],indivs[1]);
                    x.set_real_r_value(r_value);
                    dyad_dict_ptr->insert({dyad_name,matidx});
                    normalize_factor += r_value;
                    all_dyads_ptr->push_back(x);
                    subset_idx_ptr->push_back(all_dyads_ptr->size());
                    dyads_with_erroneous_r_values.push_back(all_dyads_ptr->size());
                }
            }else{
                cout << "check: Erroneous r_value ("<<r_value<<") for dyad "<<dyad_name<<" ... is currently ignored"<<endl;
            }
            dl_count += 1;
        }
        data_full_dyads.close(); 
        cout<<" ... complete pedigree dyad list loaded ..."<<endl;
    }else{
        int dl_count = 0;
        ifstream data_dyads(file_dyads+".txt");// load dyad information which belongs to pedigree with gaps -> gap r values // please ensure a complete dyad_list with all dyads of interest (from WGS)
        if (! data_dyads) {
            cout << "unable to open file '"<<file_dyads <<".txt' for reading" << endl;
        }else{
            while(data_dyads>>indiv1>>indiv2>>dyad_name>>r_value>>paths>>pathline>>kinline>>common_anc>>depth>>kinlabel>>fullhalf>>real_r){
                if(dl_count%1000000==0){
                    cout << dl_count<<endl;
                }
                dyad x(indiv1,indiv2);
                if(r_value>=0 && r_value<=1){
                    x.set_gap_r_value(r_value);
                    x.set_path_str(paths);
                }else{
                    cout << "Erroneous r_value ("<<r_value<<") for dyad "<<dyad_name<<endl;
                    dyads_with_erroneous_r_values.push_back(dl_count);
                }
                if(real_r>=0 && real_r<=1){
                    x.set_real_r_value(real_r);
                    subset_idx_ptr->push_back(dl_count);
                }else{
                    cout << "check: Erroneous real_r ("<<real_r<<") for dyad "<<dyad_name<<" ... is currently ignored"<<endl;
                }
                dyad_dict_ptr->insert({dyad_name,dl_count});
                all_dyads_ptr->push_back(x);
                dl_count += 1;
            }
        }
        data_dyads.close(); 
        cout<<" ... dyad list loaded..."<<endl;
    }
    // post_calculation of r (erroneous values)
    cout << "dyads_with_erroneous_r_values:"<<endl;
    for (int i = 0;i<dyads_with_erroneous_r_values.size();i++){
        cout << dyads_with_erroneous_r_values[i]<<": "<<all_dyads_ptr->at(dyads_with_erroneous_r_values[i]).get_dyad_name()<<endl;
    }
    // STATS
    int noof_diff = 0; // total amount of r differences in dyads
    double sum_diff = 0;
    double max_diff = 0;
    for(int i = 0;i<all_dyads_ptr->size();i++){
        //cout <<all_dyads_ptr->at(i).get_r_infos();
        double rmg = all_dyads_ptr->at(i).get_r_diff("real_minus_gap");
        if(rmg>0 && rmg != -1){
            noof_diff += 1;
            sum_diff += rmg;
        }
        if(rmg>max_diff && rmg != -1){
            max_diff = rmg;
        }
    }
    cout << noof_diff<<"/"<< subset_idx_ptr->size()<<" of total " <<all_dyads_ptr->size()<<" dyads show a discrepancy in r (mean: "<<sum_diff/noof_diff<<", max: "<<max_diff<<")"<<endl;
    return normalize_factor;
}
std::string get_next_solution_all_gaps(int count,double temperature,std::deque<gap>*gaps_ptr,int default_year,int maturation_age_m,int maturation_age_f,int gestation_length,std::deque<node>*all_nodes_ptr,std::deque<node>*new_nodes_ptr){
    //cout << "\n\n ---------------------------- ["<<count<<"|"<<temperature<<"] NEXT SOLUTION  ---------------------------- \n";
    cout << "["<<count<<"|"<<temperature<<"] ";
    string all_changed_indivs = "";
    for(int i = 0;i<gaps_ptr->size();i++){ // update gap pointer
        gaps_ptr->at(i).set_indiv(&(new_nodes_ptr->at(gaps_ptr->at(i).get_indiv()->get_matidx())));
        gaps_ptr->at(i).set_pot_candidate(&(new_nodes_ptr->at(gaps_ptr->at(i).get_pot_candidate()->get_matidx())));
    }
    int rand_gap = rand()%gaps_ptr->size(); // choose random gap from pool -> which will be changed
    // ALTERING the choosen gap
    node* indiv_to_change = gaps_ptr->at(rand_gap).get_indiv(); // get direct access to the offspring of the to-be-altered gap-parent
    string parent_sex = gaps_ptr->at(rand_gap).get_sexseq()[0]; // get parent sex
    std::deque <node*> parent_pool = get_parent_pool(parent_sex,maturation_age_m,maturation_age_f, gestation_length,all_nodes_ptr,&(all_nodes_ptr->at(indiv_to_change->get_matidx()))); // get pool of potential parent to choose from   
    if(parent_pool.size()!=0 && indiv_to_change->get_DOB().get_year()!=default_year){ // theoretically indiv_to_change cannot be a founder individual, since  gap was choosen from already altered gaps 
        int rand_parent = -1;
        if(parent_pool.size()>1){ // choose new parent if possible
            rand_parent = rand()%parent_pool.size(); // choose random parent for the gap from pool of potential parents
            if(parent_pool.size()>1 && parent_pool[rand_parent]->get_name()==gaps_ptr->at(rand_gap).get_pot_candidate()->get_name()){ // prevent choosing the same parent as before (in this case: dont choose another randomly parent (it will be the same again), instead use next one in line)
                if(rand_parent < (parent_pool.size()-1)){
                    rand_parent += 1;
                }else{
                    rand_parent = 0;
                }  
            }
        }else{
            cout <<"WARNING. no other parent is possible. check."<<endl;
        }
        //cout << "indiv_to_change ("<<new_nodes_ptr->at(indiv_to_change->get_matidx()).get_name() << ")'s new random parent: "<< parent_pool[rand_parent]->get_name();
        all_changed_indivs += (new_nodes_ptr->at(indiv_to_change->get_matidx()).get_name() + "," + parent_pool[rand_parent]->get_name() + "," + gaps_ptr->at(rand_gap).get_pot_candidate()->get_name());
        if(parent_sex=="m"){ // gap filling: set missing sire to randomly choosen parent (from pool of potential parents): in new_all_nodes 
            new_nodes_ptr->at(indiv_to_change->get_matidx()).set_sire(new_nodes_ptr->at(parent_pool[rand_parent]->get_matidx()).get_name());
            new_nodes_ptr->at(indiv_to_change->get_matidx()).set_sire_node(&(new_nodes_ptr->at(parent_pool[rand_parent]->get_matidx())));
        }else if(parent_sex=="f"){// gap filling: set missing mom to randomly choosen parent (from pool of potential parents): in new_all_nodes 
            new_nodes_ptr->at(indiv_to_change->get_matidx()).set_mom(new_nodes_ptr->at(parent_pool[rand_parent]->get_matidx()).get_name());
            new_nodes_ptr->at(indiv_to_change->get_matidx()).set_mom_node(&(new_nodes_ptr->at(parent_pool[rand_parent]->get_matidx())));
        }
        //cout << " --> new indiv info: "<<new_nodes_ptr->at(indiv_to_change->get_matidx()).get_infos();
        //cout << "change gap "<< gaps_ptr->at(rand_gap).get_gap_infos()<< " --to--> ";
        gaps_ptr->at(rand_gap).set_pot_candidate(&(new_nodes_ptr->at(parent_pool[rand_parent]->get_matidx()))); // set missing parent in the gap_infos too
        //cout << gaps_ptr->at(rand_gap).get_gap_infos()<<endl;
    }else if(indiv_to_change->get_DOB().get_year()== default_year){
        cout << "WARNING. founder individual. shouldn't be possible.\n#\n#\n#\n#\n";
    }else{
        cout << "WARNING. possible ERROR. no possible parent could be found (empty parent_pool).\n#\n#\n#\n#\n";
    }
    return all_changed_indivs;
}
double compare_f_matrix(std::deque<std::deque<double>> &f_mat_1,std::deque<std::deque<double>> & f_mat_2,std::deque<node>*all_nodes_ptr,std::deque<node>*nodes_ptr_2, std::deque<int>*subset_idx_ptr, std::deque<dyad>*all_dyads_ptr){
    double diff_sum = 0;
    if(f_mat_1.size()==f_mat_2.size()){
        for(int i = 0;i<subset_idx_ptr->size();i++){
            auto [x,y] = all_dyads_ptr->at(subset_idx_ptr->at(i)).get_dyad_idx_in_f_matrix();
            if(x >= 0 && y >= 0){
            //for(int i = 0;i<f_mat_1.size();i++){
                //for(int j = 0;j<f_mat_1[i].size();j++){
                if(f_mat_1[x][y] < 0 || f_mat_1[x][y] > 1){
                    cout << f_mat_1[x][y];
                    cout << ": ERROR. f matrix (orig) comprise NA values in relevant places (x,y: "<<x<<" ("<<all_nodes_ptr->at(x).get_name()<<"),"<<y<<" ("<<all_nodes_ptr->at(y).get_name()<<") dyad "<<all_dyads_ptr->at(subset_idx_ptr->at(i)).get_dyad_infos()<<". please check. currently ignored"<<endl;
                }else if(f_mat_2[x][y] < 0 || f_mat_2[x][y] > 1){
                    cout << f_mat_2[x][y];
                    cout << ": ERROR. f matrix (current/new) comprise NA values in relevant places (x,y: "<<x<<" ("<<nodes_ptr_2->at(x).get_name()<<"),"<<y<<" ("<<nodes_ptr_2->at(y).get_name()<<") dyad "<<all_dyads_ptr->at(subset_idx_ptr->at(i)).get_dyad_infos()<<". please check. currently ignored"<<endl;
                }else{
                    diff_sum += (sqrt((f_mat_1[x][y] - f_mat_2[x][y])*(f_mat_1[x][y] - f_mat_2[x][y])));
                    //if((f_mat_1[i][j] - f_mat_2[i][j]) !=0){
                    //   cout << all_nodes_ptr->at(i).get_name() <<"_"<<all_nodes_ptr->at(j).get_name()<<": "<< f_mat_1[i][j] <<"-"<< f_mat_2[i][j]<<"="<<f_mat_1[i][j] - f_mat_2[i][j]<<endl;
                    //}
                }
            }else{
                cout << "set_idx() error: dyad index == -1"<<endl;
            }
        }
        return diff_sum;
    }else{
        cout << "ERROR. f matrices with unequal size"<<endl;
        return -1;
    }
}
double get_init_temp_factor(std::deque<std::deque<double>> &f_mat,std::deque<node>*all_nodes_ptr){
    double init_temp_factor = 0;
    for(int i = 2;i<f_mat.size();i++){
        int i_count = 0;
        double temp_factor = 0;
        for(int j = 2;j<f_mat[i].size();j++){
            if(i>j && f_mat[i][j]<=1 && f_mat[i][j]>=0){
                temp_factor += f_mat[i][j];
                i_count += 1;
            }else if(i<j && f_mat[i][j]<=1 && f_mat[i][j]>=0){
                temp_factor += f_mat[j][i];
                i_count += 1;
            }
        }
        if (temp_factor/i_count > init_temp_factor){
            init_temp_factor = temp_factor/i_count;
        }
    }
    return init_temp_factor;
}
void fill_pure_f_matrix(std::deque<dyad>* all_dyads_ptr,std::deque<int>* subset_idx_ptr,std::deque<node>* all_nodes_ptr,std::deque<node>* current_nodes_ptr,std::deque<std::deque<double>>* f_matrix,bool full_ped,int dyads_start,int dyads_end,bool multithreading,int thread_no){
    //ofstream out; //ofstream is the class for fstream package
    //out.open("output_thread_"+to_string(thread_no)+".txt");
    //if(next_solution == false){
    if(multithreading == false){
        dyads_start = 0;
        dyads_end = subset_idx_ptr->size();
    }else if(multithreading == true && dyads_end == 0){
        cout << "WARNING. unable to fill f matrix by multiple threads (dyads_end == 0)"<<endl;
    }
    //out << dyads_start<< " to "<<dyads_end<<endl;
    for(int i = dyads_start;i<dyads_end;i++){
        //out <<"["<< thread_no<<"] "<<ceil(((i-dyads_start)*100)/(dyads_end-dyads_start))<<"% \t("<< i <<" of "<<dyads_start<<".."<<dyads_end<<")"<<endl;
        auto [x,y] = all_dyads_ptr->at(subset_idx_ptr->at(i)).get_dyad_idx_in_f_matrix();// getting index for matrix cell --> [x][y]
        //out << "auto["<<x<<","<<y<<"]"<<endl;
        if(all_dyads_ptr->at(subset_idx_ptr->at(i)).get_indiv_1_name()!=all_nodes_ptr->at(0).get_name() // check that dyad definitely consists of non-imaginary nodes (real individuals)
                &&all_dyads_ptr->at(subset_idx_ptr->at(i)).get_indiv_2_name()!=all_nodes_ptr->at(0).get_name()
                &&all_dyads_ptr->at(subset_idx_ptr->at(i)).get_indiv_1_name()!=all_nodes_ptr->at(1).get_name()
                &&all_dyads_ptr->at(subset_idx_ptr->at(i)).get_indiv_2_name()!=all_nodes_ptr->at(1).get_name() 
                && x >= 0 && y >= 0){ // if == -1: set_dyad_idx error
                //out<<"calculate_pure_f_xy: ";
            double f_xy = calculate_pure_f_xy(&(current_nodes_ptr->at(all_dyads_ptr->at(subset_idx_ptr->at(i)).get_indiv_1_idx())),&(current_nodes_ptr->at(all_dyads_ptr->at(subset_idx_ptr->at(i)).get_indiv_2_idx())),all_dyads_ptr,all_nodes_ptr,f_matrix, full_ped);
            //out<< f_xy<<endl;
            f_matrix->at(x)[y]=f_xy;
        }  
    }
    /*}else{
        if(multithreading == false){
            dyads_start = 0;
            dyads_end = subset_idx_ptr->size();
        }else if(multithreading == true && dyads_end == 0){
            cout << "WARNING. unable to fill f matrix by multiple threads (dyads_end == 0)"<<endl;
        }
        out << dyads_start<< " to "<<dyads_end<<endl;
        for(int i = dyads_start;i<dyads_end;i++){
            out <<"["<< thread_no<<"] "<<ceil(((i-dyads_start)*100)/(dyads_end-dyads_start))<<"% \t("<< i <<" of "<<dyads_start<<".."<<dyads_end<<")"<<endl;
            ####
            auto [x,y] = all_dyads_ptr->at(subset_idx_ptr->at(i)).get_dyad_idx_in_f_matrix();// getting index for matrix cell --> [x][y]
            if(x!=999 && y!=999){ // != 999 if set_dyad_idx had worked
                double f_xy = calculate_pure_f_xy(&(current_nodes_ptr->at(all_dyads_ptr->at(subset_idx_ptr->at(i)).get_indiv_1_idx())),&(current_nodes_ptr->at(all_dyads_ptr->at(subset_idx_ptr->at(i)).get_indiv_2_idx())),all_dyads_ptr,all_nodes_ptr,f_matrix, full_ped);
                f_mat_new[x][y]=f_xy;
            }
        }
    }*/
    //out.close();
}
void simulated_annealing_complete_pedigree(double init_temp,double stop_temp,double temp_decay,std::deque<node>*all_nodes_ptr,std::deque<dyad>*all_dyads_ptr,std::deque<int>*subset_idx_ptr,int gestation_length,int maturation_age_f,int maturation_age_m,int default_year,string filepath,bool visualization,bool multithreading, int n_cores){ // simulated annealing not focal-wise -> instead simulate solution for complete pedigree at once
    std::chrono::steady_clock::time_point init_annealing_start = std::chrono::steady_clock::now();    
    cout << "start simulated annealing ... "<<endl;
    double new_r_diff = -1, current_r_diff = -1, best_r_diff = -1;
    std::deque<std::deque<double>> visualization_data = {};
    std::deque<node> * best_pedigree = nullptr;
    std::deque<gap> gaps = {};
    std::deque<string> sexseq = {};
    std::deque<std::deque<double>> f_mat_orig = {},f_mat_new = {};
    bool full_ped = all_nodes_ptr->size()==subset_idx_ptr->size();
    for(int i = 0;i<all_nodes_ptr->size();i++){ // init original f matrix 
        f_mat_orig.push_back(std::deque<double>(i+1, -1.0f));
        f_mat_orig[i][i] = 1;
        f_mat_orig[i][0] = 0;
        if(i>0){
            f_mat_orig[i][1] = 0;
        }
    }
    std::deque<std::deque<double>> f_mat_current = f_mat_orig, f_mat_default = f_mat_orig; // init current & new f matrix (based on init-orig-f-matrix)
    for(int i = 0;i<all_dyads_ptr->size();i++){
        all_dyads_ptr->at(i).set_idx(all_nodes_ptr);
    }
    cout<<"set_idx..."<<endl;
    for(int i = 0;i<subset_idx_ptr->size();i++){//all_dyads_ptr->size();i++){ // fill original f matrix based on the realized relatedness values (from input data, e.g. WGS or real r, stored in "real_r_value")
        //cout << i<<": "<<all_dyads_ptr->at(subset_idx_ptr->at(i)).get_dyad_infos()<<endl;
        auto [x,y] = all_dyads_ptr->at(subset_idx_ptr->at(i)).get_dyad_idx_in_f_matrix();
        if(x >= 0 && y >= 0){
            f_mat_orig[x][y] = all_dyads_ptr->at(subset_idx_ptr->at(i)).get_real_r_value();
        }else{
            cout << "set_idx() error: dyad index == -1"<<endl;
        }
    }
    cout<<"set fmatorig with realR..."<<endl;
    if(init_temp<stop_temp){
        init_temp = get_init_temp_factor(f_mat_orig,all_nodes_ptr)*all_nodes_ptr->size()*1.5;
    }else{
        cout << "calculated init_temp (not used): "<<get_init_temp_factor(f_mat_orig,all_nodes_ptr)*all_nodes_ptr->size()*1.5<<endl;
    }
    double temperature = init_temp;
    //cout << "ORIGINAL MATRIX"<<endl;
    //print_f_matrix(f_mat_orig,all_nodes_ptr);
    cout<<"temp: "<<temperature<<endl;
    std::deque<double> parent_stats = get_all_gaps(&gaps, all_nodes_ptr, default_year); // find all gaps within pedigree
    std::deque<node> new_nodes = {}, current_nodes = *(all_nodes_ptr);
    for(int i = 0;i<current_nodes.size();i++){
        current_nodes[i].create_parent_ptr(&(current_nodes),true);
    }
    cout<<"create parent ptr for current_nodes finished..."<<endl;
    //cout<<gaps.size()<<" gaps:"<<endl;
    //for(int i = 0;i<gaps.size();i++){
    //    cout << gaps[i].get_gap_infos()<<endl;
    //}
    std::chrono::steady_clock::time_point init_annealing_end = std::chrono::steady_clock::now();
    std::cout << "init_annealing = " << std::chrono::duration_cast<std::chrono::nanoseconds> (init_annealing_end - init_annealing_start).count() << "[ns]" << std::endl;
    std::chrono::steady_clock::time_point start_solution_start = std::chrono::steady_clock::now();
    get_start_solution_all_gaps(temperature,&gaps,all_nodes_ptr,&current_nodes,maturation_age_m,maturation_age_f,gestation_length,default_year);
    int dyad_start = 0, dyad_end = subset_idx_ptr->size();
    if(n_cores>dyad_end){
        n_cores = dyad_end;
    }
    //cout << dyad_start << " to "<<dyad_end<<endl;
    if(multithreading == true){
        std::vector<std::thread> threads = {};
        for(int i = 0;i<n_cores;i++){ 
            int dyads_start = (int) (i*floor(subset_idx_ptr->size()/n_cores));
            int dyads_end = (int) ((i+1)*floor(subset_idx_ptr->size()/n_cores));
            if(i==(n_cores-1)){
                dyads_end = subset_idx_ptr->size();
            }
            //string printer = "fill_pure_f_matrix thread "+to_string(i)+" from subset_idx "+to_string(dyads_start) +" to "+to_string(dyads_end);
            //cout <<printer <<endl;
            thread th1(fill_pure_f_matrix,all_dyads_ptr,subset_idx_ptr,all_nodes_ptr,&current_nodes,&f_mat_current, full_ped, dyads_start, dyads_end,multithreading, i+1);
            threads.push_back(std::move(th1));
        }
        for(int i = 0;i<threads.size();i++){
            threads[i].join();
        }
        cout << "treads joined after fill_pure_f_matrix (start_solution)"<<endl;
    }else{
        fill_pure_f_matrix(all_dyads_ptr,subset_idx_ptr,all_nodes_ptr,&current_nodes,&f_mat_current, full_ped);
    }
    //cout << " #### FIRST SOLUTION MATRIX ####"<<endl;
    //print_f_matrix(f_mat_current,&current_nodes);
    current_r_diff = compare_f_matrix(f_mat_orig,f_mat_current,all_nodes_ptr,&current_nodes,subset_idx_ptr,all_dyads_ptr);// make sure its greater than 0, just the magnitude of its difference
    //cout << "\n #### EVALUATION R ####\t\t diff (orig-first): "<< current_r_diff <<endl;
    all_nodes_to_population_file(&current_nodes,filepath+"_start_solution_"+to_string(current_r_diff));
    best_r_diff = current_r_diff; // set start_solution as currently best solution
    best_pedigree = &current_nodes;
    std::chrono::steady_clock::time_point start_solution_end = std::chrono::steady_clock::now();
    auto start_solution = std::chrono::duration_cast<std::chrono::nanoseconds> (start_solution_end - start_solution_start).count();
    std::cout << "start_solution = " << start_solution << "[ns]" << std::endl;
    std::deque<int> next_solution_time = {};
    int count = 0;
    cout << "start annealing with start temperature "<<temperature<<endl;
    while(temperature > stop_temp){
        //cout <<"NEXT SOLUTION (new): "<<to_string(count)<<endl;
        std::chrono::steady_clock::time_point next_solution_start = std::chrono::steady_clock::now();
        count += 1;
        new_nodes = current_nodes;
        //all_nodes_to_population_file(&new_nodes,filepath+"_next_"+to_string(count));
        f_mat_new = f_mat_current; // statt f_mat_default --> f_mat_current & dann über reduce_node_space dyads rausfinden, die neuberechnet werden müssen!!! oder einfach alle, die mit indiv, indivs old parent & indivs new parent dyad bilden (und in subset_idx enthalten)
        //write_matrix(all_nodes_ptr,&f_mat_new,filepath+"_fmat_"+to_string(count));
        //cout << "f_mat_new set to f_mat_current"<<endl;
        new_r_diff = 0;
        for(int i = 0;i<new_nodes.size();i++){ // make sure pointers point to the correct individual in the correct list of nodes
            new_nodes[i].create_parent_ptr(&(new_nodes),true);
            //cout << new_nodes[i].get_name()<<"|";
        }
        string all_changed_indivs = get_next_solution_all_gaps(count,temperature,&gaps,default_year,maturation_age_m,maturation_age_f,gestation_length,all_nodes_ptr,&new_nodes); // generate new solution
        cout << "[" <<all_changed_indivs<<"]";
        std::deque <string> changed_indivs = str_split(all_changed_indivs,",");
        std::deque<int> subset_changed_indivs = {};
        for(int i = 0;i<subset_idx_ptr->size();i++){ // iterate through all dyads and calculate their new relatedness coeffiecient r (after vs. before gap filling)
            if(all_dyads_ptr->at(subset_idx_ptr->at(i)).get_indiv_1_name()==changed_indivs[0] // check that dyad definitely consists of non-imaginary nodes (real individuals)
                ||all_dyads_ptr->at(subset_idx_ptr->at(i)).get_indiv_2_name()==changed_indivs[0]
                ||all_dyads_ptr->at(subset_idx_ptr->at(i)).get_indiv_1_name()==changed_indivs[1]
                ||all_dyads_ptr->at(subset_idx_ptr->at(i)).get_indiv_2_name()==changed_indivs[1]
                ||all_dyads_ptr->at(subset_idx_ptr->at(i)).get_indiv_1_name()==changed_indivs[2]
                ||all_dyads_ptr->at(subset_idx_ptr->at(i)).get_indiv_2_name()==changed_indivs[2]){
                subset_changed_indivs.push_back(subset_idx_ptr->at(i));
            }  
        }
        //cout << "subset_changed_indivs.size(): "<<subset_changed_indivs.size()<<endl;
        if(n_cores>subset_changed_indivs.size()){
            n_cores = subset_changed_indivs.size();
        }
        //cout << dyad_start << " to "<<dyad_end<<endl;
        if(multithreading == true){
            std::vector<std::thread> threads = {};
            for(int i = 0;i<n_cores;i++){ 
                int dyads_start = (int) (i*floor(subset_changed_indivs.size()/n_cores));
                int dyads_end = (int) ((i+1)*floor(subset_changed_indivs.size()/n_cores));
                if(i==(n_cores-1)){
                    dyads_end = subset_changed_indivs.size();
                }
                //cout << "fill_pure_f_matrix thread "<<i<<" from subset_changed_indivs "<<dyads_start<<" to "<<dyads_end<<endl;
                thread th1(fill_pure_f_matrix,all_dyads_ptr,&subset_changed_indivs,all_nodes_ptr,&new_nodes,&f_mat_new, full_ped, dyads_start, dyads_end,multithreading, i+1);
                threads.push_back(std::move(th1));
            }
            for(int i = 0;i<threads.size();i++){
                threads[i].join();
            }
        }else{
            fill_pure_f_matrix(all_dyads_ptr,&subset_changed_indivs,all_nodes_ptr,&new_nodes,&f_mat_new, full_ped);
        }
        
        //cout << "\n #### NEW MATRIX #### "<<endl;
        //print_f_matrix(f_mat_new,&new_nodes);
        new_r_diff = compare_f_matrix(f_mat_orig,f_mat_new,all_nodes_ptr,&new_nodes,subset_idx_ptr,all_dyads_ptr);// make sure its greater than 0, just the magnitude of its difference
        //cout << "\n #### EVALUATION R ####\t\t diff (orig-new): "<< new_r_diff <<" ";   
        
        // COMPARISON OF CURRENT AND NEW SOLUTION
        if(new_r_diff<best_r_diff){ // swap and go on with new solution as current solution if it is the best solution so far
            best_r_diff = new_r_diff;
            best_pedigree = &new_nodes;
            cout << " --> new_solution ("<<new_r_diff<<") is the next current_solution ("<<current_r_diff<<") (swapped because it is the best solution until now)"<<endl;
            if(visualization){
                std::deque<double> vis_vec = {temperature,new_r_diff,current_r_diff,1};
                visualization_data.push_back(vis_vec);
            }
            swap(current_r_diff,new_r_diff);
            swap(current_nodes,new_nodes);
            swap(f_mat_current,f_mat_new);
        }else{ // whether swap or not is based on probability (smaller r diff, the higher probility to go on with new solution)
            double prob = exp((-(new_r_diff-current_r_diff)/temperature)); // determine probability; is >1 if new_r_diff < current_r_diff (definitely swap)
            double random_number = randomnumber(); // choose random number between 0..1
            //cout << "\nprob: "<<prob<<" (diff (new-current): "<<new_r_diff-current_r_diff<<", temp: "<<temperature<<", div by temp: "<<(-(new_r_diff-current_r_diff)/temperature)<<") > randomnr: "<<random_number<<endl;
            if(prob > random_number){
                cout << " --> new_solution ("<<new_r_diff<<") is the next current_solution ("<<current_r_diff<<") (swapped";
                if(new_r_diff<current_r_diff){
                    cout << " because it is better)"<<endl;
                }else{
                    cout << " although it is worse)"<<endl;
                }
                if(visualization){
                    std::deque<double> vis_vec = {temperature,new_r_diff,current_r_diff,1};
                    visualization_data.push_back(vis_vec);
                }
                swap(current_nodes,new_nodes);
                swap(current_r_diff,new_r_diff);
                swap(f_mat_current,f_mat_new);
            }else{
                if(visualization){
                    std::deque<double> vis_vec = {temperature,new_r_diff,current_r_diff,0};
                    visualization_data.push_back(vis_vec);
                }
                cout << " --> not swapped"<<endl;
            }
        }
        /*cout << " ####   USED GAPS  #### \t";
        for(int i = 0;i<gaps.size();i++){
            cout << gaps[i].get_gap_infos()<<" ";
        }*/
        temperature = temperature * temp_decay; // reduce temperature
        std::chrono::steady_clock::time_point next_solution_end = std::chrono::steady_clock::now();
        int time_elapsed = std::chrono::duration_cast<std::chrono::nanoseconds> (next_solution_end - next_solution_start).count();
        next_solution_time.push_back(time_elapsed);
    }
    cout <<"\n\n ----------------------------------------------------------------------------\n END ------------------------------------------------------------------------\n";
    cout << "init_temp: "<<init_temp<<", temp_decay: "<<temp_decay<<endl;
    cout << "best_r_diff: " << best_r_diff<<endl;
    cout << "end_r_diff: "<< current_r_diff<<endl;
    string indivs_with_gap = "";
    for(int i = 0;i<gaps.size();i++){
        indivs_with_gap = indivs_with_gap + "|" + gaps[i].get_indiv()->get_name();
    }
    if(fs::exists(filepath+".txt")){
        int false_count = 0;
        for(int i = 2;i<current_nodes.size();i++){
            //cout << current_nodes[i].get_name()<<" "<<current_nodes[i].compare_pedigree_infos()<<endl;
            if(current_nodes[i].compare_pedigree_infos()!="true"){
                cout<<current_nodes[i].compare_pedigree_infos()<<endl;
                false_count += 1;
                if(indivs_with_gap.find(current_nodes[i].get_name()) == std::string::npos){
                    cout << "Warning. Indiv not to compare (not in gap list): "<<current_nodes[i].get_name()<<endl;
                }
            }
        }
        cout << to_string(false_count)<<"/"<<gaps.size()<<" pedigree errors"<<endl;
    }
    cout << "\n ----------------------------------------------------------------------------\n ----------------------------------------------------------------------------\n"<<endl;          
    if(visualization){
        ofstream out; 
        out.open(filepath+"_visualization.txt"); 
        out << "temperature\tnew_r_diff\tcurrent_r_diff\tswapped\n";
        for(int i = 0;i<visualization_data.size();i++){
            for(int j = 0;j<visualization_data[i].size();j++){
                out << visualization_data[i][j];
                if(j<visualization_data[i].size()-1){
                    out << "\t";
                }else{
                    out << "\n";
                }
            }
        }
        out.close();
        cout<<"writing "<<filepath+"_visualization.txt"<<endl;
    }
    cout << "parent_stats:\nmin: "<<parent_stats[0]<<"\nmax: "<<parent_stats[1]<<"\nmean: "<<parent_stats[2]<<"\ntotal sum: "<<parent_stats[3]<<"\nn_kombination: "<<parent_stats[4]<<endl;
    cout << "writing "<<all_nodes_to_population_file(&current_nodes,filepath+"_annealed");
    std::cout << "\n--> TIMING\nstart_solution = " << start_solution << "[ns]" << std::endl;
    std::chrono::steady_clock::time_point annealing_end = std::chrono::steady_clock::now();
    cout << "simulated_annealing_complete_pedigree = " << std::chrono::duration_cast<std::chrono::nanoseconds> (annealing_end - init_annealing_start).count() << "[ns]" << std::endl;
    int next_solution_time_sum = 0;
    for(int i = 0;i<next_solution_time.size();i++){
        next_solution_time_sum += next_solution_time[i];
    }
    cout << "next_solution_time_mean: "<< next_solution_time_sum/count<<" [ns]"<<endl;
}