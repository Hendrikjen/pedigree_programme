#include "simulated_annealing.h"

double randomnumber() { // get random double number between 0 and 1
    static std::default_random_engine randomnumber;
    std::uniform_real_distribution<double> dist(0.0, 1.0); 
    return dist(randomnumber); 
}
double load_data_for_sim_annealing(string file_gaps, string file_dyads, std::deque<dyad>*all_dyads_ptr,std::deque<int>*subset_idx_ptr, std::deque<node>*all_nodes_ptr,map<string, int> *dyad_dict_ptr,int maturation_age_f,int maturation_age_m,int gestation_length,string file_full){ // Load and prepare all relevant data for simulated annealing (gap pedigree and their respectively dyadic information + calculates basic statistic (ratio, mean, max and min of gap-resulting discrepancy)
    try{
        // PEDIGREE INFORMATION
        string sex,name,mom,sire,DOB,DOD,nonsires,nondams,indiv1, indiv2,dyad_name,paths,pathline,kinline,common_anc,depth,kinlabel,fullhalf;
        int birthseason;
        double r_value,real_r;
        node imaginary_female(name="unkn_f","f","NA","NA","NA","NA","NA","NA",0,0); // create male & female imaginary nodes as dummy in case of an unknown parent
        node imaginary_male(name="unkn_m","m","NA","NA","NA","NA","NA","NA",0,1);
        all_nodes_ptr->push_back(imaginary_female);
        all_nodes_ptr->push_back(imaginary_male);
        cout << "load pedigree (with gaps)"<<endl;
        ifstream data(file_gaps); // load pedigree (with gaps)
        if (! data) {
            throw std::runtime_error("Unable to open file '"+file_gaps+"' to load pedigree");
        }
        int matidx = 2;
        while(data >> name >> sex >> birthseason >> mom >> sire >> DOB >> DOD >> nonsires >> nondams){
            if(mom=="unknown"||mom=="unkn_f"||mom=="UNK"||mom=="NA"){
                mom=="unkn_f";
            }
            if(sire=="unknown"||sire=="unkn_m"||sire=="UNK"||sire=="NA"){
                sire=="unkn_m";
            }
            node x(name,sex,mom,sire,nonsires,nondams,DOB,DOD,birthseason,matidx); // create individuals (nodes) for each row
            all_nodes_ptr->push_back(x);
            matidx += 1;
        }
        data.close(); 
        if(file_full != "NA"){ // if correctly filled pedigree version is given (load it too, to compare accuracy of simulated annealing in the end)
            cout << "load pedigree (complete)"<<endl;
            ifstream data_full(file_full); // load completely known pedigree
            if (! data_full) {
                throw std::runtime_error("Unable to open file '"+file_full+"' to load full pedigree");
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
        }
        cout << "create parent pointers"<<endl;
        for(int i = 0;i<all_nodes_ptr->size();i++){ 
            all_nodes_ptr->at(i).create_parent_ptr(all_nodes_ptr);
        }
        
        // DYAD INFOMATION
        cout << "load dyad data"<<endl;
        double normalize_factor = 0; // sum of all realized relatedness values
        std::deque<int> dyads_with_erroneous_r_values = {};
        ifstream data_dyads(file_dyads); // load dyad information which belongs to pedigree with gaps -> gap r values // please ensure a complete dyad_list with all dyads of interest (from WGS)
        if (! data_dyads) {
            throw std::runtime_error("Unable to open file '"+file_dyads+"' to load dyad information");
        }
        int dl_count = 0;
        while(data_dyads>>indiv1>>indiv2>>r_value>>real_r){
            dyad x(indiv1,indiv2);
            if(r_value>=0 && r_value<=1){
                x.set_gap_r_value(r_value);
            }else{
                throw std::runtime_error("Erroneous r_value (pedigree-derived r) ("+to_string(r_value)+") for dyad "+dyad_name);
            }
            if(real_r>=0 && real_r<=1){
                x.set_real_r_value(real_r);
                normalize_factor += r_value;
            }else{
                throw std::runtime_error("Erroneous r_value (realized r) ("+to_string(r_value)+") for dyad "+dyad_name);
            }
            dyad_dict_ptr->insert({dyad_name,dl_count});
            all_dyads_ptr->push_back(x);
            subset_idx_ptr->push_back(dl_count);
            dl_count += 1;
        }
        data_dyads.close(); 

        int total_diff = 0; // number of dyads with discrepancy between realized relatedness and pedigree-derived relatedness
        double sum_diff = 0; // total amount of r difference in dyads
        double max_diff = 0; // maximal difference
        for(int i = 0;i<all_dyads_ptr->size();i++){ // iterate through all nodes
            double rmg = all_dyads_ptr->at(i).get_r_diff("real_minus_gap"); // get difference between realized relatedness and calculated relatedness from incomplete pedigree
            if(rmg>0 && rmg != -1){
                total_diff += 1;
                sum_diff += rmg;
            }
            if(rmg>max_diff && rmg != -1){
                max_diff = rmg;
            }
        }
        cout << total_diff<<"/"<< subset_idx_ptr->size()<<" dyads show a discrepancy in r (mean: "<<sum_diff/total_diff<<", max: "<<max_diff<<")"<<endl;
        if(subset_idx_ptr->size()==0){
            throw std::runtime_error("zero dyads to analyze - potential error while loading dyad data; please check input requirements.");
        }
        return normalize_factor;
    }catch(const std::exception &ex) {
        std::cerr << "Error in load_data_for_sim_annealing(): " << ex.what() << std::endl;
        return 0.0;
    }
}
string get_next_solution(int count,double temperature,std::deque<gap>*gaps_ptr,int maturation_age_m,int maturation_age_f,int gestation_length,bool twins,std::deque<node>*all_nodes_ptr,std::deque<node>*new_nodes_ptr){ // get new pedigree solution: randomly choose one gap to exchange the parental candidate & calculate the new overall discrepancy between realized and pedigree-derived relatedness values
    try{
        cout << "["<<count<<"|"<<temperature<<"] ";
        string all_changed_indivs = "";
        for(int i = 0;i<gaps_ptr->size();i++){ // update gap pointer (to point to individuals in new_nodes)
            gaps_ptr->at(i).set_indiv(&(new_nodes_ptr->at(gaps_ptr->at(i).get_indiv()->get_matidx())));
            gaps_ptr->at(i).set_pot_candidate(&(new_nodes_ptr->at(gaps_ptr->at(i).get_pot_candidate()->get_matidx())));
        }
        int rand_gap = rand()%gaps_ptr->size(); // choose random gap from gap pool -> which will be changed
        node* indiv_to_change = gaps_ptr->at(rand_gap).get_indiv(); // get direct access to the offspring of the to-be-altered gap-parent
        string parent_sex = gaps_ptr->at(rand_gap).get_sexseq()[0]; // get parent sex
        std::deque <node*> parent_pool = get_parent_pool(parent_sex,maturation_age_m,maturation_age_f, gestation_length,all_nodes_ptr,&(all_nodes_ptr->at(indiv_to_change->get_matidx())),twins); // get pool of selected individual to change to choose new potential parent from   
        if(parent_pool.size()!=0){ // ensuring that indiv to change is not a founder (should not be)
            int rand_parent = -1;
            if(parent_pool.size()>1){ // choose new parent if possible (only one parental candidate - should also not be)
                rand_parent = rand()%parent_pool.size(); // choose random parent for the gap from pool of potential parents
                if(parent_pool[rand_parent]->get_name()==gaps_ptr->at(rand_gap).get_pot_candidate()->get_name()){ // prevent choosing the same parent as before (in this case: dont choose another randomly parent (it will be the same again), instead use next one in line)
                    if(rand_parent < (parent_pool.size()-1)){
                        rand_parent += 1;
                    }else{
                        rand_parent = 0;
                    }  
                }
                all_changed_indivs += (new_nodes_ptr->at(indiv_to_change->get_matidx()).get_name() + "," + parent_pool[rand_parent]->get_name() + "," + gaps_ptr->at(rand_gap).get_pot_candidate()->get_name()); // save the changed individual in a string representation to return later
                if(parent_sex=="m"){ // gap filling: set missing sire to randomly choosen parent (from pool of potential parents): in new_all_nodes 
                    new_nodes_ptr->at(indiv_to_change->get_matidx()).set_sire(new_nodes_ptr->at(parent_pool[rand_parent]->get_matidx()).get_name());
                    new_nodes_ptr->at(indiv_to_change->get_matidx()).set_sire_node(&(new_nodes_ptr->at(parent_pool[rand_parent]->get_matidx())));
                }else if(parent_sex=="f"){// gap filling: set missing mom to randomly choosen parent (from pool of potential parents): in new_all_nodes 
                    new_nodes_ptr->at(indiv_to_change->get_matidx()).set_mom(new_nodes_ptr->at(parent_pool[rand_parent]->get_matidx()).get_name());
                    new_nodes_ptr->at(indiv_to_change->get_matidx()).set_mom_node(&(new_nodes_ptr->at(parent_pool[rand_parent]->get_matidx())));
                }
                gaps_ptr->at(rand_gap).set_pot_candidate(&(new_nodes_ptr->at(parent_pool[rand_parent]->get_matidx()))); // set missing parent in the gap_infos too
            }else{
                throw std::runtime_error("Only one potential parent candidate in get_next_solution() for "+indiv_to_change->get_name());
            }            
        }else{
            throw std::runtime_error("Empty parent pool in get_next_solution()");
        }
        return all_changed_indivs;
    }catch(const std::exception &ex) {
        std::cerr << "Error in get_next_solution(): " << ex.what() << std::endl;
        return "unable to get next solution";
    }
}
double compare_f_matrix(std::deque<std::deque<double>> &f_mat_1,std::deque<std::deque<double>> & f_mat_2,std::deque<node>*all_nodes_ptr,std::deque<node>*nodes_ptr_2, std::deque<int>*subset_idx_ptr, std::deque<dyad>*all_dyads_ptr){ // compare two relatedness matrices (current solution vs. realized relatedness) to get the current relatedness over all difference between both 
    try{
        double diff_sum = 0; // to add up all r differences
        if(f_mat_1.size()==f_mat_2.size()){ // check if matrix dimensions match (should be)
            for(int i = 0;i<subset_idx_ptr->size();i++){ // iterate through all relevant dyads
                auto [x,y] = all_dyads_ptr->at(subset_idx_ptr->at(i)).get_dyad_idx_in_f_matrix(); // get matrix coordinates for specific dyad
                if(x >= 0 && y >= 0){
                    if(f_mat_1[x][y] < 0 || f_mat_1[x][y] > 1){ // invalid r value in matrix 1
                        throw std::runtime_error("Original relatedness matrix comprise NA values in relevant places (x,y: "+to_string(x)+" ("+all_nodes_ptr->at(x).get_name()+"),"+to_string(y)+" ("+all_nodes_ptr->at(y).get_name()+"): "+to_string(f_mat_1[x][y]));
                    }else if(f_mat_2[x][y] < 0 || f_mat_2[x][y] > 1){ // invalid r value in matrix 2
                        throw std::runtime_error("Original relatedness matrix comprise NA values in relevant places (x,y: "+to_string(x)+" ("+nodes_ptr_2->at(x).get_name()+"),"+to_string(y)+" ("+nodes_ptr_2->at(y).get_name()+"): "+to_string(f_mat_1[x][y]));
                    }else{
                        diff_sum += (sqrt((f_mat_1[x][y] - f_mat_2[x][y])*(f_mat_1[x][y] - f_mat_2[x][y]))); // get difference & add to diff_sum
                    }
                }else{
                    throw std::runtime_error("set_idx() error: dyad index == -1");
                }
            }
            return diff_sum;
        }else{
            throw std::runtime_error("Relatedness matrices of unequal size");
            return -1;
        }
    }catch(const std::exception &ex) {
        std::cerr << "Error in compare_f_matrix(): " << ex.what() << std::endl;
        return -1;
    }
}
double get_init_temp_factor(std::deque<std::deque<double>> &f_mat,std::deque<node>*all_nodes_ptr){ // if not given, start temperature for simulated annealing is calculated automatically -> here: get init_temp_factor
    try{
        double init_temp_factor = 0;
        for(int i = 2;i<f_mat.size();i++){ 
            int i_count = 0;
            double temp_factor = 0;
            for(int j = 2;j<f_mat[i].size();j++){// iterate through all cells [i][j] of the relatedness matrix (realized)
                if(i>j && f_mat[i][j]<=1 && f_mat[i][j]>=0){ // only if its a cell of a relevant dyad (filled) and not of individuals to itself (i==j)
                    temp_factor += f_mat[i][j]; // add relatedness value to temp_factor
                    i_count += 1;
                }else if(i<j && f_mat[i][j]<=1 && f_mat[i][j]>=0){
                    temp_factor += f_mat[j][i];// add relatedness value to temp_factor
                    i_count += 1;
                }
            }
            if (temp_factor/i_count > init_temp_factor){ // get the highest mean relatedness of an individual for relevant dyads
                init_temp_factor = temp_factor/i_count;
            }
        }
        return init_temp_factor;
    }catch(const std::exception &ex) {
        std::cerr << "Error in get_init_temp_factor(): " << ex.what() << std::endl;
        return 0;
    }
}
void fill_pure_f_matrix(std::deque<dyad>* all_dyads_ptr,std::deque<int>* subset_idx_ptr,std::deque<node>* all_nodes_ptr,std::deque<node>* current_nodes_ptr,std::deque<std::deque<double>>* f_matrix,bool full_ped,int dyads_start,int dyads_end,bool multithreading,int thread_no){ // help function to start calculating the relatedness coefficient (without path characteristics)
    try{
        if(multithreading == false){ // prepare multithreading if requested
            dyads_start = 0;
            dyads_end = subset_idx_ptr->size();
        }else if(multithreading == true && dyads_end == 0){
            throw std::runtime_error("Unable to fill relatedness matrix by multiple threads (no dyads -> dyads_end == 0)");
        }
        for(int i = dyads_start;i<dyads_end;i++){
            auto [x,y] = all_dyads_ptr->at(subset_idx_ptr->at(i)).get_dyad_idx_in_f_matrix();// getting index for matrix cell --> [x][y]
            if(all_dyads_ptr->at(subset_idx_ptr->at(i)).get_indiv_1_name()!=all_nodes_ptr->at(0).get_name() // check that dyad definitely consists of non-imaginary nodes (real individuals)
                    &&all_dyads_ptr->at(subset_idx_ptr->at(i)).get_indiv_2_name()!=all_nodes_ptr->at(0).get_name()
                    &&all_dyads_ptr->at(subset_idx_ptr->at(i)).get_indiv_1_name()!=all_nodes_ptr->at(1).get_name()
                    &&all_dyads_ptr->at(subset_idx_ptr->at(i)).get_indiv_2_name()!=all_nodes_ptr->at(1).get_name() 
                    && x >= 0 && y >= 0){ // if == -1: set_dyad_idx error
                double f_xy = calculate_pure_f_xy(&(current_nodes_ptr->at(all_dyads_ptr->at(subset_idx_ptr->at(i)).get_indiv_1_idx())),&(current_nodes_ptr->at(all_dyads_ptr->at(subset_idx_ptr->at(i)).get_indiv_2_idx())),all_dyads_ptr,all_nodes_ptr,f_matrix, full_ped); // calculate dyadic relatedness
                f_matrix->at(x)[y]=f_xy;
            }  
        }
    }catch(const std::exception &ex) {
        std::cerr << "Error in fill_pure_f_matrix(): " << ex.what() << std::endl;
    }
}
void simulated_annealing(string input_pedigree,string input_dyadlist,string output,double init_temp,double stop_temp,double temp_decay,int cores,int maturation_age_f,int maturation_age_m,int gestation_length,bool twins,bool visualization,string complete_pedigree){ // start simulated annealing to get pedigree solution which explains best the relatedness value differences in all dyads while various combination of potential parent candidates were tested
    try{
        // load data and initialize variables, objects, deques...
        std::deque <dyad> all_dyads={};
        std::deque<int> subset_idx = {};
        std::deque <node> all_nodes = {};
        std::map<string, int> dyad_dict = {};
        load_data_for_sim_annealing(input_pedigree,input_dyadlist,&all_dyads,&subset_idx,&all_nodes,&dyad_dict,maturation_age_f,maturation_age_m,gestation_length,complete_pedigree); // load all given pedigree and dyadic data 
        cout << "set parent pool"<<endl;
        set_parent_pool(&all_nodes,maturation_age_m,maturation_age_f,gestation_length,twins); // save potential parents as node attribute
        cout << "number of nodes: "<<all_nodes.size()<<endl;
        cout << "number of dyads: "<<all_dyads.size()<<", subset (relevant dyads): "<<subset_idx.size()<<endl;
        double new_r_diff = -1, current_r_diff = -1, best_r_diff = -1;
        std::deque<std::deque<double>> visualization_data = {};
        std::deque<node> * best_pedigree = nullptr;
        std::deque<gap> gaps = {};
        std::deque<string> sexseq = {};
        std::deque<std::deque<double>> f_mat_orig = {},f_mat_new = {};
        bool full_ped = all_nodes.size()==subset_idx.size();
        cout << "initialize and fill relatedness matrices based on dyad data"<<endl;
        for(int i = 0;i<all_nodes.size();i++){ // initialize original f matrix
            f_mat_orig.push_back(std::deque<double>(i+1, -1.0f)); // -1 as default relatedness value
            f_mat_orig[i][i] = 1; // relatedness to itself
            f_mat_orig[i][0] = 0; // not related to imaginary mom
            if(i>0){
                f_mat_orig[i][1] = 0; // not related to imaginary sire
            }
        }
        std::deque<std::deque<double>> f_mat_current = f_mat_orig; // init current f matrix (based on init-orig-f-matrix)
        for(int i = 0;i<all_dyads.size();i++){ // for each dyad, find node index and link it to the focal
            all_dyads.at(i).set_idx(&all_nodes);
        }
        for(int i = 0;i<subset_idx.size();i++){// fill original f matrix based on the realized relatedness values (from input data, e.g. WGS or real r, stored in "real_r_value")
            auto [x,y] = all_dyads.at(subset_idx.at(i)).get_dyad_idx_in_f_matrix(); // determine matrix indices
            if(x >= 0 && y >= 0){
                f_mat_orig[x][y] = all_dyads.at(subset_idx.at(i)).get_real_r_value();
            }else{
                throw std::runtime_error("Invalid allocation of matrix indices -> dyad index of "+all_dyads[subset_idx[i]].get_dyad_name()+" is -1");
            }
        }
        
        if(init_temp<stop_temp){ // if no user-defined initial temperature is given, calculate an appropiate init_temp
            init_temp = get_init_temp_factor(f_mat_orig,&all_nodes)*all_nodes.size()*1.5;
        }else{
            cout << "calculated init_temp (but not used): "<<get_init_temp_factor(f_mat_orig,&all_nodes)*all_nodes.size()*1.5<<endl;
        }
        double temperature = init_temp;
        std::deque<double> parent_stats = get_all_gaps(&gaps, &all_nodes); // find all gaps within pedigree
        
        cout << "create parent pointer for first pedigree solution (current_nodes)"<<endl;
        std::deque<node> new_nodes = {}, current_nodes = all_nodes;
        for(int i = 0;i<current_nodes.size();i++){
            current_nodes[i].create_parent_ptr(&(current_nodes),true);
        }
        
        cout << "get start solution"<<endl;
        get_start_solution(temperature,&gaps,&all_nodes,&current_nodes,maturation_age_m,maturation_age_f,gestation_length,twins);
        cout << "calculation of relatedness coefficients in start solution"<<endl;
        int dyad_start = 0, dyad_end = subset_idx.size(); // prepare splitting up for multithreading 
        if(cores>dyad_end){ //limit cores if they are overabundant
            cores = dyad_end;
        }
        if(cores > 1){ // if multithreading on multiple cores was requested
            std::vector<std::thread> threads = {};
            for(int i = 0;i<cores;i++){ // assign task to each core
                int dyads_start = (int) (i*floor(subset_idx.size()/cores)); // determine an almost uniform number of dyads to process for each core
                int dyads_end = (int) ((i+1)*floor(subset_idx.size()/cores));
                if(i==(cores-1)){ // last core
                    dyads_end = subset_idx.size();
                }
                thread th1(fill_pure_f_matrix,&all_dyads,&subset_idx,&all_nodes,&current_nodes,&f_mat_current, full_ped, dyads_start, dyads_end,cores>1, i+1); // calculation of relatedness coefficient -> saved in f matrix
                threads.push_back(std::move(th1));
            }
            for(int i = 0;i<threads.size();i++){ // wait for all threads to finish
                threads[i].join();
            }
        }else{
            fill_pure_f_matrix(&all_dyads,&subset_idx,&all_nodes,&current_nodes,&f_mat_current, full_ped); // calculation of relatedness coefficient without multithreading
        }
        current_r_diff = compare_f_matrix(f_mat_orig,f_mat_current,&all_nodes,&current_nodes,&subset_idx,&all_dyads);// evaluate the new solution (get relatedness variance)
        all_nodes_to_population_file(&current_nodes,output+"_start_solution");
        best_r_diff = current_r_diff; // set start_solution as currently best solution
        best_pedigree = &current_nodes;
        
        int count = 0;
        cout << "start simulated annealing loop"<<endl;
        cout << "start temperature: "<<temperature<<endl;
        cout << "stop temperature: "<<stop_temp<<endl;
        cout << "factor to moderately reduce temperature: "<<temp_decay<<endl;
        while(temperature > stop_temp){
            count += 1;
            new_nodes = current_nodes; // use current solution as starting point to the next step
            f_mat_new = f_mat_current; 
            new_r_diff = 0;
            for(int i = 0;i<new_nodes.size();i++){ // make sure pointers point to the correct individual in the correct list of nodes
                new_nodes[i].create_parent_ptr(&(new_nodes),true);
            }
            string all_changed_indivs = get_next_solution(count,temperature,&gaps,maturation_age_m,maturation_age_f,gestation_length,twins,&all_nodes,&new_nodes); // generate new solution (exchange one potential parent candidate)
            cout << "[" <<all_changed_indivs<<"]"; // get_next_solution returns also the IDs of the individuals, which were changed (offspring with the originally missing parent, old and new parent candidate)
            std::deque <string> changed_indivs = str_split(all_changed_indivs,",");
            std::deque<int> subset_changed_indivs = {};
            for(int i = 0;i<subset_idx.size();i++){ // iterate through all dyads and filter for all dyads that relatedness coefficient might have changed due to the new solution
                if(all_dyads.at(subset_idx.at(i)).get_indiv_1_name()==changed_indivs[0] 
                    ||all_dyads.at(subset_idx.at(i)).get_indiv_2_name()==changed_indivs[0]
                    ||all_dyads.at(subset_idx.at(i)).get_indiv_1_name()==changed_indivs[1]
                    ||all_dyads.at(subset_idx.at(i)).get_indiv_2_name()==changed_indivs[1]
                    ||all_dyads.at(subset_idx.at(i)).get_indiv_1_name()==changed_indivs[2]
                    ||all_dyads.at(subset_idx.at(i)).get_indiv_2_name()==changed_indivs[2]){
                    subset_changed_indivs.push_back(subset_idx.at(i));
                }  
            }
            if(cores>subset_changed_indivs.size()){ // prepare multiprocessing (limit cores if they are overabundant)
                cores = subset_changed_indivs.size();
            }
            if(cores > 1){ // if multiprocessing was requested
                std::vector<std::thread> threads = {};
                for(int i = 0;i<cores;i++){  // assing task to each core
                    int dyads_start = (int) (i*floor(subset_changed_indivs.size()/cores)); // determine an almost uniform number of dyads to process for each core
                    int dyads_end = (int) ((i+1)*floor(subset_changed_indivs.size()/cores));
                    if(i==(cores-1)){
                        dyads_end = subset_changed_indivs.size();
                    }
                    thread th1(fill_pure_f_matrix,&all_dyads,&subset_changed_indivs,&all_nodes,&new_nodes,&f_mat_new, full_ped, dyads_start, dyads_end,cores>1, i+1); // calculate their new relatedness coeffiecient r (after vs. before gap filling)
                    threads.push_back(std::move(th1));
                }
                for(int i = 0;i<threads.size();i++){ // wait for all threads to finish
                    threads[i].join();
                }
            }else{
                fill_pure_f_matrix(&all_dyads,&subset_changed_indivs,&all_nodes,&new_nodes,&f_mat_new, full_ped); // calculate their new relatedness coeffiecient r (after vs. before gap filling) without multiprocessing
            }
            new_r_diff = compare_f_matrix(f_mat_orig,f_mat_new,&all_nodes,&new_nodes,&subset_idx,&all_dyads);// evaluate the new solution (get relatedness variance)
            
            if(new_r_diff<best_r_diff){ // swap and definitely go on with new solution as current solution if it is the best solution so far
                best_r_diff = new_r_diff;
                best_pedigree = &new_nodes;
                cout << " --> new_solution ("<<new_r_diff<<") is the next current_solution ("<<current_r_diff<<") (swapped because it is the best solution until now)"<<endl;
                if(visualization){ // if requested save some details about the step in extra deque
                    std::deque<double> vis_vec = {temperature,new_r_diff,current_r_diff,1};
                    visualization_data.push_back(vis_vec);
                }
                swap(current_r_diff,new_r_diff);
                swap(current_nodes,new_nodes);
                swap(f_mat_current,f_mat_new);
            }else{ // whether to swap or not is based on probability (the smaller the r diff is, the higher the probability to go on with new solution)
                double prob = exp((-(new_r_diff-current_r_diff)/temperature)); // determine probability; is >1 if new_r_diff < current_r_diff (definitely swap)
                double random_number = randomnumber(); // choose random number between 0..1
                if(prob > random_number){ // go on with the new solution as next starting point for the next solution (swap)
                    cout << " --> new_solution ("<<new_r_diff<<") is the next current_solution ("<<current_r_diff<<") (swapped";
                    if(new_r_diff<current_r_diff){
                        cout << " because it is better)"<<endl;
                    }else{
                        cout << " although it is worse)"<<endl;
                    }
                    if(visualization){// if requested save some details about the step in extra deque
                        std::deque<double> vis_vec = {temperature,new_r_diff,current_r_diff,1};
                        visualization_data.push_back(vis_vec);
                    }
                    swap(current_nodes,new_nodes);
                    swap(current_r_diff,new_r_diff);
                    swap(f_mat_current,f_mat_new);
                }else{ // acceptance criterion was not fulfilled, use current solution again as starting point for next solution (no swap)
                    if(visualization){// if requested save some details about the step in extra deque
                        std::deque<double> vis_vec = {temperature,new_r_diff,current_r_diff,0};
                        visualization_data.push_back(vis_vec);
                    }
                    cout << " --> not swapped"<<endl;
                }
            }
            temperature = temperature * temp_decay; // reduce temperature
        }
        cout <<"\n\n ----------------------------------------------------------------------------\n END ------------------------------------------------------------------------\n";
        cout << "init_temp: "<<init_temp<<", temp_decay: "<<temp_decay<<endl;
        cout << "best_r_diff: " << best_r_diff<<endl;
        cout << "end_r_diff: "<< current_r_diff<<endl;
        
        if(complete_pedigree != "NA"){ // if the correctly reconstructed pedigree exist, check the accuracy of the assigned potential parents in all original gaps
            int false_count = 0; // counts the misassigned parent candidates
            for(int i = 2;i<current_nodes.size();i++){ // iterate through the end solution
                if(current_nodes[i].compare_pedigree_infos()!="true"){ // if parental mismatch -> print info
                    cout<<current_nodes[i].compare_pedigree_infos()<<endl;
                    false_count += 1;
                }
            }
            cout << to_string(false_count)<<"/"<<gaps.size()<<" pedigree errors"<<endl; // print how many parental mismatches at all
        }
        cout << "\n ----------------------------------------------------------------------------\n ----------------------------------------------------------------------------\n"<<endl;          
        if(visualization){ // if visualization requested, save simulated annealing step info in extra file
            ofstream out; 
            out.open(output+"_visualization.txt"); 
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
            cout<<"writing "<<output+"_visualization.txt"<<endl;
        }
        cout << "parent_stats:\nmin: "<<parent_stats[0]<<"\nmax: "<<parent_stats[1]<<"\nmean: "<<parent_stats[2]<<"\ntotal sum: "<<parent_stats[3]<<"\nn_kombination: "<<parent_stats[4]<<endl;
        cout << "writing "<<all_nodes_to_population_file(&current_nodes,output+"_annealed");
    }catch(const std::exception &ex) {
        std::cerr << "Error in simulated_annealing(): " << ex.what() << std::endl;
    }
}