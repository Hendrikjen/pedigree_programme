#include "general_class_functions.h"

void set_parent_pool(std::deque <node> *all_nodes_ptr, int maturation_age_m,int maturation_age_f,int gestation_length,bool twins){ //set the pool of potential parents for each node (index in deque pool, string representation & full generation)
    try{
        for(int i = 2;i<all_nodes_ptr->size();i++){ // iterate over all nodes (except imaginary nodes)
            string pot_moms = "NA";
            string pot_sires = "NA";
            if(all_nodes_ptr->at(i).get_mom()=="unknown" || all_nodes_ptr->at(i).get_mom()=="NA" || all_nodes_ptr->at(i).get_mom()=="unkn_f"){ // if mom not known
                std::deque <node*> mom_pool = get_parent_pool("f",maturation_age_m,maturation_age_f,gestation_length,all_nodes_ptr,&(all_nodes_ptr->at(i)),twins); // get all potential moms
                for(int j = 0;j<mom_pool.size();j++){ // iterate through potential moms
                    all_nodes_ptr->at(i).push_back_pot_parent(mom_pool[j]->get_sex(),mom_pool[j]->get_matidx()); //add potential parent index (to mom pool deque)
                    if(j==0){ // string them together with @ as delimiter
                        pot_moms = mom_pool[j]->get_name();
                    }else{
                        pot_moms = pot_moms + "@" + mom_pool[j]->get_name();
                    }
                }
            }
            if(all_nodes_ptr->at(i).get_sire()=="unknown" || all_nodes_ptr->at(i).get_sire()=="NA" || all_nodes_ptr->at(i).get_sire()=="unkn_m"){ // if sire not known
                std::deque <node*> sire_pool = get_parent_pool("m",maturation_age_m,maturation_age_f,gestation_length,all_nodes_ptr,&(all_nodes_ptr->at(i)),twins); // get all potential sires
                for(int j = 0;j<sire_pool.size();j++){ // iterate through potential sire
                    all_nodes_ptr->at(i).push_back_pot_parent(sire_pool[j]->get_sex(),sire_pool[j]->get_matidx()); // add potential parent index (to sire pool deque)
                    if(j==0){ // string them together with @ as delimiter
                        pot_sires = sire_pool[j]->get_name();
                    }else{
                        pot_sires = pot_sires + "@" + sire_pool[j]->get_name();
                    }
                }
            }
            all_nodes_ptr->at(i).set_pot_moms(pot_moms); // add string representation as additional attribute (mom)
            all_nodes_ptr->at(i).set_pot_sires(pot_sires); // add string representation as additional attribute (sire)
            all_nodes_ptr->at(i).set_full_generations(get_full_generation_depth(&all_nodes_ptr->at(i))); // determine & set the number of complete generations for the focal
        }
    }catch(const std::exception &ex) {
        std::cerr << "Error in set_parent_pool(): " << ex.what() << std::endl;
    }
}
void get_ancestor_of_focal_indiv(std::deque<node*>*anc_nodes_ptr,node*indiv,std::map<string,std::deque<node*>>*offspring_dict){ // get deque of all ancestor pointers of a given individual & fill parent-offspring dictionary
    try{
        int further_ancs_idx = 0; // "counter" which indiv in the deque will be the next to process (get mom/sire)
        anc_nodes_ptr->push_back(indiv); // to collect all ancestors of indiv (and indiv itself)
        while(further_ancs_idx < anc_nodes_ptr->size()){ // as long as there are unprocessed indivs (who might come along with more ancestors)
            if(anc_nodes_ptr->at(further_ancs_idx)->get_mom()!="unkn_f"){ // if mom is known 
                anc_nodes_ptr->push_back(anc_nodes_ptr->at(further_ancs_idx)->get_mom_node()); // add to anc_nodes_ptr
                if(offspring_dict->count(anc_nodes_ptr->at(further_ancs_idx)->get_mom_node()->get_name())==0){// if key (mom) not already initialized -> initialize key-value pair with empty std::deque
                    offspring_dict->insert(pair<string,std::deque<node*>>(anc_nodes_ptr->at(further_ancs_idx)->get_mom_node()->get_name(),std::deque<node*>())); 
                }
                offspring_dict->at(anc_nodes_ptr->at(further_ancs_idx)->get_mom_node()->get_name()).push_back(anc_nodes_ptr->at(further_ancs_idx)); // add current indiv to mother(key)/offspring_list(value) dictionary
            }
            if(anc_nodes_ptr->at(further_ancs_idx)->get_sire()!="unkn_m"){ // if sire is known 
                anc_nodes_ptr->push_back(anc_nodes_ptr->at(further_ancs_idx)->get_sire_node()); // add to anc_nodes_ptr
                if(offspring_dict->count(anc_nodes_ptr->at(further_ancs_idx)->get_sire_node()->get_name())==0){ // key (sire) not already initialized -> initialize key-value pair with empty std::deque 
                    offspring_dict->insert(pair<string,std::deque<node*>>(anc_nodes_ptr->at(further_ancs_idx)->get_sire_node()->get_name(),std::deque<node*>()));
                }
                offspring_dict->at(anc_nodes_ptr->at(further_ancs_idx)->get_sire_node()->get_name()).push_back(anc_nodes_ptr->at(further_ancs_idx)); // add current indiv to mother(key)/offspring_list(value) dictionary
            }
            further_ancs_idx += 1;
        }
    }catch(const std::exception &ex) {
        std::cerr << "Error in get_ancestor_of_focal_indiv(): " << ex.what() << std::endl;
    }
}
bool contains_ancestor(node* pot_offspring,node* pot_ancestor,std::deque <node> *all_nodes_ptr){ // returns boolean value whether an individual (pot_ancestor) is an ancestor of an other individual (pot_offspring) 
    try{
        bool is_ancestor = false; // default
        std::deque <node*> indiv_stack;
        indiv_stack.push_back(pot_offspring);
        if(((*pot_offspring).get_birthseason()==0)||((*pot_ancestor).get_birthseason()==0)||
            ((*pot_ancestor).get_birthseason()<(*pot_offspring).get_birthseason())){ // preselection -> exclude ancestors born after pot_offspring
            while(is_ancestor==false&& indiv_stack.size()>0){ // once ancestry is proven -> stop, else iterate through indiv_stack
                node* mom = (*indiv_stack[0]).mom_node; // get mom of first individual on stack 
                node* sire = (*indiv_stack[0]).sire_node; // get sire of first individual on stack
                if((*pot_ancestor).get_name()==(*mom).get_name()||(*pot_ancestor).get_name()==(*sire).get_name()){ // check whether mom/sire of current individual match with pot_ancestor
                    is_ancestor = true; // ancestor found
                }
                if((*mom).get_name()!=all_nodes_ptr->at(0).get_name()&&(*mom).get_name()!="NA"){ //push mom to stack if mom is a known individual (not unknown)
                    indiv_stack.push_back(mom);
                }
                if((*sire).get_name()!=all_nodes_ptr->at(1).get_name()&&(*sire).get_name()!="NA"){ // push sire to stack if sire is a known individual (not unknown)
                    indiv_stack.push_back(sire);
                } 
                indiv_stack.erase(indiv_stack.begin()); // remove first individual from stack -> next while loop with next (now first) individual
            }
        }  
        return is_ancestor;
    }catch(const std::exception &ex) {
        std::cerr << "Error in contains_ancestor(): " << ex.what() << std::endl;
        return false;
    }
}
double calculate_pure_f_xy(node*indiv_1,node*indiv_2,std::deque<dyad>*dyads_ptr,std::deque <node>*nodes_ptr,std::deque<std::deque<double>>* f_matrix,bool fill_in){ // calculates only the dyadic relatedness coefficient without path characteristics for r comparison
    try{
        int i = (*indiv_1).get_matidx();
        int j = (*indiv_2).get_matidx();
        if(i<j){ // even it is a square matrix, only the part where i<j should be filled to prevent duplicate values
            swap(i,j);
        }
        double f_xy = f_matrix->at(i)[j];// get (default) r value from f_matrix
        if(indiv_1->get_name()=="unkn_m" || indiv_1->get_name()=="unkn_f"||indiv_2->get_name()=="unkn_m" || indiv_2->get_name()=="unkn_f"){ //if one of the individuals is just a imaginary node
            f_xy = 0;
        }else if(indiv_1->get_name()==indiv_2->get_name()){ // if indiv_1 == indiv_2
            f_xy = 1;
        }else{
            if((*indiv_1).get_mom_node()->get_name()==(*indiv_2).get_name()){  // offspring|mom 
                f_xy = 0.5*(1+calculate_pure_f_xy((*indiv_1).sire_node,indiv_2,dyads_ptr,nodes_ptr,f_matrix,fill_in));
            }else if((*indiv_1).get_sire_node()->get_name()==(*indiv_2).get_name()){  // offspring|sire
                f_xy = 0.5*(1+calculate_pure_f_xy((*indiv_1).mom_node,indiv_2,dyads_ptr,nodes_ptr,f_matrix,fill_in));
            }else if((*indiv_1).get_name()==(*indiv_2).get_mom_node()->get_name()){  // mom|offspring
                f_xy = 0.5*(1+calculate_pure_f_xy(indiv_1,(*indiv_2).sire_node,dyads_ptr,nodes_ptr,f_matrix,fill_in));
            }else if((*indiv_1).get_name()==(*indiv_2).get_sire_node()->get_name()){ // sire|offspring
                f_xy = 0.5*(1+calculate_pure_f_xy(indiv_1,(*indiv_2).mom_node,dyads_ptr,nodes_ptr,f_matrix,fill_in));
            }else if(contains_ancestor(indiv_2,indiv_1,nodes_ptr)){ // indiv 1 is ancestor of indiv 2
                f_xy = 0.5*(calculate_pure_f_xy(indiv_1,(*indiv_2).sire_node,dyads_ptr,nodes_ptr,f_matrix,fill_in)
                    +calculate_pure_f_xy(indiv_1,(*indiv_2).mom_node,dyads_ptr,nodes_ptr,f_matrix,fill_in));
            }else if(contains_ancestor(indiv_1,indiv_2,nodes_ptr)){ // indiv 2 is ancestor of indiv 1
                f_xy = 0.5*(calculate_pure_f_xy((*indiv_1).sire_node,indiv_2,dyads_ptr,nodes_ptr,f_matrix,fill_in)
                    +calculate_pure_f_xy((*indiv_1).mom_node,indiv_2,dyads_ptr,nodes_ptr,f_matrix,fill_in));
            }else{ // indiv 1 and indiv 2 are not (in)direct relatives (one is the ancestor of the other), i.e. if they are related, the path needs to change direction (but only once!)
                f_xy = 0.25*(calculate_pure_f_xy((*indiv_1).mom_node,(*indiv_2).mom_node,dyads_ptr,nodes_ptr,f_matrix,fill_in)
                    +calculate_pure_f_xy((*indiv_1).mom_node,(*indiv_2).sire_node,dyads_ptr,nodes_ptr,f_matrix,fill_in)
                    +calculate_pure_f_xy((*indiv_1).sire_node,(*indiv_2).mom_node,dyads_ptr,nodes_ptr,f_matrix,fill_in)
                    +calculate_pure_f_xy((*indiv_1).sire_node,(*indiv_2).sire_node,dyads_ptr,nodes_ptr,f_matrix,fill_in));
            }
        }
        if(fill_in==true){ // fill f_matrix
            f_matrix->at(i)[j]=f_xy;
        }
        return f_xy;
    }catch(const std::exception &ex) {
        std::cerr << "Error in calculate_pure_f_xy(): " << ex.what() << std::endl;
        return 999.0;
    }
}
std::deque <node*> get_parent_pool(string sex,int mat_age_m,int mat_age_f,int gest_length, std::deque <node>*all_nodes_ptr, node* indiv,bool twins, string nonparents){ // Returns pool of potential parents for a specific individual
    std::deque<node*> potparents = {}; // to collect all potential candidates
    try{
        int mat_age = 0;
        if(sex=="f"){ // get sex specific maturation age and excluded parents
            mat_age = mat_age_f;
            nonparents = nonparents + indiv->get_nondams();
            if(twins==false){
                for(node &node: *all_nodes_ptr){ // iterate through all nodes and add all moms to nonparents who has another offspring in the birth cohort of the focal individual (no twins)
                    if(node.get_birthseason()==indiv->get_birthseason() && node.get_name()!=indiv->get_name()){
                        if(nonparents==""){
                            nonparents = node.get_mom();
                        }else{
                            nonparents = nonparents + "@" + node.get_mom();
                        }
                    }
                }
            }
        }else{
            mat_age = mat_age_m;
            nonparents = nonparents + indiv->get_nonsires();
        }
        std::deque<string> nonparents_deq = str_split(nonparents,"@");
        for(int i = 0;i<all_nodes_ptr->size();i++){ // iterate again through all nodes and rule out all parents who did not match the mom/sire criteria
            if(all_nodes_ptr->at(i).get_sex()==sex && all_nodes_ptr->at(i).get_name() != "unkn_m" && all_nodes_ptr->at(i).get_name() != "unkn_f" // individuals with the correct sex & not imaginary nodes
                &&((sex=="m" && (date_diff(all_nodes_ptr->at(i).get_DOD(),indiv->get_DOB())<gest_length||all_nodes_ptr->at(i).get_DOD().get_date()=="0-0-0")) // if male: potential sire has to be alive at the time of conception
                 ||(sex=="f" && (date_diff(all_nodes_ptr->at(i).get_DOD(),indiv->get_DOB())<0||all_nodes_ptr->at(i).get_DOD().get_date()=="0-0-0"))) // if female: potential mom has to be alive at the time of birth
                && date_diff(all_nodes_ptr->at(i).get_DOB(),indiv->get_DOB())>=(mat_age+gest_length) // if potential parent is mature enough at the time of conception
                && nonparents.find("@"+all_nodes_ptr->at(i).get_name()+"@") == string::npos && nonparents_deq[0] != all_nodes_ptr->at(i).get_name() && nonparents_deq[nonparents_deq.size()-1] != all_nodes_ptr->at(i).get_name()){ // if not excluded as parent due to STR marker (sire) or having already an child in this birthcohort (moms)
                potparents.push_back(&(all_nodes_ptr->at(i)));
            }
        }
        return potparents;
    }catch(const std::exception &ex) {
        std::cerr << "Error in get_parent_pool(): " << ex.what() << std::endl;
        return potparents;
    }
}
string all_nodes_to_population_file(std::deque <node> *all_nodes_ptr,string filename){ // save all_nodes as file
    try{
        ofstream out; 
        out.open(filename+".txt"); 
        for(int i = 0;i<all_nodes_ptr->size();i++){
            if(all_nodes_ptr->at(i).get_name()!="unkn_f"&&all_nodes_ptr->at(i).get_name()!="unkn_m"){ // skip imaginary nodes while writing
                out << all_nodes_ptr->at(i).get_name()<<"\t"<<all_nodes_ptr->at(i).get_sex()<<"\t"<<all_nodes_ptr->at(i).get_birthseason()<<"\t"<<all_nodes_ptr->at(i).get_mom("unknown")<<"\t"<<all_nodes_ptr->at(i).get_sire("unknown")<<"\t"<<all_nodes_ptr->at(i).get_DOB().get_date("NA")<<"\t"<<all_nodes_ptr->at(i).get_DOD().get_date("NA")<<"\tNA\tNA"<<endl; // can be input for next pip_forward -> instead of "_info"-file write pedigree file with nondams/nonsires instead of potsire/potdam/generation_depth
            }
        }
        out.close();
        cout<<"write "<<filename+".txt"<<endl;
        return filename;
    }catch(const std::exception &ex) {
        cerr << "Error in all_nodes_to_population_file(): " << ex.what() << std::endl;
        return "unable to save all nodes to pedigree file";
    }
}
std::deque<double> get_all_gaps(std::deque<gap>*gaps_ptr,std::deque<node>*all_nodes_ptr){ // put all gaps from the given pedigree into gaps_ptr & return some general statistics about parent pool (e.g. min/max size of parent pool)
    try{
        double ps_min = 0, ps_max = 0,n_gaps = 0, ps_sum = 0,ps_prod = 1;// n(min/max/mean/total_sum/total_product) of parent_pool_size
        for(int i = 0;i<all_nodes_ptr->size();i++){ // iterate through all nodes
            if(all_nodes_ptr->at(i).get_name()!="unkn_m" && all_nodes_ptr->at(i).get_name()!="unkn_f"){ // skip imaginary nodes
                if(all_nodes_ptr->at(i).get_mom()=="unkn_f"||all_nodes_ptr->at(i).get_mom()=="unknown"){ // if mom is unknown
                    std::deque<int> mom_pool = all_nodes_ptr->at(i).get_mom_pool(); // get potential mom candidates
                    if(mom_pool.size()>1){ // if multiple mom candidates exist
                        gap g(&(all_nodes_ptr->at(i)),{"f"},gaps_ptr->size()); // create gap object
                        gaps_ptr->push_back(g);
                        if(ps_min == 0||mom_pool.size()<ps_min){ // update minimal size of parent pool in parent stats
                            ps_min=mom_pool.size();
                        }
                        if(mom_pool.size()>ps_max){ // update maximal size of parent pool in parent stats
                            ps_max=mom_pool.size();
                        }
                        ps_sum += mom_pool.size();
                        ps_prod *= mom_pool.size();
                    }else if(mom_pool.size()==1){ // if only one mom candidate is available -> set this mom as "true" mom instead as a candidate within a gap
                        cerr << all_nodes_ptr->at(i).get_name()<< " has only one possible mom: "<<all_nodes_ptr->at(mom_pool[0]).get_name()<<endl;
                        all_nodes_ptr->at(i).set_mom(all_nodes_ptr->at(mom_pool[0]).get_name());
                        all_nodes_ptr->at(i).set_mom_node(&all_nodes_ptr->at(mom_pool[0]));
                    }else{
                        cerr << "Empty pool of potential moms for "+all_nodes_ptr->at(i).get_name()+". Further treated as founder individual/gap remains empty."<<endl;
                    }
                }
                if(all_nodes_ptr->at(i).get_sire()=="unkn_m"||all_nodes_ptr->at(i).get_sire()=="unknown" ){ // if sire is unknown
                    std::deque<int> sire_pool = all_nodes_ptr->at(i).get_sire_pool(); // get potential sire candidates
                    if(sire_pool.size()>1){ // if multiple sire candidates exist
                        gap g(&(all_nodes_ptr->at(i)),{"m"},gaps_ptr->size()); // create gap object
                        gaps_ptr->push_back(g); 
                        if(ps_min == 0||sire_pool.size()<ps_min){ // update minimal size of parent pool in parent stats
                            ps_min=sire_pool.size();
                        }
                        if(sire_pool.size()>ps_max){ // update maximal size of parent pool in parent stats
                            ps_max=sire_pool.size();
                        }
                        ps_sum += sire_pool.size();
                        ps_prod *= sire_pool.size();
                    }else if(sire_pool.size()==1){ // if only one sire candidate is available -> set this sire as "true" sire instead as a candidate within a gap
                        cerr << all_nodes_ptr->at(i).get_name()<< " has only one possible sire: "<<all_nodes_ptr->at(sire_pool[0]).get_name()<<endl;
                        all_nodes_ptr->at(i).set_sire(all_nodes_ptr->at(sire_pool[0]).get_name());
                        all_nodes_ptr->at(i).set_sire_node(&all_nodes_ptr->at(sire_pool[0]));
                    }else{
                        cerr << "Empty pool of potential sires for "+all_nodes_ptr->at(i).get_name()+". Further treated as founder individual/gap remains empty."<<endl;                   
                    }
                }
            }
        }
        std::deque<double> parent_stats = {ps_min,ps_max,ps_sum/gaps_ptr->size(),ps_sum,ps_prod};
        return parent_stats; 
    }catch(const std::exception &ex) {
        std::cerr << "Error in get_all_gaps(): " << ex.what() << std::endl;
        std::deque<double> parent_stats = {};
        return parent_stats;
    }
}
void get_start_solution(double temperature,std::deque<gap>*gaps_ptr,std::deque<node>*all_nodes_ptr,std::deque<node>*current_nodes_ptr,int maturation_age_m,int maturation_age_f,int gestation_length,bool twins){ // fills all gaps randomly in the pedigree in order to get a start solution for simulated annealing
    try{
        for(int i = 0;i<gaps_ptr->size();i++){ // iterate through all gaps
            node* offspring = &(current_nodes_ptr->at(gaps_ptr->at(i).get_indiv()->get_matidx())); // get offspring pointer from gap
            string parent_sex = gaps_ptr->at(i).get_sexseq()[0]; // determine sex of the missing parent of the gap
            std::deque <node*> parent_pool = get_parent_pool(parent_sex,maturation_age_m, maturation_age_f,gestation_length,all_nodes_ptr,&(all_nodes_ptr->at(offspring->get_matidx())),twins); // get pool of potential parents
            if(parent_pool.size()!=0){ // if no founder individual (potential parents exist) --> fill gap with randomly choosen parent from parent_pool
                int rand_parent = rand()%parent_pool.size(); // choose parent
                if(parent_sex=="m"){ // gap filling: set missing sire to newly choosen sire in the std::deque of nodes       
                    current_nodes_ptr->at(offspring->get_matidx()).set_sire(current_nodes_ptr->at(parent_pool[rand_parent]->get_matidx()).get_name());
                    current_nodes_ptr->at(offspring->get_matidx()).set_sire_node(&(current_nodes_ptr->at(parent_pool[rand_parent]->get_matidx())));
                }else if(parent_sex=="f"){// gap filling: set missing mom to newly choosen mom in the std::deque of nodes  
                    current_nodes_ptr->at(offspring->get_matidx()).set_mom(current_nodes_ptr->at(parent_pool[rand_parent]->get_matidx()).get_name());
                    current_nodes_ptr->at(offspring->get_matidx()).set_mom_node(&(current_nodes_ptr->at(parent_pool[rand_parent]->get_matidx())));
                }
                gaps_ptr->at(i).set_pot_candidate(&(current_nodes_ptr->at(parent_pool[rand_parent]->get_matidx()))); // set parent also in std::deque of gaps
                cout << " --> new indiv info: "<<current_nodes_ptr->at(offspring->get_matidx()).get_infos();   
            }else{ // no parent pool -> no parent candidates to fill the gap -> treat as founder & remove gap from gaps_ptr
                gaps_ptr->erase(gaps_ptr->begin() + i); // remove gap
                cerr << "Empty parent pool for "+offspring->get_name()+". Further treated as founder individual."<<endl;
            }
        }
    }catch(const std::exception &ex) {
        std::cerr << "Error in get_start_soluation(): " << ex.what() << std::endl;
    }
}
void get_random_full_ped(string file,std::deque<node>*all_nodes_ptr,int maturation_age_f,int maturation_age_m,int gestation_length,bool twins){ //help function to generate pedigree solution, where gaps were filled randomly with one potential parent candidate (not usable currently)
    try{
        //initialize variables
        string sex,name,mom,sire,DOB,DOD,nonsires,nondams,indiv1, indiv2,dyad_name,paths,pathline,kinline,common_anc,depth,kinlabel,fullhalf;
        int birthseason;
        double r_value;
        node imaginary_female(name="unkn_f","f","NA","NA","NA","NA","NA","NA",0,0); // create male & female imaginary nodes as dummy in case of an unknown parent
        node imaginary_male(name="unkn_m","m","NA","NA","NA","NA","NA","NA",0,1);
        all_nodes_ptr->push_back(imaginary_female);
        all_nodes_ptr->push_back(imaginary_male);
        cout<<"load gap pedigree ("<<endl;
        ifstream data(file+".txt"); // load given pedigree (with gaps)
        if (! data) {
            throw std::runtime_error("Unable to open file '"+file+".txt' for reading");
        }
        int matidx = 2;
        while(data >> name >> sex >> birthseason >> mom >> sire >> DOB >> DOD >> nonsires >> nondams){
            if(mom=="unknown"||mom=="unkn_f"||mom=="UNK"||mom=="NA"){ // if mom is unknown -> set mom name to unkn_f
                mom=="unkn_f";
            }
            if(sire=="unknown"||sire=="unkn_m"||sire=="UNK"||sire=="NA"){ // if sire in unknown -> set sire name to unkn_m
                sire=="unkn_m";
            }
            node x(name,sex,mom,sire,nonsires,nondams,DOB,DOD,birthseason,matidx); // generate for each row a node individual with the given attributes
            all_nodes_ptr->push_back(x);
            matidx += 1;
        }
        data.close(); 
        cout << matidx-2<<" individuals)\n create parent pointer and set parent pool"<<endl;
        for(int i = 0;i<all_nodes_ptr->size();i++){ // create parent pointers for each node
            all_nodes_ptr->at(i).create_parent_ptr(all_nodes_ptr);
        }
        set_parent_pool(all_nodes_ptr,maturation_age_m,maturation_age_f,gestation_length,twins);
        cout << "get all gaps"<<endl;
        std::deque<gap> gaps = {};
        std::deque<double> parent_stats = get_all_gaps(&gaps, all_nodes_ptr); // find all gaps within pedigree
        for(int i = 0;i<parent_stats.size();i++){ // print min,max,sum/gaps_ptr->size(),sum,prod of potential parents for all gaps
            cout << parent_stats[i]<<endl;
        }
        cout << "create parent pointer for current nodes"<<endl;
        std::deque<node>  current_nodes = *(all_nodes_ptr);
        for(int i = 0;i<current_nodes.size();i++){
            current_nodes[i].create_parent_ptr(&(current_nodes),true);
        }
        cout << "get start solution"<<endl;
        get_start_solution(0,&gaps,all_nodes_ptr,&current_nodes,maturation_age_m,maturation_age_f,gestation_length,twins); // fill pedigree randomly with potential parent candidates
        cout << "gaps.size(): "<<gaps.size()<<endl;
        cout << "writing "<<all_nodes_to_population_file(&current_nodes,file+"_complete");
    }catch(const std::exception &ex) {
        std::cerr << "Error in get_random_full_ped(): " << ex.what() << std::endl;
    }
}
int get_full_generation_depth(node* indiv){ // determines how many full known generations an individual has
    int full_gen = 0;
    try{
        bool complete = true;
        std::deque<node*> current_generation; 
        current_generation.push_back(indiv);// in the beginning current generation is individual itself
        while (complete == true){ // as long as current generation is completely known -> try next generation
            std::deque<node*> next_generation;
            for(int i = 0;i<current_generation.size();i++){ // for all individuals of the current generation add their parents to the next generation deque
                next_generation.push_back(current_generation[i]->get_mom_node());
                next_generation.push_back(current_generation[i]->get_sire_node());
                if(current_generation[i]->get_mom()=="unkn_f" || current_generation[i]->get_sire()=="unkn_m"){// if unknown mom/sire is in next generation, stop while loop (generation not longer complete)
                    complete = false; 
                }
            }
            swap(current_generation,next_generation); // set next_generation as current_generation and start again
            full_gen += 1; // increase full_gen counter
        }
        return full_gen;
    }catch(const std::exception &ex) {
        std::cerr << "Error in get_full_generation_depth(): " << ex.what() << std::endl;
        return full_gen;
    }
}

