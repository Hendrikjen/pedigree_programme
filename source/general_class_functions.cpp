#include "general_class_functions.h"

void set_parent_pool(std::deque <node> *all_nodes_ptr, int maturation_age_m,int maturation_age_f,int gestation_length){
    try{
            for(int i = 2;i<all_nodes_ptr->size();i++){
            string pot_moms = "NA";
            string pot_sires = "NA";
            if(all_nodes_ptr->at(i).get_mom()=="unknown" || all_nodes_ptr->at(i).get_mom()=="NA" || all_nodes_ptr->at(i).get_mom()=="unkn_f"){
                std::deque <node*> mom_pool = get_parent_pool("f",maturation_age_m,maturation_age_f,gestation_length,all_nodes_ptr,&(all_nodes_ptr->at(i)));
                for(int j = 0;j<mom_pool.size();j++){
                    all_nodes_ptr->at(i).push_back_pot_parent(mom_pool[j]->get_sex(),mom_pool[j]->get_matidx());
                    if(j==0){
                        pot_moms = mom_pool[j]->get_name();
                    }else{
                        pot_moms = pot_moms + "@" + mom_pool[j]->get_name();
                    }
                }
            }
            if(all_nodes_ptr->at(i).get_sire()=="unknown" || all_nodes_ptr->at(i).get_sire()=="NA" || all_nodes_ptr->at(i).get_sire()=="unkn_m"){
                std::deque <node*> sire_pool = get_parent_pool("m",maturation_age_m,maturation_age_f,gestation_length,all_nodes_ptr,&(all_nodes_ptr->at(i)));
                for(int j = 0;j<sire_pool.size();j++){
                    all_nodes_ptr->at(i).push_back_pot_parent(sire_pool[j]->get_sex(),sire_pool[j]->get_matidx());
                    if(j==0){
                        pot_sires = sire_pool[j]->get_name();
                    }else{
                        pot_sires = pot_sires + "@" + sire_pool[j]->get_name();
                    }
                }
            }
            all_nodes_ptr->at(i).set_pot_moms(pot_moms);
            all_nodes_ptr->at(i).set_pot_sires(pot_sires);
            all_nodes_ptr->at(i).set_full_generations(get_full_generation_depth(&all_nodes_ptr->at(i)));
        }
    }catch(const std::exception &ex) {
        std::cerr << "Error in set_parent_pool(): " << ex.what() << std::endl;
    }
}
void get_ancestor_of_focal_indiv(std::deque<node*>*anc_nodes_ptr,node*indiv,std::map<string,std::deque<node*>>*offspring_dict){
    try{
        int further_ancs_idx = 0;
        anc_nodes_ptr->push_back(indiv); // to collect all ancestors of indiv (and indiv itself)
        while(further_ancs_idx < anc_nodes_ptr->size()){ // anc_nodes_ptr additionally used as stack, position of to-do-"indivs" = further_ancs_idx
            if(anc_nodes_ptr->at(further_ancs_idx)->get_mom()!="unkn_f"){ // if mom is known -> add current node/"indiv" to mother(key)/offspring_list(value) dictionary as well as to ancestor node std::deque (anc_nodes_ptr) (== expand stack)
                anc_nodes_ptr->push_back(anc_nodes_ptr->at(further_ancs_idx)->get_mom_node());
                if(offspring_dict->count(anc_nodes_ptr->at(further_ancs_idx)->get_mom_node()->get_name())==0){// if key not already initialized -> initialize key-value pair with empty std::deque
                    offspring_dict->insert(pair<string,std::deque<node*>>(anc_nodes_ptr->at(further_ancs_idx)->get_mom_node()->get_name(),std::deque<node*>())); 
                }
                offspring_dict->at(anc_nodes_ptr->at(further_ancs_idx)->get_mom_node()->get_name()).push_back(anc_nodes_ptr->at(further_ancs_idx)); // add current node/"indiv" from stack (offspring of the known mother)
            }
            if(anc_nodes_ptr->at(further_ancs_idx)->get_sire()!="unkn_m"){ // if sire is known -> add current node/"indiv" to sire(key)/offspring_list(value) dictionary as well as to ancestor node std::deque (anc_nodes_ptr) (== expand stack)
                anc_nodes_ptr->push_back(anc_nodes_ptr->at(further_ancs_idx)->get_sire_node());
                if(offspring_dict->count(anc_nodes_ptr->at(further_ancs_idx)->get_sire_node()->get_name())==0){
                    offspring_dict->insert(pair<string,std::deque<node*>>(anc_nodes_ptr->at(further_ancs_idx)->get_sire_node()->get_name(),std::deque<node*>()));
                }
                offspring_dict->at(anc_nodes_ptr->at(further_ancs_idx)->get_sire_node()->get_name()).push_back(anc_nodes_ptr->at(further_ancs_idx)); // add current node/"indiv" from stack (offspring of the known sire)
            }
            further_ancs_idx += 1;
        }
    }catch(const std::exception &ex) {
        std::cerr << "Error in get_ancestor_of_focal_indiv(): " << ex.what() << std::endl;
    }
}
bool contains_ancestor(node* pot_offspring,node* pot_ancestor,std::deque <node> *all_nodes_ptr){ // determines whether an individual (pot_ancestor) is an ancestor of an other individual (pot_offspring) 
    try{
        bool is_ancestor = false; // default
        std::deque <node*> indiv_stack;
        indiv_stack.push_back(pot_offspring);
        if(((*pot_offspring).get_birthseason()==0)||((*pot_ancestor).get_birthseason()==0)||
            ((*pot_ancestor).get_birthseason()<(*pot_offspring).get_birthseason())){ // preselection -> exclude ancestors born after pot_offspring
            while(is_ancestor==false&& indiv_stack.size()>0){ // if ancestry is proven once -> stop
                node* mom = (*indiv_stack[0]).mom_node; // check first individual on stack 
                node* sire = (*indiv_stack[0]).sire_node;
                if((*pot_ancestor).get_name()==(*mom).get_name()||(*pot_ancestor).get_name()==(*sire).get_name()){ // check whether mom/sire of current individual match with pot_ancestor
                    is_ancestor = true;
                }
                if((*mom).get_name()!=all_nodes_ptr->at(0).get_name()&&(*mom).get_name()!="NA"){ //push mom to stack
                    indiv_stack.push_back(mom);
                }
                if((*sire).get_name()!=all_nodes_ptr->at(1).get_name()&&(*sire).get_name()!="NA"){ // push sire to stack
                    indiv_stack.push_back(sire);
                } 
                indiv_stack.erase(indiv_stack.begin()); // remove first individual -> next while loop with next (now first) individual
            }
        }  
        return is_ancestor;
    }catch(const std::exception &ex) {
        std::cerr << "Error in contains_ancestor(): " << ex.what() << std::endl;
    }
}
double calculate_pure_f_xy(node*indiv_1,node*indiv_2,std::deque<dyad>*dyads_ptr,std::deque <node>*nodes_ptr,std::deque<std::deque<double>>* f_matrix,bool fill_in){ // calculates only the dyadic relatedness coefficient without path characteristics for r comparison
    try{
        int i = (*indiv_1).get_matidx();
        int j = (*indiv_2).get_matidx();
        if(i<j){
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
            }else{ // indiv 1 and indiv 2 are not (in)direct ancestors, i.e. if they are related, the path changes direction once
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
    }
}
std::deque <node*> get_parent_pool(string sex,int mat_age_m,int mat_age_f,int gest_length, std::deque <node>*all_nodes_ptr, node* indiv, string nonparents){ // Returns pool of potential parents for a specific individual
    try{
        std::deque<node*> filtered_nodes; // output std::deque
        int mat_age = 0;
        if(sex=="f"){
            mat_age = mat_age_f;
            nonparents = nonparents + indiv->get_nondams();
            for(node &node: *all_nodes_ptr){ // iterate through all nodes and save all moms in nonparents who has another offspring in the birth cohort of the focal individual (no twins)
                if(node.get_birthseason()==indiv->get_birthseason() && node.get_name()!=indiv->get_name()){
                    if(nonparents==""){
                        nonparents = node.get_mom();
                    }else{
                        nonparents = nonparents + "@" + node.get_mom();
                    }
                }
            }
        }else{
            mat_age = mat_age_m;
            nonparents = nonparents + indiv->get_nonsires();
        }
        //cout << "non_parents: "<<nonparents<<endl;
        for(int i = 0;i<all_nodes_ptr->size();i++){ // iterate again through all nodes and rule out all parents who did not match the mom/sire criteria
            if(all_nodes_ptr->at(i).get_sex()==sex && all_nodes_ptr->at(i).get_name() != "unkn_m" && all_nodes_ptr->at(i).get_name() != "unkn_f"
                &&((sex=="m" && (date_diff(all_nodes_ptr->at(i).get_DOD(),indiv->get_DOB())<gest_length||all_nodes_ptr->at(i).get_DOD().get_date()=="0-0-0")) // if potential sire is alive at the time of conception
                    ||(sex=="f" && (date_diff(all_nodes_ptr->at(i).get_DOD(),indiv->get_DOB())<0||all_nodes_ptr->at(i).get_DOD().get_date()=="0-0-0"))) // if potential mom is alive at the time of birth
                && date_diff(all_nodes_ptr->at(i).get_DOB(),indiv->get_DOB())>=(mat_age+gest_length) // if potential parent is mature enough at the time of conception
                && nonparents.find(all_nodes_ptr->at(i).get_name()) == string::npos){ // if not excluded as parent due to STR marker (sire) or having already an child in this birthcohort (moms)
                filtered_nodes.push_back(&(all_nodes_ptr->at(i)));
            }
        }
        return filtered_nodes;
    }catch(const std::exception &ex) {
        std::cerr << "Error in get_parent_pool(): " << ex.what() << std::endl;
    }
}
string all_nodes_to_population_file(std::deque <node> *all_nodes_ptr,string filename){ // save all_nodes as file
    try{
        ofstream out; 
        out.open(filename+".txt"); 
        for(int i = 0;i<all_nodes_ptr->size();i++){
            if(all_nodes_ptr->at(i).get_name()!="unkn_f"&&all_nodes_ptr->at(i).get_name()!="unkn_m"){ // no imaginary nodes within population file
                out << all_nodes_ptr->at(i).get_name()<<"\t"<<all_nodes_ptr->at(i).get_sex()<<"\t"<<all_nodes_ptr->at(i).get_birthseason()<<"\t"<<all_nodes_ptr->at(i).get_mom()<<"\t"<<all_nodes_ptr->at(i).get_sire()<<"\t"<<all_nodes_ptr->at(i).get_DOB().get_date()<<"\t"<<all_nodes_ptr->at(i).get_DOD().get_date()<<"\tNA\tNA"<<endl; // will be input for next pip_forward -> instead of "_info"-file write pedigree file with nondams/nonsires instead of potsire/potdam/generation_depth
            }
        }
        out.close();
        cout<<"write "<<filename+".txt"<<endl;
        return filename;
    }catch(const std::exception &ex) {
        std::cerr << "Error in all_nodes_to_population_file(): " << ex.what() << std::endl;
    }
}
std::deque<double> get_all_gaps(std::deque<gap>*gaps_ptr,std::deque<node>*all_nodes_ptr,int default_year){ // puts all gaps from the complete pedigree inside gaps_ptr (for simulated_annealing of the complete pedigree)
    try{
        double ps_min = 0, ps_max = 0,n_gaps = 0, ps_sum = 0,ps_prod = 1;// n(min/max/mean/total_sum/total_product) of parent_pool_size
        for(int i = 0;i<all_nodes_ptr->size();i++){
            if(all_nodes_ptr->at(i).get_name()!="unkn_m" && all_nodes_ptr->at(i).get_name()!="unkn_f" && ((all_nodes_ptr->at(i).get_mom()=="unkn_f" && all_nodes_ptr->at(i).get_sire()=="unkn_m" && all_nodes_ptr->at(i).get_birthseason()<=default_year) == false)){//all_nodes_ptr->at(i).get_birthseason()!=default_year){
                if(all_nodes_ptr->at(i).get_mom()=="unkn_f"||all_nodes_ptr->at(i).get_mom()=="unknown"){
                    std::deque<int> mom_pool = all_nodes_ptr->at(i).get_mom_pool();
                    if(mom_pool.size()>1){
                        gap g(&(all_nodes_ptr->at(i)),{"f"},gaps_ptr->size());
                        gaps_ptr->push_back(g);
                        if(ps_min == 0||mom_pool.size()<ps_min){ps_min=mom_pool.size();}
                        if(mom_pool.size()>ps_max){ps_max=mom_pool.size();}
                        ps_sum += mom_pool.size();
                        ps_prod *= mom_pool.size();
                    }else if(mom_pool.size()==1){
                        cout << all_nodes_ptr->at(i).get_name()<< " has only one possible mom: "<<all_nodes_ptr->at(mom_pool[0]).get_name()<<endl;
                        all_nodes_ptr->at(i).set_mom(all_nodes_ptr->at(mom_pool[0]).get_name());
                        all_nodes_ptr->at(i).set_mom_node(&all_nodes_ptr->at(mom_pool[0]));
                    }else{
                        cout << "WARNING. "<<all_nodes_ptr->at(i).get_name()<<" does not have any mom candidate"<<endl;
                    }
                }
                if(all_nodes_ptr->at(i).get_sire()=="unkn_m"||all_nodes_ptr->at(i).get_sire()=="unknown" ){
                    std::deque<int> sire_pool = all_nodes_ptr->at(i).get_sire_pool();
                    if(sire_pool.size()>1){
                        gap g(&(all_nodes_ptr->at(i)),{"m"},gaps_ptr->size());
                        gaps_ptr->push_back(g);
                        if(ps_min == 0||sire_pool.size()<ps_min){ps_min=sire_pool.size();}
                        if(sire_pool.size()>ps_max){ps_max=sire_pool.size();}
                        ps_sum += sire_pool.size();
                        ps_prod *= sire_pool.size();
                    }else if(sire_pool.size()==1){
                        cout << all_nodes_ptr->at(i).get_name()<< " has only one possible sire: "<<all_nodes_ptr->at(sire_pool[0]).get_name()<<endl;
                        all_nodes_ptr->at(i).set_sire(all_nodes_ptr->at(sire_pool[0]).get_name());
                        all_nodes_ptr->at(i).set_sire_node(&all_nodes_ptr->at(sire_pool[0]));
                    }else{
                        cout << "WARNING. "<<all_nodes_ptr->at(i).get_name()<<" does not have any sire candidate"<<endl;
                    }
                }
            }
        }
        std::deque<double> parent_stats = {ps_min,ps_max,ps_sum/gaps_ptr->size(),ps_sum,ps_prod};
        return parent_stats; 
    }catch(const std::exception &ex) {
        std::cerr << "Error in get_all_gaps(): " << ex.what() << std::endl;
    }
}
void get_start_solution_all_gaps(double temperature,std::deque<gap>*gaps_ptr,std::deque<node>*all_nodes_ptr,std::deque<node>*current_nodes_ptr,int maturation_age_m,int maturation_age_f,int gestation_length,int default_year){ // fills all gaps in the pedigree in order to get a start solution for simulated annealing
    try{
        for(int i = 0;i<gaps_ptr->size();i++){
            node* offspring = &(current_nodes_ptr->at(gaps_ptr->at(i).get_indiv()->get_matidx()));
            string parent_sex = gaps_ptr->at(i).get_sexseq()[0]; // is it mom or sire of offspring which is missing in the gap
            std::deque <node*> parent_pool = get_parent_pool(parent_sex,maturation_age_m, maturation_age_f,gestation_length,all_nodes_ptr,&(all_nodes_ptr->at(offspring->get_matidx()))); // get pool of potential parents
            if(parent_pool.size()!=0 && offspring->get_DOB().get_year()!=default_year){ // if no founder individual --> fill gap with randomly choosen parent from parent_pool
                int rand_parent = rand()%parent_pool.size(); // choose parent
                //cout << "offspring ("<<current_nodes_ptr->at(offspring->get_matidx()).get_name() << "): random parent: "<< parent_pool[rand_parent]->get_name();
                if(parent_sex=="m"){ // gap filling: set missing parents to newly choosen parent in the std::deque of nodes       
                    current_nodes_ptr->at(offspring->get_matidx()).set_sire(current_nodes_ptr->at(parent_pool[rand_parent]->get_matidx()).get_name());
                    current_nodes_ptr->at(offspring->get_matidx()).set_sire_node(&(current_nodes_ptr->at(parent_pool[rand_parent]->get_matidx())));
                }else if(parent_sex=="f"){
                    current_nodes_ptr->at(offspring->get_matidx()).set_mom(current_nodes_ptr->at(parent_pool[rand_parent]->get_matidx()).get_name());
                    current_nodes_ptr->at(offspring->get_matidx()).set_mom_node(&(current_nodes_ptr->at(parent_pool[rand_parent]->get_matidx())));
                }
                gaps_ptr->at(i).set_pot_candidate(&(current_nodes_ptr->at(parent_pool[rand_parent]->get_matidx()))); // set parent also in std::deque of gaps
                //gaps_ptr->at(i).set_was_altered(true);
                cout << " --> new indiv info: "<<current_nodes_ptr->at(offspring->get_matidx()).get_infos();   
            }else if(offspring->get_DOB().get_year()== default_year){ // population start year is reached --> offspring == founder individual)
                cout << "founder individual. No parental pool."<<endl;
            }else{
                cout << "WARNING. possible ERROR. no possible parent could be found (empty parent_pool).\n#\n#\n#\n#\n";
            }
        }
    }catch(const std::exception &ex) {
        std::cerr << "Error in get_start_soluation_all_gaps(): " << ex.what() << std::endl;
    }
}
void get_random_full_ped(string file,std::deque<node>*all_nodes_ptr,int maturation_age_f,int maturation_age_m,int gestation_length,int default_year){
    try{// LOAD PEDIGREE
        string sex,name,mom,sire,DOB,DOD,nonsires,nondams,indiv1, indiv2,dyad_name,paths,pathline,kinline,common_anc,depth,kinlabel,fullhalf;
        int birthseason;
        double r_value;
        node imaginary_female(name="unkn_f","f","NA","NA","NA","NA","NA","NA",0,0); // create male & female imaginary nodes as dummy in case of an unknown parent
        node imaginary_male(name="unkn_m","m","NA","NA","NA","NA","NA","NA",0,1);
        all_nodes_ptr->push_back(imaginary_female);
        all_nodes_ptr->push_back(imaginary_male);
        ifstream data(file+".txt"); // load pedigree with gaps
        if (! data) {
            cout << "unable to open file '"<<file <<".txt' for reading" << endl;
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
        cout<<" ... load gap pedigree ("<< matidx-2<<" individuals) ..."<<endl;
        for(int i = 0;i<all_nodes_ptr->size();i++){
            all_nodes_ptr->at(i).create_parent_ptr(all_nodes_ptr);
        }
        cout << "create parent ptr"<<endl;
        set_parent_pool(all_nodes_ptr,maturation_age_m,maturation_age_f,gestation_length);
        std::deque<gap> gaps = {};
        std::deque<double> parent_stats = get_all_gaps(&gaps, all_nodes_ptr, default_year); // find all gaps within pedigree
        for(int i = 0;i<parent_stats.size();i++){ // min,max,sum/gaps_ptr->size(),sum,prod
            cout << parent_stats[i]<<endl;
        }
        /*for(int i = 0;i<gaps.size();i++){
            cout << gaps[i].get_gap_infos()<<endl;
        }*/
        std::deque<node>  current_nodes = *(all_nodes_ptr);
        for(int i = 0;i<current_nodes.size();i++){
            current_nodes[i].create_parent_ptr(&(current_nodes),true);
        }
        get_start_solution_all_gaps(0,&gaps,all_nodes_ptr,&current_nodes,maturation_age_m,maturation_age_f,gestation_length,default_year);
        cout << "gaps.size(): "<<gaps.size()<<endl;
        cout << "writing "<<all_nodes_to_population_file(&current_nodes,file+"_complete");
    }catch(const std::exception &ex) {
        std::cerr << "Error in get_random_full_ped(): " << ex.what() << std::endl;
    }
}
int get_full_generation_depth(node* indiv){
    try{
        int full_gen = 0;
        bool complete = true;
        std::deque<node*> current_generation;
        current_generation.push_back(indiv);
        while (complete == true){
            std::deque<node*> next_generation;
            for(int i = 0;i<current_generation.size();i++){
                next_generation.push_back(current_generation[i]->get_mom_node());
                next_generation.push_back(current_generation[i]->get_sire_node());
                if(current_generation[i]->get_mom()=="unkn_f" || current_generation[i]->get_sire()=="unkn_m"){
                    complete = false;
                }
            }
            swap(current_generation,next_generation);
            full_gen += 1;
        }
        return full_gen;
    }catch(const std::exception &ex) {
        std::cerr << "Error in get_full_generation_depth(): " << ex.what() << std::endl;
    }
}

