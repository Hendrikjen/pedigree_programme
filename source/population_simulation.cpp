#include "population_simulation.h"

string generate_name(int number) { 
    try {//generate names for individuals in a simulated pedigree based on index number (max. 456976)
        string name = "0000"; // initialize name (four placeholders)
        if (number < (26 * 26 * 26 * 26)) {// if the given number is less than the maximum number of names (26^4), generate a name
            name[0] = (char)((number / (26 * 26 * 26) + 65));       // 1st letter
            name[1] = (char)((number / (26 * 26) % 26 + 65));       // 2nd letter
            name[2] = (char)((number / 26 % 26 + 65));              // 3rd letter
            name[3] = (char)((number % 26 + 65));                   // 4th letter
        } else { // if the number is too large, throw an exception.
            throw std::out_of_range("The maximum number of simulated individuals has been reached.");
        }
        return name; 
    } catch (const std::exception &ex) {
        std::cerr << "Error in generate_name(): " << ex.what() << std::endl; 
        return "unable to generate indiv name";
    }
}
std::deque <int> filter_for_pot_anc(string sex,int maturation_age,int current_year,std::deque<node>*all_nodes_ptr,int max_age){ // get pool of potential ancestors (parents) to choose from for pedigree simulation
    std::deque<int>filtered_nodes = {};
    try{ // get deque of suitable parental candidates
        datefmt current_date(current_year,1,1);
        datefmt alive(0,0,0);
        for(node &node: *all_nodes_ptr){ // check all nodes for potential candidates
            if(node.get_name()!="unkn_f" && node.get_name()!="unkn_m" // if no imaginary node
            && date_diff(node.get_DOB(),current_date)<=(365*max_age) // AND not too old (else node might die this year -> therefore not sufficient as parent)
            &&node.get_sex()==sex && date_diff(node.get_DOB(),current_date)>=maturation_age // AND old enough to be fertile
            && date_diff(node.get_DOD(),alive)==0){ // AND still alive at the current date
                filtered_nodes.push_back(node.get_matidx()); // than add node to deque with potential parental candidates
            }   
        }
    }catch(const std::exception &ex) {
        std::cerr << "Error in filter_for_pot_anc(): " << ex.what() << std::endl;
        return filtered_nodes;
    }
    return filtered_nodes;
}
void expand_paths(std::deque<simpath>parent_paths,std::string parent_sex,node* offspring,std::deque<simpath>*all_dyadic_paths_ptr,std::map<string,std::deque<simpath>>*path_dict_ptr){
    try{ // add offspring to current paths and update path_dict & all_dyadic_paths
        for(int i = 0;i<parent_paths.size();i++){ // expand parent paths from dict
            if(parent_paths[i].get_first() == offspring->get_name() or parent_paths[i].get_last() == offspring->get_name()){
                continue; // if the parent path already includes the offspring (at the end/start), continue (no extension necessary)
            }else if(parent_paths[i].get_first() == offspring->get_parent(parent_sex)){ // parent is located at the path beginning -> extend frontside
                simpath sp(offspring->get_name()+ "_" + parent_paths.at(i).get_path_name(),parent_paths.at(i).get_r_val() * 0.5); // new simpath with path_name expanded with offspring's name and updated r_value (*0.5) because one edge is added to each parent_path
                all_dyadic_paths_ptr->push_back(sp); // add new simpath to all_dyadic_paths
                path_dict_ptr->at(offspring->get_name()).push_back(sp); // add new simpath_ptr to offspring-key
                path_dict_ptr->at(parent_paths[i].get_last()).push_back(sp); // add new simpath_ptr to indiv-key which is on the other end (not-offspring-end) = anc_name of the new simpath
            }else if(parent_paths[i].get_last() == offspring->get_parent(parent_sex)){ // parent is located at the end of the path -> extend the end
                simpath sp(parent_paths.at(i).get_path_name() + "_" + offspring->get_name(),parent_paths.at(i).get_r_val() * 0.5); // new simpath with path_name expanded with offspring's name and updated r_value (*0.5) because one edge is added to each parent_path
                all_dyadic_paths_ptr->push_back(sp); // add new simpath to all_dyadic_paths
                path_dict_ptr->at(offspring->get_name()).push_back(sp); // add new simpath_ptr to offspring-key
                path_dict_ptr->at(parent_paths[i].get_first()).push_back(sp); // add new simpath_ptr to indiv-key which is on the other end (not-offspring-end) = anc_name of the new simpath
            }else{
                throw std::runtime_error("Incompatible path encountered. Unable to expand path (" + parent_paths[i].get_info() + ")");
            }
        }
    }catch(const std::exception &ex) {
        std::cerr << "Error in expand_path(): " << ex.what() << std::endl;
    }
}
void write_simulated_r_values_by_dyads(int simulation_duration,int start_individuals,std::deque<simpath>*all_dyadic_paths_ptr,string output){
    try{ //writes simulated dyadic relatedness coefficients to an output file
        std::map<std::string,double> dyadic_r; //map to store dyad names and their corresponding r values
        for(int i = 0;i<all_dyadic_paths_ptr->size();i++){
            if(dyadic_r.count(all_dyadic_paths_ptr->at(i).get_dyad_name())==0){ // if there is not already an entry for this dyad name, add new entry with respective r value
                dyadic_r[all_dyadic_paths_ptr->at(i).get_dyad_name()] = all_dyadic_paths_ptr->at(i).get_r_val();
            }else{
                double old_r = dyadic_r[all_dyadic_paths_ptr->at(i).get_dyad_name()]; // else get the current r value 
                dyadic_r[all_dyadic_paths_ptr->at(i).get_dyad_name()] = old_r + all_dyadic_paths_ptr->at(i).get_r_val();// AND add the new r value to the existing value
            }
        }
        ofstream out; 
        out.open(output+"_dyadic_paths.txt");  // create output file and write map information (dyadname + r) to the file
        std::map<string,double>::iterator it = dyadic_r.begin();
        for(it;it!=dyadic_r.end();++it){
            out << it->first<<"\t"<<to_string_with_precision(it->second,15)<<endl;
        }
        out.close();
        cout << "write "<<output<<"_dyadic_paths.txt"<<endl;
    }catch(const std::exception &ex) {
        std::cerr << "Error in write_simulated_r_values_by_dyads(): " << ex.what() << std::endl;
    }
}
string population_simulation(int simulation_duration, int start_individuals,int gestation_length,int maturation_age_f,int maturation_age_m,int max_age, int default_year,string output,double birth_rate,double death_rate){ // Pipeline to simulate pedigree with given parameters (how many start_individuals and duration (years))
    try{ // simulate pedigree based on a given number of founder individuals for a given number of years
        cout << "initialize basic objects"<<endl;
        node imaginary_female("unkn_f","f","NA","NA","NA","NA","NA","NA",0,0);
        node imaginary_male("unkn_m","m","NA","NA","NA","NA","NA","NA",0,1);
        std::deque <node> all_nodes;
        std::deque <simpath> all_dyadic_paths;
        std::map <string,std::deque<simpath>> path_dict;
        all_nodes.push_back(imaginary_female);
        all_nodes.push_back(imaginary_male);
        std::deque <string> all_node_names;
        all_node_names.push_back(imaginary_female.get_name());
        all_node_names.push_back(imaginary_male.get_name());
        cout << "create founder individuals"<<endl;
        for(int i = 0;i<start_individuals;i++){ 
            // simulate start individuals with randomly choosen sex and DOB
            int rand_sex = rand() % 2;
            string sex = "m";
            if(rand_sex==0){
                sex="f";
            }
            int rand_day = (rand()%28)+1;
            int rand_month = (rand()%12)+1;
            node n = node(generate_name(i),sex,"unkn_f","unkn_m","NA","NA",to_string(default_year)+"-"+to_string(rand_month)+"-"+to_string(rand_day),"NA",default_year,i+2); // add as node
            simpath sp = simpath(n.get_name(),1); // initialize the first simpath (path with only the individual itself) for tracking the relatedness while simulation
            all_nodes.push_back(n);
            all_node_names.push_back(n.get_name());
            path_dict.insert(pair<string,std::deque<simpath>>(n.get_name(),std::deque<simpath>()));
            path_dict[n.get_name()].push_back(sp);        
        }
        cout << "start simulating descendants"<<endl;// after a time leap to ensure all start_individuals are mature & fertile
        int name_i = start_individuals;
        int time_leap = (max(maturation_age_f,maturation_age_m) + 2*365)/365;
        int current_year = default_year+time_leap;
        std::deque<int> start_females = filter_for_pot_anc("f",maturation_age_f+gestation_length,current_year,&all_nodes,max_age-5); 
        int no_births = (int)(start_females.size()*0.5); // set number of birth to 50% of all potential moms
        int no_deaths = (int) no_births-1; // adapt number of deaths
        while(current_year<(default_year+time_leap+simulation_duration)){ // simulate births and deaths for each year (as long as set by simulation_duration)
            cout << "simulate year: "<<current_year<<"\t";
            std::deque<int> males = filter_for_pot_anc("m",maturation_age_m,current_year,&all_nodes,max_age-5); // get all mature and fertile females who are able to be a parent for this cohort
            std::deque<int> females = filter_for_pot_anc("f",maturation_age_f+gestation_length,current_year,&all_nodes,max_age-5); // get all mature and fertile males who are able to be a parent for this cohort
            std::deque<int> death_pool; // pool of all individuals who can not die in the current year because they had born/sired an offspring
            if(no_births>=females.size()){ // ênsure that there are always less births than potential moms
                no_births = females.size()-1;
                no_deaths = no_births-1;
            }
            cout << no_births<< " births ... ";
            for(int i = 0;i<=no_births;i++){ // simulate offsprings
                int mom = rand() % females.size(); // choose random mom from pool
                int sire = rand() % males.size(); // choose random sire from pool
                int rand_sex = rand() % 2; // choose random sex for offspring
                string sex = "m";
                if(rand_sex==0){
                    sex="f";
                }
                int rand_day = (rand()%28)+1; // choose random DOB
                int rand_month = (rand()%12)+1;
                node offspring(generate_name(name_i),sex,all_nodes[females[mom]].get_name(),all_nodes[males[sire]].get_name(),"NA","NA",to_string(current_year)+"-"+to_string(rand_month)+"-"+to_string(rand_day),"NA",current_year,name_i+2); // add new offspring as node
                all_nodes.push_back(offspring);
                all_node_names.push_back(offspring.get_name());
                simpath sp(offspring.get_name(),1); // init offsprings own simpath 
                path_dict.insert(pair<string,std::deque<simpath>>(offspring.get_name(),std::deque<simpath>())); // create dictionary key for offspring
                path_dict[offspring.get_name()].push_back(sp); // and add its path to path_dict
                expand_paths(path_dict[offspring.get_mom()],"f",&offspring,&all_dyadic_paths,&path_dict); // expand mom_paths
                expand_paths(path_dict[offspring.get_sire()],"m",&offspring,&all_dyadic_paths,&path_dict); // expand sire_paths
                death_pool.push_back(all_nodes[females[mom]].get_matidx()); // put mom to death_pool -> to prevent her from dying this current year
                death_pool.push_back(all_nodes[males[sire]].get_matidx()); // put sire to death_pool -> to prevent him from dying this current year
                females.erase(females.begin()+mom); // make sure mom cannot be choosen a second time in the current year (no twins)
                name_i += 1;
            }
            cout << no_deaths<<" deaths"<<endl;
            // check if there are too old individuals (> max_age) & set random DOD
            int happened_deaths = 0;
            datefmt current_date(current_year,1,1);
            for(int i = 2;i<all_nodes.size();i++){
                if(date_diff(all_nodes[i].get_DOB(),current_date)>(365*max_age)){ 
                    happened_deaths = happened_deaths + 1;
                    if(all_nodes[i].get_DOB().get_month()==12){
                        all_nodes[i].set_DOD(current_date.get_date());
                    }else{
                        int rand_day = (rand()%28)+1;
                        int rand_month = (rand()%(12-all_nodes[i].get_DOB().get_month()))+(1+all_nodes[i].get_DOB().get_month());
                        all_nodes[i].set_DOD(to_string(current_year)+"-"+to_string(rand_month)+"-"+to_string(rand_day));
                    }
                }
            }
            while(happened_deaths<=no_deaths){ // choose random indiv and check whether its allowed to die or not (if yes, set random DOD, else choose another one, as long as the specified number of death in this current year is reached)
                int dying_indiv = (rand() % (all_nodes.size()-2))+2; 
                bool allowed_to_die = true;
                for(int j = 0;allowed_to_die &&j<death_pool.size();j++){ // if an individual is stacked in death_pool - means it is a parent this current year - its not allowed to die (so that it is not possible that a parent dies before the offspring is born (due to random DOB and DOD))
                    if(death_pool[j]==dying_indiv){
                        allowed_to_die = false;
                    }
                }
                if(allowed_to_die && all_nodes[dying_indiv].get_DOD().get_date() == "0-0-0"){ // only if they are no current parent and if they have not deceased earlier
                    if(all_nodes[dying_indiv].get_DOB().get_month()==12 && all_nodes[dying_indiv].get_DOB().get_year()== current_year){ // to simplify and to make sure an infant did not die before it might be born -> DOB == DOD
                        all_nodes[dying_indiv].set_DOD(all_nodes[dying_indiv].get_DOB().get_date());
                    }else{ // choose random DOD (with some limitations so that they are not dying before born)
                        int rand_month;
                        int rand_day = (rand()%28)+1;
                        if(all_nodes[dying_indiv].get_DOB().get_month()==12 && all_nodes[dying_indiv].get_DOB().get_year()!=current_year){
                            rand_month = (rand()%12)+1;
                        }else{
                            rand_month = (rand()%(12-all_nodes[dying_indiv].get_DOB().get_month()))+(1+all_nodes[dying_indiv].get_DOB().get_month());
                        }
                        all_nodes[dying_indiv].set_DOD(to_string(current_year)+"-"+to_string(rand_month)+"-"+to_string(rand_day));
                    }
                    death_pool.push_back(all_nodes[dying_indiv].get_matidx()); // put dead individual to death pool so that it doesnt die twice
                    happened_deaths += 1;
                }
            }
            no_births += birth_rate;  // rhesus: 4
            no_deaths += death_rate; // rhesus: 3
            current_year += 1;
        }
        write_simulated_r_values_by_dyads(simulation_duration,start_individuals,&all_dyadic_paths,output); // save relatedness data
        return all_nodes_to_population_file(&all_nodes,output); // save simulated population
    }catch(const std::exception &ex) {
        std::cerr << "Error in population_simulation(): " << ex.what() << std::endl;
        return "unable to simulate pedigree";
    }
}
string add_parental_gaps(string file,double mat_gaps, double pat_gaps, int start_individuals){ // load completely known pedigree and insert artificial gaps (depending on whether paternal or maternal, a different number of artificial gaps is generated)
    try{ // add artificial gaps within a pedigree
        string sex,name,mom,sire,DOB,DOD,nonsires,nondams;
        int birthseason;
        std::deque <node> all_nodes;
        ifstream data(file+".txt"); // load pedigree & save it in all_nodes
        if (! data) {
            cout << "unable to open file for reading" << endl;
        }
        int matidx = 0;
        while(data >> name >> sex >> birthseason >> mom >> sire >> DOB >> DOD >> nonsires >> nondams){
            if(mom=="unknown"||mom=="unkn_f"||mom=="UNK"||mom=="NA"){
                mom=="unkn_f";
            }
            if(sire=="unknown"||sire=="unkn_m"||sire=="UNK"||sire=="NA"){
                sire=="unkn_m";
            }
            node x(name,sex,mom,sire,nonsires,nondams,DOB,DOD,birthseason,matidx);
            all_nodes.push_back(x);
            matidx += 1;
        }
        data.close(); 
        mat_gaps = (int)(mat_gaps*(all_nodes.size()-start_individuals)); // determine amount of artificial gaps depending on total number of descendants 
        pat_gaps = (int)(pat_gaps*(all_nodes.size()-start_individuals));
        cout << "missing_moms: "<<mat_gaps<<endl;
        cout << "missing_sires: "<<pat_gaps<<endl;
        while(mat_gaps>0){ // add maternal gaps
            int idx = (rand() % (matidx+1-start_individuals))+start_individuals; // choose random individual (not start individuals) to add unknown mom
            if(all_nodes[idx].get_mom()!="unkn_f"||all_nodes[idx].get_DOD().get_date()=="1900-1-1"){
                all_nodes[idx].set_mom("unkn_f");
                mat_gaps--;
            }
        }
        while(pat_gaps>0){ // add paternal gaps
            int idx = (rand() % (matidx+1-start_individuals))+start_individuals; // choose random individual (not start individuals) to add unknown sire
            if(all_nodes[idx].get_sire()!="unkn_m"||all_nodes[idx].get_DOD().get_date()=="1900-1-1"){
                all_nodes[idx].set_sire("unkn_m");
                pat_gaps--;
            }
            
        }
        return all_nodes_to_population_file(&all_nodes,file+"_gaps"); // save population pedigree with gaps
    }catch(const std::exception &ex) {
        std::cerr << "Error in add_parental_gaps(): " << ex.what() << std::endl;
        return "unable to add parental gaps to pedigree";
    }
}