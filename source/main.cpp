//#include <sstream>
//#include <algorithm>
//#include <tuple> 
#include <chrono>
#include <time.h>
#include <cstdio>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <deque>
#include <string>
#include <filesystem>
#include <cmath>
#include <thread>
#include <getopt.h>

#include "datefmt.h"
#include "basic_functions.h"
#include "node.h"
#include "dyad.h"
#include "gap.h"
#include "general_class_functions.h"
#include "simpath.h"
#include "relatedness_calculation.h"
#include "population_simulation.h"
#include "simulated_annealing.h"

using namespace std::chrono; 
namespace fs = std::filesystem;


int main(int argc, char *argv[]) {
    auto start = high_resolution_clock::now();
    std::chrono::steady_clock::time_point init_main_start = std::chrono::steady_clock::now();
    int option;
    string functionality, input_pedigree, input_dyadlist,output;
    string output_extend = "full";
    int gestation_length = 200;
    int maturation_age_f = 1095;
    int maturation_age_m = 1250;
    int generation_limit = -1;
    int max_age = 30;
    int default_year = 1900; 
    int simulation_duration = -1;
    int start_individuals = -1;
    double birth_rate = 4.0;
    double death_rate = 3.0;
    int cores = 1;
    double init_temp = 0; 
    double stop_temp = 1; 
    double temp_decay = 0.99; 
    bool reduce_node_space = false;

    while ((option = getopt(argc, argv, "a:b:c:d:e:f:g:i:l:m:n:o:p:q:r:s:t:w:x:y:z:")) != -1) {
        switch (option) {
            case 'f':
                functionality = optarg;
                break;
            case 'p':
                input_pedigree = optarg;
                break;
            case 'd':
                input_dyadlist = optarg;
                break;
            case 'e':
                output_extend = optarg;
                break;
            case 'o':
                output = optarg;
                break;
            case 'g':
                gestation_length = stoi(optarg);
                break;
            case 'm':
                maturation_age_m = stoi(optarg); 
                break;
            case 'w':
                maturation_age_f = stoi(optarg); 
                break;
            case 'l':
                generation_limit = stoi(optarg);
                break;
            case 'y':
                default_year = stoi(optarg);
                break;
            case 's':
                simulation_duration = stoi(optarg); 
                break;
            case 'n':
                start_individuals = stoi(optarg);
                break;
            case 'b':
                birth_rate = stod(optarg);
                break;
            case 'q':
                death_rate = stod(optarg);
                break;
            case 'a':
                max_age = stoi(optarg);
                break;
            case 'c':
                cores = stoi(optarg);
                break;
            case 'i':
                init_temp = stod(optarg); 
                break;
            case 't':
                stop_temp = stod(optarg); 
                break;
            case 'x':
                temp_decay = stod(optarg);
                break;
            case 'r':
                if(optarg == "true"){
                    reduce_node_space = true;
                }else if(optarg == "false"){
                    reduce_node_space = false;
                }else{
                    cerr << "Invalid reduce_node_space argument. Please choose 'true' or 'false' (default = false)"<<endl;
                }
                break;
            default:
                std::cerr << "Invalid argument. Please check the documentation for available options, clarification and the required data types." << endl;
                return 1;
        }
    }

    if(functionality == "relatedness") {
        if (input_pedigree.empty()) {
            cerr << "Missing argument for relatedness calculation. Please add '-p [path to pedigree file.txt]'" << endl;
            return 1;
        }else{
            cout << "start relatedness calculation"<<endl;
            if(output.empty()){
                output = str_split(input_pedigree,".txt")[0];
            }
            if(input_dyadlist.empty()){
                input_dyadlist = "";
            }
            pip_forward(input_pedigree,output,input_dyadlist,maturation_age_f,maturation_age_m,gestation_length,output_extend,generation_limit,cores,reduce_node_space);
        }
    }else if(functionality == "simulation"){
        if(simulation_duration < 0 || start_individuals < 0){
            cerr << "Missing argument for population simulation. Please make sure '-s [simulation_duration]' and  '-n [number of start individuals]' are called correctly." << endl;
            return 1;
        }else{
            cout << "start population simulation"<<endl;
            if(output.empty()){
                output = "simulated_population_"+to_string(start_individuals)+"n_"+to_string(simulation_duration)+"y";
            }
            population_simulation(simulation_duration,start_individuals,gestation_length,maturation_age_f,maturation_age_m,max_age,default_year,output,birth_rate,death_rate);
        }
    }else if(functionality == "annealing"){
        if(input_pedigree.empty() || input_dyadlist.empty()){
            cerr << "Missing argument for simulated annealing. Please make sure '-p [input_pedigree]' and  '-d [dyadlist with realized relatedness values]' are called correctly." << endl;
            return 1;
        }else{
            cout << "start simulated annealing"<<endl;
        }
    }else{
        cerr << "No assigned task. Please make sure '-f [relatedness|simulation|annealing]' is called correctly." << endl;
        return 1;
    }

    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<milliseconds>(stop - start);
    cout << "\n\n####################### finished after ... #######################\n--> "<< duration_cast<microseconds>(stop - start).count()<< " microseconds\n--> "<< duration_cast<milliseconds>(stop - start).count()<< " milliseconds\n--> "<< duration_cast<seconds>(stop - start).count()<< " seconds\n--> "<< duration_cast<minutes>(stop - start).count()<< " minutes\n";
    return 0;
}

    //string working_directory = fs::current_path().generic_string() + "/";

    // PART ONE: relatedness coefficient calculation based on rhesus macaque pedigree (Cayo Santiago population) >> real pedigree with gaps
    //pip_forward(working_directory+"pedigree_since_1986_plus_parents",maturation_age_f,maturation_age_m,gestation_length,"full",generation_limit,true,n_cores);
    //pip_forward(working_directory+"pedigree_since_1986_plus_parents_complete_unk",maturation_age_f,maturation_age_m,gestation_length,"full",generation_limit,true,n_cores);
    //pip_forward(working_directory+"mini_example_git",maturation_age_f,maturation_age_m,gestation_length,"full",generation_limit);
    
    // PART TWO: simulated annealing
    /*string pedigree_full = working_directory+"pedigree_since_1986_plus_parents_complete_unk";//"popfile_simulation_3y_5s";//
    string pedigree_gaps = working_directory+"pedigree_since_1986_plus_parents";//"popfile_simulation_3y_5s_gaps";//
    string dyadlist = working_directory+"pedigree_since_1986_plus_parents_all_dyads";//"popfile_simulation_3y_5s_all_dyads";//

    std::deque <dyad> all_dyads={};
    std::deque<int> subset_idx = {};
    std::deque <node> all_nodes = {};
    std::map<string, int> dyad_dict = {};
    auto init_main = std::chrono::duration_cast<std::chrono::nanoseconds> (std::chrono::steady_clock::now() - init_main_start).count();
    cout << "init_main = " << init_main << "[ns]" << std::endl;
    std::chrono::steady_clock::time_point load_data_start = std::chrono::steady_clock::now();
    load_data_for_sim_annealing(pedigree_full,pedigree_gaps,&all_dyads,&subset_idx,&all_nodes,&dyad_dict,maturation_age_f,maturation_age_m,gestation_length,dyadlist);
    set_parent_pool(&all_nodes,maturation_age_m,maturation_age_f,gestation_length);
    auto load_data_main = std::chrono::duration_cast<std::chrono::nanoseconds> (std::chrono::steady_clock::now() - load_data_start).count();
    cout << "load_data_main = " << load_data_main << "[ns]" << std::endl;
    cout << "dyad size: "<<all_dyads.size()<<", subset: "<<subset_idx.size()<<endl;
    double init_temp = 0; // start temperature
    double stop_temp = 1; // stopping temperature
    double temp_decay = 0.9999992916945; // 0.9995factor to reduce temperature (0.945 -> 123 runs per annealing; )
    
    simulated_annealing_complete_pedigree(init_temp,stop_temp,temp_decay,&all_nodes,&all_dyads,&subset_idx,gestation_length,maturation_age_f,maturation_age_m,default_year,pedigree_full,true,true,n_cores);

    cout << "init_main = " << init_main << "[ns]" << std::endl;
    cout << "load_data_main = " << load_data_main << "[ns]" << std::endl;
    
    cout << "\nn_nodes: "<<all_nodes.size()<<endl;
    cout << "n_dyads: "<<all_dyads.size()<<endl;
    cout << "n_relevant_dyads: "<<subset_idx.size()<<endl;
    */



