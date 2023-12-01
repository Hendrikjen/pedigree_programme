#include <chrono>
#include <time.h>
#include <cstdio>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <deque>
#include <string>
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

void display_help(){ // output of -h (command line argument help)
    cout << "\nCommand Line Arguments:\nPlease make sure to select a functionality (relatedness, simulation, or annealing) \nand to specify all required arguments for the respective function --> see info in \nsquare brackets for each argument: [r|-|o]\n  o = optional\n  r = required\n  - = not necessary/used)"<<endl;
    cout << "Format: \n  -argument_identifier <argument> [data type][relatedness|simulation|annealing]\n  + further explanation (options & default values)"<<endl;
    cout << "Example (required arguments): \n  ./pedigree_programme -f relatedness -p pedigree.txt \n  ./pedigree programme -f simulation -n 20 -s 10\n  ./pedigree programme -f annealing -p pedigree.txt -d dyads.txt\n\n-----------------------------------------------------------------------------\n"<<endl;
    cout << "-a <max_age> [int][-|o|-]\n   options: age maximum in population (individuals who reach the maximum age \n            will decease in the following year)\n   default: 200"<<endl;
    cout << "-b <birth_rate> [double][-|o|-]\n   options: specifies the annual increment in the number of offsprings born \n            each year during the population simulation\n   default: 4.0"<<endl;
    cout << "-c <cores> [int][o|-|o]\n   options: number of cores for multiprocessing\n   default: 1 (no multiprocessing)"<<endl;
    cout << "-d <input_dyadlist> [string][o|-|r]\n   options: path to file with selected dyads (relatedness calculation) or \n            path to dyadlist with realized relatedness values (simulated \n            annealing)\n   default: [empty] (all dyads within the pedigree will be analysed)"<<endl;
    cout << "-e <output_extend> [string][o|-|-]\n   options:\n     - full: returns the full dyadlist output, including path characteristics\n     - reduced: returns only dyadlist with dyadic relatedness coefficients\n   default: full"<<endl;
    cout << "-f <functionality> [string][r|r|r]\n   options: \n     - relatedness: calculates the dyadic relatedness (+ path \n       characteristics) from a given pedigree\n     - simulation: simulates a pedigree\n     - annealing: starts a simulated annealing algorithm to fill the \n       parental gaps within a pedigree based on realized relatedness values\n   default: [empty] (the programme starts without task)"<<endl;
    cout << "-g <gestation_length> [int][o|o|o]\n   options: gestation length in days\n   default: 200"<<endl;
    cout << "-h (without argument): display help"<<endl;
    cout << "-i <init_temp> [double][-|-|o]\n   options: start temperature \n   default: [empty] (automatically calculated)"<<endl;
    cout << "-j <twins> [bool][o|o|o]\n   options \n     - T (true): twins are possible \n     - F (false): twins are not possible or that rare that potential mom \n       candidates can be excluded if the have already an offspring in the \n       respective birth cohort"<<endl;
    cout << "-k <visualization> [bool][-|-|o]\n   options \n     - T (true): keep track of simulated annealing steps (the respective \n       relatedness variance and if they are rejected) \n     - F (false): prior simulated annealing steps are not recorded/returned"<<endl;
    cout << "-l <generation_limit> [int][o|-|-]\n   options: restricts the distance to potential lowest common ancestors, \n            e.g. if generationlimit == 3, only paths up to the grandparent \n            generation will be returned, great-grand-parents will be considered as \n            unrelated\n   default: [empty] (no limitation; all ancestors of a focal will be \n            considered as potential lowest common ancestor)"<<endl;
    cout << "-m <maturation_age_m> [int][o|o|o]\n   options: maturation age of males in days\n   default: 1250"<<endl;
    cout << "-n <start_individual> [int][-|r|-]: number of individuals at the start of the \n   simulation"<<endl;
    cout << "-o <output> [string][o|o|o]\n   options: custom output name (prefix) e.g. if output == programmeoutput, the \n            resulting output files will be named 'programmeoutputdyadlist.txt' \n            and 'programmeoutputinfo.txt'\n   default: [empty] (the input file name will be used as prefix)"<<endl;
    cout << "-p <input_pedigree> [string][r|-|r]: path to pedigree file, e.g. pedigree.txt"<<endl;
    cout << "-q <death_rate> [double][-|o|-]\n   options: specifies the annual increment in the number of deaths each year \n            during the population simulation\n   default: 3.0"<<endl;
    cout << "-r <reduce_node_space> [bool][o|-|-]\n   options: \n     - T (true): before calculating the dyadic relatedness, the number of \n       individuals will be reduced which means that only descendants of \n       the focal's common ancestors will be considered in the analysis\n       (it effectively reduces the search space without affecting the \n       result, but might be only beneficial in almost completely reconstructed \n       pedigrees with a long history due to the extra computational cost)\n     - F (false): no prior narrowing of the search space\n   default: false"<<endl;
    cout << "-s <simulation_duration> [int][-|r|-]: number of years to restrict the duration \n   of the simulation"<<endl;
    cout << "-t <stop_temp> [double][-|-|o]\n   options: stop temperature, if current temperature falls below stop temperature, \n            the algorithm terminates\n   default: 1.0"<<endl;
    cout << "-v (without argument): version information"<<endl;
    cout << "-w <maturation_age_f> [int][o|o|o]\n   options: maturation age of females in days\n   default: 1095"<<endl;
    cout << "-x <temp_decay> [double][-|-|o]\n   options: the temperature multiplication factor to determine the number of \n            iterations \n   default: 0.99"<<endl;
    cout << "-y <default_year> [int][-|o|-]\n   options: start year for population simulation\n   default: 1900\n";
    cout << "-z <complete_pedigree> [string][-|-|o]\n   options: path to complete pedigree if full known pedigree exists (with \n            all gaps correctly filled) and if it should be used to evaluate accuracy \n            of simulated annealing output\n   default: [empty] (no comparison whether gaps are correctly filled in \n            the end of the simulated annealing)"<<endl;
}

int main(int argc, char *argv[]) {
    // initialization and setting of default values
    auto start = high_resolution_clock::now();
    std::chrono::steady_clock::time_point init_main_start = std::chrono::steady_clock::now();
    int option;
    string functionality, input_pedigree, input_dyadlist,output,complete_pedigree;
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
    bool twins = false;
    bool visualization = true;

    if (argc == 2 && std::string(argv[1]) == "-h") {
        display_help();
        return 0;
    } else if (argc == 2 && std::string(argv[1]) == "-v") {
        std::cout << "Version 1.0.0" << std::endl;
        return 0;
    }
    while ((option = getopt(argc, argv, "a:b:c:d:e:f:g:i:j:k:l:m:n:o:p:q:r:s:t:w:x:y:z:")) != -1) {
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
            case 'j':
                if(optarg == "T"){
                    twins = true;
                }else if(optarg == "F"){
                    twins = false;
                }else{
                    cerr << "Invalid twins argument. Please choose 'T' (true) or 'F' (false); default = false"<<endl;
                }
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
            case 'k':
                if(optarg == "T"){
                    visualization = true;
                }else if(optarg == "F"){
                    visualization = false;
                }else{
                    cerr << "Invalid visualization argument. Please choose 'T' (true) or 'F' (false); default = true"<<endl;
                }
                break;
            case 'r':
                if(optarg == "T"){
                    reduce_node_space = true;
                }else if(optarg == "F"){
                    reduce_node_space = false;
                }else{
                    cerr << "Invalid reduce_node_space argument. Please choose 'T' (true) or 'F' (false); default = false"<<endl;
                }
                break;
            case 'z':
                complete_pedigree = optarg;
                break;
            default:
                std::cerr << "Invalid argument. Please check the documentation for available options, clarification and the required data types." << endl;
                return 1;
        }
    }

    if(functionality == "relatedness") { // start relatedness calculation (check if all required arguments are given)
        if (input_pedigree.empty()) {
            cerr << "Missing argument for relatedness calculation. Please add '-p [path to pedigree file.txt]'" << endl;
            return 1;
        }else if(input_pedigree.length()>=4 && input_pedigree.substr(input_pedigree.length() - 4) == ".txt"){
            cout << "start relatedness calculation"<<endl;
            if(output.empty()){
                output = str_split(input_pedigree,".txt")[0];
            }
            if(input_dyadlist.empty()){
                input_dyadlist = "NA";
            }
            relatedness_calculation(input_pedigree,output,input_dyadlist,maturation_age_f,maturation_age_m,gestation_length,twins,output_extend,generation_limit,cores,reduce_node_space);
        }else{
            cerr << "please use a .txt file as input pedigree and dyad selection"<<endl;
            return 1;
        }
    }else if(functionality == "simulation"){ // start population simulation (check if all required arguments are given)
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
    }else if(functionality == "annealing"){ // start simulated annealing (check if all required arguments are given)
        if(input_pedigree.empty() || input_dyadlist.empty()){
            cerr << "Missing argument for simulated annealing. Please make sure '-p [input_pedigree]' and  '-d [dyadlist with realized relatedness values]' are called correctly." << endl;
            return 1;
        }else if((input_pedigree.length()>=4 && input_pedigree.substr(input_pedigree.length() - 4) == ".txt")==false){
            cerr << "please use a .txt file as input pedigree"<<endl;
            return 1;
        }else if((input_dyadlist.length()>=4 && input_dyadlist.substr(input_dyadlist.length() - 4) == ".txt")==false){
            cerr << "please use a .txt file as input dyad data"<<endl;
            return 1;
        }else{
            if(output.empty()){
                output = str_split(input_pedigree,".txt")[0];
            }
            if(complete_pedigree.empty()){
                complete_pedigree = "NA";
            }
            cout << "start simulated annealing"<<endl;
            simulated_annealing(input_pedigree,input_dyadlist,output,init_temp,stop_temp,temp_decay,cores,maturation_age_f,maturation_age_m,gestation_length,twins,visualization,complete_pedigree);
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
