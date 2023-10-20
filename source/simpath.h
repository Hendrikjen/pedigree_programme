#ifndef SIMPATH_H
#define SIMPATH_H

#include <string>
#include <deque>
#include "basic_functions.h"

using namespace std;

class simpath{
protected:
    string path_name;
    double r_val;
    string first;
    string last;
    string dyad_name;
public:
    simpath(string path_name, double r_val);
    string get_path_name();
    double get_r_val();
    string get_first();
    string get_last();
    string get_dyad_name();
    string get_info();
};


#endif