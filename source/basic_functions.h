#ifndef BASIC_FUNCTIONS_H
#define BASIC_FUNCTIONS_H

#include <string>
#include <deque>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>

#include "datefmt.h"

using namespace std;

std::string to_string_with_precision(double number, int n);
std::string to_string_with_precision(float number, int n);
std::string to_string_with_precision(int number, int n);
std::deque <string> str_split(string str, string pattern);
bool is_leapyear(int year);
int date_diff(datefmt date_1, datefmt date_2);
#endif