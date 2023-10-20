#include "basic_functions.h"

std::string to_string_with_precision(double number, int n) {
    std::ostringstream out;
    out.precision(n);
    out << std::fixed << number;
    return out.str();
}
std::string to_string_with_precision(float number, int n) {
    std::ostringstream out;
    out.precision(n);
    out << std::fixed << number;
    return out.str();
}
std::string to_string_with_precision(int number, int n) {
    std::ostringstream out;
    out.precision(n);
    out << std::fixed << number;
    return out.str();
}
std::deque <string> str_split(string str, string pattern){
    std::deque <string> output;
    int i = 0;
    int j = str.find(pattern);
    while (j != -1) {
        output.push_back(str.substr(i, j-i));
        i = j + pattern.size();
        j = str.find(pattern, i);
    }
    output.push_back(str.substr(i, j-i));
    return output;
}
bool is_leapyear(int year){ // return true if year is a leap year, else false
    bool is_leapyear = false;
    if(year%4==0){
        is_leapyear = true;
        if(year%100==0){
            is_leapyear = false;
            if(year%400==0){
                is_leapyear = true;
            }
        }
    }
    return is_leapyear;
}
int date_diff(datefmt date_1, datefmt date_2){ // determines difference between two dates in days
    int difference = 0;
    bool swapping = false;
    if(date_1.get_year()>=date_2.get_year()){ // swap so that date_1 < date_2
        swap(date_1,date_2);
        swapping = true;
        if((date_1.get_year()==date_2.get_year()) && (date_1.get_month()>date_2.get_month())){
            swap(date_1,date_2);
            swapping = false;
            if((date_1.get_month()==date_2.get_month())&&(date_1.get_day()>date_2.get_day())){
                swap(date_1,date_2);
                swapping = true;
            }
        }
    }
    for(int i = date_1.get_year();i<date_2.get_year();i++){ // sum up all days by date discrepancy in the time frame of years
        if(is_leapyear(i)){
            difference += 366;
        }else{
            difference += 365;
        }
    }
    for(int i = 1;i<date_1.get_month();i++){ // substract within year discrepancy (1.1. .. date_1)
        if(i==2){
            if(is_leapyear(date_1.get_year())){
                difference -= 29;
            }else{
                difference -= 28;
            }
        }else if(i==1||i==3||i==5||i==7||i==8||i==10||i==12){
            difference -= 31;
        }else{
            difference -= 30;
        }
    }
    difference -= date_1.get_day(); // subtract day discrepancy
    for(int i = 1;i<date_2.get_month();i++){ // add within year discrepancy (1.1 .. date_2)
        if(i==2){
            if(is_leapyear(date_2.get_year())){
                difference += 29;
            }else{
                difference += 28;
            }
        }else if(i==1||i==3||i==5||i==7||i==8||i==10||i==12){
            difference += 31;
        }else{
            difference += 30;
        }
    }
    difference += date_2.get_day(); // add day discrepancy
    if(swapping){
        return difference * (-1);
    }else{
        return difference; 
    }
}