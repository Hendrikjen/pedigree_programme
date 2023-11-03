#include "basic_functions.h"

string to_string_with_precision(double number, int n) { // converts a double to a string representation with a specified number of decimal places
    try{
        ostringstream out;
        out.precision(n);
        out << fixed << number;
        return out.str();
    }catch(const exception &ex) {
        cerr << "Error in to_string_with_precision(double): " << ex.what() << endl;
        return "unable to convert to string with precision";
    }
}
string to_string_with_precision(float number, int n) {// converts a float to a string representation with a specified number of decimal places
    try{
        ostringstream out;
        out.precision(n);
        out << fixed << number;
        return out.str();
    }catch(const exception &ex) {
        cerr << "Error in to_string_with_precision(float): " << ex.what() << endl;
        return "unable to convert to string with precision";
    }
}
string to_string_with_precision(int number, int n) { // converts an integer to a string representation with a specified number of decimal places
    try{
        ostringstream out;
        out.precision(n);
        out << fixed << number;
        return out.str();
    }catch(const exception &ex) {
        cerr << "Error in to_string_with_precision(int): " << ex.what() << endl;
        return "unable to convert to string with precision";
    }
}
std::deque <string> str_split(string str, string pattern){ // based on the delimiter pattern, it splits the given string into substrings
    std::deque <string> output = {};
    try{
        int i = 0;
        int j = str.find(pattern);
        while (j != -1) { // as long as the pattern is found another time
            output.push_back(str.substr(i, j-i)); // add substring to deque
            i = j + pattern.size(); // remove the former substring from the original string (to start the search again)
            j = str.find(pattern, i); // find the next occurance of the pattern
        }
        output.push_back(str.substr(i, j-i));
        return output;
    }catch(const std::exception &ex) {
        std::cerr << "Error in str_split(): " << ex.what() << std::endl;
        return output;
    }
}
bool is_leapyear(int year){ // returns true if year is a leap year, else false
    try{
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
    }catch(const std::exception &ex) {
        std::cerr << "Error in is_leapyear(): " << ex.what() << std::endl;
        return false;
    }
}
int date_diff(datefmt date_1, datefmt date_2){ // returns difference between two dates in days
    try{
        int difference = 0;
        bool swapping = false; // to keep track, if it was swapped (to invert the date difference in the end if date_2 < date_1)
        if(date_1.get_year()>=date_2.get_year()){ // if necessary, swap so that date_1 < date_2 (check first the year)
            swap(date_1,date_2);
            swapping = true;
            if((date_1.get_year()==date_2.get_year()) && (date_1.get_month()>date_2.get_month())){ // if necessary, swap so that date_1 < date_2 (then the month)
                swap(date_1,date_2);
                swapping = false;
                if((date_1.get_month()==date_2.get_month())&&(date_1.get_day()>date_2.get_day())){ // if necessary, swap so that date_1 < date_2 (at last the day)
                    swap(date_1,date_2);
                    swapping = true;
                }
            }
        }
        for(int i = date_1.get_year();i<date_2.get_year();i++){ // get the difference in days based on the discrepancy of the years
            if(is_leapyear(i)){
                difference += 366;
            }else{
                difference += 365;
            }
        }
        for(int i = 1;i<date_1.get_month();i++){ // get the amount of days of the first year before date_1 (from January 1st to day of date_1)
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
        difference -= date_1.get_day(); // substract day discrepancy
        for(int i = 1;i<date_2.get_month();i++){ // get the amount of days of the last year (not considered before) --> (from January 1st to day of date_2)
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
    }catch(const std::exception &ex) {
        std::cerr << "Error in date_diff(): " << ex.what() << std::endl;
        return 0;
    }
}