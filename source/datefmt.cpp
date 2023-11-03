#include "datefmt.h"

datefmt::datefmt(int year, int month, int day){
    datefmt::day=day;
    datefmt::month=month;
    datefmt::year=year;
}
string datefmt::get_date(string type){ // if type == "NA" return "NA" if date is not known, else (default) return 0-0-0
    string date = to_string(year)+"-"+to_string(month)+"-"+to_string(day);
    if(type == "NA" && date == "0-0-0"){
        return "NA";
    }else{
        return date;
    }
}
int datefmt::get_year(){
    return datefmt::year;
}        
int datefmt::get_month(){
    return datefmt::month;
}
int datefmt::get_day(){
    return datefmt::day;
}