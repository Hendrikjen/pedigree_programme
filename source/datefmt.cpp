#include "datefmt.h"

datefmt::datefmt(int year, int month, int day){
    datefmt::day=day;
    datefmt::month=month;
    datefmt::year=year;
}
string datefmt::get_date(){
    return to_string(year)+"-"+to_string(month)+"-"+to_string(day);
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