#ifndef DATEFMT_H
#define DATEFMT_H

#include<string>

using namespace std;

class datefmt{
protected:
    int day, month, year;
public:
    datefmt(int year=0, int month=0, int day=0);
    string get_date();
    int get_year();
    int get_month();
    int get_day();
};
#endif