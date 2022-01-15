//
// Created by Carol Bez on 20-Jun-17.
//

#ifndef UNTITLED_PSFUNCTIONS_H
#define UNTITLED_PSFUNCTIONS_H


#include <complex>
#define MAX_FREQUENCIES 10000


#define PI  3.141592653589793238463




using namespace std;



void SetMNAPS(int ,int ,std::complex<double> [MAX_NODES + 1][MAX_NODES + 2], Element [MAX_ELEM],double);

int makeNetlistPS(Element[MAX_ELEM], Element [MAX_ELEM], int , double [MAX_NODES + 1]);

int makeFrequenciesList(double , double , double , double [MAX_FREQUENCIES + 1], double *, int );

int linspace(double ,double , double ,double [MAX_FREQUENCIES+1], int = 0 );

#endif //UNTITLED_PSFUNCTIONS_H
