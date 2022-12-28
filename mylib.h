#define _USE_MATH_DEFINES

#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include "MyVector.h"
#include <stdio.h> // надо потом убрать
#include <iomanip>
#include <ios>

#define EPS 1e-8
#define DIFF 1e-10
#define gig 3.14159265358979323846 // как наладить точность?
// #define C_1_default -(pow(gig, 4) / 384. - 2) / gig// -----//-----
#define C_1_default 0.5558742601793006
#define C_2_default -0.3926990816987242
// #define C_2_default -1 * gig / 8. // начальные значения, которые надо рассчитать
#define ALPHA 0
#define K 32

using namespace std;

MyVector f(MyVector k, double alph);
void countK(MyVector *k, double h, MyVector y, double alph);
MyVector vec2Vec(vector<double> v);        

void RungeKutta(vector<vector<double> > &x, double C1, double C2, double alph, vector<double> &h);
void shooting(double *C, double alph);
void invJac(double J[2][2]);