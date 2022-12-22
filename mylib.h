#define _USE_MATH_DEFINES

#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include "MyVector.h"

#define EPS 1e-7
#define C_1_default M_PI / 4. // начальные значения, которые надо рассчитать
#define C_2_default (1 - pow(M_PI, 4) / 12.) / M_PI * 2 // -----//-----
#define ALPHA 0
#define K 32

using namespace std;

MyVector f(double t, MyVector k, double alph);
MyVector vec2Vec(vector<double> v);              

void RungeKutta(vector<vector<double> > &x, double C1, double C2, double alph);
void shooting(double *C, double alph);
void invJac(double J[2][2]);