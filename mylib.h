#define _USE_MATH_DEFINES

#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include "MyVector.h"
#include <stdio.h> // надо потом убрать
#include <iomanip>
#include <ios>

#define EPS 1e-13
#define DIFF 1e-10
// #define M_PI 3.14159265358979323846 как наладить точность?
#define C_1_default -(pow(M_PI, 4) / 384. - 2) / M_PI// -----//-----
#define C_2_default -1 * M_PI_4 / 2. // начальные значения, которые надо рассчитать
#define ALPHA 0
#define K 32

using namespace std;

MyVector f(MyVector k, double alph);

MyVector countStep(double h, MyVector y, double alph);
double L(MyVector x_curr, MyVector x_next, double alph);

void RungeKutta(vector<MyVector> &X, double C1, double C2, double alph, vector<double> &H);

void shooting(double *C, double alph);
