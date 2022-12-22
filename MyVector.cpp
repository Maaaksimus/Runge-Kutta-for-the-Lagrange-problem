#include "MyVector.h"
#include <iostream>

MyVector::MyVector()
{
	x[0] = x[1] = x[2] = x[3] = 0;
}
MyVector::MyVector(double xx, double yy, double zz, double ww)
{
	x[0] = xx; x[1] = yy; x[2] = zz, x[3] = ww;
}
MyVector::MyVector(const MyVector &v)
{
	x[0] = v.x[0]; x[1] = v.x[1]; x[2] = v.x[2];
}

void MyVector::print() {
	std::cout << "[" << x[0] << ", " << x[1] << ", " << x[2] << ", " << x[3] << "]" << std::endl;
}

double operator*(const MyVector &a, const MyVector &b)
{
	return a.x[0]*b.x[0] + a.x[1]*b.x[1] + a.x[2]*b.x[2];
}

MyVector operator+(const MyVector &a, const MyVector &b)
{
	return MyVector(a.x[0]+b.x[0], a.x[1]+b.x[1], a.x[2]+b.x[2], a.x[3]+b.x[3]);
}

MyVector operator*(double c, const MyVector &a)
{
	return MyVector(c*a.x[0], c*a.x[1], c*a.x[2], c*a.x[3]);
}
/*const MyVector & operator=(const MyVector &a, const MyVector &b)
{
	a.x = b.x; a.x[1] = b.x[1]; a.x[2] = b.x[2];
	return a;
}*/

const MyVector & MyVector::operator=(const MyVector &b)
{
	x[0] = b.x[0]; x[1] = b.x[1]; x[2] = b.x[2], x[3] = b.x[3];
	return *this;
}