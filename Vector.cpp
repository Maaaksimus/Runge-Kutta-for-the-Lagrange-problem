#include "Vector.h"

Vector::Vector()
{
	x1 = x2 = x3 = x4 = 0;
}
Vector::Vector(double xx, double yy, double zz, double ww)
{
	x1 = xx; x2 = yy; x3 = zz, x4 = ww;
}
Vector::Vector(const Vector &v)
{
	x1 = v.x1; x2 = v.x2; x3 = v.x3;
}

double operator*(const Vector &a, const Vector &b)
{
	return a.x1*b.x1 + a.x2*b.x2 + a.x3*b.x3;
}

Vector operator+(const Vector &a, const Vector &b)
{
	return Vector(a.x1+b.x1, a.x2+b.x2, a.x3+b.x3, a.x4+b.x4);
}

Vector operator*(double c, const Vector &a)
{
	return Vector(c*a.x1, c*a.x2, c*a.x3, c*a.x4);
}
/*const Vector & operator=(const Vector &a, const Vector &b)
{
	a.x = b.x; a.x2 = b.x2; a.x3 = b.x3;
	return a;
}*/

const Vector & Vector::operator=(const Vector &b)
{
	x1 = b.x1; x2 = b.x2; x3 = b.x3, x4 = b.x4;
	return *this;
}