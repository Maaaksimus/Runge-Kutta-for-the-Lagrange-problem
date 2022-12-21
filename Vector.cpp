#include "Vector.h"

Vector::Vector()
{
	x[0] = x[1] = x[2] = x[3] = 0;
}
Vector::Vector(double xx, double yy, double zz, double ww)
{
	x[0] = xx; x[1] = yy; x[2] = zz, x[3] = ww;
}
Vector::Vector(const Vector &v)
{
	x[0] = v.x[0]; x[1] = v.x[1]; x[2] = v.x[2];
}

double operator*(const Vector &a, const Vector &b)
{
	return a.x[0]*b.x[0] + a.x[1]*b.x[1] + a.x[2]*b.x[2];
}

Vector operator+(const Vector &a, const Vector &b)
{
	return Vector(a.x[0]+b.x[0], a.x[1]+b.x[1], a.x[2]+b.x[2], a.x[3]+b.x[3]);
}

Vector operator*(double c, const Vector &a)
{
	return Vector(c*a.x[0], c*a.x[1], c*a.x[2], c*a.x[3]);
}
/*const Vector & operator=(const Vector &a, const Vector &b)
{
	a.x = b.x; a.x[1] = b.x[1]; a.x[2] = b.x[2];
	return a;
}*/

const Vector & Vector::operator=(const Vector &b)
{
	x[0] = b.x[0]; x[1] = b.x[1]; x[2] = b.x[2], x[3] = b.x[3];
	return *this;
}