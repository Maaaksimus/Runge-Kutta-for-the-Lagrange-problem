#include <stdio.h>

class Vector
{
	public:

    double x[4];
	Vector();
	Vector(double xx, double yy, double zz, double ww);
	Vector(const Vector &v);
	~Vector() = default;

	friend double operator*(const Vector &a, const Vector &b);
	friend Vector operator+(const Vector &a, const Vector &b);
	friend Vector operator*(double c, const Vector &b);
    friend Vector operator-(double c, const Vector &b);
	
	//friend const Vector & operator=(Vector &a, const Vector &b);
	const Vector & operator=(const Vector &b);
};

double operator*(const Vector &a, const Vector &b);