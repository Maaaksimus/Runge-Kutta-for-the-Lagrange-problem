#include <stdio.h>

class MyVector
{
	public:

    double x[4];
	MyVector();
	MyVector(double xx, double yy, double zz, double ww);
	MyVector(const MyVector &v);
	~MyVector() = default;

	friend double operator*(const MyVector &a, const MyVector &b);
	friend MyVector operator+(const MyVector &a, const MyVector &b);
	friend MyVector operator*(double c, const MyVector &b);
    // friend MyVector operator-(double c, const MyVector &b);
	
	//friend const MyVector & operator=(MyVector &a, const MyVector &b);
	const MyVector & operator=(const MyVector &b);
};

double operator*(const MyVector &a, const MyVector &b);