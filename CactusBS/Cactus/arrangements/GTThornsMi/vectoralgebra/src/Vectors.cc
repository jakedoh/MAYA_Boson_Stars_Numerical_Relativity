#include <VectorAlgebra.hh>
#include <cmath>
using namespace VectorAlgebra;

Vector::Vector()
{
	x = 0.;
	y = 0.;
	z = 0.;
}

Vector::Vector(CCTK_REAL x_in, CCTK_REAL y_in, CCTK_REAL z_in)
{
	x = x_in;
	y = y_in;
	z = z_in;
}

Vector& VectorAlgebra::Vector::operator+=(const Vector rhs)
{
	x += rhs.x;
	y += rhs.y;
	z += rhs.z;
	return *this;
}

Vector& VectorAlgebra::Vector::operator=(const Vector rhs)
{
	x = rhs.x;
	y = rhs.y;
	z = rhs.z;
	return *this;
}

Vector& VectorAlgebra::Vector::operator-=(const Vector rhs)
{
	*this+=(-rhs);
	return *this;
}

Vector& VectorAlgebra::Vector::operator*=(const CCTK_REAL number)
{
	x *= number;
	y *= number;
	z *= number;
	return *this;
}

Vector& VectorAlgebra::Vector::operator/=(const CCTK_REAL number)
{
	*this*=(1./number);
	return *this;
}

Vector VectorAlgebra::Vector::operator-() const
{
	return Vector(-x, -y, -z);
}

Vector VectorAlgebra::operator+(const Vector lhs, const Vector rhs)
{
	Vector result = rhs;
	result += lhs;
	return result;
}

Vector VectorAlgebra::operator-(const Vector lhs, const Vector rhs)
{
	return lhs+(-rhs);
}

Vector VectorAlgebra::operator*(const Vector vec, const CCTK_REAL number)
{
	Vector result = vec;
	result*=number;
	return result;
}

Vector VectorAlgebra::operator*(const CCTK_REAL number, const Vector vec)
{
	return vec*number;
}

Vector VectorAlgebra::operator/(const Vector vec, const CCTK_REAL number)
{
	Vector result = vec;
	result/=number;
	return result;
}

Vector VectorAlgebra::crossProduct(const Vector lhs, const Vector rhs)
{
	return Vector(lhs.y*rhs.z - lhs.z*rhs.y, lhs.z*rhs.x-lhs.x*rhs.z, lhs.x*rhs.y-lhs.y*rhs.x);
}

CCTK_REAL VectorAlgebra::dotProduct(const Vector lhs, const Vector rhs)
{
	return lhs.x*rhs.x + lhs.y*rhs.y + lhs.z*rhs.z;
}

CCTK_REAL Vector::getMagnitude() const
{
	return std::sqrt(dotProduct(*this, *this));
}

Vector Vector::getUnitVector() const
{
	return Vector(x, y, z)/getMagnitude();
}

LinearOperator::LinearOperator()
{
	x = zero; y = zero; z = zero;
}

LinearOperator::LinearOperator(Vector x_in, Vector y_in, Vector z_in)
{
	x = x_in; y = y_in; z = z_in;
}

LinearOperator VectorAlgebra::LinearOperator::operator-() const
{
	return LinearOperator(-x, -y, -z);
}

CCTK_REAL LinearOperator::getTrace() const
{
	return x.x+y.y+z.z;
}

CCTK_REAL LinearOperator::getDeterminant() const
{
	return dotProduct(crossProduct(x, y), z);
}

Vector LinearOperator::transform(const Vector input) const
{
	return x*input.x + y*input.y + z*input.z;
}

LinearOperator& VectorAlgebra::LinearOperator::operator+=(const LinearOperator rhs)
{
	x += rhs.x;
	y += rhs.y;
	z += rhs.z;
	return *this;
}

LinearOperator& VectorAlgebra::LinearOperator::operator-=(const LinearOperator rhs)
{
	*this+=(-rhs);
	return *this;
}

LinearOperator& VectorAlgebra::LinearOperator::operator*=(const LinearOperator rhs)
{
	VectorAlgebra::Vector xnew, ynew, znew;
	xnew = transform(rhs.transform(ex));
	ynew = transform(rhs.transform(ey));
	znew = transform(rhs.transform(ez));

	x = xnew;
	y = ynew;
	z = znew;
	return *this;
}

LinearOperator& VectorAlgebra::LinearOperator::operator*=(const CCTK_REAL scalar)
{
	x *= scalar;
	y *= scalar;
	z *= scalar;
	return *this;
}

LinearOperator& VectorAlgebra::LinearOperator::operator=(const LinearOperator rhs)
{
	x = rhs.x;
	y = rhs.y;
	z = rhs.z;
	return *this;
}

LinearOperator VectorAlgebra::operator+(const LinearOperator lhs, const LinearOperator rhs)
{
	LinearOperator result = rhs;
	result += lhs;
	return result;
}

LinearOperator VectorAlgebra::operator-(const LinearOperator lhs, const LinearOperator rhs)
{
	return lhs+(-rhs);
}
	
LinearOperator VectorAlgebra::operator*(const LinearOperator lhs, const LinearOperator rhs)
{
	LinearOperator result = lhs;
	result*=rhs;
	return result;
}

LinearOperator VectorAlgebra::operator*(const CCTK_REAL scalar, const LinearOperator map)
{
	LinearOperator result = map;
	result*=scalar;
	return result;
}

LinearOperator VectorAlgebra::operator*(const LinearOperator map, const CCTK_REAL scalar)
{
	return scalar*map;
}

