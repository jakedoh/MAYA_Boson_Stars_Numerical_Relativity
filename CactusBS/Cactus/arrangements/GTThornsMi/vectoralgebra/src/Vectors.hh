#ifndef VECTOR_ALGEBRA_UTILS
#define VECTOR_ALGEBRA_UTILS
#include "cctk.h"

namespace VectorAlgebra
{
	struct Vector
	{
		CCTK_REAL x, y, z;
		Vector(CCTK_REAL x_in, CCTK_REAL y_in, CCTK_REAL z_in);
		Vector();
		Vector operator-() const;
		CCTK_REAL getMagnitude() const;
		Vector getUnitVector() const;
		Vector& operator+=(const Vector rhs);
		Vector& operator-=(const Vector rhs);
		Vector& operator*=(const CCTK_REAL number);
		Vector& operator/=(const CCTK_REAL number);
		Vector& operator=(const Vector rhs);
	};

	Vector operator+(const Vector lhs, const Vector rhs);
	Vector operator-(const Vector lhs, const Vector rhs);
	Vector crossProduct(const Vector lhs, const Vector rhs);
	CCTK_REAL dotProduct(const Vector lhs, const Vector rhs);
	Vector operator*(const Vector vec, const CCTK_REAL number);
	Vector operator*(const CCTK_REAL number, const Vector vec);
	Vector operator/(const Vector vec, const CCTK_REAL number);

	struct LinearOperator
	{
		Vector x, y, z;
		LinearOperator();
		LinearOperator(Vector x_in, Vector y_in, Vector z_in);
		LinearOperator operator-() const;
		CCTK_REAL getTrace() const;
		CCTK_REAL getDeterminant() const;
		Vector transform(const Vector input) const;
		LinearOperator& operator+=(const LinearOperator rhs);
		LinearOperator& operator-=(const LinearOperator rhs);
		LinearOperator& operator*=(const LinearOperator rhs);
		LinearOperator& operator*=(const CCTK_REAL scalar);
		LinearOperator& operator=(const LinearOperator rhs);
	};

	LinearOperator operator+(const LinearOperator lhs, const LinearOperator rhs);
	LinearOperator operator-(const LinearOperator lhs, const LinearOperator rhs);
	LinearOperator operator*(const LinearOperator lhs, const LinearOperator rhs);
	LinearOperator operator*(const CCTK_REAL scalar, const LinearOperator map);
	LinearOperator operator*(const LinearOperator map, const CCTK_REAL scalar);

	//not even going there with matrix inversion
	
	const Vector ex(1., 0., 0.);
	const Vector ey(0., 1., 0.);
	const Vector ez(0., 0., 1.);
	const Vector zero(0., 0., 0.);

	const LinearOperator identity(ex, ey, ez);
}
#endif
