/**
 * Demetrius Hullum Scott  2/6/25
 * Exerpt of CSE 286 Codebase (2016-2024 Eric Bachmann and Mike Zmuda)
 * All functional code is my own.
 */

#include <iostream>
#include <istream>
#include <iomanip>
#include <cstdlib>

#include "defs.h"
#include "framebuffer.h"
#include "utilities.h"
#include "ishape.h"

/*
 * NOTE: There are several dozen functions in this file that you will implement over the course of the semester.
 * Implement only the ones listed for the given assignment (read the handout); you needn't do all of them in one sitting.
 */

 /**
  * @fn	void swap(double &a, double &b)
  * @brief	Swaps that values of two doubleing point numbers, without
  * 			using std.
  * @param [in,out]	a	First double.
  * @param [in,out]	b	Second double.
  */

void swap(double& a, double& b) {
	// simple variable swap
	double temp = a;
	a = b;
	b = temp;
}

/**
 * @fn	bool approximatelyEqual(double a, double b)
 * @brief	Determines if two values are approximately equal.
 * 			That is, their values within EPSILON of each other.
 * Programming constraint: Use EPSILON defined in defs.h
 * @param	a	The first value.
 * @param	b	The second value.
 * @return	true iff (a-b) is in [-EPSILON, EPSILON].
 * @test	approximatelyEqual(3.000000, 3.0000000001) --> true
 * @test	approximatelyEqual(3.000000, 3.1) --> false
*/

bool approximatelyEqual(double a, double b) {
	return fabs(a-b) < EPSILON; // |a-b| < e
}

/**
 * @fn	bool approximatelyZero(double a)
 * @brief	Determines if a value is approximately zero.
 * 			That is, the value is within EPSILON of zero.
 * Programming constraint: Use EPSILON defined in defs.h
 * @param	a	The value.
 * @return	true iff a is in [-EPSILON, EPSILON].
 * @test	approximatelyZero(0.0000001) --> true
 * @test	approximatelyZero(0.1) --> false
 */

bool approximatelyZero(double a) {
	return approximatelyEqual(a, 0.0);
}

/**
 * @fn	double normalizeDegrees(double degrees)
 * @brief	Converts an arbitrary number of degrees to an equivalent
 * 			number of degrees in the range [0, 360). Loops should NOT
 *          be used in this function. Recursion should also not be used.
 * Programming constraint: Do not use recursion or loops. Consider using glm::mod.
 * @param	degrees	The degrees.
 * @return	Normalized degrees in the range [0, 360).
 * @test	normalizeDegrees(0) --> 0
 * @test	normalizeDegrees(1.75) --> 1.75
 * @test	normalizeDegrees(-1) --> 359
 * @test	normalizeDegrees(-721) --> 359
 */

double normalizeDegrees(double degrees) {
	if (degrees < 0)
		return 360 - fabs(degrees);
	else
		return glm::mod(degrees, 360.0); // degrees % 360
}

/**
 * @fn	double normalizeRadians(double rads)
 * @brief	Converts an arbitrary number of radians to an equivalent
 * 			number of radians in the range [0, 2*PI). Loops should NOT
 *          be used in this function.
 * Programming constraint: Do not use recursion or loops.
 * @param	rads	The radians.
 * @return	Normalized radians in the range [0, 2*PI).
 * @test	normalizeRadians(0) --> 0
 * @test	normalizeRadians(1) --> 1
 * @test	normalizeRadians(3*PI) --> 3.141.....
 * @test	normalizeRadians(-31*PI) --> 3.141.....
 */

double normalizeRadians(double rads) {
	return fabs(glm::mod(rads, 2*PI)); // | rads % 2pi |
}

/**
 * @fn	double rad2deg(double rads)
 * @brief	Converts radians to degrees.  This function behaves like glm::degrees,
 * without using glm::degrees.
 * Programming constraint: Do not glm::degrees.
 * @param	rads	The radians.
 * @return	Degrees.
 */

double rad2deg(double rads) {
	return (rads * 180) / PI; 
}

/**
 * @fn	double deg2rad(double degs)
 * @brief	Converts degrees to radians. This function behaves like glm::radians,
 * without using glm::radians.
 * Programming constraint: Do not use glm::radians.
 * @param	degs	The degrees.
 * @return	Radians.
 */

double deg2rad(double degs) {
	return (degs * PI) / 180;
}

/**
* @fn	double min(double A, double B, double C)
* @brief	Determines the minimum of three values, using glm::min.
* Programming constraint: Use glm::min, which provides the minimum of two numbers
* @param	A	First value.
* @param	B	Second value
* @param	C	Third value.
* @return	The minimum value.
*/

double min(double A, double B, double C) {
	return glm::min(A, B);
}

/**
* @fn	double max(double A, double B, double C)
* @brief	Determines the maximum of three values, using glm::max.
* Programming constraint: Use glm::max
* @param	A	First value.
* @param	B	Second value
* @param	C	Third value.
* @return	The maximum value.
*/

double max(double A, double B, double C) {
	return glm::max(glm::max(A, B), C);
}

/**
* @fn	distanceFromOrigin(double x, double y)
* @brief	Determines the distance of the point (x, y) to (0, 0).
* The distance is defined by sqrt(x^2 + y^2). Note: ^ is not how
* C++ does exponentiation; you can use glm::pow instead.
* @param	x	The x coordinate
* @param	y	The 7 coordinate.
* @return	The distance of (x, y) to the origin.
* @test	distanceFromOrigin(0, 1) --> 1.0
* @test	distanceFromOrigin(1, 0) --> 1.0
* @test	distanceFromOrigin(1, 1) --> 1.41421356237309514547
* @test	distanceFromOrigin(-10, 30) --> 31.62277660168379256334
*/

double distanceFromOrigin(double x, double y) {
	return sqrt(glm::pow(x, 2) + glm::pow(y, 2));
}

/**
* @fn	distanceBetween(double x1, double y1, double x2, double y2)
* @brief	Determines the distance of the point (x1, y1) to (x2, y2)
* The distance is defined by sqrt((x1-x2)^2 + (y1-y2)^2). Note: ^ is not how
* C++ does exponentiation; you can use glm::pow instead.
* @param	x1	The first x coordinate
* @param	y1	The first y coordinate.
* @param	x2	The second x coordinate
* @param	y2	The second y coordinate.
* @return	The distance between (x1, y1) and (x2, y2).
* @test	distanceBetween(0, 0, 1, 1) --> 1.41421356237309514547
* @test	distanceBetween(1, 1, 0, 0) --> 1.41421356237309514547
* @test	distanceBetween(10, 10, 11, 11) --> 1.41421356237309514547
* @test	distanceBetween(100, 100, 99, 99) --> 1.41421356237309514547
* @test	distanceBetween(54, 1, -34, -99) --> 133.2066064427737
*/

double distanceBetween(double x1, double y1, double x2, double y2) {
	// Application of Euclidean distance formula:
	double dx = glm::pow(x1, 2) - glm::pow(x2, 2); // x1^2-x2^2
	double dy = glm::pow(y1, 2) - glm::pow(y2, 2); // y1^2 - y2^2
	return dx + dy;
}

/**
 * @fn	double areaOfTriangle(double a, double b, double c)
 * @brief	Computes the area of triangle using Heron's formula.
 * @param	a length of first side.
 * @param	b length of second side.
 * @param	c length of third side.
 * @return	Area of triangle. Returns -1.0 if the triangle is illegal (i.e.
 * negative lengths). Legal values will yield v > 0.
 * @test	areaOfTriangle(3, 4, 5) --> 6
 * @test	areaOfTriangle(-3, 4, 5) --> -1
 * @test	areaOfTriangle(3, 4, 50) --> -1
 */

double areaOfTriangle(double a, double b, double c) {
	double s = (a + b + c) / 2; 
	double result = sqrt(s * (s - a) * (s - b) * (s - c)); // Heron's formula
	return result;
}

/**
 * @fn	double areaOfTriangle(double x1, double y1, double x2, double y2, double x3, double y3)
 * @brief	Computes the area of triangle formed by the three vertices (x1, y1), (x2, y2), and
 * (x3, y3). You can assume all vertices are distinct.
 * @param	x1 the x value of the first vertice
 * @param	y1 the y value of the first vertice
 * @param	x2 the x value of the second vertice
 * @param	y2 the y value of the second vertice
 * @param	x3 the x value of the third vertice
 * @param	y3 the y value of the third vertice
 * @return	Area of triangle.
 * @test	areaOfTriangle(0, 0, 3, 0, 0, 4) --> 6
 */

double areaOfTriangle(double x1, double y1, double x2, double y2, double x3, double y3) {
	double a = distanceBetween(x1, y1, x2, y2);
	double b = distanceBetween(x2, y2, x3, y3);
	double c = distanceBetween(x3, y3, x1, x2);
	return areaOfTriangle(a, b, c);
}

/**
 * @fn	void pointOnUnitCircle(double angleRads, double &x, double &y)
 * @brief	Determines the (x,y) position of a point on the standard
 * 			unit circle.
 * @param 		  	angleRads	The angle in radians.
 * @param [in,out]	x		 	A double to process.
 * @param [in,out]	y		 	A double to process.
 */

void pointOnUnitCircle(double angleRads, double& x, double& y) {
	x = glm::cos(angleRads); // cosine represents x-axis
	y = glm::sin(angleRads); // sine represents y-axis
}

/**
* @fn	dvec2 pointOnCircle(const dvec2 &center, double R, double angleRads)
* @brief	Computes the (x,y) value of a point on the circle centered on 'center',
* 			having radius R. The point is determined by sweeping an arc 'angleRads'.
* @param	center   	The center of the circle
* @param	R		 	Radius of circle.
* @param	angleRads	The angle in radians.
* @return	The point on the circle.
*/

dvec2 pointOnCircle(const dvec2& center, double R, double angleRads) {
	// stretch cosine/sine expression by the radius
	// perform vector addition to translate vector to correct position
	double x = center.x + (R * glm::cos(angleRads));
	double y = center.y + (R * glm::sin(angleRads));
	return dvec2(x, y); // a 2D vector
}

/**
* @fn	double directionInRadians(const dvec2 &referencePt, const dvec2 &targetPt)
* @brief	Compute the direction/heading of 'targetPt', relative
* 			to referencePt. The return angle should be in [0, 2PI)
* @param	referencePt	Reference point.
* @param	targetPt	Target point point.
* @return	A double.
* The test cases below casually use (0,0). When writing code, consider
* (0,0) to be dvec2(0,0).
* @test	directionInRadians(dvec2(0,0), dvec2(2,2)) --> 0.7853981634
* @test	directionInRadians(dvec2(2,10), dvec2(3,11)) --> 0.7853981634
* @test	directionInRadians(dvec2(2,2), dvec2(2,0)) --> 4.7123889804
* @test directionInRadians(dvec2(1,-1), dvec2(1.3420, -1.93969)) --> 5.06144
*/

double directionInRadians(const dvec2& referencePt, const dvec2& targetPt) {
	double dx = targetPt.x - referencePt.x; // centers vector on x-coordinate
	double dy = targetPt.y - referencePt.y; // centers vector on y-coordinates 
	double angleRads = glm::atan(dy, dx);   // compute angle using arctangent with resulting coordinates
	return glm::mod(angleRads + (2.0 * PI), (2.0 * PI)); // fits result onto interval [0, 2pi)
}

/**
* @fn	double directionInRadians(const dvec2 &targetPt)
* @brief	Compute the direction/heading of 'targetPt', relative
* 			to the origin. The return angle should be in [0, 2PI)
* @param	targetPt	Target point point.
* @return	A double.
* @test	directionInRadians(dvec2(2,2)) --> 0.7853981634
* @test	directionInRadians(dvec2(0,-2)) --> 4.7123889804
*/

double directionInRadians(const dvec2& targetPt) {
	double angleRads = glm::atan(targetPt.x, targetPt.y); // computes angle using arctangent given param coords
	return glm::mod(angleRads + (2.0 * PI), (2.0 * PI)); // fits result onto interval [0, 2pi)
}

/**
* @fn	double directionInRadians(double x1, double y1, double x2, double y2)
* @brief	Compute the direction/heading of (x2, y2), relative
* 			to (x1, y1). The return angle should be in [0, 2PI)
* @param	x1  x coordinate of the reference point.
* @param	y1	y coordinate of the reference point.
* @param  x2    x coordinate of the target point.
* @param  y2    y coordinate of the target point.
* @return	A double.
* @test	directionInRadians(0,0,2,2) --> 0.7853981634
* @test	directionInRadians(2,10,3,11) --> 0.7853981634
* @test	directionInRadians(2,2,2,0) --> 4.7123889804
*/

double directionInRadians(double x1, double y1, double x2,  double y2) {
	double dx = x2 - x1; // centers vector on x-coordinate
	double dy = y2 - y1; // centers vector on x-coordinate
	double angleRads = glm::atan(dy, dx); // computes angle using arctangent with resulting coordinates
	return glm::mod(angleRads + (2.0 * PI), (2.0 * PI)); // fits result onto interval [0, 2pi)
}



/**
 * @fn	double map(double x, double fromLo, double fromHi, double toLow, double toHigh)
 * @brief	Linearly map a value from one interval to another.
 * @param 		  	x	 	x value.
 * @param 		  	fromLo 	The low value of the x range.
 * @param 		  	fromHi	The high value of the x range.
 * @param 		  	toLow 	The low value of the new range.
 * @param 		  	toHigh	The high value of the new range.
 * @test	map(2, 0, 5, 10, 11) --> 10.4
 */

double map(double x, double fromLo, double fromHi, double toLow, double toHigh) {
	double term1 = (x - fromLo) / (fromHi - fromLo); 
	double term2 = (toHigh - toLow);
	return term1 * term2 + toLow;
}

/**
 * @fn	vector<double> quadratic(double A, double B, double C)
 * @brief	Solves the quadratic equation, given A, B, and C.
 * 			0, 1, or 2 roots are inserted into the vector and returned.
 * 			The roots are placed into the vector sorted in ascending order.
 *          vector is somewhat like Java's ArrayList. Do a little research on
 *          it. The length of the vector will correspond to the number of roots.
 * @param	A	A.
 * @param	B	B.
 * @param	C	C.
 * @return	Vector containing the real roots.
 * @test	quadratic(1,4,3) --> [-3,-1]
 * @test	quadratic(1,0,0) --> [0]
 * @test	quadratic(-4, -2, -1) --> []
 */

vector<double> quadratic(double A, double B, double C) {
	double roots[2] = { 0,0 };
	quadratic(A, B, C, roots);
	vector<double> result;	// put only the roots in here
	result.push_back(roots[1]);
	result.push_back(roots[0]);
	return result;
}

/**
 * @fn	int quadratic(double A, double B, double C, double roots[2])
 * @brief	Solves the quadratic equation, given A, B, and C.
 * 			0, 1, or 2 roots are inserted into the array 'roots'.
 * 			The roots are sorted in ascending order.
 * Here is an example of how this is to be used:
 *
 * 	double roots[2];
 *	int numRoots = quadratic(1, 2, -3, roots);
 *	if (numRoots == 0) {
 *		cout << "There are no real roots" << endl;
 *	} else if (numRoots == 1) {
 *		cout << "Only one root: " << roots[0] << endl;
 *	} else if (numRoots == 2) {
 *      if (roots[0] > roots[1])
 *			cout << "Something is wrong. This should not happen" << endl;
 *		else
 *			cout << "Two roots: " << roots[0] << " and " << roots[1] << endl;
 *	} else {
 *		cout << "Something is wrong. This should not happen" << endl;
 *	}
 *
 * @param	A	 	A.
 * @param	B	 	B.
 * @param	C	 	C.
 * @param	roots	The real roots.
 * @test	quadratic(1, 4, 3, ary) --> returns 2 and fills in ary with: [-3,-1]
 * @test	quadratic(1 ,0, 0, ary) --> returns 1 and fills in ary with: [0]
 * @test	quadratic(-4, -2, -1, ary) --> returns 0 and does not modify ary.
 * @return	The number of real roots put into the array 'roots'
*/

int quadratic(double A, double B, double C, double roots[2]) {
	/* CSE 386 - todo  */
	int rootCnt = 0;
	double disc = glm::pow(B, 2) - 4 * A * C;
	if (disc < 0) {
		return 0;
	} if (disc == 0) {
		roots[0] = -B / (2 * A);
		return 1;
	}
	else {
		double r1 = (-B + glm::sqrt(glm::pow(B, 2) - 4 * A * C)) / (2 * A);
		double r2 = (-B - glm::sqrt(glm::pow(B, 2) - 4 * A * C)) / (2 * A);
		roots[0] = glm::min(r1, r2);
		roots[1] = glm::max(r1, r2);
		return 2;
	}
	return -1; // error case
}

/**
* @fn	dvec2 doubleIt(const dvec2 &V)
* @brief	Computes 2*V
* @param	V	The vector.
* @return	2*V.
*/

dvec2 doubleIt(const dvec2& V) {
	return 2.0 * V;
}

/**
* @fn	dvec3 myNormalize(const dvec3 &V)
* @brief	Computes the normalization of V, without calling glm::normalize.
*           The input vector is not be the zero vector.
* Programming constraint: Do NOT use glm::normalize
* @param	V	The vector.
* @return	Normalized vector V.
*/

dvec3 myNormalize(const dvec3& V) {
	double norm = glm::sqrt(glm::dot(V, V));
	return dvec3(0, 0, 0) + (V / norm);
}

/**
* @fn	bool isOrthogonal(const dvec3 &a, const dvec3 &b)
* @brief	Determines if two vectors are orthogonal, or nearly orthogonal. The inputs are non-zero vectors.
Two vectors are nearly orthogonal if the cosine of the angle formed by these
two vectors is approximatelyZero().
* @param	a	The first vector.
* @param	b	The second vector.
* @return	True iff the two vector are orthogonal.
*/

bool isOrthogonal(const dvec3& a, const dvec3& b) {
	if (approximatelyZero(glm::dot(a, b)))
		return true;
	return false;
}

/**
* @fn	bool formAcuteAngle(const dvec3 &a, const dvec3 &b)
* @brief	Determines if two vectors form an angle that is < 90 degrees. The inputs are non-zero vectors.
* Programming constraint: Do NOT use acos, atan, or asin (you CAN use dot, cos, etc)
* @param	a	The first vector.
* @param	b	The second vector.
* @return	True iff the two vectors form an acute angle.
*/

bool formAcuteAngle(const dvec3& a, const dvec3& b) {
	double dtpd = glm::dot(a, b);
	if (!approximatelyZero(dtpd)) {
		if (dtpd > 0)
			return true;
	}
	return false;
}

/**
 * @fn	double cosBetween(const dvec2 &v1, const dvec2 &v2)
 * @brief	Cosine between v1 and v2. The inputs are non-zero vectors.
 * @param	v1	The first vector.
 * @param	v2	The second vector.
 * @test	cosBetween(dvec2(1.0, 0.0), dvec2(1.0, 0.0)) --> 1.0
 * @test	cosBetween(dvec2(1.0, 0.0), dvec2(1.0, 1.0)) --> 0.707107
 * @test	cosBetween(dvec2(-1.0, glm::sqrt(3.0)), dvec2(-1.0, 0.0)) --> 0.5
 * @test	cosBetween(dvec2(-1.0, glm::sqrt(3.0)), dvec2(1.0, glm::sqrt(3.0))) --> 0.5
 * @return	The cosine between v1 and v2.
 */

double cosBetween(const dvec2& v1, const dvec2& v2) {
	double dtpd = glm::dot(v1, v2);
	double v1Norm = glm::sqrt(glm::dot(v1, v1));
	double v2Norm = glm::sqrt(glm::dot(v2, v2));
	return dtpd / (v1Norm * v2Norm);
}

/**
 * @fn	double cosBetween(const dvec3 &v1, const dvec3 &v2)
 * @brief	Computes the cosine between v1 and v2.
 * @param	v1	The first vector.
 * @param	v2	The second vector.
 * @return	A double.
 */

double cosBetween(const dvec3& v1, const dvec3& v2) {
	double dtpd = glm::dot(v1, v2);
	double v1Norm = glm::sqrt(glm::dot(v1, v1));
	double v2Norm = glm::sqrt(glm::dot(v2, v2));
	return dtpd / (v1Norm * v2Norm);
}

/**
 * @fn	double cosBetween(const dvec4 &v1, const dvec4 &v2)
 * @brief	Computes the cosine between v1 and v2.
 * @param	v1	The first vector.
 * @param	v2	The second vector.
 * @return	A double.
 */

double cosBetween(const dvec4& v1, const dvec4& v2) {
	return glm::dot(v1, v2) / (glm::length(v1) * glm::length(v2));;
}

/**
 * @fn	double areaOfParallelogram(const dvec3 &v1, const dvec3 &v2)
 * @brief	Computes the area of parallelogram, given two vectors eminating
 * 			from the same corner of the parallelogram.
 * @param	v1	The first vector.
 * @param	v2	The second vector.
 * @test	areaOfParallelogram(dvec3(1.0, 0.0, 0.0), dvec3(0.0, 1.0, 0.0)) --> 1.0
 * @test	areaOfParallelogram(dvec3(1.0, 1.0, 1.0), dvec3(1.0, 0.0, 1.0)) --> 1.41421
 * @return	Area of parallelogram.
 */

double areaOfParallelogram(const dvec3& v1, const dvec3& v2) {
	return glm::length(glm::cross(v1, v2));
}

/**
 * @fn	double areaOfTriangle(const dvec3 &pt1, const dvec3 &pt2, const dvec3 &pt3)
 * @brief	Computes the area of triangle.
 * Programming constraint: use areaOfParalellogram to solve this one.
 * @param	pt1	The first point.
 * @param	pt2	The second point.
 * @param	pt3	The third point.
 * @test	areaOfTriangle(dvec3(0.0, 0.0, 0.0), dvec3(1.0, 0.0, 0.0), dvec3(0.0, 1.0, 0.0)) --> 0.5
 * @test	areaOfTriangle(dvec3(-10.0, -10.0, -10.0), dvec3(-11.0, -10.0, -10.0), dvec3(-10.0, -11.0, -10.0)) --> 0.5
 * @return	Area of triangle.
 */

double areaOfTriangle(const dvec3& pt1, const dvec3& pt2, const dvec3& pt3) {
	return 0.5 * areaOfParallelogram(pt2 - pt1, pt3 - pt1);
}

/**
* @fn	dvec3 pointingVector(const dvec3 &pt1, const dvec3 &pt2)
* @brief	Computes unit-length pointing vector.
* @param	pt1	The first point.
* @param	pt2	The second point.
* @return	Unit length vector that points from pt1 toward pt2.
*/

dvec3 pointingVector(const dvec3& pt1, const dvec3& pt2) {
	return glm::normalize(pt2 - pt1);
}

/**
*             ......
*             ......
*  --- Rest of code truncated ---
*             ......
*             ......
*/