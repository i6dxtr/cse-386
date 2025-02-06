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
*             ......
*             ......
*  --- Rest of code truncated ---
*             ......
*             ......
*/