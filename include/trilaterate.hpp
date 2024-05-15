#ifndef TRILATERATE_HPP
#define TRILATERATE_HPP

#include <iostream>
#include <fstream>

# define M_PI 3.14159265358979323846
#define PYTHAG(x1, x2) sqrt(pow(x1, 2.0) + pow(x2, 2.0))



typedef struct {
	double x;
	double y;
} point_t;

typedef struct {
	point_t p;
	double r;
} circle_t;

typedef struct {
	point_t p;
	double slope;
} line_t;

typedef struct {
  double lat;
  double lon;
} coord_t;

typedef struct {
  coord_t *p;
  double r;
} record_t;


// Function prototypes

int circleCircleIntersection(circle_t c0, circle_t c1, point_t *p1, point_t *p2);

int lineLineIntersect(line_t l1, line_t l2, point_t *p);

double longitudeDegreeDistane(double lon);

double latitudeDegreeDistance(double lat);

double getDistance(coord_t p1, coord_t p2);

int trilaterate(record_t r1, record_t r2, record_t r3, coord_t *trilaterationCoord, coord_t *intersectionCoord, std::ofstream &exportFile);

int pointToCoord(coord_t base, point_t p, coord_t *res);

double getDelta(double a, double b);

// #warning Disable or erase testing functions for trilateration at deployment in file: include/trilaterate.hpp!

# define LOG(x) std::cout << x << std::endl;

void printLine(const char *c, line_t line);

void printPoint(const char *c, point_t point);

void printCoord(const char *c, coord_t coord);

void printCircle(const char *c, circle_t circle);

void printPointAsCoord(const char *c, coord_t base, point_t point);

#endif // TRILATERATE_HPP