#include <math.h>
#include <trilaterate.hpp>
#include <iostream>
#include <algorithm>
#include <fstream>

// Reference: https://paulbourke.net/geometry/circlesphere/
int circleCircleIntersection(circle_t c0, circle_t c1, point_t *p1, point_t *p2)
{
  double a, dx, dy, d, h, rx, ry;
  double x0, x1, x2, y0, y1, y2, r0, r1;
  
  x0 = c0.p.x;
  x1 = c1.p.x;
  y0 = c0.p.y;
  y1 = c1.p.y;

  r0 = c0.r;
  r1 = c1.r;
  
  /* dx and dy are the vertical and horizontal distances between
   * the circle centers.
   */
  dx = x1 - x0;
  dy = y1 - y0;

  /* Determine the straight-line distance between the centers. */
  //d = sqrt((dy*dy) + (dx*dx));
  d = hypot(dx,dy); // Suggested by Keith Briggs

  /* Check for solvability. */
  if (d > (r0 + r1))
  {
    /* no solution. circles do not intersect. */
    return 0;
  }
  if (d < fabs(r0 - r1))
  {
    /* no solution. one circle is contained in the other */
    return 0;
  }

  /* 'point 2' is the point where the line through the circle
   * intersection points crosses the line between the circle
   * centers.  
   */

  /* Determine the distance from point 0 to point 2. */
  a = ((r0*r0) - (r1*r1) + (d*d)) / (2.0 * d) ;

  /* Determine the coordinates of point 2. */
  x2 = x0 + (dx * a/d);
  y2 = y0 + (dy * a/d);

  /* Determine the distance from point 2 to either of the
   * intersection points.
   */
  h = sqrt((r0*r0) - (a*a));

  /* Now determine the offsets of the intersection points from
   * point 2.
   */
  rx = -dy * (h/d);
  ry = dx * (h/d);

  /* Determine the absolute intersection points. */
  p1->x = x2 + rx;
  p1->y = y2 + ry;

  p2->x = x2 - rx;
  p2->y = y2 - ry;
  // l->slope = (x1 - x0) / (y1 - y0) * -1; // inverse slope
  
  return 1;
}


// Calculate intersection between 2 lines. Store result in point_t p
// Returns:
//  0     Fail
//  1     Succes (result in p)
int lineLineIntersect(line_t l1, line_t l2, point_t *p)
{
  if (l1.slope == l2.slope) // No cross section
    return 0;
  
  // y=-2.49981(x+-687.318)+3490.75
  double a1, c1, a2, c2, b1, b2;
  double xi, yi;

  a1 = l1.slope;
  c1 = l1.p.x *l1.slope + l1.p.y;

  a2 = l2.slope;
  c2 = l2.p.x *l2.slope + l2.p.y;

  xi = (c2-c1) / (a1-a2);
  yi = (a2*c1 - a1*c2) / ((a1)-(a2)) * -1;

  p->x = xi;
  p->y = yi;

  return 1;
}

// Calculate length of degree of longitude at latitude x in meters
// Reference: https://en.wikipedia.org/wiki/Geographic_coordinate_system
double longitudeDegreeDistane(double lon)
{
  // Degrees to radians
  lon = lon * M_PI / 180.0;

  double res = (111412.84*cos(lon)) - (93.5 * cos(3*lon)) + (0.118 * cos(5*lon));
  return res;
}

// Calculate length of degree of latitude at latitude x in meters
// Reference: https://en.wikipedia.org/wiki/Geographic_coordinate_system
double latitudeDegreeDistance(double lat)
{
	lat = lat * M_PI / 180.0;
	return (111132.92 - 559.82*cos(2*lat)) + (1.175*cos(4*lat)) + (0.0023*cos(6*lat));
}

int trilaterate(record_t r1, record_t r2, record_t r3, coord_t *trilaterationCoord, coord_t *intersectionCoord, std::ofstream &exportFile)
{
  // Convert GPS into relative positions in meters
  // r1 is base
  double lonLength = longitudeDegreeDistane(r1.p->lat); // Length of degree longitude in meters
  double latLength = latitudeDegreeDistance(r1.p->lat); // length of degree latitude in meters

  circle_t c1, c2, c3;
  c1 = {{0.0, 0.0}, r1.r};
  c2 = {{
        (r2.p->lon - r1.p->lon) * lonLength, 
        (r2.p->lat - r1.p->lat) * latLength},
         r2.r };
  c3 = {{
        (r3.p->lon - r1.p->lon) * lonLength,
        (r3.p->lat - r1.p->lat) * latLength}, 
        r3.r};
  // Calculate lines for circle intersections
  // l1 = intersection c1-c2
  // l2 = intersection c2-c3
  // l3 = intersection c1-c3
  line_t l1, l2, l3;

  point_t p1, p2, p3, p4, p5, p6;

  
  if (!circleCircleIntersection(c1, c2, &p1, &p2))
  {
    // TODO: handle error
  }
  if (!circleCircleIntersection(c2, c3, &p3, &p4))
  {
    // TODO: Handle error
  }
  if (!circleCircleIntersection(c1, c3, &p5, &p6))
  {
    // TODO: Handle error
  }

  // printPoint("p", p1);
  // printPoint("p", p2);
  // printPoint("p", p3);
  // printPoint("p", p4);
  // printPoint("p", p5);
  // printPoint("p", p6);
  
  // printPointAsCoord("its 1", *(r1.p), p1);
  // printPointAsCoord("its 2", *(r1.p), p2);

  // printPointAsCoord("its 3", *(r1.p), p3);
  // printPointAsCoord("its 4", *(r1.p), p4);

  // printPointAsCoord("its 5", *(r1.p), p5);
  // printPointAsCoord("its 6", *(r1.p), p6);

  point_t *points[6] = {&p1, &p2, &p3, &p4, &p5, &p6};
  point_t *closestPoints[3];
  double shortestDistance = 10000000;


  for (int a = 0; a < 6; a++)
  {
    for (int b = 0; b < 6; b++)
    {
      for (int c = 0; c < 6; c++)
      {
        if(a != b && b != c && a != c)
        {
          double distance = 
            PYTHAG(getDelta(points[a]->x, points[b]->x), getDelta(points[a]->y, points[b]->y)) +  // a-b
            PYTHAG(getDelta(points[b]->x, points[c]->x), getDelta(points[b]->y, points[c]->y)) +  // b-c
            PYTHAG(getDelta(points[a]->x, points[c]->x), getDelta(points[a]->y, points[c]->y));   // a-c

          if (distance < shortestDistance)
          {
            closestPoints[0] = points[a];
            closestPoints[1] = points[b];
            closestPoints[2] = points[c];
            shortestDistance = distance;
          }

        }
      }
    }
  }

  // calculate average of intersection
  // TODO: Handle average if no 3 points
  double xAvg, yAvg;
  for (size_t i = 0; i < 3; i++)
  {
    xAvg += closestPoints[i]->x;
    yAvg += closestPoints[i]->y;
  }

  xAvg /= 3;
  yAvg /= 3;  


  
  point_t trilaterationResult = {xAvg, yAvg};
  pointToCoord(*(r1.p), trilaterationResult, trilaterationCoord);

  l1 = {p1, (p2.y - p1.y) / (p2.x - p1.x) * -1};
  l2 = {p3, (p4.y - p3.y) / (p4.x - p3.x) * -1};
  l3 = {p5, (p6.y - p5.y) / (p6.x - p5.x)* -1};
  // At this point lines between intersections have been constructed

  point_t i1, i2, i3, iRelative;
  lineLineIntersect(l1, l2, &i1);
  lineLineIntersect(l2, l3, &i2);
  lineLineIntersect(l1, l3, &i3);

  point_t intersectionResult = {(i1.x + i2.x + i3.x) / 3 * -1, (i1.y + i2.y + i3.y) / 3};
  pointToCoord(*(r1.p), intersectionResult, intersectionCoord);

//
  coord_t x1, x2, x3, x4, x5, x6, x7, x8, x9;
  pointToCoord(*r1.p, p1, &x1);
  pointToCoord(*r1.p, p2, &x2);
  pointToCoord(*r1.p, p3, &x3);
  pointToCoord(*r1.p, p4, &x4);
  pointToCoord(*r1.p, p5, &x5);
  pointToCoord(*r1.p, p6, &x6);

  pointToCoord(*r1.p, *closestPoints[0], &x7);
  pointToCoord(*r1.p, *closestPoints[1], &x8);
  pointToCoord(*r1.p, *closestPoints[2], &x9);

  exportFile << "[intersections:point]\n";
  exportFile << "monitor_1-2_p1" << "=" << x1.lat << "," << x1.lon << "\n";
  exportFile << "monitor_1-2_p2" << "=" << x2.lat << "," << x2.lon << "\n";

  exportFile << "monitor_2-3_p1"  << "=" << x3.lat << "," << x3.lon << "\n";
  exportFile << "monitor_2-3_p2"  << "=" << x4.lat << "," << x4.lon << "\n";

  exportFile << "monitor_1-3_p1"  << "=" << x5.lat << "," << x5.lon << "\n";
  exportFile << "monitor_1-3_p2"  << "=" << x6.lat << "," << x6.lon << "\n";

  exportFile << "[trilaterate_intersections:point]\n";
  exportFile << "chosen_point_1"  << "=" << x7.lat << "," << x7.lon << "\n";
  exportFile << "chosen_point_2"  << "=" << x8.lat << "," << x8.lon << "\n";
  exportFile << "chosen_point_3"  << "=" << x9.lat << "," << x9.lon << "\n";

  exportFile << "[results:point]\n";
  exportFile << "trilaterate prediction"  << "=" << trilaterationCoord->lat << "," << trilaterationCoord->lon << "\n";
  exportFile << "intersect prediction"  << "=" << intersectionCoord->lat << "," << intersectionCoord->lon << "\n";

  exportFile << "[lines:line]\n";
  exportFile << "line_1_2"  << "=" << x1.lat << "," << x1.lon << " " << x2.lat << "," << x2.lon << "\n";
  exportFile << "line_2_3"  << "=" << x3.lat << "," << x3.lon << " " << x4.lat << "," << x4.lon << "\n";
  exportFile << "line_1_3"  << "=" << x5.lat << "," << x5.lon << " " << x6.lat << "," << x6.lon << "\n";

  exportFile << "[circles:circle]\n";
  circle_t circles[3] = {c1, c2, c3};
  int resolution = 50;
  for (int n = 0; n < 3; n++)
  {
    exportFile << "circle" << n+1 << "=";
    for (int i = 0; i < resolution + 1; i++)
    {
      circle_t c = circles[n];
      double angle = (2*M_PI / resolution) * i;
      double x, y;
      x = c.p.x + c.r * cos(angle);
      y = c.p.y + c.r * sin(angle);
      coord_t coordOfAngle;
      pointToCoord(*r1.p, {x, y}, &coordOfAngle);
      exportFile << coordOfAngle.lon << "," << coordOfAngle.lat << ",0 ";
    }
    exportFile << "\n";
  }

  // exportFile << "trilaterate"  << "=" << trilaterationCoord->lat << "," << trilaterationCoord->lon << "\n";
//
  return 1;
}

double getDistance(coord_t p1, coord_t p2)
{
    double res, xi, yi;

    double latLength = (latitudeDegreeDistance(p1.lat) + latitudeDegreeDistance(p2.lat)) / 2;
    double lonLength = (longitudeDegreeDistane(p1.lat) + longitudeDegreeDistane(p2.lat)) / 2;

    xi = (p1.lon - p2.lon) * lonLength;
    yi = (p1.lat - p2.lat) * latLength;
    
    return abs(PYTHAG(xi, yi));
}

int pointToCoord(coord_t base, point_t p, coord_t *res)
{
    double latLength = latitudeDegreeDistance(base.lat);
    double lonLength = longitudeDegreeDistane(base.lat);

    res->lat = p.y / latLength + base.lat; 
    res->lon = p.x / lonLength + base.lon; 

    return 1;
}

double getDelta(double a, double b)
{
  if (a < b)
    return abs(a - b);
  else
    return abs(b - a);
}

// THESE FUNCTIONS ARE FOR DEBUGGING ONLY
void printLine(const char *c, line_t line)
{
  double degree = atan(line.slope) * (180/M_PI) + 180;
	std::cout << "Line "<< c << "\nX = " << line.p.x << "\nY = " << line.p.y << "\nSlope = " << line.slope << 
    "(" << degree << ") (" << fmod(degree + 180.0, 360.0) << ")\n" <<
    "y = " << line.slope << "(x+"<< line.p.x << ")+" << line.p.y << "\n\n";
}

void printPoint(const char *c, point_t point)
{
  std::cout << "Point " << c << "\nX = " << point.x << "\nY = " << point.y << "\n\n";
}

void printCoord(const char *c, coord_t coord)
{
  std::cout << "Coord " << c << "\nLatitude = " << coord.lat << "\nLongitude = " << coord.lon << "\n\n";
}

void printCircle(const char *c, circle_t circle)
{
  std::cout << "Circle " << c << "\nX = " << circle.p.x << "\nY = " << circle.p.y << "\nRadius = " << circle.r << "\n\n";
}

void printPointAsCoord(const char *c, coord_t base, point_t point)
{
  double latLength, lonLength;
  latLength = latitudeDegreeDistance(base.lat);
  lonLength = longitudeDegreeDistane(base.lat);
  std::cout << "Point " << c << "\nLatitude = " << (point.y / latLength) + base.lat <<  "\nLongitude = " << (point.x / lonLength) + base.lon << "\n\n";
}