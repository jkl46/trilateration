#include <math.h>
#include <trilaterate.hpp>

#include <iostream>
#include <algorithm>

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
  double m, x, y;

  if (l1.slope == l2.slope) // No cross section
    return 0;

  // xi = ((l1.slope * l1.p.x + (l1.p.y * -1 )) - (l2.slope * l2.p.x + (l2.p.y * -1 ))) / (l1.slope - l2.slope);
  // yi = l1.slope*(xi + (l1.p.x * -1) + l1.p.y);


  x =  (l1.p.x * l1.slope - l2.p.x * l2.slope ) / (l1.slope - l2.slope) + (l1.p.y - l2.p.y);
  y = (x + l1.p.x) * l1.slope + l1.p.y;

  p->x = x;
  p->y = y;
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

int trilaterate(record_t r1, record_t r2, record_t r3, coord_t *p)
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
            PYTHAG(abs(abs(points[a]->x) - abs(points[b]->x)), abs(abs(points[a]->y) - abs(points[b]->y))) + // a-b
            PYTHAG(abs(abs(points[b]->x) - abs(points[c]->x)), abs(abs(points[b]->y) - abs(points[c]->y))) + // b-c 
            PYTHAG(abs(abs(points[a]->x) - abs(points[c]->x)), abs(abs(points[a]->y) - abs(points[c]->y))); // a-c
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

  point_t result = {xAvg, yAvg};
  pointToCoord(*(r1.p), result, p);

  // SUCCES
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

// THESE FUNCTIONS ARE FOR DEBUGGING ONLY
void printLine(const char *c, line_t line)
{
	std::cout << "Line "<< c << "\nX = " << line.p.x << "\nY = " << line.p.y << "\nSlope = " << line.slope << "\n\n";
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