#include <trilaterate.hpp>
#include <iostream>

int main(int argc, char **argv)
{
    double realLength1, realLength2, realLength3;
    coord_t trilaterationCoord, intersectionCoord;
    coord_t monitor1, monitor2, monitor3, hive;

    monitor1 = {52.42209520256235,6.384834648239419};
    monitor2 = {52.42401961969845,6.384311413786394};
    monitor3 = {52.41931876736916,6.382396416741736};
    hive = {52.42238050706415,6.379774525289168};

    double scale = 1.0;

    realLength1 = getDistance(monitor1, hive) * scale;
    realLength2 = getDistance(monitor2, hive) * scale;
    realLength3 = getDistance(monitor3, hive) * scale;
    
    record_t r1 = {&monitor1, realLength1};
    record_t r2 = {&monitor2, realLength2};
    record_t r3 = {&monitor3, realLength3};

    
    trilaterate(r1, r2, r3, &trilaterationCoord, &intersectionCoord);

    std::cout << trilaterationCoord.lat;
    std::cout << trilaterationCoord.lon;
}

