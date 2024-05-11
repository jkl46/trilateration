#include <trilaterate.hpp>
#include <iostream>
#include <stdlib.h>
#include <time.h>       
#include <math.h>

int main()
{
    double realLength1, realLength2, realLength3;
    coord_t result;

    coord_t monitor1 = {52.47422209267754, 6.302148212997174};
    coord_t monitor2 = {52.47694409402182, 6.29636923721738};
    coord_t monitor3 = {52.47890770097344, 6.306102581303577};

    coord_t hive = {52.4766533222011, 6.302939665868415};

    realLength1 = getDistance(monitor1, hive);
    realLength2 = getDistance(monitor2, hive);
    realLength3 = getDistance(monitor3, hive);
    
    std::cout << "Length to hive " << "\nMonitor 1: " << realLength1 << "\nMonitor 2: " << realLength2 << "\nMonitor 3: " << realLength3 << "\n\n";

    double scale = 1.0;

    record_t r1 = {&monitor1, realLength1 * 1.1};
    record_t r2 = {&monitor2, realLength2 * 1.1};
    record_t r3 = {&monitor3, realLength2 * 0.9};

    trilaterate(r1, r2, r3, &result);

    printCoord("result", result);
    
	return 0;
}

