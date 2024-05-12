#include <trilaterate.hpp>
#include <iostream>
#include <stdlib.h>
#include <time.h>       
#include <math.h>
#include <fstream>


int main()
{
    double realLength1, realLength2, realLength3;
    coord_t trilaterationCoord, intersectionCoord;

    coord_t monitor1 = {53.02315707965995, 6.051166966717045};
    coord_t monitor2 = {53.04565562998146, 6.051170158444927};
    coord_t monitor3 = {53.03045156041614, 6.087672338575723};
    coord_t hive = {53.04334815287351, 6.078294935686204};

    realLength1 = getDistance(monitor1, hive) * 1.05;
    realLength2 = getDistance(monitor2, hive) * 0.9;
    realLength3 = getDistance(monitor3, hive) * 0.9;
    
    std::cout << "Length to hive " << "\nMonitor 1: " << realLength1 << "\nMonitor 2: " << realLength2 << "\nMonitor 3: " << realLength3 << "\n\n";

    record_t r1 = {&monitor1, realLength1};
    record_t r2 = {&monitor2, realLength2};
    record_t r3 = {&monitor3, realLength3};
    trilaterate(r1, r2, r3, &trilaterationCoord, &intersectionCoord);

    printCoord("Trilaterate estimation", trilaterationCoord);
    printCoord("Intersection estimation", intersectionCoord);

    #define TEST_WITH_SCALE
    #ifdef TEST_WITH_SCALE
    double newLength1, newLength2, newLength3;
    double start = 0.8;
    double step = 0.0004;
    double n = 1000;
    double scale = start;

    std::ofstream trilaterateFile("trilaterate.txt");
    std::ofstream intersectionFile("intersection.txt");

    std::cout << "Testing " << n << " Points with length scaled from " << start << " up too " << scale + (n * step) << "\n";
    for (int i = 0; i < n; i++)
    {
        scale += step;
        r1.r = realLength1 * scale;
        r2.r = realLength2 * scale;
        r3.r = realLength3 * scale;

        trilaterate(r1, r2, r3, &trilaterationCoord, &intersectionCoord);
        trilaterateFile << trilaterationCoord.lat << " " << trilaterationCoord.lon << "\n";
        intersectionFile << intersectionCoord.lat << " " << intersectionCoord.lon << "\n";
    }
    trilaterateFile.close();
    intersectionFile.close();

    std::cout << "Coordinates written to trilaterate.txt and intersection.txt";

    #endif
	return 0;
}


