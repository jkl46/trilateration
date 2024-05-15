#include <trilaterate.hpp>
#include <iostream>
#include <stdlib.h>
#include <time.h>       
#include <math.h>
#include <fstream>
#include <string>


int main(int argc, char **argv)
{
//
    std::ofstream exportFile("trilaterate.ini");
//

    double realLength1, realLength2, realLength3;
    coord_t trilaterationCoord, intersectionCoord;
    coord_t monitor1, monitor2, monitor3, hive;

    if (argc == 9)
    {
        #define stod std::stod
        monitor1 = {stod(argv[1]),stod(argv[2])};
        monitor2 = {stod(argv[3]),stod(argv[4])};
        monitor3 = {stod(argv[5]),stod(argv[6])};
        hive = {stod(argv[7]),stod(argv[8])};

        #pragma endregion
    } else
    {
        LOG("No valid input provided! using sample coordinates")
        monitor1 = {52.42209520256235,6.384834648239419};
        monitor2 = {52.42401961969845,6.384311413786394};
        monitor3 = {52.41931876736916,6.382396416741736};
        hive = {52.42238050706415,6.379774525289168};
    }

    double scale = 1.0;

    realLength1 = getDistance(monitor1, hive) * scale;
    realLength2 = getDistance(monitor2, hive) * scale;
    realLength3 = getDistance(monitor3, hive) * scale;
    
    std::cout << "length to hive " << "\nMonitor 1: " << realLength1 << "\nMonitor 2: " << realLength2 << "\nMonitor 3: " << realLength3 << "\n\n";

    record_t r1 = {&monitor1, realLength1};
    record_t r2 = {&monitor2, realLength2};
    record_t r3 = {&monitor3, realLength3};


    
    trilaterate(r1, r2, r3, &trilaterationCoord, &intersectionCoord, exportFile);


    printCoord("Trilaterate estimation", trilaterationCoord);
    printCoord("Intersection estimation", intersectionCoord);

//
    exportFile << "[main_points:point]\n";
    exportFile << "monitor1" << "=" << monitor1.lat << "," << monitor1.lon << "\n";
    exportFile << "monitor2" << "=" << monitor2.lat << "," << monitor2.lon << "\n";
    exportFile << "monitor3" << "=" << monitor3.lat << "," << monitor3.lon << "\n";
    exportFile << "hive" << "=" << hive.lat << "," << hive.lon << "\n";
    exportFile.close();
//
    // #define TEST_WITH_SCALE
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


