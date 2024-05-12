#include <trilaterate.hpp>
#include <iostream>
#include <stdlib.h>
#include <time.h>       
#include <math.h>

int main()
{
    double realLength1, realLength2, realLength3;
    coord_t result;
    coord_t monitor1 = {52.42361031942352, 6.26869432090041};
    coord_t monitor2 = {52.44489551955696, 6.282622149172687};
    coord_t monitor3 = {52.44134150054891, 6.249388194168882};
    
    coord_t hive = {52.45497912284191, 6.2585747806824};
    realLength1 = getDistance(monitor1, hive);
    realLength2 = getDistance(monitor2, hive);
    realLength3 = getDistance(monitor3, hive);
    
    std::cout << "Length to hive " << "\nMonitor 1: " << realLength1 << "\nMonitor 2: " << realLength2 << "\nMonitor 3: " << realLength3 << "\n\n";


    record_t r1 = {&monitor1, realLength1};
    record_t r2 = {&monitor2, realLength2};
    record_t r3 = {&monitor3, realLength3};

    trilaterate(r1, r2, r3, &result);

    printCoord("result", result);
    
	return 0;
}

