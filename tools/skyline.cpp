#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <cmath>
using namespace std;
ofstream outfile;

int skyline( unsigned int i, unsigned int p)
{
    if(i==p)
    {
        //printf("%c", p+'A');
        outfile<<char(p+'A');

        return 0;
    }
    skyline(i+1, p);
   // printf("%c", i+'A');
    outfile<<char(i+'A');
    skyline(i+1, p);
    return 0;
}
int main(int argc, char** argv)
{

    if(argc==1)
    {
        std::cout << "The file size is not specified.\n";
        std::cout << "The second parameter is the size of file in GiB." << std::endl;
        return 0;
    }

    std::cout << "The file size must be 2, 4, 8,..., 2^n, and the second parameter is 2, 4, 8,...,2^n.\n";

    unsigned int n= std::log(atoi(argv[1])) / std::log(2) + 30; // produce a skyline string of 2^n characters

    std::cout << "n = " << n << std::endl;

    std::cout << "The skyline file size = " << atoi(argv[1]) << "GiB" << std::endl;

    std::string fileName("skyline_");

    fileName += std::to_string(atoi(argv[1])) + "G";

    std::cout << "fileName = " << fileName << std::endl;

    outfile.open(fileName);
    skyline(1, n);
    //printf("A");
    outfile<<char('A');

    std::cout << "The file " << fileName << " is generated.\n";
    return 0;
    //return 0;
}

