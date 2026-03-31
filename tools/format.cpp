// Author: Ling Bo Han,  Email: hanlb@mail2.sysu.edu.cn
// School of Data and Computer Science, Sun Yat-sen University,
// Guangzhou, China
// Date: December 20, 2019
//
// This is a tool to replace the alphabet of an input with a new alphabet that does not contain zero characters.


#include <iostream>
#include <fstream>
#include <vector>
#include <stdio.h>
#include<stdlib.h>
#include<string>

using uint64 = unsigned long long int;
using uint8 = unsigned char;


const uint64 K_512 = 512 * 1024;
const uint64 K_1024 = 1024 * 1024;

uint64 MAX_MEM = 3 * 1024 * K_1024;


int main(int argc, char** argv)
{


    if (argc != 2)
    {
        std::cerr << "Two parameters are required.\n";

        exit(-1);
    }

    std::cout << "Format the alphabet of an input with a non-zero alphabet.\n";

    std::string s_fname(argv[1]);

    std::cout << "The input file is " << s_fname << std::endl;


    std::ifstream fin(s_fname, std::ios_base::in | std::ios_base::binary);

    if (!fin.is_open())
    {
        std::cout << " Source file does not exist." << std::endl;
        return 0;
    }

    fin.seekg(0, fin.end);

    uint64 s_len = fin.tellg();

    std::cout << "File length = " << s_len << std::endl;

    uint8 * buf = new uint8[MAX_MEM];

    uint64 loopNum = s_len / MAX_MEM + ( (s_len % MAX_MEM) == 0 ? 0 : 1);

    uint8 * transTable = new uint8[256];

    for(uint64 i = 0; i < 256; i++)
        transTable[i] = 0;

    bool is_zero = false;

    for (uint64 i = 0, readSize = 0; i < loopNum; i++)
    {

        fin.seekg( MAX_MEM * i, fin.beg );

        readSize = (s_len - MAX_MEM * i) <= MAX_MEM ? s_len - MAX_MEM * i : MAX_MEM;

        fin.read((char*)buf,readSize);

        for (uint64 j = 0; j < readSize; j++)
        {

            if(buf[j] != 0)
                transTable[buf[j]] = buf[j];
            else
            {
                is_zero = true;
            }
        }

    }

    delete[] buf;


    for (uint64 i = 1, j = ( is_zero == true ? 1 : 0 ); i < 256; i++)
        if (transTable[i] != 0)
            transTable[i] = ++j;

    if(is_zero)
        transTable[0] = 1;

    // copy original file to one new file
    std::string s_fname_bak = s_fname + ".format";

    std::ofstream fot(s_fname_bak.c_str(), std::ios_base::out | std::ios_base::binary );

    uint64 buf_len = MAX_MEM / 2;  /// buffer length

    uint8 * buf_0 = new uint8[buf_len], * buf_1 = new uint8[buf_len];

    loopNum = s_len / buf_len + ((s_len % buf_len) == 0 ? 0 : 1);

    for (uint64 i = 0, readSize = 0; i < loopNum; i++)
    {

        fin.seekg( buf_len * i, fin.beg);

        readSize =  (s_len - buf_len * i) <= buf_len ? s_len - buf_len * i : buf_len;

        fin.read((char*)buf_0, readSize);

        for (uint64 j = 0; j < readSize; j++)
            buf_1[j] = transTable[buf_0[j]];

        fot.write((char *)buf_1,readSize);
    }

    fin.close(), fot.close();

    delete[] buf_0;
    delete[] buf_1;

    buf_0 = nullptr;
    buf_1 = nullptr;

    uint64 max_alpha = 0;// the character set number

    for (uint64 i = 255; i >= 1; --i)
    {

        if (transTable[i] != 0)
        {
            max_alpha = transTable[i];
            break;
        }

    }

    std::cout << " The maximum character in the formatted file is " << max_alpha <<  std::endl;

    delete [] transTable;


}
