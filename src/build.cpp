/******************************************************************************
 * build.cpp
 *
 * This is the demo source code for the algorithm to build SA on external memory:
 *
 ******************************************************************************
 * Ling Bo Han <hanlb (at) mail2.sysu.edu.cn>
 * Yi Wu <wu.yi.christian (at) gmail.com>
 ******************************************************************************
*/

#include "common.h"
#include "builder4.h"
#include "utility.h"


int main(int argc, char** argv)
{

    if (argc < 2)
    {
        std::cerr << "one param required: input string name. The program is returned\n";
        return 0;
    }

    Timer::prog_beg_time = Timer::get_wall_time();

    std::cout << "The problem begins: ";

    Timer::show_time();

    std::cout << "The MAX memory is " << MAX_MEM / K_1024 << "M.\n";
    std::cout << "The IO buffer size is " << VEC_BUF_RAM / K_1024 << "M.\n";
    std::cout << "The physical file size is " << PHI_VEC_EM / K_1024 << "M.\n";

    // retrieve file name for input string
    std::string s_fname(argv[1]);
    std::string prog_name(argv[0]);
    std::cout << "The program name is " << prog_name << std::endl;
    std::cout << "The input file is " << s_fname << std::endl;

    if( argc == 3)
    {
        if(std::string(argv[2]) == "-bwtsa")  // build SA and BWT
        {
            is_getBWT = true;
            is_getSA = true;
            std::cout << "Build SA and BWT." << std::endl;
        }
        else if(std::string(argv[2]) == "-sa")   // build SA only
        {
            is_getBWT = false;
            is_getSA = true;
            std::cout << "Build SA only." << std::endl;
        }
        else if(std::string(argv[2]) == "-bwt")   // build SA only
        {
            is_getBWT = true;
            is_getSA = false;
            std::cout << "Build BWT only." << std::endl;
        }
        else
        {
            std::cout << "Build SA and BWT,  using 'build_bwt_sa input.txt -bwtsa' .\n";
            std::cout << "Build SA only,  using 'build_bwt_sa input.txt -sa' .\n";
            std::cout << "Build BWT only, using 'build_bwt_sa input.txt -bwt' .\n";

            return 0;
        }

    }
    else
    {
        std::cout << "Build SA and BWT,  using 'build_bwt_sa input.txt -bwtsa' .\n";
        std::cout << "Build SA only,  using 'build_bwt_sa input.txt -sa' .\n";
        std::cout << "Build BWT only, using 'build_bwt_sa input.txt -bwt' .\n";

        return 0;
    }


    MyVector<uint40> * sa_reverse;

    /// byte-alphabet
    uint8 alpha = 255;

    alphabetSet_vec.push_back(alpha);

    /// <alphabet_type, offset_type, compress_type, relative_offset_type>
    DSAComputation<uint8, uint40, uint8, uint32> dsa(s_fname, 0, sa_reverse, alpha);

    if(!dsa.get_build_way())dsa.run();


    double beg_time = Timer::get_wall_time();

    uint64 wrt_sa_iv = Logger::cur_iv;

    uint64 wrt_sa_ov = Logger::cur_ov;


    /// if the SA has been generated, write to disk.
    if(sa_reverse != nullptr && !sa_reverse->is_empty())
    {

        /// write the final SA to disk
        std::cout << "The sa_reverse size = " << sa_reverse->size() << std::endl;

        std::string sa_name = s_fname + ".sa5";

        FILE * fp = fopen(sa_name.c_str(),"wb");

        uint64 buf_size = 128 * K_1024 * sizeof(uint40);

        uint40 * sa_buf = new uint40[ buf_size ];

        uint64 index = 0;

        sa_reverse->start_read_reverse(); // start read reversely

        sa_reverse->next_remove_reverse(); // skip the sentinel

        while (!sa_reverse->is_eof())
        {
            sa_buf[index++] = sa_reverse->get_reverse();

            if (index == buf_size)
            {

                double bg = Timer::get_wall_time();
                fwrite(sa_buf,sizeof(uint40),index,fp);
                Timer::add_write_time(Timer::get_wall_time()-bg);

                Logger::addPDU(index * sizeof(uint40));

                Logger::addOV(index * sizeof(uint40));

                index = 0;
            }

            sa_reverse->next_remove_reverse();
        }

        double bg = Timer::get_wall_time();

        fwrite(sa_buf, sizeof(uint40), index, fp);

        Timer::add_write_time(Timer::get_wall_time() - bg);

        Logger::addPDU(index * sizeof(uint40));

        Logger::addOV(index * sizeof(uint40));

        fclose(fp);

        if(sa_reverse->is_eof())
        {
            delete sa_reverse;
            sa_reverse = nullptr;
        }

        delete sa_buf;
        sa_buf = nullptr;


    }
    else
    {
        std::cout << "â€˜sa_reverseâ€? is nullptr, and it is right.\n";
    }

    Timer::prog_end_time = Timer::get_wall_time();

    Timer::prog_total_time = Timer::prog_end_time - Timer::prog_beg_time;

    Timer::add_write_sa_time(Timer::get_wall_time() - beg_time);

    Logger::add_wrt_sa_iv(Logger::cur_iv - wrt_sa_iv);

    Logger::add_wrt_sa_ov(Logger::cur_ov - wrt_sa_ov);

    Logger::report(UtilityFunctions::getFileLength(s_fname));

    std::cout << "program is closed. " << std::endl;


}
