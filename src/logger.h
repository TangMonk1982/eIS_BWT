/******************************************************************************
 * logger.h
 *
 ******************************************************************************
 * Ling Bo Han <hanlb (at) mail2.sysu.edu.cn>
 * Yi Wu <wu.yi.christian (at) gmail.com>
 * Ge Nong <issng (at) mail.sysu.edu.cn>
 ******************************************************************************
*/

#ifndef _LOGGER_H
#define _LOGGER_H

#include <iomanip>
#include "timer.h"



/// \brief a logger for recording pdu and iov
///
class Logger
{

public:
    static uint64 max_pdu; ///< maximum peak disk use

    static uint64 cur_pdu; ///< current peak disk use

    static uint64 cur_iv; ///< current input volume

    static uint64 cur_ov; ///< current output volume

    static bool is_max_pdu;

    /// 20190712
    static uint64 red_iv; /// input volume on reduce phase

    static uint64 red_ov; /// output volume on reduce phase

    static uint64 ind_iv; /// input volume on induce phase

    static uint64 ind_ov; /// output volume on induce phase

    static uint64 cpt_bwt_iv; /// input volume on computing bwt phase

    static uint64 cpt_bwt_ov; /// output volume on computing bwt phase

    static uint64 get_s1_iv; /// input volume on computing bwt phase

    static uint64 get_s1_ov; /// output volume on computing bwt phase

    static uint64 cpt_lms_pos_iv;

    static uint64 cpt_lms_pos_ov;

    static uint64 cpt_sa_iv;

    static uint64 cpt_sa_ov;

    static uint64 wrt_sa_iv;

    static uint64 wrt_sa_ov;

public:

    /// \brief increase pdu
    ///
    static void addPDU(const uint64 _delta)
    {

        cur_pdu += _delta;

        if (cur_pdu >= max_pdu)
            max_pdu = cur_pdu;


    }

    /// \brief decrease pdu
    static void delPDU(const uint64 _delta)
    {

        cur_pdu -= _delta;
    }

    /// \brief increase iv
    ///
    static void addIV(const uint64 _delta)
    {

        cur_iv += _delta;
    }

    /// \brief increase ov
    ///
    static void addOV(const uint64 _delta)
    {

        cur_ov += _delta;
    }

    static void cur_IOV()
    {

        std::cout << "cur_IV = " << cur_iv << std::endl;
        std::cout << "cur_OV = " << cur_ov << std::endl;
        std::cout << "cur_IOV = " << cur_iv + cur_ov << std::endl;

    }

    static uint64 report_IOV()
    {
        return cur_iv + cur_ov;
    }

    static void output_error(const char * _file, uint32 _line)
    {
        std::cout << "Error. File: " << _file << ", Line: " << _line << std::endl;

        std::cin.get();
    }

    static void pause(const char * _file, uint32 _line)
    {
        std::cout << " Pause. File: " << _file << ", Line: " << _line << std::endl;
        std::cin.get();
    }

    static void output_separator_line(std::string _s)
    {
        std::cout << "\n --------- ";
        std::cout << _s ;
        std::cout << "    --------- \n";
    }

    static void output_ram_use()
    {

//        FILE *fp;
//        char buffer[1024]= {0};
//        fp=popen("./ram_use.sh","r");
//        fread(buffer,1,sizeof(buffer),fp);
//        printf("%s",buffer);
//        pclose(fp);

    }

    static void report(const uint64 _corpora_size)
    {

        std::cerr << "--------------------------------------------------------------\n";

        std::cerr << "Statistics collection:\n";

        std::cerr << "peak disk use: " << max_pdu + _corpora_size << " bytes." << std::endl;

        std::cerr << "peak disk use (per char): " << double(max_pdu + _corpora_size) / double(_corpora_size) << std::endl;

        std::cerr << "--------------------------------------------------------------\n";

        std::cerr << "read volume: " << cur_iv << " bytes." << std::endl;

        std::cerr << "read volume (per char) "  << double(cur_iv) / double(_corpora_size) << std::endl;

        std::cerr << "write volume: " << cur_ov << " bytes." << std::endl;

        std::cerr << "write volume (per char)"  << double(cur_ov) / double(_corpora_size) << std::endl;

       /* std::cerr << "------------------------------------------------------------\n";

        std::cerr << "the time of computing bwt :" << Timer::bwt_time << ", about " << (Timer::bwt_time / Timer::prog_total_time) * 100 << "%." << std::endl;
        std::cerr << "the time of reduced sorting subs :" << Timer::reduce_sort_time << std::endl;
        std::cerr << "the time of induced sorting subs :" << Timer::induce_sort_time << std::endl;
        std::cerr << "the time of generating s1 :" << Timer::generate_s1_time << std::endl;
        std::cerr << "the time of computing LMS position :" << Timer::compute_lms_pos_time << std::endl;
        std::cerr << "the time of computing SA in RAM :" << Timer::sais_time << std::endl;
        std::cerr << "the time of writing SA to disk :" << Timer::write_sa_time << std::endl;

        std::cerr << "the total time of statistics time : " << (Timer::bwt_time + Timer::reduce_sort_time + Timer::induce_sort_time + Timer::generate_s1_time  + Timer::compute_lms_pos_time + Timer::sais_time + Timer::write_sa_time ) << std::endl;
        std::cerr << "------------------------------------------------------------\n";
        std::cerr << "the iv in reducing stage : " << Logger::red_iv / double(_corpora_size) << std::endl;
        std::cerr << "the ov in reducing stage : " << Logger::red_ov / double(_corpora_size) << std::endl;
        std::cerr << "the iv in inducing stage : " << Logger::ind_iv / double(_corpora_size) << std::endl;
        std::cerr << "the ov in inducing stage : " << Logger::ind_ov / double(_corpora_size) << std::endl;

        std::cerr << "the total IO volume on induced sorting stage : " << (Logger::red_iv + Logger::red_ov + Logger::ind_iv + Logger::ind_ov) / double(_corpora_size) << std::endl;
        std::cerr << "------------------------------------------------------------\n";
        std::cerr << "the iv volume for computing bwt : " << Logger::cpt_bwt_iv / double(_corpora_size) << std::endl;
        std::cerr << "the ov volume for computing bwt : " << Logger::cpt_bwt_ov / double(_corpora_size) << std::endl;
        std::cerr << "the total volume for computing bwt : " << (Logger::cpt_bwt_ov + Logger::cpt_bwt_iv)/ double(_corpora_size) << " , about " << double(Logger::cpt_bwt_ov + Logger::cpt_bwt_iv) / (cur_iv + cur_ov) * 100 << "%" << std::endl;

        std::cerr << "------------------------------------------------------------\n";
        std::cerr << "the iv volume for generating s1 : " << Logger::get_s1_iv / double(_corpora_size) << std::endl;
        std::cerr << "the ov volume for generating s1 : " << Logger::get_s1_ov / double(_corpora_size) << std::endl;
        std::cerr << "the total volume for generating s1 : " << (Logger::get_s1_iv + Logger::get_s1_ov)/ double(_corpora_size) << std::endl;
        std::cerr << "------------------------------------------------------------\n";
        std::cerr << "the iv volume for computing lms pos : " << Logger::cpt_lms_pos_iv / double(_corpora_size) << std::endl;
        std::cerr << "the ov volume for computing lms pos : " << Logger::cpt_lms_pos_ov / double(_corpora_size) << std::endl;
        std::cerr << "the total volume for computing lms pos : " << (Logger::cpt_lms_pos_iv + Logger::cpt_lms_pos_ov)/ double(_corpora_size) << std::endl;
        std::cerr << "------------------------------------------------------------\n";
        std::cerr << "the iv volume for computing sa in RAM : " << Logger::cpt_sa_iv / double(_corpora_size) << std::endl;
        std::cerr << "the ov volume for computing sa in RAM : " << Logger::cpt_sa_ov / double(_corpora_size) << std::endl;
        std::cerr << "the total volume for computing sa in RAM : " << (Logger::cpt_sa_iv + Logger::cpt_sa_ov)/ double(_corpora_size) << std::endl;
        std::cerr << "------------------------------------------------------------\n";
        std::cerr << "the iv volume for writing sa : " << Logger::wrt_sa_iv / double(_corpora_size) << std::endl;
        std::cerr << "the ov volume for writing sa : " << Logger::wrt_sa_ov / double(_corpora_size) << std::endl;
        std::cerr << "the total volume for writing sa : " << (Logger::wrt_sa_iv + Logger::wrt_sa_ov)/ double(_corpora_size) << std::endl;*/
        std::cerr << "total IOV: " << (Logger::cur_iv + Logger::cur_ov) / double(_corpora_size) << std::endl;

        std::cerr << "------------------------------------------------------------\n";

        std::cerr << "total time: " << Timer::prog_total_time << "S" << std::endl;

        std::cerr << "------------------------------------------------------------\n";


    }

    static void time_report()
    {

        std::cout << "program running time = " << Timer::get_wall_time() - Timer::prog_beg_time << std::endl;

        std::cerr << "the time of computing bwt :" << Timer::bwt_time << std::endl;
        std::cerr << "the time of reduced sorting subs :" << Timer::reduce_sort_time << std::endl;
        std::cerr << "the time of induced sorting subs :" << Timer::induce_sort_time << std::endl;
        std::cerr << "the time of generating s1 :" << Timer::generate_s1_time << std::endl;
        std::cerr << "the time of computing LMS position :" << Timer::compute_lms_pos_time << std::endl;
        std::cerr << "the time of computing SA in RAM :" << Timer::sais_time << std::endl;
        std::cerr << "the time of writing SA to disk :" << Timer::write_sa_time << std::endl;

        std::cerr << "the total time of statistics time : " << (Timer::bwt_time + Timer::reduce_sort_time + Timer::induce_sort_time + Timer::generate_s1_time  + Timer::compute_lms_pos_time + Timer::sais_time + Timer::write_sa_time ) << std::endl;



    }

    /// 20190712
    static void add_red_iv(const uint64 _delta)
    {

        red_iv += _delta;
    }

    static void add_red_ov(const uint64 _delta)
    {

        red_ov += _delta;
    }

    static void add_ind_iv(const uint64 _delta)
    {

        ind_iv += _delta;
    }

    static void add_ind_ov(const uint64 _delta)
    {

        ind_ov += _delta;
    }

    static void add_cpt_bwt_iv(const uint64 _delta)
    {

        cpt_bwt_iv += _delta;
    }

    static void add_cpt_bwt_ov(const uint64 _delta)
    {

        cpt_bwt_ov += _delta;
    }

    static void add_get_s1_iv(const uint64 _delta)
    {

        get_s1_iv += _delta;
    }

    static void add_get_s1_ov(const uint64 _delta)
    {

        get_s1_ov += _delta;
    }

    static void add_cpt_lms_pos_iv(const uint64 _delta)
    {

        cpt_lms_pos_iv += _delta;
    }

    static void add_cpt_lms_pos_ov(const uint64 _delta)
    {

        cpt_lms_pos_ov += _delta;
    }

    static void add_cpt_sa_iv(const uint64 _delta)
    {

        cpt_sa_iv += _delta;
    }

    static void add_cpt_sa_ov(const uint64 _delta)
    {

        cpt_sa_ov += _delta;
    }

    static void add_wrt_sa_iv(const uint64 _delta)
    {

        wrt_sa_iv += _delta;
    }

    static void add_wrt_sa_ov(const uint64 _delta)
    {

        wrt_sa_ov += _delta;
    }

};


uint64 Logger::max_pdu = 0;

uint64 Logger::cur_pdu = 0;

uint64 Logger::cur_iv = 0;

uint64 Logger::cur_ov = 0;

bool Logger::is_max_pdu = false;

/// 20190712
uint64 Logger::red_ov = 0;

uint64 Logger::red_iv = 0;

uint64 Logger::ind_ov = 0;

uint64 Logger::ind_iv = 0;

uint64 Logger::cpt_bwt_iv = 0;

uint64 Logger::cpt_bwt_ov = 0;

uint64 Logger::get_s1_iv = 0;

uint64 Logger::get_s1_ov = 0;

uint64 Logger::cpt_lms_pos_iv = 0;

uint64 Logger::cpt_lms_pos_ov = 0;

uint64 Logger::cpt_sa_iv = 0;

uint64 Logger::cpt_sa_ov = 0;

uint64 Logger::wrt_sa_iv = 0;

uint64 Logger::wrt_sa_ov = 0;

#endif
