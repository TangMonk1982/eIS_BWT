#ifndef TIMER_H
#define TIMER_H

#include<sys/time.h>
#include<time.h>
#include<iostream>

class Timer
{
public:

    static double prog_beg_time;

    static double prog_mid_time;

    static double prog_end_time;

    static double prog_total_time; /// 20190811

    static double time_level_0;

    static double read_time;

    static double write_time;


    static double bwt_time;

    static double induce_sort_time;

    static double reduce_sort_time;

    static double generate_s1_time;

    static double compute_lms_pos_time;

    static double sais_time;

    static double write_sa_time;

    /**< get the wall time in linux system. */
    static double get_wall_time()
    {
        struct timeval time;
        if (gettimeofday(&time, NULL))
        {
            return 0;
        }
        return (double)time.tv_sec + (double)time.tv_usec * .000001;
    }

    static void add_read_time(double _time) {

		read_time += _time;
	}

	static void add_write_time(double _time) {

		write_time += _time;
	}

	static void add_bwt_time(double _time) {

		bwt_time += _time;
	}

	static void add_induce_sort_time(double _time) {

		induce_sort_time += _time;
	}

	static void add_reduce_sort_time(double _time) {

		reduce_sort_time += _time;
	}

	static void add_generate_s1_time(double _time) {

		generate_s1_time += _time;
	}

	static void add_compute_lms_pos_time(double _time) {

		compute_lms_pos_time += _time;
	}

	static void add_sais_time(double _time) {

		sais_time += _time;
	}


	static void add_write_sa_time(double _time) {

		write_sa_time += _time;
	}

    static void outputTime(){

    	time_t timep;
    	time (&timep);
    	char tmp[64];
    	strftime(tmp, sizeof(tmp), "%Y-%m-%d %H:%M:%S",localtime(&timep) );
    	std::cout << "The current time = " << tmp << std::endl;

    }

    static void show_time(){

        outputTime();
    }
};


double Timer::prog_beg_time = 0;
double Timer::prog_mid_time = 0;
double Timer::prog_end_time = 0;
double Timer::prog_total_time = 0;

double Timer::time_level_0= 0;
double Timer::read_time= 0;
double Timer::write_time= 0;

double Timer::bwt_time= 0;

double Timer::induce_sort_time = 0;

double Timer::reduce_sort_time = 0;

double Timer::generate_s1_time = 0;

double Timer::compute_lms_pos_time = 0;

double Timer::sais_time = 0;

double Timer::write_sa_time = 0;

#endif // TIMER_H
