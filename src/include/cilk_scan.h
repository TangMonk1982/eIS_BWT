/******************************************************************************
 * include/cilk_scan.h
 *
 * Top-level code for tiled parallel scan.
 *
 ******************************************************************************
 * Copyright (c) 2012 Michael McCool, Arch Robison, and James Reinders.
 * All rights reserved.
 *
 * The files in this directory and all subdirectories are licensed to you
 * according to the file LICENSE.txt (which resides in the same directory as
 * this file).
*****************************************************************************/

#include "split.h"
#include "upsweep.h"
#include "downsweep.h"
#include "temp_space.h"

template<typename T, typename R, typename C, typename S>
void cilk_scan(const size_t n, const T & initial, const size_t tilesize, const R reduce, const C combine, const S scan) {
    if( n>0 ) {
        size_t m = (n-1)/tilesize;
        temp_space<T> r(m+1);
        upsweep(0, m+1, tilesize, r.data(),
				n - m*tilesize, reduce, combine);
        downsweep(0, m+1, tilesize, r.data(), 
                n-m*tilesize, initial, combine, scan);
    }
}
