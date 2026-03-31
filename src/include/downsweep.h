/******************************************************************************
 * include/downsweep.h
 *
 * Downsweep phase for tiled parallel scan in Cilk Plus.
 *
 ******************************************************************************
 * Copyright (c) 2012 Michael McCool, Arch Robison, and James Reinders.
 * All rights reserved.
 *
 * The files in this directory and all subdirectories are licensed to you
 * according to the file LICENSE.txt (which resides in the same directory as
 * this file).
 *****************************************************************************/

template<typename T, typename C, typename S>
void downsweep(const size_t i, const size_t m, const size_t tilesize, const T r[], const size_t lastsize, T initial, const C combine, const S scan) {
    if( m==1 ) {
        scan(i*tilesize, lastsize, initial);
    } else {
        size_t k = split(m);
        cilk_spawn downsweep(i, k, tilesize, r, tilesize, initial, combine, scan);
		combine(initial, r[k - 1]);
		downsweep(i + k, m - k, tilesize, r + k, lastsize, initial, combine, scan);
	}
}