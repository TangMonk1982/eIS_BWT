/******************************************************************************
 * include/upsweep.h
 *
 * Upsweep phase for tiled parallel scan in Cilk Plus.
 *
 ******************************************************************************
 * Copyright (c) 2012 Michael McCool, Arch Robison, and James Reinders.
 * All rights reserved.
 *
 * The files in this directory and all subdirectories are licensed to you
 * according to the file LICENSE.txt (which resides in the same directory as
 * this file).
 *****************************************************************************/

template<typename T, typename R, typename C>
void upsweep(const size_t i, const size_t m, const size_t tilesize, T r[], const size_t lastsize, const R reduce, const C combine) {
    if( m==1 ) {
		reduce(i*tilesize, lastsize, r[0]);
	} else {
        size_t k = split(m);
        cilk_spawn upsweep(i, k, tilesize, r, tilesize, reduce, combine);
        upsweep(i+k, m-k, tilesize, r+k, lastsize, reduce, combine);
        cilk_sync;
        if( m==2*k)
			combine(r[m - 1], r[k - 1]);
	}
}