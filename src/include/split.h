/******************************************************************************
 * include/split.h
 *
 * Implementation of function split mentioned on page 243 of the book:
 * Structured Parallel Programming.
 * by Michael McCool, Arch Robison, and James Reinders.
 * published by Morgan Kaufman, July 2012.
 *
 ******************************************************************************
 * Copyright (c) 2012 Michael McCool, Arch Robison, and James Reinders.
 * All rights reserved.
 *
 * The files in this directory and all subdirectories are licensed to you
 * according to the file LICENSE.txt (which resides in the same directory as
 * this file).
 *****************************************************************************/

#include <cstddef>

// Return greatest power of 2 less than m
size_t split( size_t m ) {
    size_t k=1; 
    while( 2*k<m ) k*=2;
    return k;   
}
