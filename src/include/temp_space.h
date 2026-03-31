/******************************************************************************
 * include/temp_space.h
 *
 * Type for scratch space. It is optimized for allocating short arrays of a
 * type T with a trivial constructor and destructor.
 *
 ******************************************************************************
 * Copyright (c) 2012 Michael McCool, Arch Robison, and James Reinders.
 * All rights reserved.
 *
 * The files in this directory and all subdirectories are licensed to you
 * according to the file LICENSE.txt (which resides in the same directory as
 * this file).
 *****************************************************************************/

template<typename T>
class temp_space {
    static const size_t n = 4096/sizeof(T);
    T temp[n];
    T* base;
public:
    T* data() {return base;}
    T& operator[]( size_t k ) {return base[k];}
    temp_space( size_t size ) {
        base = size<=n ? temp : new T[size];
    }
    ~temp_space() {
        if( base!=temp ) 
            delete[] base;
    }
};