#ifndef _TYPE_H
#define _TYPE_H


using int8 = char;
using uint8 = unsigned char;
using int16 = short;
using uint16 = unsigned short;
using int32 = int;
using uint32 = unsigned int;
using int64 = long long int;
using uint64 = unsigned long long int;

//using string = string;
constexpr int32 my_pointer_size = sizeof(void *);

template<int32 ptrSize>
struct choose_int_types {};

template<>
struct choose_int_types<4> { //32-bit
	using int_type = int32;
	using unsigned_type = uint32;
};

template<>
struct choose_int_types<8> { //64-bit
	using int_type = int64;
	using unsigned_type = uint64;
};

using int_type = choose_int_types<my_pointer_size>::int_type;
using uint_type = choose_int_types<my_pointer_size>::unsigned_type;




#endif // !_TYPE_H#pragma once
