
#ifndef _UINT40_H
#define _UINT40_H

#include <iostream>
#include"type.h"


namespace spt
{

#if _MSC_VER
#pragma pack(push, 1)
#endif
template <typename HighType>
class uint_pair
{
public:
    //! lower part type, always 32-bit
    typedef uint32 low_type;
    //! higher part type, currently either 8-bit or 16-bit
    typedef HighType high_type;

private:
    low_type low;

    high_type high;

    static uint_type low_max()
    {
        return std::numeric_limits<low_type>::max();
    }

    static const size_t low_bits = 8 * sizeof(low_type);

    static uint_type high_max()
    {
        return std::numeric_limits<high_type>::max();
    }

    static const size_t high_bits = 8 * sizeof(high_type);

public:
    static const size_t digits = low_bits + high_bits;

    static const size_t bytes = sizeof(low_type) + sizeof(high_type);

    inline uint_pair() {}

    inline uint_pair(const low_type& l, const high_type& h) : low(l), high(h) {}

    inline uint_pair(const uint_pair& a) : low(a.low), high(a.high) {}

    inline uint_pair(const uint32& a) : low(a), high(0) {}

    //! assumption: a >= 0
    inline uint_pair(const int32& a) : low(a), high(0) {}


    //!assumption: no overflow
    inline uint_pair(const uint64& a) : low((low_type)(a & low_max())), high((high_type)((a >> low_bits) & high_max())) {}

    //!assumption: a >= 0 and no overflow
    inline uint_pair(const int64& a) : low((low_type)(a & low_max())), high((high_type)((a >> low_bits) & high_max())) {}

    inline uint64 u64() const
    {
        return ((uint64)high) << low_bits | (uint64)low;
    }

    inline operator uint64 () const
    {
        //std::cout << "u64() executed.\n";
        return u64();
    }

    //add devision
    /*inline uint64 operator / (uint8 const &a) const
    {
        return u64() / a;
    }

    inline uint64 operator / (uint16 const &a) const
    {
        return u64() / a;
    }

    inline uint64 operator / (uint32 const &a) const
    {
        return u64() / a;
    }

    inline uint64 operator / (uint64 const &a) const
    {
        return u64() / a;
    }*/


    inline uint_pair& operator ++ ()
    {
        if (low == low_max())
            ++high, low = 0;
        else
            ++low;
        return *this;
    }

    inline uint_pair& operator -- ()
    {
        if (low == 0)
            --high, low = (low_type)low_max();
        else
            --low;
        return *this;
    }

    inline uint_pair& operator += (const uint_pair& b)
    {
        uint64 add = (uint64)low + b.low;
        low = (low_type)(add & low_max());
        high = (high_type)(high + b.high + ((add >> low_bits) & high_max()));
        return *this;
    }

    inline bool operator == (const uint_pair& b) const
    {
        return (low == b.low) && (high == b.high);
    }

    inline bool operator == (const uint32 & b) const
    {
        return !high && ( low == b);
    }

    inline bool operator == (const int32 & b) const
    {
        return !high && ( low == uint32(b) );
    }

    inline bool operator == (const uint64 & b) const
    {
        return ( u64() == b );
    }

    inline bool operator != (const uint_pair& b) const
    {
        return (low != b.low) || (high != b.high);
    }

    inline bool operator != (const uint32& b) const
    {
        //std::cout << "!= uint32 integer.\n";
        return high || (low != b);
    }

    inline bool operator != (const int32& b) const
    {
       // std::cout << "!= uint32 integer.\n";
        return high || (low != uint32(b));
    }

    inline bool operator != (const uint64 & b) const
    {
        //std::cout << "< uint64 integer.\n";
        return u64() != b;
    }

    inline bool operator < (const uint_pair& b) const
    {
        return (high < b.high) || (high == b.high && low < b.low);
    }

    inline bool operator < (const uint32 & b) const
    {
        //std::cout << "< uint32 integer.\n";
        return !high && low < b;
    }

    inline bool operator < (const uint64 & b) const
    {
        //std::cout << "< uint64 integer.\n";
        return u64() < b;
    }

    inline bool operator <= (const uint_pair& b) const
    {
        return (high < b.high) || (high == b.high && low <= b.low);
    }

    inline bool operator <= (const uint32 & b) const
    {
        //std::cout << " <= uint32 integer.\n";
        return !high && low <= b;
    }

    inline bool operator <= (const uint64 & b) const
    {
        //std::cout << " <= uint64 " << std::endl;

        return u64() <= b;
    }

    inline bool operator > (const uint_pair& b) const
    {
        return (high > b.high) || (high == b.high && low > b.low);
    }

    // !assumption: b >= 0
    inline bool operator > (const uint32 & b) const
    {
        //std::cout << "> uint32 integer.\n";
        return high || low > b;
    }

    inline bool operator > (const uint64 & b) const
    {
        //std::cout << "> uint64 integer.\n";
        return u64() > b;
    }

    // !assumption: b >= 0
    /*inline bool operator > (const int32 & b) const
    {

        std::cout << " uint40 can not compare to the int32.\n";

        std::cin.get();

        return (high > 0) || (high == 0 && low > b);
    }*/


    inline bool operator >= (const uint_pair& b) const
    {
        return (high > b.high) || (high == b.high && low >= b.low);
    }

    inline bool operator >= (const uint32 & b) const
    {
      //  std::cout << " >= uint32 " << std::endl;

        return high || low >= b;
    }

    inline bool operator >= (const uint64 & b) const
    {
    //    std::cout << " >= uint64 " << std::endl;

        return u64() >= b;
    }



    friend std::ostream& operator << (std::ostream& os, const uint_pair& a)
    {
        return os << a.u64();
    }

    static uint_pair min()
    {
        return uint_pair(std::numeric_limits<low_type>::min(),
                         std::numeric_limits<high_type>::min());
    }

    static uint_pair max()
    {
        return uint_pair(std::numeric_limits<low_type>::max(),
                         std::numeric_limits<high_type>::max());
    }
}
#if _MSC_VER
;
#pragma pack(pop)
#else
__attribute__((packed));
#endif

#if _MSC_VER
#pragma pack(push, 1)
#endif
template <typename HighType>
class uint_pair_24b
{
public:
    //! lower part type, always 32-bit
    typedef uint16 low_type;
    //! higher part type, currently either 8-bit or 16-bit
    typedef HighType high_type;

private:
    low_type low;

    high_type high;

    static uint_type low_max()
    {
        return std::numeric_limits<low_type>::max();
    }

    static const size_t low_bits = 8 * sizeof(low_type);

    static uint_type high_max()
    {
        return std::numeric_limits<high_type>::max();
    }

    static const size_t high_bits = 8 * sizeof(high_type);

public:
    static const size_t digits = low_bits + high_bits;

    static const size_t bytes = sizeof(low_type) + sizeof(high_type);

    inline uint_pair_24b() {}

    inline uint_pair_24b(const low_type& l, const high_type& h) : low(l), high(h) {}

    inline uint_pair_24b(const uint_pair_24b& a) : low(a.low), high(a.high) {}

//    inline uint_pair_24b(const uint32& a) : low(a), high(0) {}
    inline uint_pair_24b(const uint16& a) : low(a), high(0) {}

    //! assumption: a >= 0
    inline uint_pair_24b(const int16& a) : low(a), high(0) {}


    //!assumption: no overflow
    inline uint_pair_24b(const uint32& a) : low((low_type)(a & low_max())), high((high_type)((a >> low_bits) & high_max())) {}

    //!assumption: a >= 0 and no overflow
    inline uint_pair_24b(const int32& a) : low((low_type)(a & low_max())), high((high_type)((a >> low_bits) & high_max())) {}

        //!assumption: no overflow
    inline uint_pair_24b(const uint64& a) : low((low_type)(a & low_max())), high((high_type)((a >> low_bits) & high_max())) {}

    //!assumption: a >= 0 and no overflow
    inline uint_pair_24b(const int64& a) : low((low_type)(a & low_max())), high((high_type)((a >> low_bits) & high_max())) {}


    inline uint32 u32() const
    {
        return ((uint32)high) << low_bits | (uint32)low;
    }

    inline operator uint32 () const
    {
        //std::cout << "u64() executed.\n";
        return u32();
    }

    inline uint64 u64() const
    {
        return ((uint64)high) << low_bits | (uint64)low;
    }

    inline operator uint64 () const
    {
        //std::cout << "u64() executed.\n";
        return u64();
    }


    inline uint_pair_24b& operator ++ ()
    {
        if (low == low_max())
            ++high, low = 0;
        else
            ++low;
        return *this;
    }

    inline uint_pair_24b& operator -- ()
    {
        if (low == 0)
            --high, low = (low_type)low_max();
        else
            --low;
        return *this;
    }

   /* inline uint_pair_24b& operator - (const uint32 &b)
    {
        return uint_pair_24b(u32() - b);
    }*/


    inline uint_pair_24b& operator += (const uint_pair_24b& b)
    {
        uint32 add = (uint32)low + b.low;
        low = (low_type)(add & low_max());
        high = (high_type)(high + b.high + ((add >> low_bits) & high_max()));
        return *this;
    }

    inline bool operator == (const uint_pair_24b& b) const
    {
        return (low == b.low) && (high == b.high);
    }

    inline bool operator == (const uint16 & b) const
    {
        return (high == 0) && ( low == b);
    }

    inline bool operator == (const uint32 & b) const
    {
        return ( u32() == b );
    }

    inline bool operator != (const uint_pair_24b& b) const
    {
        return (low != b.low) || (high != b.high);
    }

    inline bool operator < (const uint_pair_24b& b) const
    {
        return (high < b.high) || (high == b.high && low < b.low);
    }

    inline bool operator <= (const uint_pair_24b& b) const
    {
        return (high < b.high) || (high == b.high && low <= b.low);
    }

    inline bool operator > (const uint_pair_24b& b) const
    {
        return (high > b.high) || (high == b.high && low > b.low);
    }

    inline bool operator >= (const uint_pair_24b& b) const
    {
        return (high > b.high) || (high == b.high && low >= b.low);
    }

    friend std::ostream& operator << (std::ostream& os, const uint_pair_24b& a)
    {
        return os << a.u32();
    }

    static uint_pair_24b min()
    {
        return uint_pair_24b(std::numeric_limits<low_type>::min(),
                         std::numeric_limits<high_type>::min());
    }

    static uint_pair_24b max()
    {
        return uint_pair_24b(std::numeric_limits<low_type>::max(),
                         std::numeric_limits<high_type>::max());
    }
}

#if _MSC_VER
;
#pragma pack(pop)
#else
__attribute__((packed));
#endif

}
using uint24 = spt::uint_pair_24b<uint8>;
using uint40 = spt::uint_pair<uint8>;
using uint48 = spt::uint_pair<uint16>;

namespace std
{

//! template class providing some numeric_limits fields for uint_pair types.
template <typename HighType>
class numeric_limits<spt::uint_pair<HighType> >
{
public:
    //! yes we have information about uint_pair
    static const bool is_specialized = true;

    //! return an uint_pair instance containing the smallest value possible
    static spt::uint_pair<HighType> min()
    {
        return spt::uint_pair<HighType>::min();
    }

    //! return an uint_pair instance containing the largest value possible
    static spt::uint_pair<HighType> max()
    {
        return spt::uint_pair<HighType>::max();
    }

    //! return an uint_pair instance containing the smallest value possible
    static spt::uint_pair<HighType> lowest()
    {
        return min();
    }

    //! unit_pair types are unsigned
    static const bool is_signed = false;

    //! uint_pair types are integers
    static const bool is_integer = true;

    //! unit_pair types contain exact integers
    static const bool is_exact = true;

    //! unit_pair radix is binary
    static const int radix = 2;

    //! number of binary digits (bits) in uint_pair
    static const int digits = spt::uint_pair<HighType>::digits;

    //! epsilon is zero
    static const spt::uint_pair<HighType> epsilon()
    {
        return spt::uint_pair<HighType>(0, 0);
    }

    //! rounding error is zero
    static const spt::uint_pair<HighType> round_error()
    {
        return spt::uint_pair<HighType>(0, 0);
    }

    //! no exponent
    static const int min_exponent = 0;

    //! no exponent
    static const int min_exponent10 = 0;

    //! no exponent
    static const int max_exponent = 0;

    //! no exponent
    static const int max_exponent10 = 0;

    //! no infinity
    static const bool has_infinity = false;
};

template <typename HighType>
class numeric_limits<spt::uint_pair_24b<HighType> >
{
public:
    //! yes we have information about uint_pair
    static const bool is_specialized = true;

    //! return an uint_pair instance containing the smallest value possible
    static spt::uint_pair_24b<HighType> min()
    {
        return spt::uint_pair_24b<HighType>::min();
    }

    //! return an uint_pair instance containing the largest value possible
    static spt::uint_pair_24b<HighType> max()
    {
        return spt::uint_pair_24b<HighType>::max();
    }

    //! return an uint_pair instance containing the smallest value possible
    static spt::uint_pair_24b<HighType> lowest()
    {
        return min();
    }

    //! unit_pair types are unsigned
    static const bool is_signed = false;

    //! uint_pair types are integers
    static const bool is_integer = true;

    //! unit_pair types contain exact integers
    static const bool is_exact = true;

    //! unit_pair radix is binary
    static const int radix = 2;

    //! number of binary digits (bits) in uint_pair
    static const int digits = spt::uint_pair_24b<HighType>::digits;

    //! epsilon is zero
    static const spt::uint_pair_24b<HighType> epsilon()
    {
        return spt::uint_pair_24b<HighType>(0, 0);
    }

    //! rounding error is zero
    static const spt::uint_pair_24b<HighType> round_error()
    {
        return spt::uint_pair_24b<HighType>(0, 0);
    }

    //! no exponent
    static const int min_exponent = 0;

    //! no exponent
    static const int min_exponent10 = 0;

    //! no exponent
    static const int max_exponent = 0;

    //! no exponent
    static const int max_exponent10 = 0;

    //! no infinity
    static const bool has_infinity = false;
};

} // namespace std

#endif
