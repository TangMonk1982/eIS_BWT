
#ifndef _UINT24_H
#define _UINT24_H

#include <iostream>
#include"type.h"


namespace spt_24b
{

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

    inline uint32 u32() const
    {
        return ((uint32)high) << low_bits | (uint32)low;
    }

    inline operator uint32 () const
    {
        //std::cout << "u64() executed.\n";
        return u32();
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
//using uint24 = spt_24b::uint_pair_24b<uint8>;

namespace std
{

//! template class providing some numeric_limits fields for uint_pair types.
template <typename HighType>
class numeric_limits<spt_24b::uint_pair_24b<HighType> >
{
public:
    //! yes we have information about uint_pair
    static const bool is_specialized = true;

    //! return an uint_pair instance containing the smallest value possible
    static spt_24b::uint_pair_24b<HighType> min()
    {
        return spt_24b::uint_pair_24b<HighType>::min();
    }

    //! return an uint_pair instance containing the largest value possible
    static spt_24b::uint_pair_24b<HighType> max()
    {
        return spt_24b::uint_pair_24b<HighType>::max();
    }

    //! return an uint_pair instance containing the smallest value possible
    static spt_24b::uint_pair_24b<HighType> lowest()
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
    static const int digits = spt_24b::uint_pair_24b<HighType>::digits;

    //! epsilon is zero
    static const spt_24b::uint_pair_24b<HighType> epsilon()
    {
        return spt_24b::uint_pair_24b<HighType>(0, 0);
    }

    //! rounding error is zero
    static const spt_24b::uint_pair_24b<HighType> round_error()
    {
        return spt_24b::uint_pair_24b<HighType>(0, 0);
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
