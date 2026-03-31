/******************************************************************************
 * utility.h
 *
 ******************************************************************************
 * Ling Bo Han <hanlb (at) mail2.sysu.edu.cn>
 * Yi Wu <wu.yi.christian (at) gmail.com>
 * Ge Nong <issng (at) mail.sysu.edu.cn>
 ******************************************************************************
*/

#ifndef MY_UTILITY_H
#define MY_UTILITY_H

#include <algorithm>    // std::copy
#include <exception>

#include "common.h"
#include "tuple.h"
#include "tuple_sorter.h"
#include "logger.h"
#include "include/cilk_scan.h"


#define COPRESS_RATION

#ifdef CILK
#define parallel_for cilk_for
#else
#define parallel_for for
#endif


// set parallel
const bool is_parallel = true;

//#define DEBUG
static const uint8 MASK[] = { 0x80, 0x40, 0x20, 0x10, 0x08, 0x04, 0x02, 0x01 };

static const uint8 MASK2[] = { 0xC0, 0x30, 0x0C, 0x03 };

static const uint8 SHIFT2[] = { 6,4,2,0 };

const unsigned int cs_buf = 100 * 1024, cs_grain = 1024;

//using namespace backward;

// Pipeline for range [begin,end) in either direction.
template<typename S, typename F>
void cilk_pipeline(size_t begin, size_t end, size_t tilesize, S source, F function)
{
    bool rev = (begin > end);
    size_t n = rev ? (begin - end) : (end - begin);
    if (n == 0) return;
    size_t m = (n - 1) / tilesize, first, last;
    if (rev)
    {
        for (size_t k = 0; k <= m; k++)
        {
            first = begin - k * tilesize;
            last = (k == m) ? end : (first - tilesize);
            source(first, last, k);
            cilk_sync;
            cilk_spawn function(first, last, k);
        }
    }
    else
    {
        for (size_t k = 0; k <= m; k++)
        {
            first = begin + k * tilesize;
            last = (k == m) ? end : (first + tilesize);
            source(first, last, k);
            cilk_sync;
            cilk_spawn function(first, last, k);
        }
    }
}



/// \brief data wrapper for a bit array in RAM
///
class BitWrapper
{

private:

    char * data;

public:

    BitWrapper(char * _data = nullptr) : data(_data) {}

    bool get(size_t _idx)
    {

        return  (data[_idx / 8] & MASK[_idx % 8]) ? 1 : 0;
    }

    void set(size_t _idx, bool _val)
    {

        data[_idx / 8] = _val ? (MASK[_idx % 8] | data[_idx / 8]) : ((~MASK[_idx % 8]) & data[_idx / 8]);
    }
    ~BitWrapper()
    {
        data = nullptr;
    }
};

class BitWrapper_uint8
{

private:

    uint8 * data;

public:

    BitWrapper_uint8(uint8 * _data = nullptr) : data(_data) {}

    bool get(uint32 _idx)
    {

        return  (data[_idx / 8] & MASK[_idx % 8]) ? 1 : 0;
    }

    void set(uint32 _idx, bool _val)
    {

        data[_idx / 8] = _val ? (MASK[_idx % 8] | data[_idx / 8]) : ((~MASK[_idx % 8]) & data[_idx / 8]);
    }
    ~BitWrapper_uint8()
    {
        data = nullptr;
    }
};


#if _MSC_VER
#pragma pack(push, 1)
#endif
template<typename offset_type>
struct Alpha_Block
{

public:

    offset_type m_beg_alpha;

    offset_type m_end_alpha;

    offset_type m_size;

    Alpha_Block():m_beg_alpha(0),m_end_alpha(0),m_size(0) {}

    Alpha_Block(offset_type _beg, offset_type _end, offset_type _size):m_beg_alpha(_beg),m_end_alpha(_end),m_size(_size) {}

    void display()
    {

        std::cout << " Alpha_Block (" << m_beg_alpha << "," << m_end_alpha << "," << m_size << ").\n";
    }

}
#if _MSC_VER
;
#pragma pack(pop)
#else
__attribute__((packed));
#endif



/// \brief data warpper for a 2-bit array in RAM
///
class Bit2Wrapper
{
private:

    char * data;

public:

    Bit2Wrapper(char* _data = nullptr) : data(_data) {}

    uint8 get(size_t _idx)
    {

        return (data[_idx / 4] >> SHIFT2[_idx % 4]) & MASK2[3];
    }

    void set(size_t _idx, uint8 _val)
    {

        data[_idx / 4] = (_val << SHIFT2[_idx % 4]) | (data[_idx / 4] & ~MASK2[_idx % 4]);
    }
};



//template<typename alphabet_type>
class UtilityFunctions
{

public:

    /// \brief increase pdu
    ///
    static uint64 getFileLength(const std::string & _fname)
    {

        FILE * fp = fopen(_fname.c_str(), "rb");
        if (!fp)
        {
            std::cout << _fname << " is not exist!" << std::endl;
            return 0;
        }
        fseek(fp, 0, SEEK_END);
        uint64 size = ftell(fp);
        fclose(fp);
        return size;
    }

    template<typename alphabet_type>
    static void fileRead(alphabet_type * &_buf, const std::string & _fname, const uint64 & _offset, const uint64 & _size)
    {

        FILE * fp = fopen(_fname.c_str(), "rb");
        if (!fp)
        {
            std::cout << _fname << " is not exist!" << std::endl;
            return;
        }
        /**< _offset should be multiplied by the sizeof(alphabet_type) */
        fseek(fp, _offset * sizeof(alphabet_type), SEEK_SET);
        double bg = Timer::get_wall_time();
        fread(_buf, sizeof(alphabet_type), _size, fp);
        Timer::add_read_time(Timer::get_wall_time() - bg);
        fclose(fp);

        Logger::addIV(_size * sizeof(alphabet_type));

        return;
    }


    /// \note reduction phase
    /// alphabet <= 255
    template<typename alphabet_type, typename compress_type, typename relative_offset_type>
    static void th_reduceSort_bwt_1B( std::string  _fname
                                      , const uint64  _beg
                                      , const uint64  _end
                                      , const alphabet_type  _alpha
                                      , const uint32  _level
                                      , const uint32 _id
                                      , std::vector< std::vector<uint64> *> & _LMS_counter
                                      , MyVector< Pair<alphabet_type, compress_type> > *& _sub_l_bwt_vector
                                      , BitVectorWapper<uint64>                        *& _sub_l_bit_vector
                                      , MyVector< Pair<alphabet_type, compress_type> > *& _sub_s_bwt_vector
                                      , BitVectorWapper<uint64>                        *& _sub_s_bit_vector
                                      , MyVector<relative_offset_type>                 *& _lms_rPos_vector
                                      , BitVectorWapper<uint64>                        *& _lms_rPos_diff_bit_vector
                                      , bool _is_single
                                      , bool _is_sentinel
                                    )
    {

        if(_is_single)
        {
            std::cout << "this block is single block.\n";

            compute_singleton_bwt<alphabet_type, compress_type, relative_offset_type>(
                _fname
                , _beg
                , _end
                , _alpha
                , _level
                , _id
                , _LMS_counter
                , nullptr
                , _sub_l_bwt_vector
                , _sub_l_bit_vector
                , _sub_s_bwt_vector
                , _sub_s_bit_vector
                , _lms_rPos_vector
                , _lms_rPos_diff_bit_vector
                , _is_sentinel
            );


            return ;
        }

        uint64 i, j;

        uint64 _sLen = _end - _beg + 1;

        alphabet_type * _s = nullptr;

        if (_is_sentinel)
        {
            _sLen++;
            _s = new alphabet_type[_sLen];///\note memory leak: Invalid write of size 1
            _s[_sLen - 1] = 0;
        }
        else _s = new alphabet_type[_sLen];


        UtilityFunctions::fileRead(_s, _fname, _beg, (_is_sentinel ? _sLen - 1 : _sLen));


        /// compute sa, use the relative position
        relative_offset_type *sa = new relative_offset_type[_sLen];

        char * t_buf = new char[_sLen / 8 + 1];
        BitWrapper t(t_buf);


        t.set(_sLen - 1, S_TYPE);
        t.set(_sLen - 2, L_TYPE);
        i = _sLen - 3;

        for (; ; --i)
        {
            t.set(i, (_s[i] < _s[i + 1] || (_s[i] == _s[i + 1] && t.get(i + 1) == S_TYPE)) ? S_TYPE : L_TYPE);
            if (i == 0) break;
        }

        // sort all the S-substrings
        relative_offset_type *bkt = new relative_offset_type[_alpha + 1];

        /// set the maximum value
        relative_offset_type MAX_VALUE(std::numeric_limits<relative_offset_type>::max());

        for (i = 0; i < _sLen; ++i) sa[i] = MAX_VALUE;//init sa

/*******************************************************************************************************/
#ifdef CILK
        /// step 1: initialize the S* tuple
        putSubstr0<alphabet_type,relative_offset_type>(sa, _s, bkt, _sLen, _alpha + 1, _is_sentinel);
        std::cout << __LINE__ << " putSubstr on Level " << _level << " done. \n";
        /// step 2: L-type inducing
        induceSAl0<alphabet_type,relative_offset_type>(sa,_s,bkt,_sLen,_alpha+1,false,_is_sentinel);
        std::cout << __LINE__ << ", induceSAL on Level " << _level << " done. \n";
        /// step 3: S-type inducing
        induceSAs0<alphabet_type,relative_offset_type>(sa,_s,bkt,_sLen,_alpha+1,false);
        std::cout << __LINE__ << ", induceSAS on Level " << _level << " done. \n";

/*******************************************************************************************************/
#else
        /// step 1: initialize the S* tuple
        getBuckets<alphabet_type,relative_offset_type>(_s, _sLen, bkt, _alpha + 1, true); // find end of buckets

        /// mind the virtual sentinel
        if( _is_sentinel ) i = _sLen - 2; // not put the sentinel into SA
        else i = _sLen - 1;

        for (; i >= 1; --i)
        {
            if (t.get(i) && !t.get(i - 1))
            {
                sa[bkt[_s[i]]] = i, --bkt[_s[i]];
            }
        }

        /// step 2: L-type inducing

        getBuckets<alphabet_type,relative_offset_type>(_s, _sLen, bkt, _alpha + 1, false); // find heads of buckets

        if(_is_sentinel) sa[0] = _sLen - 1;

        for (i = 0; i < _sLen; ++i)
        {
            if (sa[i] != MAX_VALUE)// && sa[i] != 0) // sa[i] is never equal to 0, because the first item in _s is not L-type.
            {
                j = sa[i] - 1;

                if(!t.get(j))
                {
                    sa[bkt[_s[j]]] = j, ++bkt[_s[j]];
                    sa[i] = MAX_VALUE;
                }
            }
        }

        /// step 3: S-type inducing
        getBuckets<alphabet_type, relative_offset_type>(_s, _sLen, bkt, _alpha + 1, true); /// find ends of buckets
        for (i = _sLen - 1; ; --i)
        {
            if (sa[i] && sa[i] != MAX_VALUE)
            {
                j = sa[i] - 1;

                if (t.get(j))
                {
                    sa[bkt[_s[j]]] = j, --bkt[_s[j]];
                    sa[i] = MAX_VALUE;
                }
            }
            if (i == 0) break;
        }
#endif // CILK
/*******************************************************************************************************/


        /// step 4: compare the adjacent S* substrings in SA
        j = _sLen / 8 + 1;
        char *  t_buf1 = new char[j];
        for( i = 0; i < j; ++i) t_buf1[i] = 0;
        BitWrapper is_cLMS(t_buf1);

        relative_offset_type d1(1),d2(1), lg(0); /// the offset of LMS1 and LMS2


        //MyVector< relative_offset_type > * lms_rPos = new MyVector<relative_offset_type>();
        //BitVectorWapper<uint64> * lms_bit_diff = new BitVectorWapper<uint64>();

//        lms_rPos->start_write();
//        lms_bit_diff->start_write();

        /// scan SA from right to left, for each item being not equal to MAX_VALUE, compare it to its right neighboring item to determine whether they are equal
        for(uint64 i = _sLen - 1; i >= 0; --i)
        {

            if(sa[i] != MAX_VALUE)  /// the first smallest LMS substring
            {
                j = i;

                /*lms_rPos->push_back(sa[i]);
                lms_bit_diff->push_back(true);*/ /// must be different

                _lms_rPos_vector->push_back(sa[i]);
                _lms_rPos_diff_bit_vector->push_back(true);

                while(true)
                {
                    if( !(t.get(sa[i] + d1 - 1)) && t.get(sa[i] + d1) )break;
                    else ++d1;
                }

                is_cLMS.set(sa[i] + d1,1); /// set(id, value)

                if( i == 0 ) break;

                /// scan from right to left
                for(i = i - 1; i >= 0; --i)
                {

                    if(sa[i] != MAX_VALUE)
                    {
//                        lms_rPos->push_back(sa[i]);
                        _lms_rPos_vector->push_back(sa[i]);

                        while(true)
                        {
                            if( !t.get(sa[i] + d2 - 1) && t.get( sa[i] + d2 ) )break;
                            else ++d2;
                        }

                        if (d1 != d2)
                        {
//                            lms_bit_diff->push_back(true);
                            _lms_rPos_diff_bit_vector->push_back(true);

                            is_cLMS.set(sa[i]+d2,1);
                        }
                        else
                        {
                            lg = 0;

                            while(lg <= d1)
                            {
                                if( _s[sa[i]+lg] != _s[sa[j]+lg] )
                                {
                                    break;
                                }
                                ++lg;
                            }

                            if(lg <= d1) /// LMS1 != LMS2
                            {
//                                lms_bit_diff->push_back(true);
                                _lms_rPos_diff_bit_vector->push_back(true);

                                is_cLMS.set(sa[i]+d2,1);
                            }
                            else /// LMS1 == LMS2
                            {
//                                lms_bit_diff->push_back(false);
                                _lms_rPos_diff_bit_vector->push_back(false);
                            }

                        }

                        j = i;
                        d1 = d2; /// the relative offset between the first and last S* character of a S* substring
                        d2 = 1;

                    }

                    if(i == 0) break;
                }
            }

            if(i == 0) break;
        }


        getBuckets<alphabet_type,relative_offset_type>(_s, _sLen, bkt, _alpha + 1, true); // find end of buckets

        for(uint64 i = 0; i < _sLen; ++i)sa[i] = MAX_VALUE;

        //for(uint64 i = 0; i <= _alpha; ++i) _LMS_counter[i]->push_back(0);

        /// step 5.1 initialize
        for(uint64 i = ( ( _is_sentinel && is_cLMS.get(_sLen - 1)) ? _sLen - 2 : _sLen - 1); i >= 0; --i) /// \note error-prone
        {
            if(is_cLMS.get(i))
            {
                sa[bkt[_s[i]]] = i, --bkt[_s[i]];
            }
            if(i == 0) break;
        }


        for(uint64 i = 0; i < _sLen; ++i)
        {
            if(sa[i] != MAX_VALUE) (*(_LMS_counter[_s[sa[i]]]))[_id]++;//_LMS_counter[_s[sa[i]]]->back()++;
        }

        if( _is_sentinel && is_cLMS.get(_sLen - 1)) /// \note error-prone
        {
            sa[0] = _sLen - 1;
        }


        /// step 5.2 L-type inducing
        typedef Pair<alphabet_type, compress_type> pair_type;

        alphabet_type cur_bwt(0);

        compress_type cur_cnt(0), max_count(std::numeric_limits<compress_type>::max());

        getBuckets<alphabet_type,relative_offset_type>(_s, _sLen, bkt, _alpha + 1, false); // find heads of buckets

        for (uint64 i = 0; i < _sLen; ++i)
        {
            if (sa[i] != MAX_VALUE)// && sa[i] != 0) // sa[i] is never equal to 0, because the first item in _s is not L-type.
            {
                j = sa[i] - 1;

                if(!t.get(j))
                {
                    /// j - 1 is greater than or equal to zero
                    /*if (!t.get(j-1))L_bwt_type->push_back(false);   ///  t[j] = L-type
                    else L_bwt_type->push_back(true); ///  t[j] = S-type*/
                    if (!t.get(j-1))_sub_l_bit_vector->push_back(false);   ///  t[j] = L-type
                    else _sub_l_bit_vector->push_back(true); ///  t[j] = S-type

                    cur_bwt = _s[j],cur_cnt = 1; /// first item's BWT and counter

                    sa[bkt[_s[j]]] = j, ++bkt[_s[j]];

                    sa[i] = MAX_VALUE;

                    for(i = i+1 ; i < _sLen; ++i)
                    {

                        if (sa[i] != MAX_VALUE)// && sa[i] != 0) // sa[i] is never equal to 0, because the first item in _s is not L-type.
                        {
                            j = sa[i] - 1;

                            if(!t.get(j))
                            {

                                /*if (!t.get(j-1))L_bwt_type->push_back(false);   ///  t[j] = L-type
                                else L_bwt_type->push_back(true); ///  t[j] = S-type*/
                                if (!t.get(j-1))_sub_l_bit_vector->push_back(false);   ///  t[j] = L-type
                                else _sub_l_bit_vector->push_back(true);

                                if(cur_bwt != _s[j])
                                {
                                    /*L_cBWT->push_back(pair_type(cur_bwt,cur_cnt));*/
                                    _sub_l_bwt_vector->push_back(pair_type(cur_bwt,cur_cnt));

                                    cur_bwt = _s[j], cur_cnt = 1;
                                }
                                else
                                {
                                    if(cur_cnt == max_count)
                                    {
                                        /*L_cBWT->push_back(pair_type(cur_bwt,cur_cnt));*/
                                        _sub_l_bwt_vector->push_back(pair_type(cur_bwt,cur_cnt));

                                        cur_cnt = 1;
                                    }
                                    else
                                    {
                                        cur_cnt++;
                                    }

                                }

                                sa[bkt[_s[j]]] = j, ++bkt[_s[j]];

                                sa[i] = MAX_VALUE;


                            }

                        }
                    }

                }

            }
        }

        _sub_l_bwt_vector->push_back(pair_type(cur_bwt,cur_cnt));

        cur_bwt = 0, cur_cnt = 0;

        getBuckets<alphabet_type,relative_offset_type>(_s, _sLen, bkt, _alpha + 1, true); // find ends of buckets

        for (i = _sLen - 1; ; --i) /// i > 0, because sa[0] is the current smallest item, whose preceding must be L-type
        {
            if (sa[i] != MAX_VALUE && sa[i])
            {
                j = sa[i] - 1;

                if (t.get(j))   /// 'j' is greater than or equal to zero
                {

//                    if(!j || !t.get(j-1)) s_bwt_type->push_back(true); /// if j==0, t[j]'s preceding is L-type, i.e. t[j] is S*
//                    else s_bwt_type->push_back(false); /// t[j] is not S*
                    if(!j || !t.get(j-1)) _sub_s_bit_vector->push_back(true); /// if j==0, t[j]'s preceding is L-type, i.e. t[j] is S*
                    else _sub_s_bit_vector->push_back(false);

                    cur_bwt = _s[j],cur_cnt = 1; /// first largest sBWT

                    sa[bkt[_s[j]]] = j, --bkt[_s[j]];

                    if(i == 0) break;

                    for(i = i-1; ; --i)
                    {

                        if (sa[i] != MAX_VALUE && sa[i])
                        {
                            j = sa[i] - 1;

                            if (t.get(j))   /// t[j] == S-type
                            {

                                /*if(!j || !t.get(j-1)) s_bwt_type->push_back(true); /// t[j] is S*
                                else s_bwt_type->push_back(false); /// t[j] is not S**/

                                if(!j || !t.get(j-1)) _sub_s_bit_vector->push_back(true); /// t[j] is S*
                                else _sub_s_bit_vector->push_back(false); /// t[j] is not S*

                                if(cur_bwt != _s[j])
                                {
//                                    s_cBwt->push_back(pair_type(cur_bwt,cur_cnt));
                                    _sub_s_bwt_vector->push_back(pair_type(cur_bwt,cur_cnt));

                                    cur_bwt = _s[j], cur_cnt = 1;
                                }
                                else
                                {
                                    if(cur_cnt == max_count)
                                    {
                                        /*s_cBwt->push_back(pair_type(cur_bwt,cur_cnt));*/
                                        _sub_s_bwt_vector->push_back(pair_type(cur_bwt,cur_cnt));

                                        cur_cnt = 1;
                                    }
                                    else
                                    {
                                        cur_cnt++;
                                    }

                                }

                                sa[bkt[_s[j]]] = j, --bkt[_s[j]];
                            }
                        }

                        if(i == 0) break;
                    }

                }
            }

            if (i == 0) break;
        }

        _sub_s_bwt_vector->push_back(pair_type(cur_bwt,cur_cnt));



        t.~BitWrapper();
        is_cLMS.~BitWrapper();
        delete[] sa;
        sa = nullptr;
        delete[] t_buf;
        t_buf = nullptr;
        delete[] t_buf1;
        t_buf1 = nullptr;
        delete[] bkt;
        bkt = nullptr;

        delete[] _s;
        _s = nullptr;


        return;
    }


    /// \note reduction phase
    /// alphabet > 255
    template<typename alphabet_type, typename compress_type, typename relative_offset_type>
    static void th_reduceSort_bwt_4B(std::string  _fname
                                     , const uint64  _beg
                                     , const uint64  _end
                                     , const alphabet_type  _alpha
                                     , const uint32  _level
                                     , const uint32 _id
                                     , MyVector< Pair<alphabet_type, compress_type> > *& _cLMS_vector
                                     , MyVector< Pair<alphabet_type, compress_type> > *& _sub_l_bwt_vector
                                     , BitVectorWapper<uint64>                        *& _sub_l_bit_vector
                                     , MyVector< Pair<alphabet_type, compress_type> > *& _sub_s_bwt_vector
                                     , BitVectorWapper<uint64>                        *& _sub_s_bit_vector
                                     , MyVector<relative_offset_type>                 *& _lms_rPos_vector
                                     , BitVectorWapper<uint64>                        *& _lms_rPos_diff_bit_vector
                                     , bool _is_single
                                     , bool _is_sentinel
                                    )
    {
        if(_is_single)
        {


            std::vector< std::vector<uint64> *> tp_vector;

            compute_singleton_bwt<alphabet_type, compress_type, relative_offset_type>(
                _fname
                , _beg
                , _end
                , _alpha
                , _level
                , _id
                , tp_vector
                , _cLMS_vector
                , _sub_l_bwt_vector
                , _sub_l_bit_vector
                , _sub_s_bwt_vector
                , _sub_s_bit_vector
                , _lms_rPos_vector
                , _lms_rPos_diff_bit_vector
                , _is_sentinel
            );


            return ;
        }


        uint64 i, j;

        const uint64 _sLen = _end - _beg + 1;

        alphabet_type * _s = new alphabet_type[_sLen];

        UtilityFunctions::fileRead(_s,_fname,_beg,_sLen);

        /// Prepare: compress alphabet_type to small relative_offset_type

        typedef Pair<alphabet_type,relative_offset_type> alphabet_pair_type;

        typedef TupleAscCmp1< alphabet_pair_type > pair_Asc_Compatator_type;

        alphabet_pair_type * t_pair_ary = new alphabet_pair_type[_sLen]; /// convert pair

        for(i = 0; i < _sLen; i++)
            t_pair_ary[i] = alphabet_pair_type(_s[i], i);


        std::sort(t_pair_ary, t_pair_ary + _sLen, pair_Asc_Compatator_type());

        alphabet_type * s = new alphabet_type [_sLen];

        relative_offset_type cur_name(1);

        alphabet_type t_value = t_pair_ary[0].first; /// first pair

        t_pair_ary[0].first = cur_name; /// replace the original name with new name

        for(i = 1; i < _sLen; i++) /// start from the second item
        {

            if(t_pair_ary[i].first != t_value) ++cur_name;

            t_value = t_pair_ary[i].first;

            t_pair_ary[i].first = cur_name;

        }

        for(i = 0; i < _sLen; i++)
            s[ t_pair_ary[i].second ] = t_pair_ary[i].first;


        delete[] t_pair_ary;


        /// compute sa, use the relative position
        relative_offset_type *sa = new relative_offset_type[_sLen];

        char * t_buf = new char[_sLen / 8 + 1];
        BitWrapper t(t_buf);


        t.set(_sLen - 1, S_TYPE);
        t.set(_sLen - 2, L_TYPE);
        i = _sLen - 3;

        for (; ; --i)
        {
            t.set(i, (_s[i] < _s[i + 1] || (_s[i] == _s[i + 1] && t.get(i + 1) == S_TYPE)) ? S_TYPE : L_TYPE);
            if (i == 0) break;
        }

        // sort all the S-substrings
        relative_offset_type *bkt = new relative_offset_type[cur_name + 1];

        /// set the maximum value
        relative_offset_type MAX_VALUE(std::numeric_limits<relative_offset_type>::max());

        for (i = 0; i < _sLen; ++i) sa[i] = MAX_VALUE;//init sa

/*******************************************************************************************************/
#ifdef CILK
        /// step 1: initialize the S* tuple
        if(sizeof(alphabet_type) == 1)
        {
            putSubstr0<alphabet_type,relative_offset_type>(sa, s, bkt, _sLen, cur_name + 1, _is_sentinel);
        }
        else
        {
            getBuckets<alphabet_type,relative_offset_type>(s, _sLen, bkt, cur_name + 1, true); // find end of buckets

            /// mind the virtual sentinel
            if(_is_sentinel ) i = _sLen - 2; // not put the sentinel into SA
            else i = _sLen - 1;

            for (; i >= 1; --i)
            {
                if (t.get(i) && !t.get(i - 1))
                {
                    sa[bkt[s[i]]] = i, --bkt[s[i]];
                }
            }
        }
        std::cout << __LINE__ << ", pubSubstr done at Level " << _level << std::endl;

        /// step 2: L-type inducing
        induceSAl0<alphabet_type,relative_offset_type>(sa, s, bkt, _sLen, cur_name + 1, false, _is_sentinel);
        std::cout << __LINE__ << ", induceSA_L_sub done at Level " << _level << std::endl;
        induceSAs0<alphabet_type,relative_offset_type>(sa, s, bkt, _sLen, cur_name + 1, false);
        std::cout << __LINE__ << ", induceSA_S_sub done at Level " << _level << std::endl;
/*******************************************************************************************************/
#else
        /// step 1: initialize the S* tuple
        getBuckets<alphabet_type,relative_offset_type>(s, _sLen, bkt, cur_name + 1, true); // find end of buckets

        /// mind the virtual sentinel
        if(_is_sentinel ) i = _sLen - 2; // not put the sentinel into SA
        else i = _sLen - 1;

        for (; i >= 1; --i)
        {
            if (t.get(i) && !t.get(i - 1))
            {
                sa[bkt[s[i]]] = i, --bkt[s[i]];
            }
        }
        /// step 2: L-type inducing
        getBuckets<alphabet_type,relative_offset_type>(s, _sLen, bkt, cur_name + 1, false); // find heads of buckets

        for (i = 0; i < _sLen; ++i)
        {
            if (sa[i] != MAX_VALUE)// && sa[i] != 0) // sa[i] is never equal to 0, because the first item in _s is not L-type.
            {
                j = sa[i] - 1;

                if(!t.get(j))
                {
                    sa[bkt[s[j]]] = j, ++bkt[s[j]]; /// the substitute string
                    sa[i] = MAX_VALUE;
                }
            }
        }
        /// step 3: S-type inducing
        getBuckets<alphabet_type, relative_offset_type>(s, _sLen, bkt, cur_name + 1, true); /// find ends of buckets

        for (i = _sLen - 1; ; --i)
        {
            if (sa[i] && sa[i] != MAX_VALUE)
            {
                j = sa[i] - 1;

                if (t.get(j))
                {
                    sa[bkt[s[j]]] = j, --bkt[s[j]]; /// the substitute string
                    sa[i] = MAX_VALUE;
                }
            }
            if (i == 0) break;
        }
#endif // CILK
/*******************************************************************************************************/

        /// step 4: compare the adjacent S* substrings in SA
        j = _sLen / 8 + 1;
        char *  t_buf1 = new char[j];
        for( i = 0; i < j; i++) t_buf1[i] = 0;
        BitWrapper is_cLMS(t_buf1);

        relative_offset_type d1(1),d2(1), lg(0); /// the offset of LMS1 and LMS2


        /// scan SA from right to left, for each item being not equal to MAX_VALUE, compare it to its right neighboring item to determine whether they are equal
        for(i = _sLen - 1; i >= 0; i--)
        {

            if(sa[i] != MAX_VALUE)  /// the first smallest LMS substring
            {
                j = i;

                _lms_rPos_vector->push_back(sa[i]);
                _lms_rPos_diff_bit_vector->push_back(true);

                while(true)
                {
                    if( !(t.get(sa[i] + d1 - 1)) && t.get(sa[i] + d1) )break;
                    else ++d1;
                }

                is_cLMS.set(sa[i] + d1,1); /// set(id, value)

                if(i==0)break;

                ///scan from right to left
                for(i = i - 1; i >= 0; i--)
                {

                    if(sa[i] != MAX_VALUE)
                    {
                        /*lms_rPos->push_back(sa[i]);*/
                        _lms_rPos_vector->push_back(sa[i]);

                        while(true)
                        {
                            if( !t.get(sa[i] + d2 - 1) && t.get( sa[i] + d2 ) )break;
                            else ++d2;
                        }

                        if (d1 != d2)
                        {
                            /*lms_bit_diff->push_back(true);*/
                            _lms_rPos_diff_bit_vector->push_back(true);

                            is_cLMS.set(sa[i]+d2,1);
                        }
                        else
                        {
                            lg = 0;

                            while(lg <= d1)
                            {
                                if( _s[sa[i]+lg] != _s[sa[j]+lg] )
                                {
                                    break;
                                }
                                ++lg;
                            }

                            if(lg <= d1) /// LMS1 != LMS2
                            {
                                /*lms_bit_diff->push_back(true);*/
                                _lms_rPos_diff_bit_vector->push_back(true);

                                is_cLMS.set(sa[i]+d2,1);
                            }
                            else /// LMS1 == LMS2
                            {
                                /*lms_bit_diff->push_back(false);*/
                                _lms_rPos_diff_bit_vector->push_back(false);
                            }

                        }

                        j = i;
                        d1 = d2; /// the relative offset between the first and last S* character of a S* substring
                        d2 = 1;

                    }

                    if(i==0)break;
                }
            }

            if(i == 0) break;
        }

        /*_lms_seq_vector.push_back(lms_rPos);
        _lms_diff_bit_vector.push_back(lms_bit_diff);*/

        /// step 5: compute bwt
        for(i = 0; i < _sLen; ++i)sa[i] = MAX_VALUE;

        /// step 5.1 initialize
        getBuckets<alphabet_type,relative_offset_type>(s, _sLen, bkt, cur_name + 1, true); // find end of buckets

        for( i = ( ( _is_sentinel && is_cLMS.get(_sLen - 1)) ? _sLen - 2 : _sLen - 1); i >= 0; --i)
        {
            if(is_cLMS.get(i))
            {
                sa[bkt[s[i]]] = i, --bkt[s[i]]; /// the substitute string
            }
            if(i == 0) break;
        }

        /// generate the cLMS
        typedef Pair<alphabet_type,compress_type> lms_pair_type, bwt_pair_type;

        alphabet_type cur_lms(0), cur_bwt(0);

        compress_type cur_cnt(0), max_count(std::numeric_limits<compress_type>::max());


        for(i = 0; i < _sLen; ++i)
        {
            if(sa[i] != MAX_VALUE) /// first LMS char
            {
                cur_lms = _s[sa[i]], cur_cnt = 1; /// original string

                for(i = i + 1; i < _sLen; ++i)  /// from the second item
                {

                    if(sa[i] != MAX_VALUE)
                    {

                        if(cur_lms != _s[sa[i]])  /// original string
                        {

                            /*lms_cChar->push_back(lms_pair_type(cur_lms,cur_cnt));*/
                            _cLMS_vector->push_back(lms_pair_type(cur_lms,cur_cnt));

                            cur_lms = _s[sa[i]], cur_cnt = 1;

                        }
                        else
                        {

                            if(cur_cnt == max_count)
                            {

                                /*lms_cChar->push_back(lms_pair_type(cur_lms,cur_cnt));*/
                                _cLMS_vector->push_back(lms_pair_type(cur_lms,cur_cnt));

                                cur_cnt = 1;

                            }
                            else
                            {

                                cur_cnt++;
                            }

                        }
                    }

                }
            }
        }

        /*lms_cChar->push_back(lms_pair_type(cur_lms,cur_cnt));*/ /// last LMS pair
        _cLMS_vector->push_back(lms_pair_type(cur_lms,cur_cnt)); /// last LMS pair

        /*_cLms_vector.push_back(lms_cChar);*/


        /// step 5.2 L-type inducing

        /// put the sentinel position
        if( _is_sentinel && is_cLMS.get(_sLen - 1))
        {
            sa[0] = _sLen - 1;
        }


        getBuckets<alphabet_type,relative_offset_type>(s, _sLen, bkt, cur_name + 1, false); // find heads of buckets

        for (i = 0; i < _sLen; ++i)
        {
            if (sa[i] != MAX_VALUE)// && sa[i] != 0) // sa[i] is never equal to 0, because the first item in _s is not L-type.
            {
                j = sa[i] - 1;

                if(!t.get(j))
                {
                    /// j - 1 is greater than or equal to zero
                    /*if (!t.get(j-1))l_bwt_type->push_back(false);   ///  t[j] = L-type
                    else l_bwt_type->push_back(true);*/ ///  t[j] = S-type
                    if (!t.get(j-1))_sub_l_bit_vector->push_back(false);   ///  t[j] = L-type
                    else _sub_l_bit_vector->push_back(true);

                    cur_bwt = _s[j],cur_cnt = 1; /// first item's BWT and counter; cur_bwt is the original char

                    sa[bkt[s[j]]] = j, ++bkt[s[j]]; /// the substitute string

                    sa[i] = MAX_VALUE;

                    for(i = i+1 ; i < _sLen; ++i)
                    {

                        if (sa[i] != MAX_VALUE)// && sa[i] != 0) // sa[i] is never equal to 0, because the first item in _s is not L-type.
                        {
                            j = sa[i] - 1;

                            if(!t.get(j))
                            {

                                /*if (!t.get(j-1))l_bwt_type->push_back(false);   ///  t[j] = L-type
                                else l_bwt_type->push_back(true);*/ ///  t[j] = S-type
                                if (!t.get(j-1))_sub_l_bit_vector->push_back(false);   ///  t[j] = L-type
                                else _sub_l_bit_vector->push_back(true);

                                if(cur_bwt != _s[j])
                                {
                                    /*l_cBwt->push_back(bwt_pair_type(cur_bwt,cur_cnt));*/
                                    _sub_l_bwt_vector->push_back(bwt_pair_type(cur_bwt,cur_cnt));

                                    cur_bwt = _s[j], cur_cnt = 1;
                                }
                                else
                                {
                                    if(cur_cnt == max_count)
                                    {
                                        /*l_cBwt->push_back(bwt_pair_type(cur_bwt,cur_cnt));*/
                                        _sub_l_bwt_vector->push_back(bwt_pair_type(cur_bwt,cur_cnt));

                                        cur_cnt = 1;
                                    }
                                    else
                                    {
                                        cur_cnt++;
                                    }

                                }

                                sa[bkt[s[j]]] = j, ++bkt[s[j]]; /// the substitute string

                                sa[i] = MAX_VALUE;


                            }
                        }
                    }
                }
            }
        }

        /*l_cBwt->push_back(bwt_pair_type(cur_bwt,cur_cnt));*//// last cBWT pair
        _sub_l_bwt_vector->push_back(bwt_pair_type(cur_bwt,cur_cnt));/// last cBWT pair

        cur_bwt = 0, cur_cnt = 0;

        getBuckets<alphabet_type,relative_offset_type>(s, _sLen, bkt, cur_name + 1, true); // find ends of buckets

        for (i = _sLen - 1; ; --i) /// i > 0, because sa[0] is the current smallest item, whose preceding must be L-type
        {
            if (sa[i] && sa[i] != MAX_VALUE)
            {
                j = sa[i] - 1;

                if (t.get(j))   /// 'j' is greater than or equal to zero
                {

                    /*if(!j || !t.get(j-1)) s_bwt_type->push_back(true); /// if j==0, t[j]'s preceding is L-type, i.e. t[j] is S*
                    else s_bwt_type->push_back(false); /// t[j] is not S**/
                    if(!j || !t.get(j-1)) _sub_s_bit_vector->push_back(true); /// if j==0, t[j]'s preceding is L-type, i.e. t[j] is S*
                    else _sub_s_bit_vector->push_back(false); /// t[j] is not S*

                    cur_bwt = _s[j],cur_cnt = 1; /// first largest sBWT

                    sa[bkt[s[j]]] = j, --bkt[s[j]]; /// the substitute string

                    if(i == 0) break;

                    for(i = i-1; ; i--)
                    {

                        if (sa[i] && sa[i] != MAX_VALUE)
                        {
                            j = sa[i] - 1;

                            if (t.get(j))   /// t[j] == S-type
                            {

                                /*if(!j || !t.get(j-1)) s_bwt_type->push_back(true); /// t[j] is S*
                                else s_bwt_type->push_back(false); /// t[j] is not S**/
                                if(!j || !t.get(j-1)) _sub_s_bit_vector->push_back(true); /// t[j] is S*
                                else _sub_s_bit_vector->push_back(false); /// t[j] is not S*

                                if(cur_bwt != _s[j])
                                {
                                    /*s_cBwt->push_back(bwt_pair_type(cur_bwt,cur_cnt));*/
                                    _sub_s_bwt_vector->push_back(bwt_pair_type(cur_bwt,cur_cnt));

                                    cur_bwt = _s[j], cur_cnt = 1;
                                }
                                else
                                {
                                    if(cur_cnt == max_count)
                                    {
                                        /*s_cBwt->push_back(bwt_pair_type(cur_bwt,cur_cnt));*/
                                        _sub_s_bwt_vector->push_back(bwt_pair_type(cur_bwt,cur_cnt));

                                        cur_cnt = 1;
                                    }
                                    else
                                    {
                                        cur_cnt++;
                                    }

                                }

                                sa[bkt[s[j]]] = j, --bkt[s[j]]; /// the substitute string
                            }
                        }

                        if(i == 0) break;
                    }

                }
            }

            if (i == 0) break;
        }

        /*s_cBwt->push_back(bwt_pair_type(cur_bwt,cur_cnt));*//// last sBWT
        _sub_s_bwt_vector->push_back(bwt_pair_type(cur_bwt,cur_cnt));/// last sBWT


        t.~BitWrapper();
        is_cLMS.~BitWrapper();
        delete[] sa;
        sa = nullptr;
        delete[] t_buf;
        t_buf = nullptr;
        delete[] t_buf1;
        t_buf1 = nullptr;
        delete[] bkt;
        bkt = nullptr;

        delete[] s; /// \note error-prone
        s = nullptr;

        delete[] _s; /// \note error-prone
        _s = nullptr;

        return;
    }


    /// \note alphabet <= 255 for is_unique = false
    template<typename alphabet_type, typename compress_type, typename relative_offset_type>
    static void th_induceSort_bwt_1B(    std::string                                      _fname
                                         , const uint64                                     _beg
                                         , const uint64                                     _end
                                         , const alphabet_type                              _alpha
                                         , const uint32                                     _level
                                         , const uint32                                     _id
                                         , std::vector< std::vector<uint64> *>            & _LMS_counter
                                         //, MyVector< Pair<alphabet_type, compress_type> > * _cLMS_vector
                                         , MyVector< Pair<alphabet_type, compress_type> > *& _suf_l_bwt_vector
                                         , BitVectorWapper<uint64>                        *& _suf_l_bit_vector
                                         , MyVector< Pair<alphabet_type, compress_type> > *& _suf_s_bwt_vector
                                         , BitVectorWapper<uint64>                        *& _suf_s_bit_vector
                                         , MyVector<relative_offset_type>                 *& _lms_rPos_vector
                                         , bool _is_single
                                         , bool _is_sentinel
                                         , MyVector< relative_offset_type >               *& _lms_order /// save the sorted rank of LMS from big to small
                                    )
    {
        if( _is_single )
        {
            std::cout << "The current block is single block." << std::endl;

            _lms_order->start_read();

            delete _lms_order;

            _lms_order = nullptr;

            MyVector< Pair<alphabet_type, compress_type> > * cLMS_vector =  nullptr;

            compute_singleton_bwt_induce(
                _fname
                , _beg
                , _end
                , _alpha
                , _level
                , _id
                , _LMS_counter
                , cLMS_vector
                , _suf_l_bwt_vector
                , _suf_l_bit_vector
                , _suf_s_bwt_vector
                , _suf_s_bit_vector
                , _lms_rPos_vector
                , _is_sentinel
            );

            return;
        }

        uint64 i, j, sLen = _end - _beg + 1;

        bool is_virtual_sentinel = _is_sentinel;

        alphabet_type * _s = nullptr;;

        if (_is_sentinel)
        {
            sLen++;
            _s = new alphabet_type[sLen];///\note memory leak: Invalid write of size 1
            _s[sLen - 1] = 0;
        }
        else _s = new alphabet_type[sLen];

        //if(!_level && _lms_rPos_vector.empty()) is_virtual_sentinel = true;

        UtilityFunctions::fileRead(_s,_fname,_beg,(_is_sentinel ? sLen - 1 : sLen));

        /// compute sa, use the relative position
        relative_offset_type *sa = new relative_offset_type[sLen];

        char * t_buf = new char[sLen / 8 + 1];
        BitWrapper t(t_buf);


        t.set(sLen - 1, S_TYPE);
        t.set(sLen - 2, L_TYPE);
        i = sLen - 3;

        for (; ; --i)
        {
            t.set(i, (_s[i] < _s[i + 1] || (_s[i] == _s[i + 1] && t.get(i + 1) == S_TYPE)) ? S_TYPE : L_TYPE);
            if (i == 0) break;
        }

        /// scan the s from left to right, to statistic the LMS position
        relative_offset_type * _LMS_pos = new relative_offset_type [_lms_order->size()];

        _LMS_pos[_lms_order->size() - 1] = sLen - 1;

        j = _lms_order->size() - 2;

        /// find lms-chars in s[1, sLen - 3]
        for (i = sLen - 3; i >= 1; --i)
        {
            if (t.get(i) && !t.get(i - 1))
            {
                _LMS_pos[j] = i;

                if(j == 0)break;
                else --j;
            }
        }


        MyVector<relative_offset_type> * sorted_lms_pos = new MyVector<relative_offset_type>();

        sorted_lms_pos->start_write();

        _lms_order->start_read_reverse(); /// from large to small, include the sentinel

        while( !_lms_order->is_eof() )
        {
            sorted_lms_pos->push_back( _LMS_pos[ _lms_order->get_reverse() - 1] );

            _lms_order->next_remove_reverse();
        }

        if(_lms_order->is_eof())
        {
            delete _lms_order;
            _lms_order = nullptr;
        }

        delete [] _LMS_pos;
        _LMS_pos = nullptr;


        // sort all the S-substrings
        relative_offset_type *bkt = new relative_offset_type[_alpha + 1];

//        std::cout << "alpha + 1 = " << _alpha + 1 - 0 << std::endl;
        /// set the maximum value
        relative_offset_type MAX_VALUE(std::numeric_limits<relative_offset_type>::max());

        try
        {
            for (i = 0; i < sLen; ++i) sa[i] = MAX_VALUE;//init sa
        }
        catch(std::exception & e)
        {
            std::cout << "Standard exception: " << e.what() << std::endl;
        }


        /// step 1: initialize the S* tuple
        getBuckets<alphabet_type,relative_offset_type>(_s, sLen, bkt, _alpha + 1, true); // find end of buckets

//        std::cout << "compute the end of buckets.\n";

        sorted_lms_pos->start_read();

        while(!sorted_lms_pos->is_eof())
        {

            j = sorted_lms_pos->get();

            sa[ bkt[_s[j]] ] = j, --bkt[_s[j]];

            sorted_lms_pos->next_remove();
        }

        /// delete the sorted_lms_pos
        if(sorted_lms_pos->is_eof())
        {
            delete sorted_lms_pos;

            sorted_lms_pos = nullptr;
        }

        if(is_virtual_sentinel)
        {

            // ++bkt[0], sa[ bkt[0] ] = MAX_VALUE
            ++bkt[_s[j]], sa[ bkt[_s[j]] ] = MAX_VALUE;  /// \note error-prone

            sa[0] = sLen - 1;
        }



        // debug

        bool debug_on = false;

        std::vector<relative_offset_type> init_lms_pos;
        std::vector<relative_offset_type> cpt_lms_pos;


        for(i = 0; i < sLen; ++i)
        {

            if(sa[i] != MAX_VALUE)
            {

                /*_LMS_counter[_s[sa[i]]]->back()++;*/
                (*_LMS_counter[_s[sa[i]]])[_id]++;

                if(debug_on)init_lms_pos.push_back(sa[i]);

                /*lms_rPos->push_back(sa[i]);*/
                _lms_rPos_vector->push_back(sa[i]);

            }

        }

        /// step 2 L-type inducing
        typedef Pair<alphabet_type, compress_type> pair_type;


        alphabet_type cur_bwt(0);

        compress_type cur_cnt(0), max_count(std::numeric_limits<compress_type>::max());

//        std::cout << "The max count of 'compress_type' = " << max_count - 0 << std::endl;

        getBuckets<alphabet_type,relative_offset_type>(_s, sLen, bkt, _alpha + 1, false); // find heads of buckets

        for (uint64 i = 0; i < sLen; ++i)
        {
            if (sa[i] != MAX_VALUE)// && sa[i] != 0) // sa[i] is never equal to 0, because the first item in _s is not L-type.
            {
                j = sa[i] - 1;

                // \note retrieve the preceding char of L* type suffix
                if( !_level && is_getBWT )
                {
                    // if _s[sa[i]] is L* type, push _s[sa[i]-1] into external BWT sequence.
                    if(!t.get(j+1) && t.get(j))
                        ((MyVector<uint8> *)LStarPreChar_seqs[_id])->push_back(_s[j]);
                }


                if(!t.get(j))
                {
                    /// j - 1 is greater than or equal to zero
                    /*if (!t.get(j-1))L_bwt_type->push_back(false);   ///  t[j] = L-type
                    else L_bwt_type->push_back(true);*/ ///  t[j] = S-type

                    if (!t.get(j-1))_suf_l_bit_vector->push_back(false);   ///  t[j] = L-type
                    else _suf_l_bit_vector->push_back(true);

                    cur_bwt = _s[j],cur_cnt = 1; /// first item's BWT and counter

                    sa[bkt[_s[j]]] = j, ++bkt[_s[j]];

                    sa[i] = MAX_VALUE;

                    for(i = i+1 ; i < sLen; ++i)
                    {

                        if (sa[i] != MAX_VALUE)// && sa[i] != 0) // sa[i] is never equal to 0, because the first item in _s is not L-type.
                        {
                            j = sa[i] - 1;


                            if( !_level && is_getBWT )
                            {
                                // if _s[sa[i]] is L* type, push _s[sa[i]-1] into external BWT sequence.
                                if(!t.get(j+1) && t.get(j))
                                    ((MyVector<uint8> *)LStarPreChar_seqs[_id])->push_back(_s[j]);
                            }


                            if(!t.get(j))
                            {

                                /*if (!t.get(j-1))L_bwt_type->push_back(false);   ///  t[j] = L-type
                                else L_bwt_type->push_back(true);*/ ///  t[j] = S-type

                                if (!t.get(j-1))_suf_l_bit_vector->push_back(false);   ///  t[j] = L-type
                                else _suf_l_bit_vector->push_back(true);

                                if(cur_bwt != _s[j])
                                {
                                    /*L_cBWT->push_back(pair_type(cur_bwt,cur_cnt));*/
                                    _suf_l_bwt_vector->push_back(pair_type(cur_bwt,cur_cnt));

                                    cur_bwt = _s[j], cur_cnt = 1;
                                }
                                else
                                {
                                    if(cur_cnt == max_count)
                                    {
                                        /*L_cBWT->push_back(pair_type(cur_bwt,cur_cnt));*/
                                        _suf_l_bwt_vector->push_back(pair_type(cur_bwt,cur_cnt));

                                        cur_cnt = 1;
                                    }
                                    else
                                    {
                                        cur_cnt++;
                                    }

                                }

                                sa[bkt[_s[j]]] = j, ++bkt[_s[j]];

                                sa[i] = MAX_VALUE;


                            }

                        }
                    }

                }

            }
        }

        /*L_cBWT->push_back(pair_type(cur_bwt,cur_cnt));*//// last cBWT pair
        _suf_l_bwt_vector->push_back(pair_type(cur_bwt,cur_cnt));


        // \note retrieve the preceding char of the L-type suffix
//        if(!_level && is_getBWT)
//        {
//            // the writing is completed.
//            L_eBWT_seqs[_id]->end_write();
//        }


        cur_bwt = 0, cur_cnt = 0;

        getBuckets<alphabet_type,relative_offset_type>(_s, sLen, bkt, _alpha + 1, true); // find ends of buckets

        for (i = sLen - 1; ; --i) /// i > 0, because sa[0] is the current smallest item, whose preceding must be L-type
        {

            if(sa[i] == 0)
            {
                if( !_level && is_getBWT )
                {
                    uint8 * tp_char = new uint8[1];

                    if(_beg)
                    {
                        UtilityFunctions::fileRead<uint8>(tp_char,_fname,_beg-1,1);
                    }
                    else
                    {
                        Logger::output_error(__FILE__,__LINE__);
                    }

                    ((MyVector<uint8> *)SStarPreChar_seqs[_id])->push_back(tp_char[0]);

                    delete tp_char;
                }
            }


            if (sa[i] && sa[i] != MAX_VALUE)
            {
                j = sa[i] - 1;


                // \note retrieve the preceding char of S* type suffix
                if( !_level && is_getBWT )
                {
                    // if _s[sa[i]] is S* type, push _s[sa[i]-1] into external BWT sequence.
                    if(t.get(j+1) && !t.get(j))
                        ((MyVector<uint8> *)SStarPreChar_seqs[_id])->push_back(_s[j]);

                    if(debug_on)
                    {
                        if(!t.get(j+1) && t.get(j))cpt_lms_pos.push_back(j+1);
                    }
                }

                if (t.get(j))   /// 'j' is greater than or equal to zero
                {

                    /*if(!j || !t.get(j-1)) s_bwt_type->push_back(true); /// if j==0, t[j]'s preceding is L-type, i.e. t[j] is S*
                    else s_bwt_type->push_back(false);*/ /// t[j] is not S*

                    if(!j || !t.get(j-1)) _suf_s_bit_vector->push_back(true); /// if j==0, t[j]'s preceding is L-type, i.e. t[j] is S*
                    else _suf_s_bit_vector->push_back(false);

                    cur_bwt = _s[j],cur_cnt = 1; /// first largest sBWT

                    sa[bkt[_s[j]]] = j, --bkt[_s[j]];

                    if(i == 0) break;

                    for(i = i-1; ; --i)
                    {

                        if(sa[i] == 0)
                        {
                            if( !_level && is_getBWT )
                            {
                                uint8 * tp_char = new uint8[1];

                                if(_beg)
                                {
                                    UtilityFunctions::fileRead<uint8>(tp_char,_fname,_beg-1,1);
                                }
                                else
                                {
                                    Logger::output_error(__FILE__,__LINE__);
                                }

                                ((MyVector<uint8> *)SStarPreChar_seqs[_id])->push_back(tp_char[0]);

                                delete tp_char;
                            }
                        }

                        if (sa[i] && sa[i] != MAX_VALUE)
                        {
                            j = sa[i] - 1;

                            if( !_level && is_getBWT )
                            {
                                // if _s[sa[i]] is S* type, push _s[sa[i]-1] into external BWT sequence.
                                if(t.get(j+1) && !t.get(j))
                                    ((MyVector<uint8> *)SStarPreChar_seqs[_id])->push_back(_s[j]);

                                if(debug_on)
                                {
                                    if(!t.get(j+1) && t.get(j))cpt_lms_pos.push_back(j+1);
                                }
                            }


                            if (t.get(j))   /// t[j] == S-type
                            {

                                /*if(!j || !t.get(j-1)) s_bwt_type->push_back(true); /// t[j] is S*
                                else s_bwt_type->push_back(false);*/ /// t[j] is not S*

                                if(!j || !t.get(j-1)) _suf_s_bit_vector->push_back(true); /// t[j] is S*
                                else _suf_s_bit_vector->push_back(false); /// t[j] is not S*

                                if(cur_bwt != _s[j])
                                {
                                    /*s_cBwt->push_back(pair_type(cur_bwt,cur_cnt));*/
                                    _suf_s_bwt_vector->push_back(pair_type(cur_bwt,cur_cnt));

                                    cur_bwt = _s[j], cur_cnt = 1;
                                }
                                else
                                {
                                    if(cur_cnt == max_count)
                                    {
                                        /*s_cBwt->push_back(pair_type(cur_bwt,cur_cnt));*/
                                        _suf_s_bwt_vector->push_back(pair_type(cur_bwt,cur_cnt));

                                        cur_cnt = 1;
                                    }
                                    else
                                    {
                                        cur_cnt++;
                                    }

                                }

                                sa[bkt[_s[j]]] = j, --bkt[_s[j]];
                            }
                        }

                        if(i == 0) break;
                    }

                }
            }

            if (i == 0) break;
        }

        /*s_cBwt->push_back(pair_type(cur_bwt,cur_cnt));*//// last sBWT
        _suf_s_bwt_vector->push_back(pair_type(cur_bwt,cur_cnt));/// last sBWT

        // \note retrieve the preceding char of the S-type suffix
//        if(is_getBWT)
//        {
//            // the writing is completed.
//            S_eBWT_seqs[_id]->end_write();
//        }

        if(debug_on)
        {
            int n = init_lms_pos.size();

            std::cout << "init_lms_pos.size() = " << init_lms_pos.size() << std::endl;
            std::cout << "cpt_lms_pos.size() = " << cpt_lms_pos.size() << std::endl;

            for(int i = 0; i < n; i++)
            {
                if(init_lms_pos[i] != cpt_lms_pos[i])
                {
                    std::cout << "init_lms_pos[" << i << "] = " << init_lms_pos[i] << std::endl;
                    std::cout << "cpt_lms_pos[" << i << "] = " << cpt_lms_pos[i] << std::endl;
                    std::cin.get();
                }
            }
        }
//        std::cout << "Line = " << __LINE__ << std::endl;
//        std::cin.get();
        if( !_level && is_getBWT )
        {
            std::cout << "block id = " << _id << ", lms number is " << ((MyVector<uint8> *)SStarPreChar_seqs[_id])->size() << std::endl;
            std::cout << "block id = " << _id << ", lml number is " << ((MyVector<uint8> *)LStarPreChar_seqs[_id])->size() << std::endl;
        }

        t.~BitWrapper();
        delete[] sa;
        sa = nullptr;
        delete[] t_buf;
        t_buf = nullptr;
        delete[] bkt;
        bkt = nullptr;

        delete[] _s;
        _s = nullptr;


        return;
    }

    /// \note alphabet > 255 for is_unique = true
    template<typename alphabet_type, typename compress_type, typename relative_offset_type>
    static void th_induceSort_bwt_4B(std::string                                          _fname
                                     , const uint64                                     _beg
                                     , const uint64                                     _end
                                     , const alphabet_type                              _alpha
                                     , const uint32                                     _level
                                     , const uint32                                     _id
//                                       , std::vector< std::vector<uint64> *>            & _LMS_counter
                                     , MyVector< Pair<alphabet_type, compress_type> > *& _cLMS_vector
                                     , MyVector< Pair<alphabet_type, compress_type> > *& _suf_l_bwt_vector
                                     , BitVectorWapper<uint64>                        *& _suf_l_bit_vector
                                     , MyVector< Pair<alphabet_type, compress_type> > *& _suf_s_bwt_vector
                                     , BitVectorWapper<uint64>                        *& _suf_s_bit_vector
                                     , MyVector<relative_offset_type>                 *& _lms_rPos_vector
                                     , bool _is_single
                                     , bool _is_sentinel
                                     , MyVector< relative_offset_type >               *& _lms_order /// save the sorted rank of LMS from big to small
                                    )
    {
        if( _is_single )
        {

            _lms_order->start_read();

            delete _lms_order;

            _lms_order = nullptr;

            std::vector< std::vector<uint64> * > LMS_counter;

            compute_singleton_bwt_induce(
                _fname
                , _beg
                , _end
                , _alpha
                , _level
                , _id
                , LMS_counter
                , _cLMS_vector
                , _suf_l_bwt_vector
                , _suf_l_bit_vector
                , _suf_s_bwt_vector
                , _suf_s_bit_vector
                , _lms_rPos_vector
                , _is_sentinel
            );

            return;
        }


        uint64 i, j, _sLen = _end - _beg + 1;

        bool is_virtual_sentinel = _is_sentinel;

        if(_is_sentinel)
        {
            std::cout << __LINE__ << "error.\n";
            std::cin.get();
        }

//        if(!_level && _lms_seq_vector.empty()) is_virtual_sentinel = true;
        alphabet_type * _s = new alphabet_type[_sLen];

        UtilityFunctions::fileRead(_s,_fname,_beg,_sLen);


        /// Prepare: compress alphabet_type to small relative_offset_type

        typedef Pair<alphabet_type,relative_offset_type> alphabet_pair_type;

        typedef TupleAscCmp1< alphabet_pair_type > pair_Asc_Compatator_type;

        alphabet_pair_type * t_pair_ary = new alphabet_pair_type[_sLen]; /// convert pair

        for(i = 0; i < _sLen; i++)
            t_pair_ary[i] = alphabet_pair_type(_s[i], i);

        /// use STL parallel sort Library, mind the RAM capacity

        std::sort(t_pair_ary, t_pair_ary + _sLen, pair_Asc_Compatator_type());

        alphabet_type * s = new alphabet_type [_sLen];

        relative_offset_type cur_name(1);

        alphabet_type t_value = t_pair_ary[0].first; /// first pair

        t_pair_ary[0].first = cur_name; /// replace the original name with new name

        for(i = 1; i < _sLen; i++) /// start from the second item
        {

            if(t_pair_ary[i].first != t_value) ++cur_name;

            t_value = t_pair_ary[i].first;

            t_pair_ary[i].first = cur_name;

        }

        for(i = 0; i < _sLen; i++)
            s[ t_pair_ary[i].second ] = t_pair_ary[i].first;

        delete[] t_pair_ary;


        /// compute sa, use the relative position
        relative_offset_type *sa = new relative_offset_type[_sLen];

        char * t_buf = new char[_sLen / 8 + 1];
        BitWrapper t(t_buf);


        t.set(_sLen - 1, S_TYPE);
        t.set(_sLen - 2, L_TYPE);
        i = _sLen - 3;

        for (; ; --i)
        {
            t.set(i, (_s[i] < _s[i + 1] || (_s[i] == _s[i + 1] && t.get(i + 1) == S_TYPE)) ? S_TYPE : L_TYPE);
            if (i == 0) break;
        }

        /// scan the s from left to right, to compute the LMS position
        relative_offset_type * _LMS_pos = new relative_offset_type [_lms_order->size()];

        _LMS_pos[_lms_order->size() - 1] = _sLen - 1;

        j = _lms_order->size() - 2;

        /// find lms-chars in s[1, _sLen - 3]
        for (i = _sLen - 3; i >= 1; --i)
        {
            if (t.get(i) && !t.get(i - 1))
            {
                _LMS_pos[j] = i;

                if(j == 0)break;
                else --j;
            }
        }

        /// store the LMS position from big to small
        MyVector<relative_offset_type> * sorted_lms_pos = new MyVector<relative_offset_type>();

        sorted_lms_pos->start_write();

        _lms_order->start_read_reverse();

        while( !_lms_order->is_eof() )
        {
            sorted_lms_pos->push_back( _LMS_pos[ _lms_order->get_reverse() - 1] );

            _lms_order->next_remove_reverse();
        }

        if(_lms_order->is_eof())
        {

            delete _lms_order;

            _lms_order = nullptr;
        }

        delete [] _LMS_pos;
        _LMS_pos = nullptr;



        // sort all the S-substrings
        relative_offset_type *bkt = new relative_offset_type[cur_name + 1];

        /// set the maximum value
        relative_offset_type MAX_VALUE(std::numeric_limits<relative_offset_type>::max());

        for (i = 0; i < _sLen; ++i) sa[i] = MAX_VALUE;//init sa

        /// step 1: initialize the S* tuple
        getBuckets<alphabet_type,relative_offset_type>(s, _sLen, bkt, cur_name + 1, true); // find end of buckets

        sorted_lms_pos->start_read();

        while(!sorted_lms_pos->is_eof())
        {

            j = sorted_lms_pos->get();

            sa[ bkt[s[j]] ] = j, --bkt[s[j]];

            sorted_lms_pos->next_remove();
        }

        /// delete the sorted_lms_pos
        if(sorted_lms_pos->is_eof())
        {
            delete sorted_lms_pos;

            sorted_lms_pos = nullptr;
        }

        if(is_virtual_sentinel)
        {

            ++bkt[s[j]], sa[ bkt[s[j]] ] = MAX_VALUE;  /// \note error-prone

            sa[0] = _sLen - 1;
        }


        /// generate the cLMS
        typedef Pair<alphabet_type,compress_type> lms_pair_type, bwt_pair_type;

        alphabet_type cur_lms(0),cur_bwt;

        compress_type cur_cnt(0), max_count(std::numeric_limits<compress_type>::max());


        for(i = 0; i < _sLen; ++i)
        {

            if(sa[i] != MAX_VALUE)
            {

                cur_lms = _s[sa[i]]; /// \note first item

                cur_cnt = 1;

                /*lms_rPos->push_back(sa[i]);*/
                _lms_rPos_vector->push_back(sa[i]);
            }

            for(++i; i < _sLen; ++i)
            {

                if(sa[i] != MAX_VALUE)
                {

                    if(cur_lms != _s[sa[i]])  /// original string
                    {

                        /*lms_cChar->push_back(lms_pair_type(cur_lms,cur_cnt));*/
                        _cLMS_vector->push_back(lms_pair_type(cur_lms,cur_cnt));


                        cur_lms = _s[sa[i]], cur_cnt = 1;

                    }
                    else
                    {

                        if(cur_cnt == max_count)
                        {

                            /*lms_cChar->push_back(lms_pair_type(cur_lms,cur_cnt));*/
                            _cLMS_vector->push_back(lms_pair_type(cur_lms,cur_cnt));

                            cur_cnt = 1;

                        }
                        else
                        {
                            cur_cnt++;
                        }

                    }

                    /*lms_rPos->push_back( sa[i] );*/
                    _lms_rPos_vector->push_back( sa[i] );

                }

            }

        }

        /*lms_cChar->push_back(lms_pair_type(cur_lms,cur_cnt));*/
        _cLMS_vector->push_back(lms_pair_type(cur_lms,cur_cnt));

        getBuckets<alphabet_type,relative_offset_type>(s, _sLen, bkt, cur_name + 1, false); // find heads of buckets

        for (i = 0; i < _sLen; ++i)
        {
            if (sa[i] != MAX_VALUE)// && sa[i] != 0) // sa[i] is never equal to 0, because the first item in _s is not L-type.
            {
                j = sa[i] - 1;

                if(!t.get(j))
                {
                    /// j - 1 is greater than or equal to zero
                    /*if (!t.get(j-1))l_bwt_type->push_back(false);   ///  t[j] = L-type
                    else l_bwt_type->push_back(true);*/ ///  t[j] = S-type

                    if (!t.get(j-1))_suf_l_bit_vector->push_back(false);   ///  t[j] = L-type
                    else _suf_l_bit_vector->push_back(true);

                    cur_bwt = _s[j],cur_cnt = 1; /// first item's BWT and counter; cur_bwt is the original char

                    sa[bkt[s[j]]] = j, ++bkt[s[j]]; /// the substitute string

                    sa[i] = MAX_VALUE;

                    for(i = i+1 ; i < _sLen; ++i)
                    {

                        if (sa[i] != MAX_VALUE)// && sa[i] != 0) // sa[i] is never equal to 0, because the first item in _s is not L-type.
                        {
                            j = sa[i] - 1;

                            if(!t.get(j))
                            {

                                /*if (!t.get(j-1))l_bwt_type->push_back(false);   ///  t[j] = L-type
                                else l_bwt_type->push_back(true);*/ ///  t[j] = S-type

                                if (!t.get(j-1))_suf_l_bit_vector->push_back(false);   ///  t[j] = L-type
                                else _suf_l_bit_vector->push_back(true); ///  t[j] = S-type

                                if(cur_bwt != _s[j])
                                {
                                    /*l_cBwt->push_back(bwt_pair_type(cur_bwt,cur_cnt));*/
                                    _suf_l_bwt_vector->push_back(bwt_pair_type(cur_bwt,cur_cnt));

                                    cur_bwt = _s[j], cur_cnt = 1;
                                }
                                else
                                {
                                    if(cur_cnt == max_count)
                                    {
                                        /*l_cBwt->push_back(bwt_pair_type(cur_bwt,cur_cnt));*/
                                        _suf_l_bwt_vector->push_back(bwt_pair_type(cur_bwt,cur_cnt));

                                        cur_cnt = 1;
                                    }
                                    else
                                    {
                                        cur_cnt++;
                                    }

                                }

                                sa[bkt[s[j]]] = j, ++bkt[s[j]]; /// the substitute string

                                sa[i] = MAX_VALUE;


                            }
                        }
                    }
                }
            }
        }

        /*l_cBwt->push_back(bwt_pair_type(cur_bwt,cur_cnt));*//// last cBWT pair
        _suf_l_bwt_vector->push_back(bwt_pair_type(cur_bwt,cur_cnt));/// last cBWT pair

        cur_bwt = 0, cur_cnt = 0;

        getBuckets<alphabet_type,relative_offset_type>(s, _sLen, bkt, cur_name + 1, true); // find ends of buckets

        for (i = _sLen - 1; ; --i) /// i > 0, because sa[0] is the current smallest item, whose preceding must be L-type
        {
            if (sa[i] && sa[i] != MAX_VALUE )
            {
                j = sa[i] - 1;

                if (t.get(j))   /// 'j' is greater than or equal to zero
                {
                    /// if j==0, t[j]'s preceding is L-type, i.e. t[j] is S* (the block starting and ending by S* char, exclude the leftmost block)
                    /*if(!j || !t.get(j-1)) s_bwt_type->push_back(true);
                    else s_bwt_type->push_back(false);*/ /// t[j] is not S*

                    if(!j || !t.get(j-1)) _suf_s_bit_vector->push_back(true);
                    else _suf_s_bit_vector->push_back(false);

                    cur_bwt = _s[j],cur_cnt = 1; /// first largest sBWT

                    sa[bkt[s[j]]] = j, --bkt[s[j]]; /// the substitute string

                    if(i == 0) break;

                    for(i = i-1; ; i--)
                    {

                        if (sa[i] && sa[i] != MAX_VALUE)
                        {
                            j = sa[i] - 1;

                            if (t.get(j))   /// t[j] == S-type
                            {

                                /*if(!j || !t.get(j-1)) s_bwt_type->push_back(true); /// t[j] is S*
                                else s_bwt_type->push_back(false);*/ /// t[j] is not S*

                                if(!j || !t.get(j-1)) _suf_s_bit_vector->push_back(true); /// t[j] is S*
                                else _suf_s_bit_vector->push_back(false); /// t[j] is not S*

                                if(cur_bwt != _s[j])
                                {
                                    /*s_cBwt->push_back(bwt_pair_type(cur_bwt,cur_cnt));*/
                                    _suf_s_bwt_vector->push_back(bwt_pair_type(cur_bwt,cur_cnt));

                                    cur_bwt = _s[j], cur_cnt = 1;
                                }
                                else
                                {
                                    if(cur_cnt == max_count)
                                    {
                                        /*s_cBwt->push_back(bwt_pair_type(cur_bwt,cur_cnt));*/
                                        _suf_s_bwt_vector->push_back(bwt_pair_type(cur_bwt,cur_cnt));

                                        cur_cnt = 1;
                                    }
                                    else
                                    {
                                        cur_cnt++;
                                    }

                                }

                                sa[bkt[s[j]]] = j, --bkt[s[j]]; /// the substitute string
                            }
                        }

                        if(i == 0) break;
                    }

                }
            }

            if (i == 0) break;
        }

        /*s_cBwt->push_back(bwt_pair_type(cur_bwt,cur_cnt));*//// last sBWT
        _suf_s_bwt_vector->push_back(bwt_pair_type(cur_bwt,cur_cnt));/// last sBWT



        t.~BitWrapper();
        delete[] sa;
        sa = nullptr;
        delete[] t_buf;
        t_buf = nullptr;
        delete[] bkt;
        bkt = nullptr;
        delete [] s; ///\note error-prone
        s = nullptr;

        delete[] _s;
        _s = nullptr;

        return;
    }

/// \note Reduction stage
/// \note alphabet <= 255 or > 255
/// \brief compute the singleton block BWT in reduction phase, and the computation process dose not like the operation in induction phase
    template<typename alphabet_type, typename compress_type, typename relative_offset_type>
    static void compute_singleton_bwt(
        std::string _fname
        , const uint64 _start_pos
        , const uint64 _end_pos
        , const alphabet_type & _alpha
        , const uint32 & _level
        , const uint32 _id
        , std::vector< std::vector<uint64> *>               & _LMS_counter
        , MyVector< Pair<alphabet_type, compress_type> >    * _cLMS_seq
        , MyVector< Pair<alphabet_type, compress_type> >    * _sub_l_bwt_vector
        , BitVectorWapper<uint64>                           * _sub_l_bit_vector
        , MyVector< Pair<alphabet_type, compress_type> >    * _sub_s_bwt_vector
        , BitVectorWapper<uint64>                           * _sub_s_bit_vector
        , MyVector<relative_offset_type>                    * _lms_rPos_vector
        , BitVectorWapper<uint64>                           * _lms_diff_bit_vector
        , bool _is_sentinel
    )
    {
        /// avoid the virtual sentinel
        bool is_virtual_sentinel = _is_sentinel;

        /*if(_level == 0 && _lms_rPos_vector.empty()) is_virtual_sentinel = true;*/

        uint64 sLen = 0;

        uint64 read_end_pos = 0; /// position of reading last char or last two char

        if(is_virtual_sentinel)
        {
            sLen = _end_pos - _start_pos + 1 - 1; ///avoid the virtual sentinel
            read_end_pos = _end_pos;
        }
        else
        {
            sLen = _end_pos - _start_pos + 1 - 2; /// don't include the last two char
            read_end_pos = _end_pos - 1;
        }

        alphabet_type cur_bwt(0);

        compress_type cur_cnt(1), max_count(std::numeric_limits<compress_type>::max());

        uint64 buf_size = MAX_MEM / 3 * 2 / sizeof(alphabet_type);

        alphabet_type * buf = new alphabet_type[buf_size];

        uint32 loopNum = sLen / buf_size;

        uint32 remainder = ( ( (sLen % buf_size) == 0 ) ? 0 : 1 );

        if(is_virtual_sentinel)
            UtilityFunctions::fileRead<alphabet_type>(buf,_fname,read_end_pos,1);/// avoid the virtual sentinel to read the last char
        else UtilityFunctions::fileRead<alphabet_type>(buf,_fname,read_end_pos,2); /// read the last two chars

        cur_bwt = buf[0], cur_cnt = 1;

        if(_alpha <= alphabet_type(255) )
        {

            ///step 1: initialize S*
            //for(uint32 i = 0; i <= _alpha; ++i) _LMS_counter[i]->push_back(0);

            /// if not sentinel, increase the corresponding counter by one; otherwise, omit to increase
            if(!is_virtual_sentinel)(*_LMS_counter[buf[1]])[_id]++;

        }
        else
        {

            ///step 1: initialize S*
            /*MyVector< Pair<alphabet_type, compress_type> > * lms_cChar = new MyVector< Pair<alphabet_type, compress_type> >();

            lms_cChar->start_write();

            if(!is_virtual_sentinel)lms_cChar->push_back( Pair<alphabet_type, compress_type>(buf[1],1) );
            else lms_cChar->push_back( Pair<alphabet_type, compress_type>(0,0) ); /// \note reference to sortLMS.h -> c_LMSHeapSorter()

            _cLMS_seq.push_back(lms_cChar);*/

            if(!is_virtual_sentinel)_cLMS_seq->push_back( Pair<alphabet_type, compress_type>(buf[1],1) );
            else _cLMS_seq->push_back( Pair<alphabet_type, compress_type>(0,0) );

        }
        ///step 2: write the LMS relative position
        /*MyVector<relative_offset_type> * lms_rPos = new MyVector<relative_offset_type>();
        lms_rPos->start_write();*/

        /// whether the current block includes a sentinel or not, the rLMS position must be 0 for a singleton block
        /*lms_rPos->push_back(0);*/

        /*_lms_rPos_vector.push_back(lms_rPos);*/
        _lms_rPos_vector->push_back(0);


        /*BitVectorWapper<uint64> * lms_bit_diff = new BitVectorWapper<uint64>();
        lms_bit_diff->start_write();*/

        /// the diff of a unique LMS substring must be 1 for a singleton block
        /*lms_bit_diff->push_back(1);*/

        /*_lms_diff_bit_vector.push_back(lms_bit_diff);*/
        _lms_diff_bit_vector->push_back(1);

        /// step 3: L- and S-type inducing
        typedef Pair<alphabet_type,compress_type> pair_type;

        /*MyVector<pair_type> * l_cBwt = new MyVector<pair_type>();
        l_cBwt->start_write();

        BitVectorWapper<uint64> * l_bwt_type = new BitVectorWapper<uint64>();
        l_bwt_type->start_write();

        MyVector<pair_type> * s_cBwt = new MyVector<pair_type>();
        s_cBwt->start_write();

        BitVectorWapper<uint64> * s_bwt_type = new BitVectorWapper<uint64>();
        s_bwt_type->start_write();*/



        uint8 is_s_type = false; /// the first S-type char occurs

        /// for read_size == buf_size
        for(uint32 i = 1; i <= loopNum; i++)
        {

            uint64 read_size = buf_size;

            uint64 read_beg_pos = (read_end_pos - 1) - read_size * i + 1; /// avoid the last two char


            UtilityFunctions::fileRead(buf,_fname,read_beg_pos,read_size);


            for(uint64 j(read_size - 1); j >= 0; --j)
            {

                if( !is_s_type ) /// L-BWT
                {

                    if( buf[j] < cur_bwt)
                    {

                        is_s_type = true;

                        /*l_cBwt->push_back(pair_type(cur_bwt,cur_cnt));
                        l_bwt_type->push_back(true);*/

                        _sub_l_bwt_vector->push_back(pair_type(cur_bwt,cur_cnt));
                        _sub_l_bit_vector->push_back(true);


                        cur_bwt = buf[j], cur_cnt = 1;


                    }
                    else
                    {

                        if( buf[j] == cur_bwt )
                        {

                            if(cur_cnt == max_count)
                            {

                                /*l_cBwt->push_back(pair_type(cur_bwt,cur_cnt));*/
                                _sub_l_bwt_vector->push_back(pair_type(cur_bwt,cur_cnt));

                                cur_cnt = 1;
                            }
                            else
                            {

                                cur_cnt++;
                            }

                        }
                        else
                        {

                            /*l_cBwt->push_back(pair_type(cur_bwt,cur_cnt));*/
                            _sub_l_bwt_vector->push_back(pair_type(cur_bwt,cur_cnt));

                            cur_bwt = buf[j], cur_cnt = 1;

                        }

                        /*l_bwt_type->push_back(false);*/
                        _sub_l_bit_vector->push_back(false);
                    }

                }
                else
                {

                    if( buf[j] == cur_bwt )
                    {

                        if(cur_cnt == max_count)
                        {

                            /*s_cBwt->push_back(pair_type(cur_bwt,cur_cnt));*/
                            _sub_s_bwt_vector->push_back(pair_type(cur_bwt,cur_cnt));

                            cur_cnt = 1;
                        }
                        else
                        {

                            cur_cnt++;
                        }

                    }
                    else
                    {

                        /*s_cBwt->push_back(pair_type(cur_bwt,cur_cnt));*/
                        _sub_s_bwt_vector->push_back(pair_type(cur_bwt,cur_cnt));

                        cur_bwt = buf[j], cur_cnt = 1;

                    }

                    /*s_bwt_type->push_back(false);*/
                    _sub_s_bit_vector->push_back(false);


                }

                if(j == 0)break;
            }


        }


        /// for read_size == remainder
        if(remainder)
        {

            uint64 read_size = sLen - loopNum * buf_size;

            UtilityFunctions::fileRead(buf,_fname,_start_pos,read_size);

            for(uint64 j(read_size - 1); j >= 0; --j)
            {

                if( !is_s_type ) /// LBWT
                {

                    if( buf[j] < cur_bwt)
                    {

                        is_s_type = true;

                        /*l_cBwt->push_back(pair_type(cur_bwt,cur_cnt));
                        l_bwt_type->push_back(true);*/

                        _sub_l_bwt_vector->push_back(pair_type(cur_bwt,cur_cnt));
                        _sub_l_bit_vector->push_back(true);

                        cur_bwt = buf[j], cur_cnt = 1;


                    }
                    else
                    {

                        if( buf[j] == cur_bwt )
                        {

                            if(cur_cnt == max_count)
                            {

                                /*l_cBwt->push_back(pair_type(cur_bwt,cur_cnt));*/
                                _sub_l_bwt_vector->push_back(pair_type(cur_bwt,cur_cnt));

                                cur_cnt = 1;
                            }
                            else
                            {

                                cur_cnt++;
                            }

                        }
                        else
                        {

                            /*l_cBwt->push_back(pair_type(cur_bwt,cur_cnt));*/
                            _sub_l_bwt_vector->push_back(pair_type(cur_bwt,cur_cnt));

                            cur_bwt = buf[j], cur_cnt = 1;

                        }

                        /*l_bwt_type->push_back(false);*/
                        _sub_l_bit_vector->push_back(false);
                    }

                }
                else
                {

                    if( buf[j] == cur_bwt )
                    {

                        if(cur_cnt == max_count)
                        {

                            /*s_cBwt->push_back(pair_type(cur_bwt,cur_cnt));*/
                            _sub_s_bwt_vector->push_back(pair_type(cur_bwt,cur_cnt));

                            cur_cnt = 1;
                        }
                        else
                        {

                            cur_cnt++;
                        }

                    }
                    else
                    {

                        /*s_cBwt->push_back(pair_type(cur_bwt,cur_cnt));*/
                        _sub_s_bwt_vector->push_back(pair_type(cur_bwt,cur_cnt));

                        cur_bwt = buf[j], cur_cnt = 1;

                    }

                    /*s_bwt_type->push_back(false);*/
                    _sub_s_bit_vector->push_back(false);


                }

                if(j == 0) break;
            }


        }///

        /*s_cBwt->push_back(pair_type(cur_bwt,cur_cnt));
        s_bwt_type->push_back(true);*/

        _sub_s_bwt_vector->push_back(pair_type(cur_bwt,cur_cnt));
        _sub_s_bit_vector->push_back(true);

        /*_sub_l_bwt_vector.push_back(l_cBwt);
        _sub_l_bit_vector.push_back(l_bwt_type);

        _sub_s_bwt_vector.push_back(s_cBwt);
        _sub_s_bit_vector.push_back(s_bwt_type);*/

        delete [] buf;

    }


/// \note Induction stage
/// \note alphabet <= 255 or > 255
/// \brief compute the singleton-block/leftmost-block BWT in induction phase
    template<typename alphabet_type, typename compress_type, typename relative_offset_type>
    static void compute_singleton_bwt_induce(
        std::string                                                _fname
        , const uint64                                               _start_pos
        , const uint64                                               _end_pos
        , const alphabet_type                                        _alpha
        , const uint32                                               _level
        , const uint32                                               _id
        , std::vector<std::vector<uint64> *>                       & _LMS_counter
        , MyVector< Pair<alphabet_type, compress_type> >           * _cLMS_seq
        , MyVector< Pair<alphabet_type, compress_type> >           * _suf_l_cbwt_vector
        , BitVectorWapper<uint64>                                  * _suf_l_bit_vector
        , MyVector< Pair<alphabet_type, compress_type> >           * _suf_s_cbwt_vector
        , BitVectorWapper<uint64>                                  * _suf_s_bit_vector
        , MyVector<relative_offset_type>                           * _lms_rPos_vector
        , bool _is_sentinel
    )
    {
        /// avoid the virtual sentinel
        bool is_virtual_sentinel = _is_sentinel;

        std::cout << "the leftmost block's _is_sentinel = " << _is_sentinel << std::endl;


        //if(_level == 0 && _lms_rPos_vector.empty()) is_virtual_sentinel = true;

        uint64 sLen = 0;

        uint64 read_end_pos = 0; /// position of reading last char or last two char

        if(is_virtual_sentinel)
        {
            sLen = _end_pos - _start_pos + 1 - 1; ///avoid the virtual sentinel
            read_end_pos = _end_pos;
        }
        else
        {
            sLen = _end_pos - _start_pos + 1 - 2; /// don't include the last two char
            read_end_pos = _end_pos - 1;
        }

//        std::cout <<" sLen = " << sLen << std::endl;

        alphabet_type cur_bwt(0);

        compress_type cur_cnt(1), max_count(std::numeric_limits<compress_type>::max());

        uint64 buf_size = MAX_MEM / 3 * 2 / sizeof(alphabet_type);
        //uint64 buf_size = 2;

        alphabet_type * buf = new alphabet_type[buf_size];

        uint32 loopNum = sLen / buf_size;

//        std::cout << "loopNum = " << loopNum << std::endl;

        uint32 remainder = ( ( (sLen % buf_size) == 0 ) ? 0 : 1 );

        if(is_virtual_sentinel)
            UtilityFunctions::fileRead<alphabet_type>(buf,_fname,read_end_pos,1); /// avoid the virtual sentinel to read the last char
        else UtilityFunctions::fileRead<alphabet_type>(buf,_fname,read_end_pos,2); /// read the last two chars

        cur_bwt = buf[0], cur_cnt = 1;

//        // \note retrieve BWT
//        if(!_level && is_getBWT && !is_virtual_sentinel){
//
//            ((MyVector<uint8> *)L_eBWT_seqs[_id])->push_back(buf[0]);
//
//        }


//        std::cout << "buf[0] = " << buf[0] - 0 << ", buf[1] = " << buf[1] - 0 << std::endl;

        if(_alpha <= alphabet_type(255) )
        {
            std::cout << "_alpha = " << _alpha - 0 << std::endl;

            ///step 1: initialize S*
            //for(uint32 i = 0; i <= _alpha; ++i) _LMS_counter[i]->push_back(0);

            /// if not sentinel, increase the corresponding counter by one; otherwise, omit to increase
            /*if(!is_virtual_sentinel)_LMS_counter[buf[1]]->back()++;*/
            if(!is_virtual_sentinel)(*_LMS_counter[buf[1]])[_id]++;



        }
        else
        {

            ///step 1: initialize S*
            /*MyVector< Pair<alphabet_type, compress_type> > * lms_cChar = new MyVector< Pair<alphabet_type, compress_type> >();

            lms_cChar->start_write();

            if(!is_virtual_sentinel)lms_cChar->push_back( Pair<alphabet_type, compress_type>(buf[1],1) );
            else lms_cChar->push_back( Pair<alphabet_type, compress_type>(0,1) ); /// \note error-prone reference to sortLMS.h -> c_LMSHeapSorter()

            _cLMS_seq.push_back(lms_cChar);*/

            if(!is_virtual_sentinel)_cLMS_seq->push_back(Pair<alphabet_type, compress_type>(buf[1],1));
            else _cLMS_seq->push_back(Pair<alphabet_type, compress_type>(0,1)); /// \note error-prone reference to sortLMS.h -> c_LMSHeapSorter()

        }
        ///step 2: write the LMS relative position
        /*MyVector<relative_offset_type> * lms_rPos = new MyVector<relative_offset_type>();
        lms_rPos->start_write();

        ///
        if(!is_virtual_sentinel)lms_rPos->push_back( _end_pos - _start_pos );
        else lms_rPos->push_back( _end_pos - _start_pos  + 1); /// for sentinel position

        _lms_rPos_vector.push_back(lms_rPos);*/

        if(!is_virtual_sentinel)_lms_rPos_vector->push_back(_end_pos - _start_pos);
        else _lms_rPos_vector->push_back(_end_pos - _start_pos  + 1);

        /// step 3: L- and S-type inducing
        typedef Pair<alphabet_type,compress_type> pair_type;

        /*MyVector<pair_type> * l_cBwt = new MyVector<pair_type>();
        l_cBwt->start_write();

        BitVectorWapper<uint64> * l_bwt_type = new BitVectorWapper<uint64>();
        l_bwt_type->start_write();

        MyVector<pair_type> * s_cBwt = new MyVector<pair_type>();
        s_cBwt->start_write();

        BitVectorWapper<uint64> * s_bwt_type = new BitVectorWapper<uint64>();
        s_bwt_type->start_write();*/



        uint8 is_s_type = false; /// the first S-type char occurs

        /// for read_size == buf_size
        for(uint32 i = 1; i <= loopNum; i++)
        {

            uint64 read_size = buf_size;

            uint64 read_beg_pos = (read_end_pos - 1) - read_size * i + 1; /// avoid the last two char


            UtilityFunctions::fileRead(buf,_fname,read_beg_pos,read_size);


            for(uint64 j(read_size - 1); j >= 0; --j)
            {

                if( !is_s_type ) /// L-BWT
                {

                    if( buf[j] < cur_bwt)
                    {

                        is_s_type = true;

                        /*l_cBwt->push_back(pair_type(cur_bwt,cur_cnt));
                        l_bwt_type->push_back(true);*/

                        _suf_l_cbwt_vector->push_back(pair_type(cur_bwt,cur_cnt));
                        _suf_l_bit_vector->push_back(true);


                        cur_bwt = buf[j], cur_cnt = 1;

                        /// \note get LStar preceding char
                        if(!_level && is_getBWT)
                            ((MyVector<uint8> *)LStarPreChar_seqs[_id])->push_back(buf[j]);


                    }
                    else
                    {

                        if( buf[j] == cur_bwt )
                        {

                            if(cur_cnt == max_count)
                            {

                                /*l_cBwt->push_back(pair_type(cur_bwt,cur_cnt));*/
                                _suf_l_cbwt_vector->push_back(pair_type(cur_bwt,cur_cnt));

                                cur_cnt = 1;
                            }
                            else
                            {

                                cur_cnt++;
                            }

                        }
                        else
                        {
                            /*l_cBwt->push_back(pair_type(cur_bwt,cur_cnt));*/
                            _suf_l_cbwt_vector->push_back(pair_type(cur_bwt,cur_cnt));

                            cur_bwt = buf[j], cur_cnt = 1;

                        }

                        /*l_bwt_type->push_back(false);*/
                        _suf_l_bit_vector->push_back(false);

                    }

                }
                else
                {

                    if( buf[j] == cur_bwt )
                    {

                        if(cur_cnt == max_count)
                        {

                            /*s_cBwt->push_back(pair_type(cur_bwt,cur_cnt));*/
                            _suf_s_cbwt_vector->push_back(pair_type(cur_bwt,cur_cnt));

                            cur_cnt = 1;
                        }
                        else
                        {

                            cur_cnt++;
                        }

                    }
                    else
                    {

                        /*s_cBwt->push_back(pair_type(cur_bwt,cur_cnt));*/
                        _suf_s_cbwt_vector->push_back(pair_type(cur_bwt,cur_cnt));

                        cur_bwt = buf[j], cur_cnt = 1;

                    }

                    /*s_bwt_type->push_back(false);*/
                    _suf_s_bit_vector->push_back(false);


                }

                if(j == 0)break;
            }


        }


        /// for read_size == remainder
        if(remainder)
        {

            //Logger::output_separator_line("remainder");

            uint64 read_size = sLen - loopNum * buf_size;

            UtilityFunctions::fileRead(buf,_fname,_start_pos,read_size);

            for(uint64 j(read_size - 1); j >= 0; --j)
            {

                if( !is_s_type ) /// LBWT
                {



                    if( buf[j] < cur_bwt)
                    {

                        is_s_type = true;

                        /*l_cBwt->push_back(pair_type(cur_bwt,cur_cnt));
                        l_bwt_type->push_back(true);*/

                        _suf_l_cbwt_vector->push_back(pair_type(cur_bwt,cur_cnt));
                        _suf_l_bit_vector->push_back(true);


                        cur_bwt = buf[j], cur_cnt = 1;

                        /// \note get LStar preceding char
                        if(!_level && is_getBWT)
                            ((MyVector<uint8> *)LStarPreChar_seqs[_id])->push_back(buf[j]);


                    }
                    else
                    {

                        if( buf[j] == cur_bwt )
                        {

                            if(cur_cnt == max_count)
                            {

                                /*l_cBwt->push_back(pair_type(cur_bwt,cur_cnt));*/
                                _suf_l_cbwt_vector->push_back(pair_type(cur_bwt,cur_cnt));

                                cur_cnt = 1;
                            }
                            else
                            {

                                cur_cnt++;
                            }

                        }
                        else
                        {

                            /*l_cBwt->push_back(pair_type(cur_bwt,cur_cnt));*/
                            _suf_l_cbwt_vector->push_back(pair_type(cur_bwt,cur_cnt));

                            cur_bwt = buf[j], cur_cnt = 1;

                        }

                        /*l_bwt_type->push_back(false);*/
                        _suf_l_bit_vector->push_back(false);


                    }

                }
                else
                {

                    if( buf[j] == cur_bwt )
                    {

                        if(cur_cnt == max_count)
                        {

                            /*s_cBwt->push_back(pair_type(cur_bwt,cur_cnt));*/
                            _suf_s_cbwt_vector->push_back(pair_type(cur_bwt,cur_cnt));

                            cur_cnt = 1;
                        }
                        else
                        {

                            cur_cnt++;
                        }

                    }
                    else
                    {

                        /*s_cBwt->push_back(pair_type(cur_bwt,cur_cnt));*/
                        _suf_s_cbwt_vector->push_back(pair_type(cur_bwt,cur_cnt));

                        cur_bwt = buf[j], cur_cnt = 1;

                    }

                    /*s_bwt_type->push_back(false);*/
                    _suf_s_bit_vector->push_back(false);


                }

                if(j == 0) break;
            }


        }///


        if(_start_pos != 0) ///  a singleton block
        {


            if( !_level && is_getBWT )
            {
                uint8 * tp_char = new uint8[1];

                UtilityFunctions::fileRead<uint8>(tp_char,_fname,_start_pos-1,1);

                ((MyVector<uint8> *)SStarPreChar_seqs[_id])->push_back(tp_char[0]);

                delete tp_char;
            }

            /*s_cBwt->push_back(pair_type(cur_bwt,cur_cnt));
            s_bwt_type->push_back(true);*/

            _suf_s_cbwt_vector->push_back(pair_type(cur_bwt,cur_cnt));
            _suf_s_bit_vector->push_back(true);


        }
        else   /// a leftmost block
        {

            if(is_s_type) /// include S- and L-type BWT
            {

                /*s_cBwt->push_back(pair_type(cur_bwt,cur_cnt));
                s_bwt_type->push_back(true);*/ /// when inducing s-type suffixes, s[0]'s preceding type is set as LMS, for not inducing its preceding

                _suf_s_cbwt_vector->push_back(pair_type(cur_bwt,cur_cnt));
                _suf_s_bit_vector->push_back(true); /// when inducing s-type suffixes, s[0]'s preceding type is set as LMS, for not inducing its preceding
                //_suf_s_bit_vector->push_back(false); /// revised on 2021-0411
            }
            else  /// only include L-type BWT
            {

                /*l_cBwt->push_back(pair_type(cur_bwt,cur_cnt));
                l_bwt_type->push_back(false);*/ /// if it is set to true, it will be used to induce its preceding during inducing s-type suffixes

                _suf_l_cbwt_vector->push_back(pair_type(cur_bwt,cur_cnt));
                _suf_l_bit_vector->push_back(false); /// if it is set to true, it will be used to induce its preceding during inducing s-type suffixes

            }


        }

        /*_suf_l_cbwt_vector.push_back(l_cBwt);
        _suf_l_bit_vector.push_back(l_bwt_type);

        _suf_s_cbwt_vector.push_back(s_cBwt);
        _suf_s_bit_vector.push_back(s_bwt_type);*/
        if( !_level && is_getBWT )
        {
            std::cout << "block id = " << _id << ", lms number is " << ((MyVector<uint8> *)SStarPreChar_seqs[_id])->size() << std::endl;
            std::cout << "block id = " << _id << ", lml number is " << ((MyVector<uint8> *)LStarPreChar_seqs[_id])->size() << std::endl;
        }
        delete [] buf;


    }


    //@usage: compute bucket size.
    template<typename alphabet_type,typename relative_offset_type>
    static void getBuckets_seq(const alphabet_type * const & _s, const relative_offset_type _sLen, relative_offset_type * &_bkt, const relative_offset_type _bktNum, const bool _end)
    {
        relative_offset_type i;
        relative_offset_type sum = 0;
        for (i = 0; i < _bktNum; ++i) _bkt[i] = 0; // clear all buckets
        for (i = 0; i < _sLen; ++i) ++_bkt[_s[i]]; // compute the size of each bucket
        for (i = 0; i < _bktNum; ++i)
        {
            sum += _bkt[i];
            _bkt[i] = _end ? sum - 1 : sum - _bkt[i];
        }
    }

    template<typename alphabet_type,typename relative_offset_type>
    static void getBuckets(const alphabet_type * const &s, const relative_offset_type n,
                           relative_offset_type * & bkt,
                           relative_offset_type K, const bool end)
    {

        if(sizeof(alphabet_type) > 1)
        {
            getBuckets_seq(s,n,bkt,K,end);
            return ;
        }

#ifndef CILK
        getBuckets_seq(s,n,bkt,K,end);
        return ;
#endif // CILK


        relative_offset_type i, sum=0;

        // clear all buckets
        for (i = 0; i < K; i++) bkt[i] = 0;

        std::vector<relative_offset_type> bkt_temp(K, 0);

        cilk_scan(n, bkt_temp, cs_buf,
                  [ = ](const size_t begin, const size_t size, std::vector<relative_offset_type> & initial)
        {
            initial.assign(K, 0);
            size_t end = begin + size;
            for (unsigned int i = begin; i < end; i++) initial[s[i]]++;
        },
        [ = ](std::vector<relative_offset_type> & x, const std::vector<relative_offset_type> & y)
        {
            for (size_t i = 0; i < K; i++) x[i] += y[i];
        },
        [ =, &bkt_temp ](const size_t begin, const size_t size, std::vector<relative_offset_type> & initial)
        {
            size_t end = begin + size;
            if (end == n)
            {
                for (unsigned int i = begin; i < end; i++) initial[s[i]]++;
                // Save result from last tile
                bkt_temp.swap(initial);
            }
        }
                 );

        for (i = 0; i < K; i++)
        {
            sum += bkt_temp[i];
            bkt[i] = end ? sum - 1 : sum - bkt_temp[i];
        }

    }

    template<typename alphabet_type,typename relative_offset_type>
    static void putSubstr0(relative_offset_type *SA,
                           alphabet_type *s, relative_offset_type *bkt,
                           relative_offset_type n, relative_offset_type K, bool _is_sentinel)
    {

        //#ifdef _verbose
        //	otherTimer.start();
        //#endif

        // find the end of each bucket.
        getBuckets<alphabet_type,relative_offset_type>(s, n, bkt, K, true);

        // set each item in SA as empty.
        //for (relative_offset_type i = 0; i < n; i++) SA[i] = 0;

        // set the last LMS char
        if(!_is_sentinel)
        {
            SA[bkt[s[n-1]]] = n - 1;
            bkt[s[n-1]]--;
        }
        else
        {
            if( sizeof(alphabet_type) == 1 )SA[0] = n - 1;
        }

        //std::vector<size_t>
        std::vector<relative_offset_type> bkt_temp(K, 0);
        std::cout << "bkt_temp.size() = " << bkt_temp.size() << std::endl;

        cilk_scan( n - 2, bkt_temp, cs_buf,
                   [ = ](const size_t begin, const size_t size,  std::vector<relative_offset_type> & initial)
        {
            initial.assign(K, 0);
            size_t end = begin + 1 + size;
            relative_offset_type i;
            alphabet_type ch;
            bool hasL = false;
            for (i = begin + 1; i < end; i++)
            {
                if (hasL)
                {
                    if (s[i - 1] < s[i])
                    {
                        hasL = false;
                        initial[ch]++;
                    }
                    else if (s[i - 1] > s[i]) ch = s[i];
                }
                else
                {
                    if (s[i - 1] > s[i])
                    {
                        hasL = true;
                        ch = s[i];
                    }
                }
            }
            while (hasL && (i < n - 1))
            {
                if (s[i - 1] == s[i]) i++;
                else
                {
                    if (s[i - 1] < s[i]) initial[ch]++;
                    hasL = false;
                }
            }
        },
        [ = ]( std::vector<relative_offset_type> & x, const  std::vector<relative_offset_type> & y)
        {
            for (size_t i = 0; i < K; i++) x[i] += y[i];
        },
        [ = ](const size_t begin, const size_t size,  std::vector<relative_offset_type> initial)
        {
            size_t end = begin + 1 + size;
            relative_offset_type i, j;
            alphabet_type ch;
            bool hasL = false;
            for (i = begin + 1; i < end; i++)
            {
                if (hasL)
                {
                    if (s[i - 1] < s[i])
                    {
                        hasL = false;
                        SA[bkt[ch] - initial[ch]++] = j;
                    }
                    else if (s[i - 1] > s[i])
                    {
                        ch = s[i];
                        j = i;
                    }
                }
                else
                {
                    if (s[i - 1] > s[i])
                    {
                        hasL = true;
                        ch = s[i];
                        j = i;
                    }
                }
            }
            while (hasL && (i < n - 1))
            {
                if (s[i - 1] == s[i]) i++;
                else
                {
                    if (s[i - 1] < s[i]) SA[bkt[ch] - initial[ch]++] = j;
                    hasL = false;
                }
            }
        }
                 );

        // set the single sentinel LMS-substring.
        //if(_is_sentinel)SA[0] = n - 1;
        std::cout << "Line is " << __LINE__ << std::endl;
    }

    template<typename alphabet_type,typename relative_offset_type>
    static void induceSAl0(relative_offset_type *SA,
                           alphabet_type *s, relative_offset_type *bkt,
                           relative_offset_type n, relative_offset_type K, bool suffix,  bool _is_sentinel)
    {

        relative_offset_type MAX_VALUE(std::numeric_limits<relative_offset_type>::max());

        // find the head of each bucket.
        getBuckets<alphabet_type,relative_offset_type>(s, n, bkt, K, false);
        if(_is_sentinel)bkt[0]++; // skip the virtual sentinel.

        size_t len = cs_buf * 2 * 4;
        relative_offset_type *buf = new relative_offset_type[len];
        for (size_t i = 0; i < len; i++) buf[i] = 0;

        cilk_pipeline(0, n, cs_buf,
                      [ = ](const size_t begin, const size_t end, const size_t id)
        {

            relative_offset_type *buf_cur = buf + 2 * cs_buf * (id % 4);

            if (id > 1)
            {
                relative_offset_type *buf_cur1 = buf + 2 * cs_buf * ((id - 2) % 4);
                parallel_for(unsigned int i = 0; i < cs_buf; i++)
                {
                    unsigned int offset = 2 * i;
                    if (buf_cur1[offset + 1] > 0)
                    {
                        SA[buf_cur1[offset]] = buf_cur1[offset + 1] - 1;
                        if (suffix)
                        {
                            buf_cur1[offset] = 0;
                            buf_cur1[offset + 1] = 0;
                        }
                    }
                }
                if (!suffix)
                {
                    size_t begin1 = 0 + (id - 2) * cs_buf;
                    size_t end1 = begin1 + cs_buf;
                    parallel_for(unsigned int i = begin1; i < end1; i++)
                    {
                        unsigned int offset = 2 * (i - begin1);
                        if (buf_cur1[offset + 1] > 0)
                        {
                            buf_cur1[offset] = 0;
                            buf_cur1[offset + 1] = 0;
                            if (i > 0) SA[i] = MAX_VALUE; // set to be empty
                        }
                    }
                }
            }

            parallel_for(unsigned int i = begin; i < end; i++)
            {
                relative_offset_type j = SA[i], offset = 2 * (i - begin);
                alphabet_type c, c1;
                if (j > 0 and j != MAX_VALUE) // avoid MAX_VALUE
                {
                    c = s[j - 1];
                    c1 = s[j];
                    if (c >= c1)
                    {
                        buf_cur[offset] = c;
                        buf_cur[offset + 1] = j;

                        SA[i] = MAX_VALUE; // set to MAX_VALUE
                    }

                }
            }
        },
        [ = ](const size_t begin, const size_t end, const size_t id)
        {

            unsigned int i, j, offset, chr, chr1, pos;
            size_t begin1 = end;
            size_t end1 = ((n - begin1) > cs_buf) ? (begin1 + cs_buf) : n;
            relative_offset_type *buf_cur = buf + 2 * cs_buf * (id % 4);
            relative_offset_type *buf_cur1 = buf + 2 * cs_buf * ((id + 1) % 4);

            for (i = begin; i < end; i++)
            {
                offset = 2 * (i - begin);
                if (buf_cur[offset + 1] > 0)
                {
                    j = buf_cur[offset + 1] - 1;
                    chr = buf_cur[offset];
                    pos = bkt[chr]++;
                    buf_cur[offset] = pos;

                    if (pos < end1 && j > 0 && ((chr1 = s[j - 1]) >= chr))
                    {
                        if (pos < end)
                        {
                            offset = 2 * (pos - begin);
                            buf_cur[offset] = chr1;
                            buf_cur[offset + 1] = j;
                        }
                        else
                        {
                            offset = 2 * (pos - begin1);
                            buf_cur1[offset] = chr1;
                            buf_cur1[offset + 1] = j;
                        }
                    }
                }
            }

            if ((n - begin) <= cs_buf)
            {
                if (id > 0)
                {
                    relative_offset_type *buf_cur2 = buf + 2 * cs_buf * ((id - 1) % 4);
                    parallel_for(unsigned int i = 0; i < cs_buf; i++)
                    {
                        unsigned int offset = 2 * i;
                        if (buf_cur2[offset + 1] > 0)
                        {
                            SA[buf_cur2[offset]] = buf_cur2[offset + 1] - 1;
                            if (suffix)
                            {
                                buf_cur2[offset] = 0;
                                buf_cur2[offset + 1] = 0;
                            }
                        }
                    }
                    if (!suffix)
                    {
                        size_t begin2 = 0 + (id - 1) * cs_buf;
                        size_t end2 = begin2 + cs_buf;
                        parallel_for(unsigned int i = begin2; i < end2; i++)
                        {
                            unsigned int offset = 2 * (i - begin2);
                            if (buf_cur2[offset + 1] > 0)
                            {
                                buf_cur2[offset] = 0;
                                buf_cur2[offset + 1] = 0;
                                if (i > 0) SA[i] = MAX_VALUE; // set to be empty
                            }
                        }
                    }
                }

                parallel_for(unsigned int i = begin; i < end; i++)
                {
                    unsigned int offset = 2 * (i - begin);
                    if (buf_cur[offset + 1] > 0)
                    {
                        SA[buf_cur[offset]] = buf_cur[offset + 1] - 1;
                        if (suffix)
                        {
                            buf_cur[offset] = 0;
                            buf_cur[offset + 1] = 0;
                        }
                    }
                }
                if (!suffix)
                {
                    parallel_for(unsigned int i = begin; i < end; i++)
                    {
                        unsigned int offset = 2 * (i - begin);
                        if (buf_cur[offset + 1] > 0)
                        {
                            buf_cur[offset] = 0;
                            buf_cur[offset + 1] = 0;
                            if (i > 0) SA[i] = MAX_VALUE; // set to be empty
                        }
                    }
                }
            }
        }
                     );

        delete[] buf; // read and write buffer
        buf = nullptr;
    }


    template<typename alphabet_type,typename relative_offset_type>
    static void induceSAs0(relative_offset_type *SA,
                           alphabet_type *s, relative_offset_type *bkt,
                           relative_offset_type n, relative_offset_type K, bool suffix)
    {

        relative_offset_type MAX_VALUE(std::numeric_limits<relative_offset_type>::max());
        // find the end of each bucket.
        getBuckets<alphabet_type,relative_offset_type>(s, n, bkt, K, true);

        size_t len = cs_buf * 2 * 4;
        relative_offset_type *buf = new relative_offset_type[len];
        for (size_t i = 0; i < len; i++) buf[i] = 0;

        cilk_pipeline(n - 1, 0, cs_buf,
                      [ = ](const size_t begin, const size_t end, const size_t id)
        {

            relative_offset_type *buf_cur = buf + 2 * cs_buf * (id % 4);

            if (id > 1)
            {
                relative_offset_type *buf_cur1 = buf + 2 * cs_buf * ((id - 2) % 4);
                parallel_for(unsigned int i = 0; i < cs_buf; i++)
                {
                    unsigned int offset = 2 * i;
                    if (buf_cur1[offset + 1] > 0)
                    {
                        SA[buf_cur1[offset]] = buf_cur1[offset + 1] - 1;
                        if (suffix)
                        {
                            buf_cur1[offset] = 0;
                            buf_cur1[offset + 1] = 0;
                        }
                    }
                }
                if (!suffix)
                {
                    size_t begin1 = (n - 1) - (id - 2) * cs_buf;
                    size_t end1 = begin1 - cs_buf;
                    parallel_for(unsigned int i = begin1; i > end1; i--)
                    {
                        unsigned int offset = 2 * (begin1 - i);
                        if (buf_cur1[offset + 1] > 0)
                        {
                            buf_cur1[offset] = 0;
                            buf_cur1[offset + 1] = 0;
                            SA[i] = MAX_VALUE;
                        }
                    }
                }
            }

            parallel_for(unsigned int i = begin; i > end; i--)
            {
                relative_offset_type j = SA[i], offset = 2 * (begin - i);
                alphabet_type c, c1;
                if (j > 0 && j != MAX_VALUE) // avoid MAX_VALUE
                {
                    c = s[j - 1];
                    c1 = s[j];
                    if (c < c1 || (c == c1 && bkt[c] < i))
                    {
                        buf_cur[offset] = c;
                        buf_cur[offset + 1] = j;

                        SA[i] = MAX_VALUE; // set to MAX_VALUE
                    }
                }
            }
        },
        [ = ](const size_t begin, const size_t end, const size_t id)
        {

            unsigned int i, j, offset, chr, chr1, pos;
            int begin1 = end;
            int end1 = (begin1 > cs_buf) ? (begin1 - cs_buf) : 0;
            relative_offset_type *buf_cur = buf + 2 * cs_buf * (id % 4);
            relative_offset_type *buf_cur1 = buf + 2 * cs_buf * ((id + 1) % 4);

            for (i = begin; i > end; i--)
            {
                offset = 2 * (begin - i);
                if (buf_cur[offset + 1] > 0)
                {
                    j = buf_cur[offset + 1] - 1;
                    chr = buf_cur[offset];
                    pos = bkt[chr]--;
                    buf_cur[offset] = pos;

                    if (pos > end1 && j > 0 && ((chr1 = s[j - 1]) < chr || (chr1 == chr && bkt[chr1] < pos)))
                    {
                        if (pos > end)
                        {
                            offset = 2 * (begin - pos);
                            buf_cur[offset] = chr1;
                            buf_cur[offset + 1] = j;
                        }
                        else
                        {
                            offset = 2 * (begin1 - pos);
                            buf_cur1[offset] = chr1;
                            buf_cur1[offset + 1] = j;
                        }
                    }
                }
            }

            if (begin <= cs_buf)
            {
                if (id > 0)
                {
                    relative_offset_type *buf_cur2 = buf + 2 * cs_buf * ((id - 1) % 4);
                    parallel_for(unsigned int i = 0; i < cs_buf; i++)
                    {
                        unsigned int offset = 2 * i;
                        if (buf_cur2[offset + 1] > 0)
                        {
                            SA[buf_cur2[offset]] = buf_cur2[offset + 1] - 1;
                            if (suffix)
                            {
                                buf_cur2[offset] = 0;
                                buf_cur2[offset + 1] = 0;
                            }
                        }
                    }
                    if (!suffix)
                    {
                        size_t begin2 = (n - 1) - (id - 1) * cs_buf;
                        size_t end2 = begin2 - cs_buf;
                        parallel_for(unsigned int i = begin2; i > end2; i--)
                        {
                            unsigned int offset = 2 * (begin2 - i);
                            if (buf_cur2[offset + 1] > 0)
                            {
                                buf_cur2[offset] = 0;
                                buf_cur2[offset + 1] = 0;
                                SA[i] = MAX_VALUE;
                            }
                        }
                    }
                }

                parallel_for(unsigned int i = begin; i > end; i--)
                {
                    unsigned int offset = 2 * (begin - i);
                    if (buf_cur[offset + 1] > 0)
                    {
                        SA[buf_cur[offset]] = buf_cur[offset + 1] - 1;
                        if (suffix)
                        {
                            buf_cur[offset] = 0;
                            buf_cur[offset + 1] = 0;
                        }
                    }
                }
                if (!suffix)
                {
                    parallel_for(unsigned int i = begin; i > end; i--)
                    {
                        unsigned int offset = 2 * (begin - i);
                        if (buf_cur[offset + 1] > 0)
                        {
                            buf_cur[offset] = 0;
                            buf_cur[offset + 1] = 0;
                            SA[i] = MAX_VALUE;
                        }
                    }
                }
            }
        }
                     );

        delete [] buf;
        buf = nullptr;

    }

};

#endif //MY_UTILITY_H
