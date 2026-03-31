/******************************************************************************
 * builder.h
 *
 *
 ******************************************************************************
 * Ling Bo Han <hanlb (at) mail2.sysu.edu.cn>
 * Yi Wu <wu.yi.christian (at) gmail.com>
 * Ge Nong <issng (at) mail.sysu.edu.cn>
 ******************************************************************************
*/

#ifndef _BUILDER_H
#define _BUILDER_H

#include <string>
#include <fstream>
#include <cassert>
#include <deque>
#include <bitset>
#include <thread>

#include "common.h"
#include "vector.h"
#include "tuple.h"
#include "tuple_sorter.h"
#include "utility.h"
#include "sortLMS.h"
#include "sorter.h"
#include "sais.h"


/// \brief preprocess

template<typename alphabet_type, typename offset_type, typename compress_type, typename relative_offset_type>
class DSAComputation
{

private:

    typedef MyVector<alphabet_type>                 alphabet_vector_type;

    typedef MyVector<offset_type>                   offset_vector_type;

    typedef MyVector<relative_offset_type>          relative_offset_vector_type;

    typedef BitVectorWapper<uint64>                 bit_vector_type;


    typedef  Pair<alphabet_type, compress_type>     compressed_pair_type;

    /// The second parameter records the occurrences times of the first parameter
    typedef MyVector< compressed_pair_type >        compressed_pair_vector_type;

    typedef MyVector<uint16>                        blockID_vector_type;

    typedef MyVector<uint32>                        rPos_vector_type; /// The relative position of LMS char in each text block

    /// Pair< position, name >
    typedef Pair<offset_type, offset_type>          sstar_pair_type;

    typedef MyVector<sstar_pair_type>               sstar_sub_vector_type;

    /// < relative_position, name>
    typedef Pair<relative_offset_type, offset_type> s1_pair_type;

    typedef MyVector<s1_pair_type>                  s1_blocks_vector_type;




private:

    const alphabet_type ALPHA_MAX;

    const alphabet_type ALPHA_MIN;

    const offset_type OFFSET_MAX;

    const offset_type OFFSET_MIN;

private:

    struct BlockInfo
    {

    public:

        uint64 m_capacity; ///< capacity for multi-block

        uint64 m_beg_pos; ///< starting position, global value

        uint64 m_end_pos; ///< ending position, global value

        uint64 m_size; ///< number of characters in the block, non-multi-block may exceed the capacity

        uint64 m_lms_num; ///< number of S*-substrs in the block.

        uint16 m_id; ///< block id, 2Bytes for 1T input data and 4G RAM. {3G / (5B input + 5B SA) = 300M of each block, 1000GB/0.3GB = 3333.33 blocks}

    public:

        /// \brief ctor
        ///
        BlockInfo(const uint64 & _end_pos, const uint16 _id) : m_end_pos(_end_pos), m_id(_id)
        {

            m_size = 0;

            m_lms_num = 0;
        }

        /// \brief assign copy
        ///
        BlockInfo& operator=(const BlockInfo & _bi)
        {

            m_capacity = _bi.m_capacity;

            m_beg_pos = _bi.m_beg_pos;

            m_end_pos = _bi.m_end_pos;

            m_size = _bi.m_size;

            m_lms_num = _bi.m_lms_num;

            m_id = _bi.m_id;

            return *this;
        }


        /// \brief check if the block is empty
        ///
        bool is_empty() const
        {

            return m_size == 0;
        }

        /// \brief check if the block contains more than one S*-substr
        ///
        bool is_multi() const
        {

            return m_lms_num > 1;
        }

        /// \brief check if the block contains only single S*-substr
        ///
        bool is_single() const
        {

            return (m_lms_num == 1);
            //return (m_lms_num == 1 && m_beg_pos != 0);
        }

        /// \brief check if the block is the leftmost block
        ///
        bool is_left_block() const
        {

            return (m_lms_num == 1 && m_beg_pos == 0);
        }

        /// \brief check if the block can afford the S*-substr of a certain size
        ///
        /// \note if the block is empty, it can afford an S*-substr of any size
        bool try_fill(const uint64 _lms_size)
        {

            if (is_empty() || m_size + _lms_size - 1 <= m_capacity)
            {

                return true;
            }

            return false;
        }

        /// \brief fill the S*-substr into the block
        ///
        /// \note If empty, insert all the characters; otherwise, insert all but the last character
        void fill(const uint64 _lms_size)
        {

            if (is_empty())
            {

                m_size = _lms_size;
            }
            else
            {

                m_size = m_size + _lms_size - 1;
            }

            ++m_lms_num;

            return;
        }

        /// \brief close the block
        ///
        void close()
        {

            m_beg_pos = m_end_pos - m_size + 1;

            return;
        }

        /// output the block informations
        void outputBlockInfo()
        {

            std::cout << "------------------------------------------------------" << std::endl;
            std::cout << "The blockID = " << int(m_id) << std::endl;

            std::cout << "size = " << m_size << std::endl;
            std::cout << "the number of LMS = " << m_lms_num << std::endl;
            std::cout << "beg_pos = " << m_beg_pos << " , end_pos = " << m_end_pos << std::endl;
            std::cout << "------------------------------------------------------" << std::endl;

        }
    };

private:

    std::vector<BlockInfo *>                m_blocks_info;

    std::vector<uint32>                     m_block_id_of_samplings;

    std::string                             m_s; /// The input string name

    uint32                                  m_level;

    alphabet_type                           m_alpha;

    uint64                                  m_s_len;

    uint64                                  m_lms_size; /// the LMS character number

    uint64                                  m_block_size; /// the block size of computing in RAM

    std::string                             m_s1; ///< the reduced string name

    offset_vector_type *&                   sa1_reverse;

    Alpha_Block<offset_type> *              m_alpha_blocks_info;

    uint16                                  m_alpha_blocks_number; /// the alpha block number

    std::vector<alphabet_type>              m_alpha_blocks_id_of_samplings;

    //std::vector<offset_type>                m_alpha_blocks_sampling_vector;

    bool                                    is_ram_build;


    /// The compressed S-type BWT sequences
    std::vector<compressed_pair_vector_type * >         m_sub_s_cbwt_seqs;
    std::vector<bit_vector_type * >                     m_sub_s_bit_seqs;/// S*-substring = 1, or otherwise = 0;

    std::vector<compressed_pair_vector_type * >         m_sub_l_cbwt_seqs;
    std::vector<bit_vector_type * >                     m_sub_l_bit_seqs;/// L*-substring = 1, or otherwise = 0;

    std::vector<compressed_pair_vector_type * >         m_suf_s_cbwt_seqs;
    std::vector<bit_vector_type * >                     m_suf_s_bit_seqs;/// S*-substring = 1, or otherwise = 0;

    std::vector<compressed_pair_vector_type * >         m_suf_l_cbwt_seqs;
    std::vector<bit_vector_type * >                     m_suf_l_bit_seqs;/// L*-substring = 1, or otherwise = 0;

    /// for alphabet > 255
    std::vector<compressed_pair_vector_type * >         m_cLMS_seqs; /// the compressed LMS character sequence


    std::vector< MyVector<relative_offset_type> *>      m_LMS_rPos_seqs; /// the relative position of LMS character in each block
    std::vector<bit_vector_type * >                     m_LMS_rPos_repeat_bit_seqs; /// mark the difference of adjacent substrings

    /// for alphabet <= 255
    std::vector<std::vector<uint64> * >                 m_LMS_counter; /// the RAM counter


    std::vector<offset_type>                            m_lms_blocks_of_samplings;

    std::vector<offset_type>                            m_lms_blocks_index;

    /// generate m_block_info.size() - 1 blocks
    std::vector<s1_blocks_vector_type *>                s1_blocks;

    MyVector<uint16> *                                  s1_blockID_seq; /// get the LMS char block id by scanning sa1_reverse from small to large



    std::vector<MyVector<uint24> *> m_sstar_char_seqs;


    std::vector< uint32 * > alpha_bkt_size_in_each_block;


    std::vector< uint32 > alpha_bkt_noEmpty_id;


public:

    DSAComputation(std::string _s, const uint32 _level, MyVector<offset_type> *& _sa_reverse, alphabet_type _alpha);

    /// separately calculate L-type BWTs
    void compute_bwt();

    void compute_blocks_info(const uint64 _end_position);

    bool mergeSortedSStarGlobal();

    void mergeSortedSuffixGlobal();

    void run();

    void compute_alphabet_blocks_id_of_samplings();

    void compute_block_id_of_samplings();

    uint32 getBlockId(const offset_type & _pos); /// return type of uint8 is set to uint32 on 20191027

    uint16 getAlphaBlockId(const alphabet_type & _pos);

    void compute_lms_blocks_id_of_samplings();

    uint16 getLmsBlockId(const offset_type & _pos);


    void getAlphaBktID(uint32, uint32 * );

    /// get alphabet size in sorted LStar subsequence
    uint32 getLStarAlphaBktID(uint32, uint64 *);

    /// return is_ram_build value
    bool get_build_way();

    /// sort S* substring 4 Small alphabet (<=1B)
    offset_type sortSub_4S();

    /// sort S* substring 4 Large alphabet (>1B)
    offset_type sortSub_4L();

    /// sort suffixes 4 Small alphabet (<=1B)
    void sortSuf_4S();

    /// sort suffixes 4 Large alphabet (>1B)
    void sortSuf_4L();

    void compute_alphabet_blocks(std::vector< Alpha_Block<offset_type> > &, offset_type &, offset_type &, offset_type &, offset_type &, offset_type, bool);

    /// compute the alpha_block id by binary_search
    uint16 get_alpha_blockID(const alphabet_type & _beg, const alphabet_type & _key);

    uint16 get_induce_alpha_blockID(const alphabet_type & _end_block_id, const alphabet_type & _key);

    /// generate s1 of 1 bytes
    void generate_s1_1B(offset_type & name);

    /// generate s1 of 2 bytes
    void generate_s1_2B(offset_type & name);

    /// generate s1 of 3 bytes
    void generate_s1_3B(offset_type & name);

    /// generate s1 of 4 bytes
    void generate_s1_4B(offset_type & name);

    /// generate s1 of 5 bytes
    void generate_s1_5B(offset_type & name);

};

/// \brief
///
template<typename alphabet_type, typename offset_type, typename compress_type, typename relative_offset_type>
DSAComputation<alphabet_type, offset_type, compress_type, relative_offset_type>::DSAComputation(std::string _s, const uint32 _level, MyVector<offset_type> *& _sa_reverse, alphabet_type _alpha)
    : ALPHA_MAX(std::numeric_limits<alphabet_type>::max())
    , ALPHA_MIN(std::numeric_limits<alphabet_type>::min())
    , OFFSET_MAX(std::numeric_limits<offset_type>::max())
    , OFFSET_MIN(std::numeric_limits<offset_type>::min())
    , m_s(_s)
    , m_level(_level)
    , sa1_reverse(_sa_reverse)
    , m_alpha(_alpha)
    , m_alpha_blocks_number(0)
    , m_s1("s1_level_" + std::to_string(_level))
    , m_lms_size(0)
    , is_ram_build(false)
{

    m_s_len = UtilityFunctions::getFileLength(m_s) / sizeof(alphabet_type);

    std::cout << "The input size = " << m_s_len << std::endl;

    /// if the computation of sorting SA can be completed in RAM, directly build SA in RAM
    {
        offset_type sa_size = m_s_len * sizeof(offset_type);
        offset_type input_size = m_s_len * sizeof(alphabet_type);
        offset_type bucket_size = offset_type(m_alpha) * sizeof(offset_type);
        offset_type type_size = m_s_len / 8;

        if(MAX_MEM > (sa_size + input_size + bucket_size + type_size))
        {
            std::cout << "Building suffix array in RAM.\n";

            FILE * fp = fopen(m_s.c_str(),"rb");

            if(!fp)
            {
                Logger::output_error(__FILE__,__LINE__);
            }

            uint64 file_length = UtilityFunctions::getFileLength(m_s);


            /// alphabet_type buffer
            alphabet_type * original_buf;

            if(m_level == 0)
            {
                original_buf = new alphabet_type [file_length + 1];

                fread(original_buf,sizeof(alphabet_type),file_length,fp);

                fclose(fp);

                Logger::addIV( (file_length + 1) * sizeof(alphabet_type) );

                original_buf[file_length] = 0;/// append the sentinel

                offset_type *sa = new offset_type[file_length + 1];

                std::cout << "start to compute the SA.\n";

                SAComputation<alphabet_type, offset_type>(original_buf, file_length + 1, _alpha, sa);

                std::cout << "The process of computing SA is over.\n";

                if(is_getBWT)
                {
                    std::cout << "start to comoute BWT." << std::endl;

                    alphabet_type * bwt_buf = new alphabet_type[file_length + 1];

                    for(uint64 i = 0; i < file_length + 1; i++)
                    {
                        if(sa[i]!=0)
                            bwt_buf[i] = original_buf[sa[i]-1];
                        else
                            bwt_buf[i] = original_buf[file_length];
                    }


                    std::string sa_name = _s + ".bwt";

                    FILE * fot = fopen(sa_name.c_str(),"wb");

                    fwrite(bwt_buf,sizeof(alphabet_type),file_length+1, fot);

                    fclose(fot);

                    fot = nullptr;


                    std::cout << "Write the BWT to disk.\n";


                    delete [] bwt_buf;
                    bwt_buf = nullptr;

                }


                if(is_getSA)
                {

                    std::string sa_name = _s + ".sa5";

                    FILE * fot = fopen(sa_name.c_str(),"wb");

                    fwrite(sa + 1,5,file_length, fot);

                    fclose(fot);

                    fot = nullptr;

                    std::cout << __LINE__ << std::endl;

                    std::cout << "Write the SA to disk.\n";

                }

                _sa_reverse = nullptr;


                delete [] original_buf;
                original_buf = nullptr;

                delete [] sa;
                sa = nullptr;


                fp = nullptr;


                Logger::addPDU( sizeof(offset_type) * (file_length) );

                Logger::addOV( sizeof(offset_type) * (file_length) );

                std::cout << "The building procedure is completed.\n";

                is_ram_build = true;

            }
            else
            {
                std::cout << "In DSAComputation(), error occur.\n";

                Logger::output_error(__FILE__,__LINE__);

                SAIS<alphabet_type, offset_type>(m_s,_sa_reverse);

            }

            return;
        }


    }


    std::cout << "Building SA on external memory.\n";


    /// computing SA on external memory
    if (m_alpha <= alphabet_type(255) )
    {

        //double div = sizeof(alphabet_type) + double((double)1 / (double)8) + sizeof(uint32) + sizeof(uint32);

        m_block_size = (MAX_MEM / (sizeof(alphabet_type) + sizeof(offset_type) + double((double)1 / (double)8) ) ) / 4; /// four thread; read into memory by the size of alphabet

        //m_block_size = 1500;/// for test sample

        std::cout << " m_block_size : " << m_block_size / K_1024 << std::endl;


        /// initialize the m_LMS_counter
        for(uint64 i(0); i <= m_alpha; ++i)
            m_LMS_counter.push_back(new std::vector<uint64>());


    }
    else
    {

        uint64 div = std::max(sizeof(alphabet_type) + double((double)1 / (double)8) + sizeof(offset_type) + sizeof(offset_type) + sizeof(offset_type),
                              double(sizeof(alphabet_type)) + sizeof(alphabet_type) + sizeof(offset_type) + sizeof(offset_type));

        // double div = std::max(sizeof(alphabet_type) + double((double)1 / (double)8) + sizeof(uint32) + sizeof(uint32) + sizeof(uint32),
        //                     double(sizeof(alphabet_type)) + sizeof(alphabet_type) + sizeof(uint32) + sizeof(uint32));

        m_block_size = MAX_MEM / div / 4; /// four thread

        std::cout << " m_block_size : " << m_block_size / K_1024 << std::endl;


        std::string alpha_blocks_file_name = "alpha_blocks_info_" + std::to_string(m_level - 1) + ".dat";

        m_alpha_blocks_number = UtilityFunctions::getFileLength(alpha_blocks_file_name) / sizeof(Alpha_Block<offset_type>); /// \note 20260225

        std::cout << "m_alpha_blocks_number = "<< m_alpha_blocks_number << std::endl;

        m_alpha_blocks_info = new Alpha_Block<offset_type>[m_alpha_blocks_number];

        FILE * fin = fopen(alpha_blocks_file_name.c_str(), "rb");

        if(!fin)
            Logger::output_error(__FILE__,__LINE__);

        fread(m_alpha_blocks_info, sizeof(Alpha_Block<offset_type>), m_alpha_blocks_number, fin);

        fclose(fin);

        std::cout << "remove alphabet block file " << alpha_blocks_file_name << "\n";

        /// remove the alpha_blocks_info
        std::remove(alpha_blocks_file_name.c_str());


        std::cout << "set the m_alpha_blocks_info[0].m_beg_alpha = 0. \n";
        m_alpha_blocks_info[0].m_beg_alpha = 0;


    }


    std::cout << "m_s_len = " << m_s_len << std::endl;

    alpha_bkt_size_in_each_block.clear();


}

/**
@function: compute the BWTs of L-type substrings only.
@para: _end_position = the beginning position of last block, it must be LMS character position
@para: _length = the length of the leftmost S*-substring
*/
template<typename alphabet_type, typename offset_type, typename compress_type, typename relative_offset_type>
void DSAComputation<alphabet_type, offset_type, compress_type, relative_offset_type>::compute_bwt()
{
    /// don't include the leftmost block
    uint32 loop = (m_blocks_info.size() - 1) / 4; /// \note four thread

    std::cout << "loop number = " << loop << std::endl;

    int is_remainder = ((m_blocks_info.size() - 1) % 4) ;/// four thread

    std::cout << "((m_blocks_info.size() - 1) % 4) = " << ((m_blocks_info.size() - 1) % 4) << std::endl;

    for(uint32 i = 0; i < loop; ++i)
    {

//        std::cout << "start to compute the BWTs of loop number = " << i << std::endl;

        if(sizeof(alphabet_type) == 1)
        {

            /// thread 1;
            bool is_sentinel_1 = (( (m_level == 0) && (m_LMS_rPos_seqs.size() == 0) ) ? true : false);

            for(uint64 j = 0; j <= m_alpha; ++j) m_LMS_counter[j]->push_back(0);

            m_sub_l_cbwt_seqs.push_back( new compressed_pair_vector_type() );
            m_sub_l_cbwt_seqs.back()->start_write();

            m_sub_l_bit_seqs.push_back( new bit_vector_type() );
            m_sub_l_bit_seqs.back()->start_write();


            m_sub_s_cbwt_seqs.push_back( new compressed_pair_vector_type() );
            m_sub_s_cbwt_seqs.back()->start_write();

            m_sub_s_bit_seqs.push_back( new bit_vector_type() );
            m_sub_s_bit_seqs.back()->start_write();

            m_LMS_rPos_seqs.push_back( new relative_offset_vector_type() );
            m_LMS_rPos_seqs.back()->start_write();

            m_LMS_rPos_repeat_bit_seqs.push_back( new bit_vector_type() );
            m_LMS_rPos_repeat_bit_seqs.back()->start_write();

            uint32 id_1 = i * 4;

            std::cout << "id = " << id_1 << std::endl;

            uint64 beg_pos_1 = m_blocks_info[id_1]->m_beg_pos;
            uint64 end_pos_1 = m_blocks_info[id_1]->m_end_pos;

            bool is_single_1 = m_blocks_info[id_1]->is_single();

            /// thread 2;

            bool is_sentinel_2 = false;

            for(uint64 j = 0; j <= m_alpha; ++j) m_LMS_counter[j]->push_back(0);

            m_sub_l_cbwt_seqs.push_back( new compressed_pair_vector_type() );
            m_sub_l_cbwt_seqs.back()->start_write();

            m_sub_l_bit_seqs.push_back( new bit_vector_type() );
            m_sub_l_bit_seqs.back()->start_write();

            m_sub_s_cbwt_seqs.push_back( new compressed_pair_vector_type() );
            m_sub_s_cbwt_seqs.back()->start_write();

            m_sub_s_bit_seqs.push_back( new bit_vector_type() );
            m_sub_s_bit_seqs.back()->start_write();

            m_LMS_rPos_seqs.push_back( new relative_offset_vector_type() );
            m_LMS_rPos_seqs.back()->start_write();

            m_LMS_rPos_repeat_bit_seqs.push_back( new bit_vector_type() );
            m_LMS_rPos_repeat_bit_seqs.back()->start_write();

            uint32 id_2 = i * 4 + 1;

            std::cout << "id = " << id_2 << std::endl;

            uint64 beg_pos_2 = m_blocks_info[id_2]->m_beg_pos;
            uint64 end_pos_2 = m_blocks_info[id_2]->m_end_pos;

            bool is_single_2 = m_blocks_info[id_2]->is_single();


            /// thread 3;

            bool is_sentinel_3 = false;

            for(uint64 j = 0; j <= m_alpha; ++j) m_LMS_counter[j]->push_back(0);

            m_sub_l_cbwt_seqs.push_back( new compressed_pair_vector_type() );
            m_sub_l_cbwt_seqs.back()->start_write();

            m_sub_l_bit_seqs.push_back( new bit_vector_type() );
            m_sub_l_bit_seqs.back()->start_write();

            m_sub_s_cbwt_seqs.push_back( new compressed_pair_vector_type() );
            m_sub_s_cbwt_seqs.back()->start_write();

            m_sub_s_bit_seqs.push_back( new bit_vector_type() );
            m_sub_s_bit_seqs.back()->start_write();

            m_LMS_rPos_seqs.push_back( new relative_offset_vector_type() );
            m_LMS_rPos_seqs.back()->start_write();

            m_LMS_rPos_repeat_bit_seqs.push_back( new bit_vector_type() );
            m_LMS_rPos_repeat_bit_seqs.back()->start_write();

            uint32 id_3 = i * 4 + 2;

            std::cout << "id = " << id_3 << std::endl;

            uint64 beg_pos_3 = m_blocks_info[id_3]->m_beg_pos;
            uint64 end_pos_3 = m_blocks_info[id_3]->m_end_pos;

            bool is_single_3 = m_blocks_info[id_3]->is_single();


            /// thread 4;

            bool is_sentinel_4 = false;

            for(uint64 j = 0; j <= m_alpha; ++j) m_LMS_counter[j]->push_back(0);

            m_sub_l_cbwt_seqs.push_back( new compressed_pair_vector_type() );
            m_sub_l_cbwt_seqs.back()->start_write();

            m_sub_l_bit_seqs.push_back( new bit_vector_type() );
            m_sub_l_bit_seqs.back()->start_write();

            m_sub_s_cbwt_seqs.push_back( new compressed_pair_vector_type() );
            m_sub_s_cbwt_seqs.back()->start_write();

            m_sub_s_bit_seqs.push_back( new bit_vector_type() );
            m_sub_s_bit_seqs.back()->start_write();

            m_LMS_rPos_seqs.push_back( new relative_offset_vector_type() );
            m_LMS_rPos_seqs.back()->start_write();

            m_LMS_rPos_repeat_bit_seqs.push_back( new bit_vector_type() );
            m_LMS_rPos_repeat_bit_seqs.back()->start_write();

            uint32 id_4 = i * 4 + 3;

            std::cout << "id = " << id_4 << std::endl;

            uint64 beg_pos_4 = m_blocks_info[id_4]->m_beg_pos;
            uint64 end_pos_4 = m_blocks_info[id_4]->m_end_pos;

            bool is_single_4 = m_blocks_info[id_4]->is_single();


            /// thread 1
            std::thread th1(UtilityFunctions::th_reduceSort_bwt_1B<alphabet_type,compress_type,relative_offset_type>
                            , m_s
                            , beg_pos_1
                            , end_pos_1
                            , m_alpha
                            , m_level
                            , id_1
                            , std::ref(m_LMS_counter)
                            , std::ref(m_sub_l_cbwt_seqs[id_1])
                            , std::ref(m_sub_l_bit_seqs[id_1])
                            , std::ref(m_sub_s_cbwt_seqs[id_1])
                            , std::ref(m_sub_s_bit_seqs[id_1])
                            , std::ref(m_LMS_rPos_seqs[id_1])
                            , std::ref(m_LMS_rPos_repeat_bit_seqs[id_1])
                            , is_single_1
                            , is_sentinel_1
                           );

            //std::cin.get();
            //th1.join(); //2024-03-04

            /// thread 2
            std::thread th2(UtilityFunctions::th_reduceSort_bwt_1B<alphabet_type,compress_type,relative_offset_type>
                            , m_s
                            , beg_pos_2
                            , end_pos_2
                            , m_alpha
                            , m_level
                            , id_2
                            , std::ref(m_LMS_counter)
                            , std::ref(m_sub_l_cbwt_seqs[id_2])
                            , std::ref(m_sub_l_bit_seqs[id_2])
                            , std::ref(m_sub_s_cbwt_seqs[id_2])
                            , std::ref(m_sub_s_bit_seqs[id_2])
                            , std::ref(m_LMS_rPos_seqs[id_2])
                            , std::ref(m_LMS_rPos_repeat_bit_seqs[id_2])
                            , is_single_2
                            , is_sentinel_2
                           );


            /// thread 3
            std::thread th3(UtilityFunctions::th_reduceSort_bwt_1B<alphabet_type,compress_type,relative_offset_type>
                            , m_s
                            , beg_pos_3
                            , end_pos_3
                            , m_alpha
                            , m_level
                            , id_3
                            , std::ref(m_LMS_counter)
                            , std::ref(m_sub_l_cbwt_seqs[id_3])
                            , std::ref(m_sub_l_bit_seqs[id_3])
                            , std::ref(m_sub_s_cbwt_seqs[id_3])
                            , std::ref(m_sub_s_bit_seqs[id_3])
                            , std::ref(m_LMS_rPos_seqs[id_3])
                            , std::ref(m_LMS_rPos_repeat_bit_seqs[id_3])
                            , is_single_3
                            , is_sentinel_3
                           );

            /// thread 4
            std::thread th4(UtilityFunctions::th_reduceSort_bwt_1B<alphabet_type,compress_type,relative_offset_type>
                            , m_s
                            , beg_pos_4
                            , end_pos_4
                            , m_alpha
                            , m_level
                            , id_4
                            , std::ref(m_LMS_counter)
                            , std::ref(m_sub_l_cbwt_seqs[id_4])
                            , std::ref(m_sub_l_bit_seqs[id_4])
                            , std::ref(m_sub_s_cbwt_seqs[id_4])
                            , std::ref(m_sub_s_bit_seqs[id_4])
                            , std::ref(m_LMS_rPos_seqs[id_4])
                            , std::ref(m_LMS_rPos_repeat_bit_seqs[id_4])
                            , is_single_4
                            , is_sentinel_4
                           );

            th1.join();
            th2.join();
            th3.join();
            th4.join();
        }
        else
        {


            /// thread 1;
            bool is_sentinel_1 = false;

            m_cLMS_seqs.push_back(new compressed_pair_vector_type());
            m_cLMS_seqs.back()->start_write();

            m_sub_l_cbwt_seqs.push_back( new compressed_pair_vector_type() );
            m_sub_l_cbwt_seqs.back()->start_write();

            m_sub_l_bit_seqs.push_back( new bit_vector_type() );
            m_sub_l_bit_seqs.back()->start_write();


            m_sub_s_cbwt_seqs.push_back( new compressed_pair_vector_type() );
            m_sub_s_cbwt_seqs.back()->start_write();

            m_sub_s_bit_seqs.push_back( new bit_vector_type() );
            m_sub_s_bit_seqs.back()->start_write();

            m_LMS_rPos_seqs.push_back( new relative_offset_vector_type() );
            m_LMS_rPos_seqs.back()->start_write();

            m_LMS_rPos_repeat_bit_seqs.push_back( new bit_vector_type() );
            m_LMS_rPos_repeat_bit_seqs.back()->start_write();

            uint32 id_1 = i * 4;

            std::cout << "id = " << id_1 << std::endl;

            uint64 beg_pos_1 = m_blocks_info[id_1]->m_beg_pos;
            uint64 end_pos_1 = m_blocks_info[id_1]->m_end_pos;

            bool is_single_1 = m_blocks_info[id_1]->is_single();



            /// thread 2;

            bool is_sentinel_2 = false;

            m_cLMS_seqs.push_back(new compressed_pair_vector_type());
            m_cLMS_seqs.back()->start_write();

            m_sub_l_cbwt_seqs.push_back( new compressed_pair_vector_type() );
            m_sub_l_cbwt_seqs.back()->start_write();

            m_sub_l_bit_seqs.push_back( new bit_vector_type() );
            m_sub_l_bit_seqs.back()->start_write();

            m_sub_s_cbwt_seqs.push_back( new compressed_pair_vector_type() );
            m_sub_s_cbwt_seqs.back()->start_write();

            m_sub_s_bit_seqs.push_back( new bit_vector_type() );
            m_sub_s_bit_seqs.back()->start_write();

            m_LMS_rPos_seqs.push_back( new relative_offset_vector_type() );
            m_LMS_rPos_seqs.back()->start_write();

            m_LMS_rPos_repeat_bit_seqs.push_back( new bit_vector_type() );
            m_LMS_rPos_repeat_bit_seqs.back()->start_write();

            uint32 id_2 = i * 4 + 1;

            std::cout << "id = " << id_2 << std::endl;

            uint64 beg_pos_2 = m_blocks_info[id_2]->m_beg_pos;
            uint64 end_pos_2 = m_blocks_info[id_2]->m_end_pos;

            bool is_single_2 = m_blocks_info[id_2]->is_single();


            /// thread 3;
            bool is_sentinel_3 = false;

            m_cLMS_seqs.push_back(new compressed_pair_vector_type());
            m_cLMS_seqs.back()->start_write();

            m_sub_l_cbwt_seqs.push_back( new compressed_pair_vector_type() );
            m_sub_l_cbwt_seqs.back()->start_write();

            m_sub_l_bit_seqs.push_back( new bit_vector_type() );
            m_sub_l_bit_seqs.back()->start_write();

            m_sub_s_cbwt_seqs.push_back( new compressed_pair_vector_type() );
            m_sub_s_cbwt_seqs.back()->start_write();

            m_sub_s_bit_seqs.push_back( new bit_vector_type() );
            m_sub_s_bit_seqs.back()->start_write();

            m_LMS_rPos_seqs.push_back( new relative_offset_vector_type() );
            m_LMS_rPos_seqs.back()->start_write();

            m_LMS_rPos_repeat_bit_seqs.push_back( new bit_vector_type() );
            m_LMS_rPos_repeat_bit_seqs.back()->start_write();

            uint32 id_3 = i * 4 + 2;

            std::cout << "id_3 = " << id_3 << std::endl;

            uint64 beg_pos_3 = m_blocks_info[id_3]->m_beg_pos;
            uint64 end_pos_3 = m_blocks_info[id_3]->m_end_pos;

            bool is_single_3 = m_blocks_info[id_3]->is_single();


            /// thread 4
            bool is_sentinel_4 = false;

            m_cLMS_seqs.push_back(new compressed_pair_vector_type());
            m_cLMS_seqs.back()->start_write();

            m_sub_l_cbwt_seqs.push_back( new compressed_pair_vector_type() );
            m_sub_l_cbwt_seqs.back()->start_write();

            m_sub_l_bit_seqs.push_back( new bit_vector_type() );
            m_sub_l_bit_seqs.back()->start_write();

            m_sub_s_cbwt_seqs.push_back( new compressed_pair_vector_type() );
            m_sub_s_cbwt_seqs.back()->start_write();

            m_sub_s_bit_seqs.push_back( new bit_vector_type() );
            m_sub_s_bit_seqs.back()->start_write();

            m_LMS_rPos_seqs.push_back( new relative_offset_vector_type() );
            m_LMS_rPos_seqs.back()->start_write();

            m_LMS_rPos_repeat_bit_seqs.push_back( new bit_vector_type() );
            m_LMS_rPos_repeat_bit_seqs.back()->start_write();

            uint32 id_4 = i * 4 + 3;

            std::cout << "id_4 = " << id_4 << std::endl;

            uint64 beg_pos_4 = m_blocks_info[id_4]->m_beg_pos;
            uint64 end_pos_4 = m_blocks_info[id_4]->m_end_pos;

            bool is_single_4 = m_blocks_info[id_4]->is_single();

            /// thread 1
            std::thread th1(UtilityFunctions::th_reduceSort_bwt_4B<alphabet_type,compress_type,relative_offset_type>
                            , m_s
                            , beg_pos_1
                            , end_pos_1
                            , m_alpha
                            , m_level
                            , id_1
                            , std::ref(m_cLMS_seqs[id_1])
                            , std::ref(m_sub_l_cbwt_seqs[id_1])
                            , std::ref(m_sub_l_bit_seqs[id_1])
                            , std::ref(m_sub_s_cbwt_seqs[id_1])
                            , std::ref(m_sub_s_bit_seqs[id_1])
                            , std::ref(m_LMS_rPos_seqs[id_1])
                            , std::ref(m_LMS_rPos_repeat_bit_seqs[id_1])
                            , is_single_1
                            , is_sentinel_1
                           );

            //std::cin.get();
            //th1.join(); //2024-03-04

            /// thread 2
            std::thread th2(UtilityFunctions::th_reduceSort_bwt_4B<alphabet_type,compress_type,relative_offset_type>
                            , m_s
                            , beg_pos_2
                            , end_pos_2
                            , m_alpha
                            , m_level
                            , id_2
                            , std::ref(m_cLMS_seqs[id_2])
                            , std::ref(m_sub_l_cbwt_seqs[id_2])
                            , std::ref(m_sub_l_bit_seqs[id_2])
                            , std::ref(m_sub_s_cbwt_seqs[id_2])
                            , std::ref(m_sub_s_bit_seqs[id_2])
                            , std::ref(m_LMS_rPos_seqs[id_2])
                            , std::ref(m_LMS_rPos_repeat_bit_seqs[id_2])
                            , is_single_2
                            , is_sentinel_2
                           );
            /// thread3
            std::thread th3(UtilityFunctions::th_reduceSort_bwt_4B<alphabet_type,compress_type,relative_offset_type>
                            , m_s
                            , beg_pos_3
                            , end_pos_3
                            , m_alpha
                            , m_level
                            , id_3
                            , std::ref(m_cLMS_seqs[id_3])
                            , std::ref(m_sub_l_cbwt_seqs[id_3])
                            , std::ref(m_sub_l_bit_seqs[id_3])
                            , std::ref(m_sub_s_cbwt_seqs[id_3])
                            , std::ref(m_sub_s_bit_seqs[id_3])
                            , std::ref(m_LMS_rPos_seqs[id_3])
                            , std::ref(m_LMS_rPos_repeat_bit_seqs[id_3])
                            , is_single_3
                            , is_sentinel_3
                           );

            /// thread4
            std::thread th4(UtilityFunctions::th_reduceSort_bwt_4B<alphabet_type,compress_type,relative_offset_type>
                            , m_s
                            , beg_pos_4
                            , end_pos_4
                            , m_alpha
                            , m_level
                            , id_4
                            , std::ref(m_cLMS_seqs[id_4])
                            , std::ref(m_sub_l_cbwt_seqs[id_4])
                            , std::ref(m_sub_l_bit_seqs[id_4])
                            , std::ref(m_sub_s_cbwt_seqs[id_4])
                            , std::ref(m_sub_s_bit_seqs[id_4])
                            , std::ref(m_LMS_rPos_seqs[id_4])
                            , std::ref(m_LMS_rPos_repeat_bit_seqs[id_4])
                            , is_single_4
                            , is_sentinel_4
                           );

            th1.join();
            th2.join();
            th3.join();
            th4.join();

//            std::cout << "th1 and th2 is closed.\n";

        }

    }

    if(is_remainder == 1)
    {
        std::cout << "************************************ is_remainder = 1 ************************************\n";

        if(sizeof(alphabet_type) == 1)
        {
            bool is_sentinel_1 = (( (m_level == 0) && (m_LMS_rPos_seqs.size() == 0) ) ? true : false);

            for(uint64 j = 0; j <= m_alpha; ++j) m_LMS_counter[j]->push_back(0);

            m_sub_l_cbwt_seqs.push_back( new compressed_pair_vector_type() );
            m_sub_l_cbwt_seqs.back()->start_write();

            m_sub_l_bit_seqs.push_back( new bit_vector_type() );
            m_sub_l_bit_seqs.back()->start_write();


            m_sub_s_cbwt_seqs.push_back( new compressed_pair_vector_type() );
            m_sub_s_cbwt_seqs.back()->start_write();

            m_sub_s_bit_seqs.push_back( new bit_vector_type() );
            m_sub_s_bit_seqs.back()->start_write();

            m_LMS_rPos_seqs.push_back( new relative_offset_vector_type() );
            m_LMS_rPos_seqs.back()->start_write();

            m_LMS_rPos_repeat_bit_seqs.push_back( new bit_vector_type() );
            m_LMS_rPos_repeat_bit_seqs.back()->start_write();

            uint32 id = m_blocks_info.size()-2;

            std::cout << "id = " << id << std::endl;

            uint64 beg_pos_1 = m_blocks_info[id]->m_beg_pos;
            uint64 end_pos_1 = m_blocks_info[id]->m_end_pos;

            bool is_single_1 = m_blocks_info[id]->is_single();

            UtilityFunctions::th_reduceSort_bwt_1B<alphabet_type,compress_type,relative_offset_type>(
                m_s
                , beg_pos_1
                , end_pos_1
                , m_alpha
                , m_level
                , id
                , std::ref(m_LMS_counter)
                , std::ref(m_sub_l_cbwt_seqs[id])
                , std::ref(m_sub_l_bit_seqs[id])
                , std::ref(m_sub_s_cbwt_seqs[id])
                , std::ref(m_sub_s_bit_seqs[id])
                , std::ref(m_LMS_rPos_seqs[id])
                , std::ref(m_LMS_rPos_repeat_bit_seqs[id])
                , is_single_1
                , is_sentinel_1
            );

        }
        else
        {
            bool is_sentinel_1 = false;

            m_cLMS_seqs.push_back(new compressed_pair_vector_type());
            m_cLMS_seqs.back()->start_write();

            m_sub_l_cbwt_seqs.push_back( new compressed_pair_vector_type() );
            m_sub_l_cbwt_seqs.back()->start_write();

            m_sub_l_bit_seqs.push_back( new bit_vector_type() );
            m_sub_l_bit_seqs.back()->start_write();


            m_sub_s_cbwt_seqs.push_back( new compressed_pair_vector_type() );
            m_sub_s_cbwt_seqs.back()->start_write();

            m_sub_s_bit_seqs.push_back( new bit_vector_type() );
            m_sub_s_bit_seqs.back()->start_write();

            m_LMS_rPos_seqs.push_back( new relative_offset_vector_type() );
            m_LMS_rPos_seqs.back()->start_write();

            m_LMS_rPos_repeat_bit_seqs.push_back( new bit_vector_type() );
            m_LMS_rPos_repeat_bit_seqs.back()->start_write();

            uint32 id = m_blocks_info.size()-2;

            std::cout << "id = " << id << std::endl;

            uint64 beg_pos_1 = m_blocks_info[id]->m_beg_pos;
            uint64 end_pos_1 = m_blocks_info[id]->m_end_pos;

            bool is_single_1 = m_blocks_info[id]->is_single();

            UtilityFunctions::th_reduceSort_bwt_4B<alphabet_type,compress_type,relative_offset_type>(
                m_s
                , beg_pos_1
                , end_pos_1
                , m_alpha
                , m_level
                , id
                , std::ref(m_cLMS_seqs[id])
                , std::ref(m_sub_l_cbwt_seqs[id])
                , std::ref(m_sub_l_bit_seqs[id])
                , std::ref(m_sub_s_cbwt_seqs[id])
                , std::ref(m_sub_s_bit_seqs[id])
                , std::ref(m_LMS_rPos_seqs[id])
                , std::ref(m_LMS_rPos_repeat_bit_seqs[id])
                , is_single_1
                , is_sentinel_1
            );
        }

    }

    if(is_remainder == 2)
    {
        std::cout << "************************************ is_remainder = 2 ************************************\n";

        if(sizeof(alphabet_type) == 1)
        {
            /// thread 1

            bool is_sentinel_1 = (( (m_level == 0) && (m_LMS_rPos_seqs.size() == 0) ) ? true : false);

            for(uint64 j = 0; j <= m_alpha; ++j) m_LMS_counter[j]->push_back(0);

            m_sub_l_cbwt_seqs.push_back( new compressed_pair_vector_type() );
            m_sub_l_cbwt_seqs.back()->start_write();

            m_sub_l_bit_seqs.push_back( new bit_vector_type() );
            m_sub_l_bit_seqs.back()->start_write();


            m_sub_s_cbwt_seqs.push_back( new compressed_pair_vector_type() );
            m_sub_s_cbwt_seqs.back()->start_write();

            m_sub_s_bit_seqs.push_back( new bit_vector_type() );
            m_sub_s_bit_seqs.back()->start_write();

            m_LMS_rPos_seqs.push_back( new relative_offset_vector_type() );
            m_LMS_rPos_seqs.back()->start_write();

            m_LMS_rPos_repeat_bit_seqs.push_back( new bit_vector_type() );
            m_LMS_rPos_repeat_bit_seqs.back()->start_write();

            uint32 id_1 = m_blocks_info.size()-3;

            std::cout << "id_1 = " << id_1 << std::endl;

            uint64 beg_pos_1 = m_blocks_info[id_1]->m_beg_pos;
            uint64 end_pos_1 = m_blocks_info[id_1]->m_end_pos;

            bool is_single_1 = m_blocks_info[id_1]->is_single();


            /// thread 2
            bool is_sentinel_2 =  false;

            for(uint64 j = 0; j <= m_alpha; ++j) m_LMS_counter[j]->push_back(0);

            m_sub_l_cbwt_seqs.push_back( new compressed_pair_vector_type() );
            m_sub_l_cbwt_seqs.back()->start_write();

            m_sub_l_bit_seqs.push_back( new bit_vector_type() );
            m_sub_l_bit_seqs.back()->start_write();


            m_sub_s_cbwt_seqs.push_back( new compressed_pair_vector_type() );
            m_sub_s_cbwt_seqs.back()->start_write();

            m_sub_s_bit_seqs.push_back( new bit_vector_type() );
            m_sub_s_bit_seqs.back()->start_write();

            m_LMS_rPos_seqs.push_back( new relative_offset_vector_type() );
            m_LMS_rPos_seqs.back()->start_write();

            m_LMS_rPos_repeat_bit_seqs.push_back( new bit_vector_type() );
            m_LMS_rPos_repeat_bit_seqs.back()->start_write();

            uint32 id_2 = id_1 + 1;

            std::cout << "id_2 = " << id_2 << std::endl;

            uint64 beg_pos_2 = m_blocks_info[id_2]->m_beg_pos;
            uint64 end_pos_2 = m_blocks_info[id_2]->m_end_pos;

            bool is_single_2= m_blocks_info[id_2]->is_single();

            /// thread 1
            std::thread th1(UtilityFunctions::th_reduceSort_bwt_1B<alphabet_type,compress_type,relative_offset_type>
                            , m_s
                            , beg_pos_1
                            , end_pos_1
                            , m_alpha
                            , m_level
                            , id_1
                            , std::ref(m_LMS_counter)
                            , std::ref(m_sub_l_cbwt_seqs[id_1])
                            , std::ref(m_sub_l_bit_seqs[id_1])
                            , std::ref(m_sub_s_cbwt_seqs[id_1])
                            , std::ref(m_sub_s_bit_seqs[id_1])
                            , std::ref(m_LMS_rPos_seqs[id_1])
                            , std::ref(m_LMS_rPos_repeat_bit_seqs[id_1])
                            , is_single_1
                            , is_sentinel_1
                           );


            /// thread 2
            std::thread th2(UtilityFunctions::th_reduceSort_bwt_1B<alphabet_type,compress_type,relative_offset_type>
                            , m_s
                            , beg_pos_2
                            , end_pos_2
                            , m_alpha
                            , m_level
                            , id_2
                            , std::ref(m_LMS_counter)
                            , std::ref(m_sub_l_cbwt_seqs[id_2])
                            , std::ref(m_sub_l_bit_seqs[id_2])
                            , std::ref(m_sub_s_cbwt_seqs[id_2])
                            , std::ref(m_sub_s_bit_seqs[id_2])
                            , std::ref(m_LMS_rPos_seqs[id_2])
                            , std::ref(m_LMS_rPos_repeat_bit_seqs[id_2])
                            , is_single_2
                            , is_sentinel_2
                           );

            th1.join();
            th2.join();

        }
        else
        {
            /// thread 1

            bool is_sentinel_1 = false;

            m_cLMS_seqs.push_back(new compressed_pair_vector_type());
            m_cLMS_seqs.back()->start_write();

            m_sub_l_cbwt_seqs.push_back( new compressed_pair_vector_type() );
            m_sub_l_cbwt_seqs.back()->start_write();

            m_sub_l_bit_seqs.push_back( new bit_vector_type() );
            m_sub_l_bit_seqs.back()->start_write();


            m_sub_s_cbwt_seqs.push_back( new compressed_pair_vector_type() );
            m_sub_s_cbwt_seqs.back()->start_write();

            m_sub_s_bit_seqs.push_back( new bit_vector_type() );
            m_sub_s_bit_seqs.back()->start_write();

            m_LMS_rPos_seqs.push_back( new relative_offset_vector_type() );
            m_LMS_rPos_seqs.back()->start_write();

            m_LMS_rPos_repeat_bit_seqs.push_back( new bit_vector_type() );
            m_LMS_rPos_repeat_bit_seqs.back()->start_write();

            uint32 id_1 = m_blocks_info.size() - 3;

            uint64 beg_pos_1 = m_blocks_info[id_1]->m_beg_pos;
            uint64 end_pos_1 = m_blocks_info[id_1]->m_end_pos;

            bool is_single_1 = m_blocks_info[id_1]->is_single();


            /// thread 2

            bool is_sentinel_2 = false;

            m_cLMS_seqs.push_back(new compressed_pair_vector_type());
            m_cLMS_seqs.back()->start_write();

            m_sub_l_cbwt_seqs.push_back( new compressed_pair_vector_type() );
            m_sub_l_cbwt_seqs.back()->start_write();

            m_sub_l_bit_seqs.push_back( new bit_vector_type() );
            m_sub_l_bit_seqs.back()->start_write();


            m_sub_s_cbwt_seqs.push_back( new compressed_pair_vector_type() );
            m_sub_s_cbwt_seqs.back()->start_write();

            m_sub_s_bit_seqs.push_back( new bit_vector_type() );
            m_sub_s_bit_seqs.back()->start_write();

            m_LMS_rPos_seqs.push_back( new relative_offset_vector_type() );
            m_LMS_rPos_seqs.back()->start_write();

            m_LMS_rPos_repeat_bit_seqs.push_back( new bit_vector_type() );
            m_LMS_rPos_repeat_bit_seqs.back()->start_write();

            uint32 id_2 = id_1 + 1;

            uint64 beg_pos_2 = m_blocks_info[id_2]->m_beg_pos;
            uint64 end_pos_2 = m_blocks_info[id_2]->m_end_pos;

            bool is_single_2 = m_blocks_info[id_2]->is_single();

            /// thread 1
            std::thread th1(UtilityFunctions::th_reduceSort_bwt_4B<alphabet_type,compress_type,relative_offset_type>
                            , m_s
                            , beg_pos_1
                            , end_pos_1
                            , m_alpha
                            , m_level
                            , id_1
                            , std::ref(m_cLMS_seqs[id_1])
                            , std::ref(m_sub_l_cbwt_seqs[id_1])
                            , std::ref(m_sub_l_bit_seqs[id_1])
                            , std::ref(m_sub_s_cbwt_seqs[id_1])
                            , std::ref(m_sub_s_bit_seqs[id_1])
                            , std::ref(m_LMS_rPos_seqs[id_1])
                            , std::ref(m_LMS_rPos_repeat_bit_seqs[id_1])
                            , is_single_1
                            , is_sentinel_1
                           );


            /// thread 2
            std::thread th2(UtilityFunctions::th_reduceSort_bwt_4B<alphabet_type,compress_type,relative_offset_type>
                            , m_s
                            , beg_pos_2
                            , end_pos_2
                            , m_alpha
                            , m_level
                            , id_2
                            , std::ref(m_cLMS_seqs[id_2])
                            , std::ref(m_sub_l_cbwt_seqs[id_2])
                            , std::ref(m_sub_l_bit_seqs[id_2])
                            , std::ref(m_sub_s_cbwt_seqs[id_2])
                            , std::ref(m_sub_s_bit_seqs[id_2])
                            , std::ref(m_LMS_rPos_seqs[id_2])
                            , std::ref(m_LMS_rPos_repeat_bit_seqs[id_2])
                            , is_single_2
                            , is_sentinel_2
                           );

            th1.join();
            th2.join();
        }

    }

    if(is_remainder == 3)
    {
        std::cout << "************************************ is_remainder = 3 ************************************\n";

        if(sizeof(alphabet_type) == 1)
        {
            /// thread 1

            bool is_sentinel_1 = (( (m_level == 0) && (m_LMS_rPos_seqs.size() == 0) ) ? true : false);

            for(uint64 j = 0; j <= m_alpha; ++j) m_LMS_counter[j]->push_back(0);

            m_sub_l_cbwt_seqs.push_back( new compressed_pair_vector_type() );
            m_sub_l_cbwt_seqs.back()->start_write();

            m_sub_l_bit_seqs.push_back( new bit_vector_type() );
            m_sub_l_bit_seqs.back()->start_write();


            m_sub_s_cbwt_seqs.push_back( new compressed_pair_vector_type() );
            m_sub_s_cbwt_seqs.back()->start_write();

            m_sub_s_bit_seqs.push_back( new bit_vector_type() );
            m_sub_s_bit_seqs.back()->start_write();

            m_LMS_rPos_seqs.push_back( new relative_offset_vector_type() );
            m_LMS_rPos_seqs.back()->start_write();

            m_LMS_rPos_repeat_bit_seqs.push_back( new bit_vector_type() );
            m_LMS_rPos_repeat_bit_seqs.back()->start_write();

            uint32 id_1 = m_blocks_info.size()-4;

            std::cout << "id_1 = " << id_1 << std::endl;

            uint64 beg_pos_1 = m_blocks_info[id_1]->m_beg_pos;
            uint64 end_pos_1 = m_blocks_info[id_1]->m_end_pos;

            bool is_single_1 = m_blocks_info[id_1]->is_single();


            /// thread 2
            bool is_sentinel_2 =  false;

            for(uint64 j = 0; j <= m_alpha; ++j) m_LMS_counter[j]->push_back(0);

            m_sub_l_cbwt_seqs.push_back( new compressed_pair_vector_type() );
            m_sub_l_cbwt_seqs.back()->start_write();

            m_sub_l_bit_seqs.push_back( new bit_vector_type() );
            m_sub_l_bit_seqs.back()->start_write();


            m_sub_s_cbwt_seqs.push_back( new compressed_pair_vector_type() );
            m_sub_s_cbwt_seqs.back()->start_write();

            m_sub_s_bit_seqs.push_back( new bit_vector_type() );
            m_sub_s_bit_seqs.back()->start_write();

            m_LMS_rPos_seqs.push_back( new relative_offset_vector_type() );
            m_LMS_rPos_seqs.back()->start_write();

            m_LMS_rPos_repeat_bit_seqs.push_back( new bit_vector_type() );
            m_LMS_rPos_repeat_bit_seqs.back()->start_write();

            uint32 id_2 = id_1 + 1;

            std::cout << "id_2 = " << id_2 << std::endl;

            uint64 beg_pos_2 = m_blocks_info[id_2]->m_beg_pos;
            uint64 end_pos_2 = m_blocks_info[id_2]->m_end_pos;

            bool is_single_2= m_blocks_info[id_2]->is_single();


            /// thread 3
            bool is_sentinel_3 =  false;

            for(uint64 j = 0; j <= m_alpha; ++j) m_LMS_counter[j]->push_back(0);

            m_sub_l_cbwt_seqs.push_back( new compressed_pair_vector_type() );
            m_sub_l_cbwt_seqs.back()->start_write();

            m_sub_l_bit_seqs.push_back( new bit_vector_type() );
            m_sub_l_bit_seqs.back()->start_write();


            m_sub_s_cbwt_seqs.push_back( new compressed_pair_vector_type() );
            m_sub_s_cbwt_seqs.back()->start_write();

            m_sub_s_bit_seqs.push_back( new bit_vector_type() );
            m_sub_s_bit_seqs.back()->start_write();

            m_LMS_rPos_seqs.push_back( new relative_offset_vector_type() );
            m_LMS_rPos_seqs.back()->start_write();

            m_LMS_rPos_repeat_bit_seqs.push_back( new bit_vector_type() );
            m_LMS_rPos_repeat_bit_seqs.back()->start_write();

            uint32 id_3 = id_1 + 2;

            std::cout << "id_3 = " << id_3 << std::endl;

            uint64 beg_pos_3 = m_blocks_info[id_3]->m_beg_pos;
            uint64 end_pos_3 = m_blocks_info[id_3]->m_end_pos;

            bool is_single_3= m_blocks_info[id_3]->is_single();

            /// thread 1
            std::thread th1(UtilityFunctions::th_reduceSort_bwt_1B<alphabet_type,compress_type,relative_offset_type>
                            , m_s
                            , beg_pos_1
                            , end_pos_1
                            , m_alpha
                            , m_level
                            , id_1
                            , std::ref(m_LMS_counter)
                            , std::ref(m_sub_l_cbwt_seqs[id_1])
                            , std::ref(m_sub_l_bit_seqs[id_1])
                            , std::ref(m_sub_s_cbwt_seqs[id_1])
                            , std::ref(m_sub_s_bit_seqs[id_1])
                            , std::ref(m_LMS_rPos_seqs[id_1])
                            , std::ref(m_LMS_rPos_repeat_bit_seqs[id_1])
                            , is_single_1
                            , is_sentinel_1
                           );


            /// thread 2
            std::thread th2(UtilityFunctions::th_reduceSort_bwt_1B<alphabet_type,compress_type,relative_offset_type>
                            , m_s
                            , beg_pos_2
                            , end_pos_2
                            , m_alpha
                            , m_level
                            , id_2
                            , std::ref(m_LMS_counter)
                            , std::ref(m_sub_l_cbwt_seqs[id_2])
                            , std::ref(m_sub_l_bit_seqs[id_2])
                            , std::ref(m_sub_s_cbwt_seqs[id_2])
                            , std::ref(m_sub_s_bit_seqs[id_2])
                            , std::ref(m_LMS_rPos_seqs[id_2])
                            , std::ref(m_LMS_rPos_repeat_bit_seqs[id_2])
                            , is_single_2
                            , is_sentinel_2
                           );

            /// thread 3
            std::thread th3(UtilityFunctions::th_reduceSort_bwt_1B<alphabet_type,compress_type,relative_offset_type>
                            , m_s
                            , beg_pos_3
                            , end_pos_3
                            , m_alpha
                            , m_level
                            , id_3
                            , std::ref(m_LMS_counter)
                            , std::ref(m_sub_l_cbwt_seqs[id_3])
                            , std::ref(m_sub_l_bit_seqs[id_3])
                            , std::ref(m_sub_s_cbwt_seqs[id_3])
                            , std::ref(m_sub_s_bit_seqs[id_3])
                            , std::ref(m_LMS_rPos_seqs[id_3])
                            , std::ref(m_LMS_rPos_repeat_bit_seqs[id_3])
                            , is_single_3
                            , is_sentinel_3
                           );

            th1.join();
            th2.join();
            th3.join();

        }
        else
        {
            /// thread 1

            bool is_sentinel_1 = false;

            m_cLMS_seqs.push_back(new compressed_pair_vector_type());
            m_cLMS_seqs.back()->start_write();

            m_sub_l_cbwt_seqs.push_back( new compressed_pair_vector_type() );
            m_sub_l_cbwt_seqs.back()->start_write();

            m_sub_l_bit_seqs.push_back( new bit_vector_type() );
            m_sub_l_bit_seqs.back()->start_write();


            m_sub_s_cbwt_seqs.push_back( new compressed_pair_vector_type() );
            m_sub_s_cbwt_seqs.back()->start_write();

            m_sub_s_bit_seqs.push_back( new bit_vector_type() );
            m_sub_s_bit_seqs.back()->start_write();

            m_LMS_rPos_seqs.push_back( new relative_offset_vector_type() );
            m_LMS_rPos_seqs.back()->start_write();

            m_LMS_rPos_repeat_bit_seqs.push_back( new bit_vector_type() );
            m_LMS_rPos_repeat_bit_seqs.back()->start_write();

            uint32 id_1 = m_blocks_info.size() - 4;

            uint64 beg_pos_1 = m_blocks_info[id_1]->m_beg_pos;
            uint64 end_pos_1 = m_blocks_info[id_1]->m_end_pos;

            bool is_single_1 = m_blocks_info[id_1]->is_single();


            /// thread 2
            bool is_sentinel_2 = false;

            m_cLMS_seqs.push_back(new compressed_pair_vector_type());
            m_cLMS_seqs.back()->start_write();

            m_sub_l_cbwt_seqs.push_back( new compressed_pair_vector_type() );
            m_sub_l_cbwt_seqs.back()->start_write();

            m_sub_l_bit_seqs.push_back( new bit_vector_type() );
            m_sub_l_bit_seqs.back()->start_write();


            m_sub_s_cbwt_seqs.push_back( new compressed_pair_vector_type() );
            m_sub_s_cbwt_seqs.back()->start_write();

            m_sub_s_bit_seqs.push_back( new bit_vector_type() );
            m_sub_s_bit_seqs.back()->start_write();

            m_LMS_rPos_seqs.push_back( new relative_offset_vector_type() );
            m_LMS_rPos_seqs.back()->start_write();

            m_LMS_rPos_repeat_bit_seqs.push_back( new bit_vector_type() );
            m_LMS_rPos_repeat_bit_seqs.back()->start_write();

            uint32 id_2 = id_1 + 1;

            uint64 beg_pos_2 = m_blocks_info[id_2]->m_beg_pos;
            uint64 end_pos_2 = m_blocks_info[id_2]->m_end_pos;

            bool is_single_2 = m_blocks_info[id_2]->is_single();

            /// thread 3
            bool is_sentinel_3 = false;

            m_cLMS_seqs.push_back(new compressed_pair_vector_type());
            m_cLMS_seqs.back()->start_write();

            m_sub_l_cbwt_seqs.push_back( new compressed_pair_vector_type() );
            m_sub_l_cbwt_seqs.back()->start_write();

            m_sub_l_bit_seqs.push_back( new bit_vector_type() );
            m_sub_l_bit_seqs.back()->start_write();


            m_sub_s_cbwt_seqs.push_back( new compressed_pair_vector_type() );
            m_sub_s_cbwt_seqs.back()->start_write();

            m_sub_s_bit_seqs.push_back( new bit_vector_type() );
            m_sub_s_bit_seqs.back()->start_write();

            m_LMS_rPos_seqs.push_back( new relative_offset_vector_type() );
            m_LMS_rPos_seqs.back()->start_write();

            m_LMS_rPos_repeat_bit_seqs.push_back( new bit_vector_type() );
            m_LMS_rPos_repeat_bit_seqs.back()->start_write();

            uint32 id_3 = id_1 + 2;

            uint64 beg_pos_3 = m_blocks_info[id_3]->m_beg_pos;
            uint64 end_pos_3 = m_blocks_info[id_3]->m_end_pos;

            bool is_single_3 = m_blocks_info[id_3]->is_single();

            /// thread 1
            std::thread th1(UtilityFunctions::th_reduceSort_bwt_4B<alphabet_type,compress_type,relative_offset_type>
                            , m_s
                            , beg_pos_1
                            , end_pos_1
                            , m_alpha
                            , m_level
                            , id_1
                            , std::ref(m_cLMS_seqs[id_1])
                            , std::ref(m_sub_l_cbwt_seqs[id_1])
                            , std::ref(m_sub_l_bit_seqs[id_1])
                            , std::ref(m_sub_s_cbwt_seqs[id_1])
                            , std::ref(m_sub_s_bit_seqs[id_1])
                            , std::ref(m_LMS_rPos_seqs[id_1])
                            , std::ref(m_LMS_rPos_repeat_bit_seqs[id_1])
                            , is_single_1
                            , is_sentinel_1
                           );


            /// thread 2
            std::thread th2(UtilityFunctions::th_reduceSort_bwt_4B<alphabet_type,compress_type,relative_offset_type>
                            , m_s
                            , beg_pos_2
                            , end_pos_2
                            , m_alpha
                            , m_level
                            , id_2
                            , std::ref(m_cLMS_seqs[id_2])
                            , std::ref(m_sub_l_cbwt_seqs[id_2])
                            , std::ref(m_sub_l_bit_seqs[id_2])
                            , std::ref(m_sub_s_cbwt_seqs[id_2])
                            , std::ref(m_sub_s_bit_seqs[id_2])
                            , std::ref(m_LMS_rPos_seqs[id_2])
                            , std::ref(m_LMS_rPos_repeat_bit_seqs[id_2])
                            , is_single_2
                            , is_sentinel_2
                           );

            /// thread 3
            std::thread th3(UtilityFunctions::th_reduceSort_bwt_4B<alphabet_type,compress_type,relative_offset_type>
                            , m_s
                            , beg_pos_3
                            , end_pos_3
                            , m_alpha
                            , m_level
                            , id_3
                            , std::ref(m_cLMS_seqs[id_3])
                            , std::ref(m_sub_l_cbwt_seqs[id_3])
                            , std::ref(m_sub_l_bit_seqs[id_3])
                            , std::ref(m_sub_s_cbwt_seqs[id_3])
                            , std::ref(m_sub_s_bit_seqs[id_3])
                            , std::ref(m_LMS_rPos_seqs[id_3])
                            , std::ref(m_LMS_rPos_repeat_bit_seqs[id_3])
                            , is_single_3
                            , is_sentinel_3
                           );

            th1.join();
            th2.join();
            th3.join();
        }

    }

}



template<typename alphabet_type, typename offset_type, typename compress_type, typename relative_offset_type>
void DSAComputation<alphabet_type, offset_type, compress_type, relative_offset_type>::compute_blocks_info( const uint64 _end_position)
{
    if (_end_position <= 2)
    {
        std::cout << "The parameter '_end_position' is ' <= 2'.\n";

        std::cout << "append the leftmost blocks information.\n";

        BlockInfo * block_info = new BlockInfo(_end_position, m_blocks_info.size());

        block_info->m_beg_pos = 0;

        /// need to count the LMS number
        block_info->m_lms_num = 1;


        m_blocks_info.push_back(block_info);

        return;
    }

    uint64 read_size = 0, read_beg_pos = 0, lms_num = 0;

    uint8 cur_type, last_type;


    if (_end_position + 1 >= m_block_size)
        read_size = m_block_size, read_beg_pos = _end_position + 1 - m_block_size;
    else
        read_size = _end_position + 1, read_beg_pos = 0;

    bool is_sentinel_on_level_0 = ((m_level == 0 && m_LMS_rPos_seqs.size() == 0) ? true : false);

    if(is_sentinel_on_level_0) read_size++;

    alphabet_type * read_buf = new alphabet_type[read_size];

    if(is_sentinel_on_level_0)
    {
        UtilityFunctions::fileRead<alphabet_type>(read_buf, m_s, read_beg_pos, read_size-1);
        read_buf[read_size - 1] = 0;/// append a virtual sentinel
    }
    else UtilityFunctions::fileRead<alphabet_type>(read_buf, m_s, read_beg_pos, read_size);

    ///scan the buffer from right to left to find the leftmost LMS character
    last_type = L_TYPE;

    uint64 i = read_size - 3;

    //std::cout << "is_sentinel_on_level_0 = " << is_sentinel_on_level_0 << std::endl;

    uint64 start_pos = read_size - 1; /// may greater than 2^32

    for (; i >= 0; --i)
    {
        if (read_buf[i] > read_buf[i + 1] || (read_buf[i] == read_buf[i + 1] && last_type == L_TYPE))
            cur_type = L_TYPE;
        else
            cur_type = S_TYPE;

        if (cur_type == L_TYPE && last_type == S_TYPE)
        {
            start_pos = i + 1;

            lms_num = lms_num + 1;
        }

        last_type = cur_type;

        if (i == 0) break;
    }

    //compute bwt, blockInfo
    if (start_pos == read_size - 1)
    {
        ///left most block, abandon
        if (read_beg_pos == 0)
        {

            std::cout << "This block is the leftmost block at line " << __LINE__ << " in function " << __FUNCTION__ << "." << std::endl;

            std::cout << "append the leftmost blocks information.\n";

            BlockInfo * block_info = new BlockInfo(_end_position, m_blocks_info.size());

            block_info->m_beg_pos = 0;

            /// need to count the LMS number
            block_info->m_lms_num = 1;

            m_blocks_info.push_back(block_info);

            delete[] read_buf;
            read_buf = nullptr;

            return ;
        }

        /// not leftmost block:
        /// (1) find the leftmost S* substring; (2) call function to compute BWTs
        {

            bool is_find_LMS = false;/// the leftmost LMS of a singleton block


            uint64 sblock_end_pos = _end_position, sblock_start_pos = 0; /// 'sblock' is short for single block

            delete[] read_buf;

            ///step 1:
            while(read_beg_pos != 0)
            {

                /// overlap a char (the last char)
                if (read_beg_pos + 1 >= m_block_size)
                    read_size = m_block_size, read_beg_pos = read_beg_pos - m_block_size + 1;
                else
                    read_size = read_beg_pos + 1, read_beg_pos = 0;

                read_buf = new alphabet_type[read_size];

                UtilityFunctions::fileRead<alphabet_type>(read_buf, m_s, read_beg_pos, read_size);

                alphabet_type last_char = read_buf[read_size - 1];/// the overlap char

                for (i = read_size - 2; i >= 0; --i)
                {
                    if (read_buf[i] > last_char || (read_buf[i] == last_char && last_type == L_TYPE))
                        cur_type = L_TYPE;
                    else
                        cur_type = S_TYPE;

                    if (cur_type == L_TYPE && last_type == S_TYPE)
                    {

                        is_find_LMS = true;

                        sblock_start_pos = read_beg_pos + i + 1; /// 'last_char' is LMS, its position is i + 1.

                        break;

                    }

                    last_type = cur_type;

                    if (i == 0) break;
                }

                delete[] read_buf;

                if(is_find_LMS) break;

            }

            /// step 2:
            if(is_find_LMS)
            {

                /// compute the block infor
                BlockInfo * block_info = new BlockInfo(_end_position, m_blocks_info.size());
                block_info->m_beg_pos = sblock_start_pos;
                block_info->m_lms_num = 1;
                m_blocks_info.push_back(block_info);

                /// set the beginning position of a new recursion
                read_beg_pos = sblock_start_pos;

            }
            else //return; /// it is the left block, direct return
            {

                std::cout << "append the leftmost blocks information.\n";

                BlockInfo * block_info = new BlockInfo(_end_position, m_blocks_info.size());

                block_info->m_beg_pos = 0;

                /// need to count the LMS number
                block_info->m_lms_num = 1;

                m_blocks_info.push_back(block_info);

                return ;
            }
        }


    }
    else
    {

        BlockInfo * block_info = new BlockInfo(_end_position, m_blocks_info.size());

        block_info->m_beg_pos = read_beg_pos + start_pos;

        /// need to count the LMS number
        block_info->m_lms_num = lms_num;

        m_blocks_info.push_back(block_info);

        read_beg_pos = read_beg_pos + start_pos;

        delete[] read_buf;
        read_buf = nullptr;
    }


//    std::cout << "Level " << m_level << ", compute blocks information recursively.\n" << std::endl;

//    std::cout << " read_beg_pos = " << read_beg_pos << std::endl;
    //recursion
    compute_blocks_info(read_beg_pos);

}



template<typename alphabet_type, typename offset_type, typename compress_type, typename relative_offset_type>
bool DSAComputation<alphabet_type, offset_type, compress_type, relative_offset_type>::mergeSortedSStarGlobal()
{

    /// global name for each substring
    offset_type name(0);

    if(m_alpha <= alphabet_type(255))
    {
        /// sort S* substrings For small alphabet ( <= 1B)
        name = sortSub_4S();
    }
    else
    {
        /// sort S* substrings For Large alphabet ( > 1B)
        name = sortSub_4L();
    }

    /// set next level alphabet set size.
    alphabetSet_vec.push_back(name);

    offset_type s1_size = 0;

    /// compute s1 size
    for(uint32 i = 0; i < s1_blocks.size(); ++i) s1_size += s1_blocks[i]->size();

    bool is_unique = true;

    /// compute is_unique
    if( name != s1_size) is_unique = false;
    else is_unique = true;


    /// generate s1
    if(!is_unique)
    {
        //generate_s1_5B(name);
        double beg_time = Timer::get_wall_time();

        uint64 get_s1_iv = Logger::cur_iv;

        uint64 get_s1_ov = Logger::cur_ov;

        if(name > offset_type(uint32_MAX - 1) ) generate_s1_5B(name);
        else if(name > offset_type(uint24_MAX - 1) ) generate_s1_4B(name);
        else if(name > offset_type(uint16_MAX - 1) ) generate_s1_3B(name);
        else if(name > offset_type(uint8_MAX - 1) ) generate_s1_2B(name);
        else generate_s1_1B(name);

        std::cout << "Generate the new string s1 of level_" + std::to_string(m_level) << std::endl;

        Timer::add_generate_s1_time(Timer::get_wall_time() - beg_time);

        Logger::add_get_s1_iv(Logger::cur_iv - get_s1_iv);

        Logger::add_get_s1_ov(Logger::cur_ov - get_s1_ov);

    }
    else
    {

        std::cout << "is_unique is true, no need to write the s1 to the disk.\n";
    }


    return is_unique;
}



/*
 * Description: sort S* substring 4 Small alphabet ( <= 1B )
 * Note:
 * (1) the left first bit of a variable 'block_id' set as 1 for different succeeding rank, or else for equal succeeding rank;
 * (2) the second bit of a variable 'block_id' is set as 1 for L* substring, or else 0 for no-L* substring;
 * (3) all ranks are denoted by uint64_MAX
 */
template<typename alphabet_type, typename offset_type, typename compress_type, typename relative_offset_type>
offset_type DSAComputation<alphabet_type, offset_type, compress_type, relative_offset_type>::sortSub_4S()
{
    std::cout << "Current Function: sortSub_4S().\n";

    Timer::show_time();

    uint64 red_iv = Logger::cur_iv;

    uint64 red_ov = Logger::cur_ov;

    double beg_time = Timer::get_wall_time();

    offset_type name(0);

    alphabet_type cur_char(0),pre_char(0);

    uint16 block_id = 0;

    bool is_LStar = false;

    uint64 g_rank = 0;/// global rank of each scanned substring

    /// define L* substring and bit sequence, the rank of last item
    std::vector<uint64> lstar_char_vector(m_alpha + 1, 0);
    MyVector<uint16> * lstar_blockID_seq = new MyVector<uint16>();
    lstar_blockID_seq->start_write();

    uint64  lstar_rank = uint64_MAX;

    /// define bucket-section type
    std::vector< MyVector<uint16> * > bks_seqs(m_alpha + 1, nullptr); /// bks_seqs(m_alpha + 1, nullptr);

    std::vector< uint64 > bks_ranks(m_alpha + 1, uint64_MAX); /// vector(number, value)

    std::cout << "m_alpha = " << m_alpha - 0<< std::endl;

    for(uint64 i = 0; i <= m_alpha; ++i)
    {
        bks_seqs[i] = new MyVector<uint16>();

        bks_seqs[i]->start_write();

        //std::cout <<  " i = " << i << std::endl;

    }


    /// define the variable of sorting L substring in same bucket. Note that the leftmost bit is the type of a substring
    InduceSortSameBkt< uint16 > * t_sortSameBkt = new InduceSortSameBkt< uint16 >(MAX_MEM / 3);
    uint64 t_name = uint64_MAX;


    /// initialize cPair_seq
    typedef Pair<alphabet_type, compress_type> cPair_type;

    /// allocate a temporary value for each compressed BWT pair
    std::vector<cPair_type> cPair_seq;

    for(uint64 i = 0; i < m_blocks_info.size()-1; i++) /// exclude the leftmost block
    {
        //std::cout << "i = " << i << std::endl;

        m_sub_l_cbwt_seqs[i]->start_read();

        /// assign the value to each cPair
        if(!m_sub_l_cbwt_seqs[i]->is_empty())
        {
            cPair_seq.push_back( m_sub_l_cbwt_seqs[i]->get() );
            m_sub_l_cbwt_seqs[i]->next_remove();
        }
        else
        {
            std::cout << "m_sub_l_cbwt_seqs[" << i << "]->is_empty() is true. Line: " << __LINE__ << std::endl;
        }

        m_sub_l_bit_seqs[i]->start_read_remove();


    }

    /// first put the preceding substring of the sentinel
    if(m_level==0)
    {

        g_rank++;

        /// the first block, i.e. the rightmost block
        block_id = 0;

        /// check whether to get next pair;
        if( !cPair_seq[block_id].second ) // if equal to zero, get next pair
        {
            cPair_seq[block_id] = m_sub_l_cbwt_seqs[block_id]->get();
            m_sub_l_cbwt_seqs[block_id]->next_remove();
        }

        is_LStar = m_sub_l_bit_seqs[block_id]->get();
        m_sub_l_bit_seqs[block_id]->next_remove();

        pre_char = cPair_seq[block_id].first;
        cPair_seq[block_id].second--; /// the number of first member decreased by one

        if( bks_ranks[pre_char] != g_rank ) block_id = (block_id | MASK_U16[1]); /// set the first bit as 1

        if(is_LStar) block_id = (block_id | MASK_U16[2]); /// set the second bit as 1

        bks_seqs[pre_char]->push_back( block_id );

        bks_ranks[pre_char] = g_rank;

    }

    std::cout << "Start to sort L-type substrings on level " << m_level << ".\n";

    std::cout << "The current PDU = " << Logger::max_pdu << ".\n";

    /// sort from 0 to m_alpha
    while( cur_char <= m_alpha)
    {

        t_sortSameBkt->reset();

        if (!bks_seqs[cur_char]->is_empty())
        {
            //std::cout << "current bks_seqs size = " << bks_seqs[cur_char]->size() << std::endl;

            /// step 1. induce sort the L-bucket
            bks_seqs[cur_char]->start_read();

            while (!bks_seqs[cur_char]->is_eof())
            {
                /// 1.1 get block_id
                block_id = bks_seqs[cur_char]->get();

                /// 1.2: get is_diff to compute name
                if (block_id & MASK_U16[1]) ++g_rank; /// the left first bit

                /// 1.3: get is_LStar
                is_LStar = block_id & MASK_U16[2]; /// the left second bit

                block_id = block_id & MASK_U16[3]; /// MASK_U16[3] = 0x3fff

                /// 1.4: if true, put this substring to 'lstar_sub_seq'
                if(is_LStar)
                {
                    if( g_rank != lstar_rank ) lstar_blockID_seq->push_back(block_id | MASK_U16[1]);
                    else
                        lstar_blockID_seq->push_back(block_id);

                    lstar_char_vector[cur_char]++;

                    lstar_rank = g_rank;
                }
                else   /// compute BWT
                {

                    /// check whether to get next pair;
                    if( !cPair_seq[block_id].second )
                    {
                        cPair_seq[block_id] = m_sub_l_cbwt_seqs[block_id]->get();
                        m_sub_l_cbwt_seqs[block_id]->next_remove();
                    }

                    pre_char = cPair_seq[block_id].first;
                    cPair_seq[block_id].second--;

                    is_LStar = m_sub_l_bit_seqs[block_id]->get();
                    m_sub_l_bit_seqs[block_id]->next_remove();

                    if(is_LStar) block_id = ( block_id | MASK_U16[2] );

                    if (pre_char > cur_char)
                    {
                        if (g_rank != bks_ranks[pre_char]) bks_seqs[pre_char]->push_back(block_id | MASK_U16[1]);
                        else bks_seqs[pre_char]->push_back( block_id );
                        bks_ranks[pre_char] = g_rank;
                    }
                    else if (pre_char == cur_char)
                    {
                        if (g_rank != t_name) t_sortSameBkt->push(block_id | MASK_U16[1]);
                        else t_sortSameBkt->push(block_id);
                        t_name = g_rank;
                    }
                    else
                    {
                        Logger::output_error(__FILE__,__LINE__);
                    }

                }

                bks_seqs[cur_char]->next_remove();

            }

            /// free disk space
            if (bks_seqs[cur_char]->is_eof())
            {
                delete bks_seqs[cur_char];
                bks_seqs[cur_char] = nullptr;
            }
            else
            {
                Logger::output_error(__FILE__,__LINE__);

            }


            /// step2: sort the same bkt
            if (t_sortSameBkt->write_size())
            {
                t_sortSameBkt->swapBuf();

                while (true)
                {
                    t_name = uint64_MAX; /// the g_rank must be different, t_name should be set as MAX_VALUE

                    while (t_sortSameBkt->read_size())
                    {

                        block_id = t_sortSameBkt->pop();

                        if (block_id & MASK_U16[1])
                            ++g_rank;

                        is_LStar = block_id & MASK_U16[2];

                        block_id = block_id & MASK_U16[3];

                        if(is_LStar)
                        {
                            if( g_rank != lstar_rank ) lstar_blockID_seq->push_back(block_id | MASK_U16[1]);
                            else lstar_blockID_seq->push_back(block_id);

                            lstar_char_vector[cur_char]++;

                            lstar_rank = g_rank;
                        }
                        else   /// cur_item is non-LStar-substring
                        {
                            /// check whether to get next pair;
                            if( !cPair_seq[block_id].second )
                            {
                                assert(!m_sub_l_cbwt_seqs[block_id]->is_eof());

                                cPair_seq[block_id] = m_sub_l_cbwt_seqs[block_id]->get();
                                m_sub_l_cbwt_seqs[block_id]->next_remove();
                            }

                            pre_char = cPair_seq[block_id].first;
                            cPair_seq[block_id].second--;

                            is_LStar = m_sub_l_bit_seqs[block_id]->get();
                            m_sub_l_bit_seqs[block_id]->next_remove();

                            if(is_LStar) block_id = ( block_id | MASK_U16[2] );

                            if (pre_char > cur_char)
                            {
                                if (g_rank != bks_ranks[pre_char]) bks_seqs[pre_char]->push_back(block_id | MASK_U16[1]);
                                else bks_seqs[pre_char]->push_back( block_id );
                                bks_ranks[pre_char] = g_rank;
                            }
                            else if (pre_char == cur_char)
                            {
                                if (g_rank != t_name) t_sortSameBkt->push(block_id | MASK_U16[1]);
                                else t_sortSameBkt->push(block_id);
                                t_name = g_rank;
                            }
                            else
                            {
                                Logger::output_error(__FILE__,__LINE__);
                            }

                        }

                    }

                    if (t_sortSameBkt->write_size())
                        t_sortSameBkt->swapBuf();
                    else
                        break;
                }
            }

        }
        else   /// if the current bks_seqs is empty, free it
        {
            delete bks_seqs[cur_char];
            bks_seqs[cur_char] = nullptr;
        }

        /// 3.induce sort the LMS-bucket
        ++g_rank;

        std::vector<uint64>::iterator it_beg = m_LMS_counter[cur_char]->begin(), it_end = m_LMS_counter[cur_char]->end();

        block_id = 0;

        while( it_beg != it_end)
        {

            uint64 cur_counter = *it_beg;

            while(cur_counter)
            {

                /// check whether to get next pair;
                if( !cPair_seq[block_id].second )
                {
                    cPair_seq[block_id] = m_sub_l_cbwt_seqs[block_id]->get();
                    m_sub_l_cbwt_seqs[block_id]->next_remove();
                }

                pre_char = cPair_seq[block_id].first;
                cPair_seq[block_id].second--;

                is_LStar = m_sub_l_bit_seqs[block_id]->get();
                m_sub_l_bit_seqs[block_id]->next_remove();


                if(is_LStar) block_id = ( block_id | MASK_U16[2] );

                if (pre_char > cur_char)
                {
                    if (g_rank != bks_ranks[pre_char]) bks_seqs[pre_char]->push_back(block_id | MASK_U16[1]);
                    else bks_seqs[pre_char]->push_back( block_id );
                    bks_ranks[pre_char] = g_rank;
                }
                else
                {
                    Logger::output_error(__FILE__,__LINE__);
                }

                block_id = block_id & MASK_U16[3]; /// error-prone

                cur_counter--;
            }

            it_beg++;

            block_id = block_id & MASK_U16[3]; /// The left two bit must be set to 0; CAUTION
            block_id++;
        }

        /// next bucket
        if(cur_char == m_alpha)break;
        else ++cur_char;
    }

    /// destructor
    {
        while (!m_sub_l_cbwt_seqs.empty())
        {
            if (!m_sub_l_cbwt_seqs.back()->is_eof())
            {
                Logger::output_error(__FILE__,__LINE__);
            }
            else
            {
                delete m_sub_l_cbwt_seqs.back();
            }

            m_sub_l_cbwt_seqs.pop_back();
        }

        while (!m_sub_l_bit_seqs.empty())
        {
            if (!m_sub_l_bit_seqs.back()->is_eof())
            {
                Logger::output_error(__FILE__,__LINE__);
            }
            else
            {
                delete m_sub_l_bit_seqs.back();
            }

            m_sub_l_bit_seqs.pop_back();
        }

    }

    std::cout << "------------------------------------------------------------------------------------------------------------------------------.\n";

    std::cout << "L-type substrings sorting on level_" << m_level << " is over.\n";

    std::cout << "The current PDU = " << Logger::max_pdu << ".\n";

    std::cout << "------------------------------------------------------------------------------------------------------------------------------.\n";

    std::cout << "Start to sort S-type substrings on level " << m_level << ".\n";

    Logger::output_ram_use();


    /// define the rank of the last S* substring
    uint64  sstar_rank = 0;

    g_rank = uint64_MAX; /// from biggest to smallest

    bool is_SStar = false;

    /// redefine bucket-section type
    bks_seqs.clear();
    bks_ranks.clear();

    for(uint64 i = 0; i <= m_alpha; ++i)
    {
        bks_seqs.push_back(new MyVector<uint16>());
        bks_seqs.back()->start_write();

        bks_ranks.push_back(0);
    }

    /// allocate a temporary value for each compressed BWT pair
    cPair_seq.clear();;

    for(uint64 i = 0; i < m_blocks_info.size()-1; i++) /// exclude the leftmost block
    {

        m_sub_s_cbwt_seqs[i]->start_read();

        /// assign the value to each cPair
        if(!m_sub_s_cbwt_seqs[i]->is_empty())
        {
            cPair_seq.push_back( m_sub_s_cbwt_seqs[i]->get() );
            m_sub_s_cbwt_seqs[i]->next_remove();
        }
        else
        {
            std::cout << "m_sub_s_cbwt_seqs[" << i << "]->is_empty() is true. Line: " << __LINE__ << std::endl;
            //Logger::output_error(__FILE__,__LINE__);
        }

        m_sub_s_bit_seqs[i]->start_read_remove();

    }

    /// start to read the positions of LMS and different bit sequence
    for(uint64 i = 0; i < m_blocks_info.size() - 1; i++) /// exclude the leftmost block
    {
        m_LMS_rPos_seqs[i]->start_read();
        m_LMS_rPos_repeat_bit_seqs[i]->start_read_remove();
    }

    /// initialize the s1_blocks
    for(uint64 i = 0; i < m_blocks_info.size() - 1; ++i)
    {
        s1_blocks.push_back( new s1_blocks_vector_type());
        s1_blocks.back()->start_write();
    }


    /// alphabet block vector
    std::vector< Alpha_Block<offset_type> > alpha_block_vec;

    /// the current alphabet block size and current alphabet size
    offset_type total_alpha_block_size = 0, cur_alpha_size = 0, pre_name(1), cur_name(1);

    /// < ch, rank, block_id> = < 4B, sizeof(relative_offset_type), 2B >
    /// bug, at inducing stage, using double working space
    offset_type alpha_block_max_size = (MAX_MEM / 4) / (4 + sizeof(relative_offset_type) + 2) / 2;
    //offset_type alpha_block_max_size = K_1024 * 1024 * 2; /// \note for test

    std::cout << "alpha_block_max_size = " << alpha_block_max_size << std::endl;

    /// start to read the L*-substrings reversely
    lstar_blockID_seq->start_read_reverse();

    cur_char = m_alpha;

    while(cur_char >= alphabet_type(0))
    {

        t_sortSameBkt->reset();

        if (!bks_seqs[cur_char]->is_empty())
        {
            bks_seqs[cur_char]->start_read();

            while (!bks_seqs[cur_char]->is_eof())
            {

                block_id = bks_seqs[cur_char]->get();

                if (block_id & MASK_U16[1]) --g_rank;

                is_SStar = block_id & MASK_U16[2];

                block_id = block_id & MASK_U16[3];

                if(is_SStar)
                {
                    if(sstar_rank != g_rank) ++name;

                    /// the first bit item must be true
                    if( !m_LMS_rPos_repeat_bit_seqs[block_id]->get() )
                    {
                        Logger::output_error(__FILE__,__LINE__);
                    }
                    m_LMS_rPos_repeat_bit_seqs[block_id]->next_remove();

                    /// push the first LMS position to its text block
                    s1_blocks[block_id]->push_back(s1_pair_type(m_LMS_rPos_seqs[block_id]->get(), name));
                    m_LMS_rPos_seqs[block_id]->next_remove();

                    compute_alphabet_blocks(
                        alpha_block_vec,
                        total_alpha_block_size,
                        alpha_block_max_size,
                        cur_alpha_size,
                        pre_name,
                        name,
                        ((sstar_rank != g_rank)? true : false)
                    );

                    /// push the position of the same LMS substrings
                    while(!m_LMS_rPos_repeat_bit_seqs[block_id]->is_eof())
                    {

                        if(m_LMS_rPos_repeat_bit_seqs[block_id]->get())break;
                        else
                        {
                            s1_blocks[block_id]->push_back(s1_pair_type(m_LMS_rPos_seqs[block_id]->get(), name));
                            m_LMS_rPos_seqs[block_id]->next_remove();

                            compute_alphabet_blocks(
                                alpha_block_vec,
                                total_alpha_block_size,
                                alpha_block_max_size,
                                cur_alpha_size,
                                pre_name,
                                name,
                                false
                            );
                        }
                        m_LMS_rPos_repeat_bit_seqs[block_id]->next_remove();
                    }

                    sstar_rank = g_rank;
                }
                else
                {
                    /// check whether to get next pair;
                    if( !cPair_seq[block_id].second )
                    {
                        cPair_seq[block_id] = m_sub_s_cbwt_seqs[block_id]->get();
                        m_sub_s_cbwt_seqs[block_id]->next_remove();
                    }

                    pre_char = cPair_seq[block_id].first;
                    cPair_seq[block_id].second--;

                    is_SStar = m_sub_s_bit_seqs[block_id]->get();
                    m_sub_s_bit_seqs[block_id]->next_remove();

                    if(is_SStar) block_id = ( block_id | MASK_U16[2] );

                    if (pre_char < cur_char)
                    {
                        if (g_rank != bks_ranks[pre_char]) block_id = (block_id | MASK_U16[1]);
                        bks_seqs[pre_char]->push_back( block_id );
                        bks_ranks[pre_char] = g_rank;
                    }
                    else if (pre_char == cur_char)
                    {
                        if (g_rank != t_name) block_id = (block_id | MASK_U16[1]);
                        t_sortSameBkt->push(block_id);
                        t_name = g_rank;
                    }
                    else
                    {
                        Logger::output_error(__FILE__,__LINE__);
                    }

                }

                bks_seqs[cur_char]->next_remove();
            }

            if (bks_seqs[cur_char]->is_eof())
            {
                delete bks_seqs[cur_char];
                bks_seqs[cur_char] = nullptr;
            }

            /// step2: Induce sort the sameBkt
            if (t_sortSameBkt->write_size())
            {
                t_sortSameBkt->swapBuf();

                while (true)
                {
                    t_name = 0;

                    while (t_sortSameBkt->read_size())
                    {
                        block_id = t_sortSameBkt->pop();

                        if (block_id & MASK_U16[1])
                            --g_rank;

                        is_SStar = block_id & MASK_U16[2];

                        block_id = block_id & MASK_U16[3];

                        if(is_SStar)
                        {
                            if(sstar_rank != g_rank) ++name;

                            /// the first bit item must be true
                            if( !m_LMS_rPos_repeat_bit_seqs[block_id]->get() )
                            {
                                Logger::output_error(__FILE__,__LINE__);
                            }
                            m_LMS_rPos_repeat_bit_seqs[block_id]->next_remove();

                            /// push the first LMS position
                            s1_blocks[block_id]->push_back(s1_pair_type(m_LMS_rPos_seqs[block_id]->get(), name));
                            m_LMS_rPos_seqs[block_id]->next_remove();

                            compute_alphabet_blocks(
                                alpha_block_vec,
                                total_alpha_block_size,
                                alpha_block_max_size,
                                cur_alpha_size,
                                pre_name,
                                name,
                                ((sstar_rank != g_rank)? true : false)
                            );

                            /// push the position of the same LMS substrings
                            while(!m_LMS_rPos_repeat_bit_seqs[block_id]->is_eof())
                            {

                                if(m_LMS_rPos_repeat_bit_seqs[block_id]->get())break;
                                else
                                {
                                    s1_blocks[block_id]->push_back(s1_pair_type(m_LMS_rPos_seqs[block_id]->get(), name));
                                    m_LMS_rPos_seqs[block_id]->next_remove();

                                    compute_alphabet_blocks(
                                        alpha_block_vec,
                                        total_alpha_block_size,
                                        alpha_block_max_size,
                                        cur_alpha_size,
                                        pre_name,
                                        name,
                                        false
                                    );

                                }
                                m_LMS_rPos_repeat_bit_seqs[block_id]->next_remove();
                            }

                            sstar_rank = g_rank;

                        }
                        else
                        {
                            /// check whether to get next pair;
                            if( !cPair_seq[block_id].second )
                            {
                                cPair_seq[block_id] = m_sub_s_cbwt_seqs[block_id]->get();
                                m_sub_s_cbwt_seqs[block_id]->next_remove();
                            }

                            pre_char = cPair_seq[block_id].first;
                            cPair_seq[block_id].second--;

                            is_SStar = m_sub_s_bit_seqs[block_id]->get();
                            m_sub_s_bit_seqs[block_id]->next_remove();

                            if(is_SStar) block_id = ( block_id | MASK_U16[2] );

                            if (pre_char < cur_char)
                            {
                                if (g_rank != bks_ranks[pre_char]) block_id = (block_id | MASK_U16[1]);
                                bks_seqs[pre_char]->push_back( block_id );
                                bks_ranks[pre_char] = g_rank;
                            }
                            else if (pre_char == cur_char)
                            {
                                if (g_rank != t_name) block_id = (block_id | MASK_U16[1]);
                                t_sortSameBkt->push(block_id);
                                t_name = g_rank;
                            }
                            else
                            {
                                Logger::output_error(__FILE__,__LINE__);
                            }

                        }

                    }

                    if (t_sortSameBkt->write_size())
                        t_sortSameBkt->swapBuf();
                    else
                        break;
                }

            }
        }
        else
        {
            delete bks_seqs[cur_char];
            bks_seqs[cur_char] = nullptr;
        }

        --g_rank;

        /// step3: Induce sort the LStar-type
        while (lstar_char_vector[cur_char])
        {

            /// set the leftmost bit of 'block_id' as 1
            block_id = lstar_blockID_seq->get_reverse() & MASK_U16[4]; /// MASK_U16[4] = 0x7fff

            /// check whether to get next pair;
            if( !cPair_seq[block_id].second )
            {
                cPair_seq[block_id] = m_sub_s_cbwt_seqs[block_id]->get();
                m_sub_s_cbwt_seqs[block_id]->next_remove();
            }

            pre_char = cPair_seq[block_id].first;
            cPair_seq[block_id].second--;

            is_SStar = m_sub_s_bit_seqs[block_id]->get();
            m_sub_s_bit_seqs[block_id]->next_remove();

            if(is_SStar) block_id = ( block_id | MASK_U16[2] );

            if (pre_char < cur_char)
            {
                if (g_rank != bks_ranks[pre_char]) bks_seqs[pre_char]->push_back(block_id | MASK_U16[1]);
                else bks_seqs[pre_char]->push_back( block_id );
                bks_ranks[pre_char] = g_rank;
            }
            else
            {
                Logger::output_error(__FILE__,__LINE__);
            }

            if (lstar_blockID_seq->get_reverse() & MASK_U16[1])
                --g_rank;

            lstar_blockID_seq->next_remove_reverse();

            lstar_char_vector[cur_char]--;

        }

        if(!cur_char) break;
        else --cur_char;
    }

    std::cout << "------------------------------------------------------------------------------------------------------------------------------.\n";

    std::cout << "S-type substrings sorting on level_" << m_level << " is over.\n";

    std::cout << "The current PDU = " << Logger::max_pdu << ".\n";

    std::cout << "------------------------------------------------------------------------------------------------------------------------------.\n";

    /// destruct
    {
        while(bks_seqs.empty())
        {
            if(bks_seqs.back())
            {
                delete bks_seqs.back();
            }
            bks_seqs.pop_back();
        }

        delete t_sortSameBkt;

        if(lstar_blockID_seq->is_eof())
        {
            delete lstar_blockID_seq;
        }
        else
        {
            Logger::output_error(__FILE__,__LINE__);

        }
        while(!m_LMS_rPos_repeat_bit_seqs.empty())
        {
            if(m_LMS_rPos_repeat_bit_seqs.back()->is_eof())delete m_LMS_rPos_repeat_bit_seqs.back();
            else
            {
                Logger::output_error(__FILE__,__LINE__);
            }
            m_LMS_rPos_repeat_bit_seqs.pop_back();
        }

        while(!m_LMS_rPos_seqs.empty())
        {
            if(m_LMS_rPos_seqs.back()->is_eof())delete m_LMS_rPos_seqs.back();
            else
            {
                Logger::output_error(__FILE__,__LINE__);
            }
            m_LMS_rPos_seqs.pop_back();
        }


        while(!m_sub_s_cbwt_seqs.empty())
        {
            if(m_sub_s_cbwt_seqs.back()->is_eof())delete m_sub_s_cbwt_seqs.back();
            else
            {
                Logger::output_error(__FILE__,__LINE__);
            }
            m_sub_s_cbwt_seqs.pop_back();
        }

        while(!m_sub_s_bit_seqs.empty())
        {
            if(m_sub_s_bit_seqs.back()->is_eof())delete m_sub_s_bit_seqs.back();
            else
            {
                Logger::output_error(__FILE__,__LINE__);
            }
            m_sub_s_bit_seqs.pop_back();
        }

    }

    /**< compute the tail elements */
    /// Notice that the size of last alphabet plus the total_alpha_block_size exceeds the alpha_block_max_size,
    cur_name = name;

    if(pre_name == cur_name)
    {
        alpha_block_vec.push_back(Alpha_Block<offset_type>(pre_name, cur_name, cur_alpha_size));
    }
    else
    {
        if(cur_alpha_size >= alpha_block_max_size)
        {
            alpha_block_vec.push_back(Alpha_Block<offset_type>(pre_name, cur_name - 1, total_alpha_block_size));
            alpha_block_vec.push_back(Alpha_Block<offset_type>(cur_name, cur_name, cur_alpha_size));
        }
        else
        {
            alpha_block_vec.push_back(Alpha_Block<offset_type>(pre_name, cur_name, cur_alpha_size + total_alpha_block_size));
        }
    }

    /*	check whether the beginning alphabet is less or equal to the ending alphabet in each bucket-section	*/
    while(true)
    {
        bool is_error(false);
        for(uint64 i = 0; i < alpha_block_vec.size(); ++i)
        {

            std::cout << i << "th alpha_block: " << alpha_block_vec[i].m_beg_alpha << " , " << alpha_block_vec[i].m_end_alpha << ".\n";

            if(alpha_block_vec[i].m_beg_alpha > alpha_block_vec[i].m_end_alpha)
            {
                std::cout << "alpha_block_vec[" << i <<"] is error.\n";
                is_error = true;
            }
        }
        if(is_error)
        {
            std::cout << "Alphabet blocks information array is error.\n";

            Logger::output_error(__FILE__,__LINE__);

        }

        break;

    }

    {
        /// compute the correct alpha_block_info

        std::vector< Alpha_Block<offset_type> > alpha_block_vec_copy;

        int index = 0;

        while(!alpha_block_vec.empty())
        {

            alpha_block_vec_copy.push_back( Alpha_Block<offset_type>( name - alpha_block_vec.back().m_end_alpha + offset_type(1), name - alpha_block_vec.back().m_beg_alpha + offset_type(1), alpha_block_vec.back().m_size) );

            std::cout << index++ << " : (" << alpha_block_vec_copy.back().m_beg_alpha << ", " << alpha_block_vec_copy.back().m_end_alpha << "," << alpha_block_vec_copy.back().m_size << ")" << std::endl;

            alpha_block_vec.pop_back();
        }

        /**< write the LMS alphabet block infor to disk, filename is alpha_blocks_info_Level_$level.dat */
        FILE * fot = fopen(("alpha_blocks_info_" + std::to_string(m_level) + ".dat").c_str(), "wb");

        fwrite(alpha_block_vec_copy.data(), sizeof(alpha_block_vec_copy[0]), alpha_block_vec_copy.size(), fot);

        fclose(fot);
    }

    Timer::show_time();

    Timer::add_reduce_sort_time(Timer::get_wall_time() - beg_time);

    Logger::add_red_iv( Logger::cur_iv - red_iv);

    Logger::add_red_ov( Logger::cur_ov - red_ov);

    return name;
}


/// \brief sort S* substring 4 Large alphabet ( >= 1B )
template<typename alphabet_type, typename offset_type, typename compress_type, typename relative_offset_type>
offset_type DSAComputation<alphabet_type, offset_type, compress_type, relative_offset_type>::sortSub_4L()
{
    std::cout << "the current function sortSub_4L()" << std::endl;

    Timer::show_time();

    uint64 red_iv = Logger::cur_iv;

    uint64 red_ov = Logger::cur_ov;

    double beg_time = Timer::get_wall_time();

    offset_type name(0);

    uint16 block_id = 0;

    bool is_LStar = false;

    uint64 g_rank = 0;/// global rank

    /// compress max count
    compress_type c_max_count(std::numeric_limits<compress_type>::max()), lstar_cnt(0);

    alphabet_type pre_lstar(std::numeric_limits<alphabet_type>::max());

    /// compress L* char
    MyVector< compressed_pair_type > * lstar_sub_char_seq = new MyVector< compressed_pair_type >();
    lstar_sub_char_seq->start_write();

    MyVector<uint16> * lstar_sub_blockId_seq = new MyVector<uint16>();
    lstar_sub_blockId_seq->start_write();

    uint64  lstar_rank = uint64_MAX;


    /// define bucket-section type
    typedef Pair< relative_offset_type, uint16> bks_pair_type;
    typedef MyVector< bks_pair_type > bks_pair_type_vector;


    std::vector< bks_pair_type_vector * > bks_seqs; /// bks_seqs(m_alpha + 1, nullptr);

    std::vector< uint64 > bks_ranks( m_alpha_blocks_number, uint64_MAX); /// vector(number, value)

    for(uint64 i = 0; i < m_alpha_blocks_number; ++i)
    {
        bks_seqs.push_back(new bks_pair_type_vector());
        bks_seqs.back()->start_write();
    }

    ///compress LMS sorter
    c_LMSHeapSorter<alphabet_type,compress_type> * cLMS_sorter = new c_LMSHeapSorter<alphabet_type, compress_type>(m_cLMS_seqs); // construct lmsSorter

    cLMS_sorter->initBkt(); // get the bucket

    /// initialize cPair_seq for compressed bwt
    //typedef Pair< alphabet_type, compress_type> cPair_type;

    /// allocate a temporary value for each compressed BWT pair
    std::vector<compressed_pair_type> cPair_seq(m_blocks_info.size());

    for(uint32 i = 0; i < m_blocks_info.size()-1; i++) /// exclude the leftmost block
    {

        m_sub_l_cbwt_seqs[i]->start_read();

        /// assign the value to each cPair
        if(!m_sub_l_cbwt_seqs[i]->is_empty())
        {
            cPair_seq[i] = m_sub_l_cbwt_seqs[i]->get();
            m_sub_l_cbwt_seqs[i]->next_remove();
        }
        else
        {
            std::cout << "m_sub_l_cbwt_seqs[" << i << "]->is_empty() is true. Line: " << __LINE__ << std::endl;
        }

        m_sub_l_bit_seqs[i]->start_read_remove();

    }


    std::cout << "Start to sort L-type substrings on level " << m_level << ".\n";

    std::cout << "The current PDU = " << Logger::max_pdu << ".\n";

    uint16 cur_alpha_block_id(0);

    /// < alphabet, s_rank, blockID >
    typedef Triple<relative_offset_type, relative_offset_type, uint16> pa_item_type;

    /// < r_alphabet, s_rank, count, blockID>
    typedef quadruple< relative_offset_type, relative_offset_type, relative_offset_type, uint16> pq_item_type;

    typedef std::priority_queue< pq_item_type, std::vector<pq_item_type>, TupleDscCmp3< pq_item_type > > min_pq_type;

    /// sort from 0 to m_alpha
    while( cur_alpha_block_id < m_alpha_blocks_number)
    {

        std::cout << "cur_alpha_block_id = " << cur_alpha_block_id << std::endl;

        offset_type cur_alpha_block_beg(m_alpha_blocks_info[cur_alpha_block_id].m_beg_alpha);

        offset_type cur_alpha_block_end(m_alpha_blocks_info[cur_alpha_block_id].m_end_alpha);

        relative_offset_type alpha_number = cur_alpha_block_end - cur_alpha_block_beg + 1;

        std::cout << "cur_alpha_block_beg = " << cur_alpha_block_beg << std::endl;

        std::cout << "cur_alpha_block_end = " << cur_alpha_block_end << std::endl;

        std::cout << "alpha_number = " << alpha_number << std::endl;

        if(alpha_number > 1)
        {

            std::cout << "the current alpha block is not single block.\n";

            std::vector< pa_item_type > pa_item_ary;

            typename std::vector< pa_item_type >::iterator it_beg, it_end;

            if(!bks_seqs[cur_alpha_block_id]->is_empty())
            {

                //pa_item_ary.resize( bks_seqs[cur_alpha_block_id]->size() ); /// \note

                pa_item_ary.reserve( bks_seqs[cur_alpha_block_id]->size() ); /// \note
                pa_item_ary.clear();

                bks_seqs[cur_alpha_block_id]->start_read();

                //relative_offset_type index(0), s_rank(1);
                relative_offset_type s_rank(1);

                while( !bks_seqs[cur_alpha_block_id]->is_eof() )
                {

                    if( bks_seqs[cur_alpha_block_id]->get().second & MASK_U16[1] ) ++s_rank;

                    pa_item_ary.push_back(pa_item_type(bks_seqs[cur_alpha_block_id]->get().first, s_rank, bks_seqs[cur_alpha_block_id]->get().second));

                    //++index;

                    bks_seqs[cur_alpha_block_id]->next_remove();
                }

                delete bks_seqs[cur_alpha_block_id];

                bks_seqs[cur_alpha_block_id] = nullptr;

                std::stable_sort(pa_item_ary.begin(), pa_item_ary.end(), TupleAscCmp1<pa_item_type>());

            }
            else
            {

                delete bks_seqs[cur_alpha_block_id];
                bks_seqs[cur_alpha_block_id] = nullptr;
            }

            it_beg = pa_item_ary.begin(), it_end = pa_item_ary.end();

            relative_offset_type cur_rel_bkt(0); /// relative alpha

            relative_offset_type t_rel_rank(0), t_rank(std::numeric_limits<relative_offset_type>::max());

            relative_offset_type pq_input_cnt(0);

            min_pq_type * m_pq = new min_pq_type();

            alphabet_type cur_char(cur_alpha_block_beg), pre_char(0);

            uint16 pre_blockID(0);

            while(cur_rel_bkt < alpha_number)
            {
                cur_char = cur_alpha_block_beg + cur_rel_bkt;

                t_rank = std::numeric_limits<relative_offset_type>::max();

                while ( it_beg != it_end && (*it_beg).first == cur_rel_bkt )
                {
                    if( (*it_beg).second != t_rank )
                    {

                        ++t_rel_rank;

                        ++g_rank;

                        t_rank = (*it_beg).second;

                    }

                    is_LStar = (*it_beg).third & MASK_U16[2]; /// the left second bit

                    /// the left two bit
                    block_id = (*it_beg).third & MASK_U16[3]; /// MASK_U16[3] = 0x3fff

                    /// 1.4: if true, put this substring to 'lstar_sub_seq'
                    if(is_LStar)
                    {
                        if( g_rank != lstar_rank ) lstar_sub_blockId_seq->push_back(block_id | MASK_U16[1]);
                        else lstar_sub_blockId_seq->push_back(block_id);

                        if(cur_char != pre_lstar)
                        {

                            lstar_sub_char_seq->push_back(compressed_pair_type(pre_lstar,lstar_cnt));

                            pre_lstar = cur_char, lstar_cnt = 1;

                        }
                        else
                        {

                            if(lstar_cnt != c_max_count)
                            {
                                lstar_cnt++;
                            }
                            else
                            {
                                lstar_sub_char_seq->push_back(compressed_pair_type(cur_char,lstar_cnt));
                                lstar_cnt = 1;
                            }

                        }

                        lstar_rank = g_rank;
                    }
                    else   /// compute BWT
                    {

                        /// check whether to get next pair;
                        if( !cPair_seq[block_id].second )
                        {
                            cPair_seq[block_id] = m_sub_l_cbwt_seqs[block_id]->get();
                            m_sub_l_cbwt_seqs[block_id]->next_remove();
                        }

                        pre_char = cPair_seq[block_id].first;
                        --cPair_seq[block_id].second;

                        is_LStar = m_sub_l_bit_seqs[block_id]->get();
                        m_sub_l_bit_seqs[block_id]->next_remove();

                        if(is_LStar) block_id = ( block_id | MASK_U16[2] );

                        pre_blockID = get_alpha_blockID(cur_alpha_block_id,pre_char);

                        if ( pre_blockID > cur_alpha_block_id )
                        {
                            if (g_rank != bks_ranks[pre_blockID]) block_id = (block_id | MASK_U16[1]);
                            bks_seqs[pre_blockID]->push_back(bks_pair_type(pre_char - m_alpha_blocks_info[pre_blockID].m_beg_alpha,block_id));
                            bks_ranks[pre_blockID] = g_rank;
                        }
                        else if (pre_blockID == cur_alpha_block_id)
                        {
                            m_pq->push( pq_item_type(relative_offset_type(pre_char - cur_alpha_block_beg), t_rel_rank, ++pq_input_cnt, block_id) );
                        }
                        else
                        {
                            Logger::output_error(__FILE__,__LINE__);
                        }

                    }

                    it_beg++;


                }

                ++g_rank;

                t_rank = std::numeric_limits<relative_offset_type>::max();

                /// step2: sort PQ
                while(!m_pq->empty() && m_pq->top().first == cur_rel_bkt)
                {

                    if (m_pq->top().second != t_rank)
                    {
                        ++g_rank;

                        ++t_rel_rank;

                        t_rank = m_pq->top().second;
                    }

                    is_LStar = m_pq->top().forth & MASK_U16[2]; /// the left second bit

                    /// the left two bit
                    block_id = m_pq->top().forth & MASK_U16[3]; /// MASK_U16[3] = 0x3fff

                    /// 1.4: if true, put this substring to 'lstar_sub_seq'
                    if(is_LStar)
                    {
                        if( g_rank != lstar_rank ) lstar_sub_blockId_seq->push_back(block_id | MASK_U16[1]);
                        else lstar_sub_blockId_seq->push_back(block_id);

                        if(cur_char != pre_lstar)
                        {
                            lstar_sub_char_seq->push_back(compressed_pair_type(pre_lstar,lstar_cnt));
                            pre_lstar = cur_char, lstar_cnt = 1;
                        }
                        else
                        {
                            if(lstar_cnt != c_max_count)
                            {
                                ++lstar_cnt;
                            }
                            else
                            {
                                lstar_sub_char_seq->push_back(compressed_pair_type(cur_char,lstar_cnt));
                                lstar_cnt = 1;
                            }

                        }
                        lstar_rank = g_rank;
                    }
                    else   /// compute BWT
                    {

                        /// check whether to get next pair;
                        if( !cPair_seq[block_id].second )
                        {
                            cPair_seq[block_id] = m_sub_l_cbwt_seqs[block_id]->get();
                            m_sub_l_cbwt_seqs[block_id]->next_remove();
                        }

                        pre_char = cPair_seq[block_id].first;
                        --cPair_seq[block_id].second;

                        is_LStar = m_sub_l_bit_seqs[block_id]->get();
                        m_sub_l_bit_seqs[block_id]->next_remove();

                        if(is_LStar) block_id = ( block_id | MASK_U16[2] );

                        pre_blockID = get_alpha_blockID(cur_alpha_block_id,pre_char);

                        if ( pre_blockID > cur_alpha_block_id )
                        {
                            if (g_rank != bks_ranks[pre_blockID]) block_id = (block_id | MASK_U16[1]);
                            bks_seqs[pre_blockID]->push_back(bks_pair_type(pre_char - m_alpha_blocks_info[pre_blockID].m_beg_alpha,block_id));
                            bks_ranks[pre_blockID] = g_rank;
                        }
                        else if (pre_blockID == cur_alpha_block_id)
                        {
                            m_pq->push( pq_item_type(relative_offset_type(pre_char - cur_alpha_block_beg), t_rel_rank, ++pq_input_cnt, block_id) );
                        }
                        else
                        {
                            Logger::output_error(__FILE__,__LINE__);
                        }

                    }

                    m_pq->pop();
                }


                /// step3: induce sort the LMS-bucket
                ++g_rank;

                ++t_rel_rank;

                while(!cLMS_sorter->is_empty() && cLMS_sorter->get_bkt() == cur_char)
                {

                    block_id = cLMS_sorter->get_block_id();

                    /// check whether to get next pair;
                    if( !cPair_seq[block_id].second )
                    {
                        cPair_seq[block_id] = m_sub_l_cbwt_seqs[block_id]->get();
                        m_sub_l_cbwt_seqs[block_id]->next_remove();
                    }

                    pre_char = cPair_seq[block_id].first;
                    --cPair_seq[block_id].second;

                    is_LStar = m_sub_l_bit_seqs[block_id]->get();
                    m_sub_l_bit_seqs[block_id]->next_remove();

                    if(is_LStar) block_id = ( block_id | MASK_U16[2] );

                    pre_blockID = get_alpha_blockID(cur_alpha_block_id,pre_char);

                    if ( pre_blockID > cur_alpha_block_id )
                    {
                        if (g_rank != bks_ranks[pre_blockID]) block_id = (block_id | MASK_U16[1]);
                        bks_seqs[pre_blockID]->push_back(bks_pair_type(pre_char - m_alpha_blocks_info[pre_blockID].m_beg_alpha,block_id));
                        bks_ranks[pre_blockID] = g_rank;
                    }
                    else if (pre_blockID == cur_alpha_block_id)
                    {
                        m_pq->push( pq_item_type(relative_offset_type(pre_char - cur_alpha_block_beg), t_rel_rank, ++pq_input_cnt, block_id) );
                    }
                    else
                    {
                        std::cout << "block_id = " << (block_id & MASK_U16[3]) << std::endl;
                        std::cout << "pre_blockID = " << pre_blockID << std::endl;
                        std::cout << "pre_char = " << pre_char << std::endl;
                        std::cout << "cur_char = " << cur_char << std::endl;
                        Logger::output_error(__FILE__,__LINE__);
                    }

                    cLMS_sorter->next();
                }

                /// next r_bucket
                ++cur_rel_bkt;
            }

            delete m_pq;
        }
        else if( alpha_number == 1)
        {

            //relative_offset_type cur_rel_bkt(0); /// relative alpha

            //relative_offset_type t_rel_rank(0);

            //relative_offset_type pq_input_cnt(0);

            alphabet_type cur_char(cur_alpha_block_beg), pre_char(0);

            uint16 pre_blockID(0);

            InduceSortSameBkt< uint16 > * t_sortSameBkt = new InduceSortSameBkt< uint16 >( K_1024 * 1024 * 1 );
            uint64 t_rank = uint64_MAX;

            /// step 1: sort the bks_seqs item
            if(!bks_seqs[cur_alpha_block_id]->is_empty())
            {

                bks_seqs[cur_alpha_block_id]->start_read();

                while(!bks_seqs[cur_alpha_block_id]->is_eof())
                {

                    if(bks_seqs[cur_alpha_block_id]->get().second & MASK_U16[1]) ++g_rank;

                    is_LStar = bks_seqs[cur_alpha_block_id]->get().second & MASK_U16[2];

                    /// the left two bit
                    block_id = bks_seqs[cur_alpha_block_id]->get().second & MASK_U16[3]; /// MASK_U16[3] = 0x3fff

                    /// if true, put this substring to 'lstar_sub_seq'
                    if(is_LStar)
                    {
                        if( g_rank != lstar_rank ) lstar_sub_blockId_seq->push_back(block_id | MASK_U16[1]);
                        else lstar_sub_blockId_seq->push_back(block_id);

                        if(cur_char != pre_lstar)
                        {

                            lstar_sub_char_seq->push_back(compressed_pair_type(pre_lstar,lstar_cnt));

                            pre_lstar = cur_char, lstar_cnt = 1;

                        }
                        else
                        {

                            if(lstar_cnt != c_max_count)
                            {
                                lstar_cnt++;
                            }
                            else
                            {
                                lstar_sub_char_seq->push_back(compressed_pair_type(cur_char,lstar_cnt));
                                lstar_cnt = 1;
                            }

                        }

                        lstar_rank = g_rank;
                    }
                    else   /// compute BWT
                    {
                        /// check whether to get next pair;
                        if( !cPair_seq[block_id].second )
                        {
                            cPair_seq[block_id] = m_sub_l_cbwt_seqs[block_id]->get();
                            m_sub_l_cbwt_seqs[block_id]->next_remove();
                        }

                        pre_char = cPair_seq[block_id].first;
                        --cPair_seq[block_id].second;

                        is_LStar = m_sub_l_bit_seqs[block_id]->get();
                        m_sub_l_bit_seqs[block_id]->next_remove();

                        if(is_LStar) block_id = ( block_id | MASK_U16[2] );

                        pre_blockID = get_alpha_blockID(cur_alpha_block_id,pre_char);

                        if ( pre_blockID > cur_alpha_block_id )
                        {
                            if (g_rank != bks_ranks[pre_blockID]) block_id = (block_id | MASK_U16[1]);
                            bks_seqs[pre_blockID]->push_back(bks_pair_type(pre_char - m_alpha_blocks_info[pre_blockID].m_beg_alpha,block_id));
                            bks_ranks[pre_blockID] = g_rank;
                        }
                        else if (pre_blockID == cur_alpha_block_id)
                        {
                            t_sortSameBkt->push( ( (g_rank != t_rank) ? (block_id | MASK_U16[1]) : block_id )   );
                            t_rank = g_rank;
                        }
                        else
                        {
                            Logger::output_error(__FILE__,__LINE__);
                        }

                    }

                    /// next item
                    bks_seqs[cur_alpha_block_id]->next_remove();

                }

                if(bks_seqs[cur_alpha_block_id]->is_eof())
                {
                    delete bks_seqs[cur_alpha_block_id];
                    bks_seqs[cur_alpha_block_id] = nullptr;

                }
                else
                {
                    Logger::output_error(__FILE__,__LINE__);
                }

            }
            else
            {

                delete bks_seqs[cur_alpha_block_id];

                bks_seqs[cur_alpha_block_id] = nullptr;
            }

            /// step 2: sort the same bucket item

            if (t_sortSameBkt->write_size())
            {

                t_sortSameBkt->swapBuf();

                while (true)
                {

                    t_rank = uint64_MAX;

                    while (t_sortSameBkt->read_size())
                    {
                        block_id = t_sortSameBkt->pop();

                        if (block_id & MASK_U16[1]) ++g_rank;

                        is_LStar = block_id & MASK_U16[2];

                        block_id = block_id & MASK_U16[3];

                        /// if true, put this substring to 'lstar_sub_seq'
                        if(is_LStar)
                        {
                            if( g_rank != lstar_rank ) lstar_sub_blockId_seq->push_back(block_id | MASK_U16[1]);
                            else lstar_sub_blockId_seq->push_back(block_id);

                            if(cur_char != pre_lstar)
                            {

                                lstar_sub_char_seq->push_back(compressed_pair_type(pre_lstar,lstar_cnt));

                                pre_lstar = cur_char, lstar_cnt = 1;

                            }
                            else
                            {

                                if(lstar_cnt != c_max_count)
                                {
                                    lstar_cnt++;
                                }
                                else
                                {
                                    lstar_sub_char_seq->push_back(compressed_pair_type(cur_char,lstar_cnt));
                                    lstar_cnt = 1;
                                }

                            }

                            lstar_rank = g_rank;
                        }
                        else   /// compute BWT
                        {
                            /// check whether to get next pair;
                            if( !cPair_seq[block_id].second )
                            {
                                cPair_seq[block_id] = m_sub_l_cbwt_seqs[block_id]->get();
                                m_sub_l_cbwt_seqs[block_id]->next_remove();
                            }

                            pre_char = cPair_seq[block_id].first;
                            --cPair_seq[block_id].second;

                            is_LStar = m_sub_l_bit_seqs[block_id]->get();
                            m_sub_l_bit_seqs[block_id]->next_remove();

                            if(is_LStar) block_id = ( block_id | MASK_U16[2] );

                            pre_blockID = get_alpha_blockID(cur_alpha_block_id,pre_char);

                            if ( pre_blockID > cur_alpha_block_id )
                            {
                                if (g_rank != bks_ranks[pre_blockID]) block_id = (block_id | MASK_U16[1]);
                                bks_seqs[pre_blockID]->push_back(bks_pair_type(pre_char - m_alpha_blocks_info[pre_blockID].m_beg_alpha,block_id));
                                bks_ranks[pre_blockID] = g_rank;
                            }
                            else if (pre_blockID == cur_alpha_block_id)
                            {
                                t_sortSameBkt->push( ( (g_rank != t_rank) ? (block_id | MASK_U16[1]) : block_id )   );
                                t_rank = g_rank;
                            }
                            else
                            {
                                Logger::output_error(__FILE__,__LINE__);
                            }

                        }

                    }

                    if (t_sortSameBkt->write_size())
                        t_sortSameBkt->swapBuf();
                    else
                    {
                        std::cout << "break.\n";
                        break;
                    }

                }
            }

            delete t_sortSameBkt;

            /// step 3: sort the LMS item

            ++g_rank;

            while(!cLMS_sorter->is_empty() && cLMS_sorter->get_bkt() == cur_char)
            {

                block_id = cLMS_sorter->get_block_id();

                /// check whether to get next pair;
                if( !cPair_seq[block_id].second )
                {
                    cPair_seq[block_id] = m_sub_l_cbwt_seqs[block_id]->get();
                    m_sub_l_cbwt_seqs[block_id]->next_remove();
                }

                pre_char = cPair_seq[block_id].first;
                --cPair_seq[block_id].second;

                is_LStar = m_sub_l_bit_seqs[block_id]->get();
                m_sub_l_bit_seqs[block_id]->next_remove();

                if(is_LStar) block_id = ( block_id | MASK_U16[2] );

                pre_blockID = get_alpha_blockID(cur_alpha_block_id,pre_char);

                if ( pre_blockID > cur_alpha_block_id )
                {
                    if (g_rank != bks_ranks[pre_blockID]) block_id = (block_id | MASK_U16[1]);
                    bks_seqs[pre_blockID]->push_back(bks_pair_type(pre_char - m_alpha_blocks_info[pre_blockID].m_beg_alpha,block_id));
                    bks_ranks[pre_blockID] = g_rank;

                }
                else
                {

                    Logger::output_error(__FILE__,__LINE__);
                }

                cLMS_sorter->next();

            }/// ending of step_3



        }
        else
        {
            Logger::output_error(__FILE__,__LINE__);
        }
        /// next bucket
        ++cur_alpha_block_id;
    }

    /// push the last item of lstar_sub
    lstar_sub_char_seq->push_back(compressed_pair_type(pre_lstar,lstar_cnt)); /// \note error-prone


    /// destructor
    {
        while (!m_sub_l_cbwt_seqs.empty())
        {
            if (!m_sub_l_cbwt_seqs.back()->is_eof())
            {
                Logger::output_error(__FILE__,__LINE__);
            }
            else
            {
                delete m_sub_l_cbwt_seqs.back();
            }

            m_sub_l_cbwt_seqs.pop_back();
        }

        while (!m_sub_l_bit_seqs.empty())
        {
            if (!m_sub_l_bit_seqs.back()->is_eof())
            {
                Logger::output_error(__FILE__,__LINE__);
            }
            else
            {
                delete m_sub_l_bit_seqs.back();
            }

            m_sub_l_bit_seqs.pop_back();
        }

        delete cLMS_sorter;

    }

    std::cout << "------------------------------------------------------------------------------------------------------------------------------.\n";

    std::cout << "L-type substrings sorting on level_" << m_level << " is over.\n";

    std::cout << "The current PDU = " << Logger::max_pdu << ".\n";

    std::cout << "------------------------------------------------------------------------------------------------------------------------------.\n";

    std::cout << "Start to sort S-type substrings on level " << m_level << ".\n";

    Logger::output_ram_use();

    // Logger::pause(__FILE__,__LINE__);


    /// define the rank of the last S* substring
    uint64  sstar_rank = 0;

    g_rank = uint64_MAX; /// from biggest to smallest

    bool is_SStar = false;

    /// redefine bucket-section type
    bks_seqs.clear();
    bks_ranks.clear();

    for(uint64 i = 0; i < m_alpha_blocks_number; ++i)
    {
        bks_seqs.push_back(new bks_pair_type_vector());
        bks_seqs.back()->start_write();

        bks_ranks.push_back(0);
    }

    /// allocate a temporary value for each compressed BWT pair
    cPair_seq.clear();;

    for(uint64 i = 0; i < m_blocks_info.size()-1; i++) /// exclude the leftmost block
    {

        m_sub_s_cbwt_seqs[i]->start_read();

        /// assign the value to each cPair
        if(!m_sub_s_cbwt_seqs[i]->is_empty())
        {
            cPair_seq.push_back( m_sub_s_cbwt_seqs[i]->get() );
            m_sub_s_cbwt_seqs[i]->next_remove();
        }
        else
        {
            std::cout << "m_sub_s_cbwt_seqs[" << i << "]->is_empty() is true. Line: " << __LINE__ << std::endl;
        }

        m_sub_s_bit_seqs[i]->start_read_remove();

    }

    /// start to read the positions of LMS and different bit sequence
    for(uint64 i = 0; i < m_blocks_info.size() - 1; i++) /// exclude the leftmost block
    {
        m_LMS_rPos_seqs[i]->start_read();
        m_LMS_rPos_repeat_bit_seqs[i]->start_read_remove();
    }

    /// initialize the s1_blocks
    for(uint64 i = 0; i < m_blocks_info.size() - 1; ++i)
    {
        s1_blocks.push_back( new s1_blocks_vector_type());
        s1_blocks.back()->start_write();
    }


    /// alphabet block vector
    std::vector< Alpha_Block<offset_type> > alpha_block_vec;

    /// the current alphabet block size and current alphabet size
    offset_type total_alpha_block_size = 0, cur_alpha_size = 0, pre_name(0), cur_name(1);

    /// < ch, rank, block_id> = < 4B, sizeof(relative_offset_type), 2B >
    /// bug, the std::stable_sort at inducing stage use double working space.
    offset_type alpha_block_max_size = (MAX_MEM / 4) / (4 + sizeof(relative_offset_type) + 2) / 2; /// \note



    /// start to read the L*-substrings reversely
    lstar_sub_blockId_seq->start_read_reverse();
    lstar_sub_char_seq->start_read_reverse(); /// \note the first L* substring of pair(alphabet_max, 0)is useless, it should be omitted.


    /// < alphabet, s_rank, blockID >
    typedef Triple<relative_offset_type, relative_offset_type, uint16> pa_item_type;

    /// < r_alphabet, s_rank, count, blockID>
    typedef quadruple< relative_offset_type, relative_offset_type, relative_offset_type, uint16> pq_item_type;

    typedef std::priority_queue< pq_item_type, std::vector<pq_item_type>, TupleAscCmp3< pq_item_type > > max_pq_type;

    cur_alpha_block_id = m_alpha_blocks_number - 1;

    std::cout << "start to sort the S-type substrings.\n";

    while(cur_alpha_block_id >= 0)
    {

        offset_type cur_alpha_block_beg(m_alpha_blocks_info[cur_alpha_block_id].m_beg_alpha);

        offset_type cur_alpha_block_end(m_alpha_blocks_info[cur_alpha_block_id].m_end_alpha);

        relative_offset_type alpha_number = cur_alpha_block_end - cur_alpha_block_beg + 1;

        if(alpha_number > 1)
        {

            std::vector< pa_item_type > pa_item_ary;

            typename std::vector< pa_item_type >::iterator it_beg, it_end;

            if(!bks_seqs[cur_alpha_block_id]->is_empty())
            {

                pa_item_ary.resize( bks_seqs[cur_alpha_block_id]->size() );

                bks_seqs[cur_alpha_block_id]->start_read();

                relative_offset_type index(0), s_rank(std::numeric_limits<relative_offset_type>::max());

                while( !bks_seqs[cur_alpha_block_id]->is_eof() )
                {

                    if( bks_seqs[cur_alpha_block_id]->get().second & MASK_U16[1] ) --s_rank;

                    pa_item_ary[index] = pa_item_type(bks_seqs[cur_alpha_block_id]->get().first, s_rank, bks_seqs[cur_alpha_block_id]->get().second);

                    bks_seqs[cur_alpha_block_id]->next_remove();

                    ++index;
                }

                delete bks_seqs[cur_alpha_block_id];

                bks_seqs[cur_alpha_block_id] = nullptr;

                std::stable_sort(pa_item_ary.begin(), pa_item_ary.end(), TupleDscCmp1<pa_item_type>());

            }

            it_beg = pa_item_ary.begin(), it_end = pa_item_ary.end();

            relative_offset_type cur_rel_bkt(alpha_number - 1); /// relative alpha

            relative_offset_type t_rel_rank(std::numeric_limits<relative_offset_type>::max());

            relative_offset_type t_rank(0);

            relative_offset_type pq_input_cnt(std::numeric_limits<relative_offset_type>::max());

            max_pq_type * m_pq = new max_pq_type();

            alphabet_type cur_char(0), pre_char(0);

            uint16 pre_blockID(0);

            while(cur_rel_bkt >= 0)
            {
                cur_char = cur_alpha_block_beg + cur_rel_bkt;

                t_rank = 0;

                while ( it_beg != it_end && (*it_beg).first == cur_rel_bkt )
                {
                    if( (*it_beg).second != t_rank )
                    {

                        --t_rel_rank;

                        --g_rank;

                        t_rank = (*it_beg).second;

                    }

                    is_SStar = (*it_beg).third & MASK_U16[2]; /// the left second bit

                    /// the left two bit
                    block_id = (*it_beg).third & MASK_U16[3]; /// MASK_U16[3] = 0x3fff

                    /// 1.4: if true, put this substring to 'lstar_sub_seq'
                    if(is_SStar)
                    {
                        if(sstar_rank != g_rank) ++name;

                        /// the first bit item must be true
                        if( !m_LMS_rPos_repeat_bit_seqs[block_id]->get() )
                        {
                            Logger::output_error(__FILE__,__LINE__);
                        }

                        m_LMS_rPos_repeat_bit_seqs[block_id]->next_remove();

                        /// push the first LMS position
                        s1_blocks[block_id]->push_back(s1_pair_type(m_LMS_rPos_seqs[block_id]->get(), name));
                        m_LMS_rPos_seqs[block_id]->next_remove();

                        compute_alphabet_blocks(
                            alpha_block_vec,
                            total_alpha_block_size,
                            alpha_block_max_size,
                            cur_alpha_size,
                            pre_name,
                            name,
                            ((sstar_rank != g_rank)? true : false)
                        );

                        /// push the position of the same LMS substrings
                        while(!m_LMS_rPos_repeat_bit_seqs[block_id]->is_eof())
                        {

                            if(m_LMS_rPos_repeat_bit_seqs[block_id]->get())break;
                            else
                            {
                                s1_blocks[block_id]->push_back(s1_pair_type(m_LMS_rPos_seqs[block_id]->get(), name));
                                m_LMS_rPos_seqs[block_id]->next_remove();

                                compute_alphabet_blocks(
                                    alpha_block_vec,
                                    total_alpha_block_size,
                                    alpha_block_max_size,
                                    cur_alpha_size,
                                    pre_name,
                                    name,
                                    false
                                );

                            }
                            m_LMS_rPos_repeat_bit_seqs[block_id]->next_remove();
                        }

                        sstar_rank = g_rank;

                    }
                    else   /// compute BWT
                    {

                        /// check whether to get next pair;
                        if( !cPair_seq[block_id].second )
                        {
                            cPair_seq[block_id] = m_sub_s_cbwt_seqs[block_id]->get();
                            m_sub_s_cbwt_seqs[block_id]->next_remove();
                        }

                        pre_char = cPair_seq[block_id].first;
                        --cPair_seq[block_id].second;

                        is_SStar = m_sub_s_bit_seqs[block_id]->get();
                        m_sub_s_bit_seqs[block_id]->next_remove();

                        if(is_SStar) block_id = ( block_id | MASK_U16[2] );

                        pre_blockID = get_induce_alpha_blockID(cur_alpha_block_id,pre_char);

                        if ( pre_blockID < cur_alpha_block_id )
                        {
                            if (g_rank != bks_ranks[pre_blockID]) block_id = (block_id | MASK_U16[1]);
                            bks_seqs[pre_blockID]->push_back(bks_pair_type(pre_char - m_alpha_blocks_info[pre_blockID].m_beg_alpha,block_id));
                            bks_ranks[pre_blockID] = g_rank;
                        }
                        else if (pre_blockID == cur_alpha_block_id)
                        {
                            m_pq->push( pq_item_type(relative_offset_type(pre_char - cur_alpha_block_beg), t_rel_rank, --pq_input_cnt, block_id) );
                        }
                        else
                        {
                            std::cout << "block_id = " << (block_id & MASK_U16[3]) << std::endl;
                            std::cout << "pre_blockID = " << pre_blockID << std::endl;
                            std::cout << "pre_char = " << pre_char << std::endl;
                            std::cout << "cur_char = " << cur_char << std::endl;
                            Logger::output_error(__FILE__,__LINE__);
                        }

                    }

                    it_beg++;


                }

                --g_rank;

                t_rank = 0;

                --t_rel_rank;

                /// step2: sort PQ
                while(!m_pq->empty() && m_pq->top().first == cur_rel_bkt)
                {

                    if (m_pq->top().second != t_rank)
                    {
                        --g_rank;

                        --t_rel_rank;

                        t_rank = m_pq->top().second;
                    }

                    is_SStar = m_pq->top().forth & MASK_U16[2]; /// the left second bit

                    /// the left two bit
                    block_id = m_pq->top().forth & MASK_U16[3]; /// MASK_U16[3] = 0x3fff

                    /// 1.4: if true, put this substring to 'lstar_sub_seq'
                    if(is_SStar)
                    {
                        if(sstar_rank != g_rank) ++name;

                        /// the first bit item must be true
                        if( !m_LMS_rPos_repeat_bit_seqs[block_id]->get() )
                        {
                            Logger::output_error(__FILE__,__LINE__);
                        }

                        m_LMS_rPos_repeat_bit_seqs[block_id]->next_remove();

                        /// push the first LMS position
                        s1_blocks[block_id]->push_back(s1_pair_type(m_LMS_rPos_seqs[block_id]->get(), name));
                        m_LMS_rPos_seqs[block_id]->next_remove();

                        compute_alphabet_blocks(
                            alpha_block_vec,
                            total_alpha_block_size,
                            alpha_block_max_size,
                            cur_alpha_size,
                            pre_name,
                            name,
                            ((sstar_rank != g_rank)? true : false)
                        );

                        /// push the position of the same LMS substrings
                        while(!m_LMS_rPos_repeat_bit_seqs[block_id]->is_eof())
                        {

                            if(m_LMS_rPos_repeat_bit_seqs[block_id]->get())break;
                            else
                            {
                                s1_blocks[block_id]->push_back(s1_pair_type(m_LMS_rPos_seqs[block_id]->get(), name));
                                m_LMS_rPos_seqs[block_id]->next_remove();

                                compute_alphabet_blocks(
                                    alpha_block_vec,
                                    total_alpha_block_size,
                                    alpha_block_max_size,
                                    cur_alpha_size,
                                    pre_name,
                                    name,
                                    false
                                );

                            }
                            m_LMS_rPos_repeat_bit_seqs[block_id]->next_remove();
                        }

                        sstar_rank = g_rank;

                    }
                    else   /// compute BWT
                    {

                        /// check whether to get next pair;
                        if( !cPair_seq[block_id].second )
                        {
                            cPair_seq[block_id] = m_sub_s_cbwt_seqs[block_id]->get();
                            m_sub_s_cbwt_seqs[block_id]->next_remove();
                        }

                        pre_char = cPair_seq[block_id].first;
                        --cPair_seq[block_id].second;

                        is_SStar = m_sub_s_bit_seqs[block_id]->get();
                        m_sub_s_bit_seqs[block_id]->next_remove();

                        if(is_SStar) block_id = ( block_id | MASK_U16[2] );

                        pre_blockID = get_induce_alpha_blockID(cur_alpha_block_id,pre_char);

                        if ( pre_blockID < cur_alpha_block_id )
                        {
                            if (g_rank != bks_ranks[pre_blockID]) block_id = (block_id | MASK_U16[1]);
                            bks_seqs[pre_blockID]->push_back(bks_pair_type(pre_char - m_alpha_blocks_info[pre_blockID].m_beg_alpha,block_id));
                            bks_ranks[pre_blockID] = g_rank;
                        }
                        else if (pre_blockID == cur_alpha_block_id)
                        {
                            m_pq->push( pq_item_type(relative_offset_type(pre_char - cur_alpha_block_beg), t_rel_rank, --pq_input_cnt, block_id) );
                        }
                        else
                        {
                            Logger::output_error(__FILE__,__LINE__);
                        }

                    }

                    m_pq->pop();
                }


                /// step3: induce sort the LMS-bucket
                --g_rank;

                --t_rel_rank;

                while(!lstar_sub_char_seq->is_eof() && lstar_sub_char_seq->get_reverse().first == cur_char)
                {

                    //compressed_pair_type t_cPair(lstar_sub_char_seq->get_reverse());
                    compress_type counter(lstar_sub_char_seq->get_reverse().second);

                    while(counter)
                    {

                        block_id = lstar_sub_blockId_seq->get_reverse() & MASK_U16[3];


                        /// check whether to get next pair;
                        if( !cPair_seq[block_id].second )
                        {
                            cPair_seq[block_id] = m_sub_s_cbwt_seqs[block_id]->get();
                            m_sub_s_cbwt_seqs[block_id]->next_remove();
                        }

                        pre_char = cPair_seq[block_id].first;
                        --cPair_seq[block_id].second;

                        is_SStar = m_sub_s_bit_seqs[block_id]->get();
                        m_sub_s_bit_seqs[block_id]->next_remove();

                        if(is_SStar) block_id = ( block_id | MASK_U16[2] );

                        pre_blockID = get_induce_alpha_blockID(cur_alpha_block_id,pre_char);

                        if ( pre_blockID < cur_alpha_block_id )
                        {
                            if (g_rank != bks_ranks[pre_blockID]) block_id = (block_id | MASK_U16[1]);
                            bks_seqs[pre_blockID]->push_back(bks_pair_type(pre_char - m_alpha_blocks_info[pre_blockID].m_beg_alpha, block_id));
                            bks_ranks[pre_blockID] = g_rank;
                        }
                        else if (pre_blockID == cur_alpha_block_id)
                        {
                            m_pq->push( pq_item_type(relative_offset_type(pre_char - cur_alpha_block_beg), t_rel_rank, --pq_input_cnt, block_id) );
                        }
                        else
                        {
                            std::cout << "block_id = " << (block_id & MASK_U16[3]) << std::endl;
                            std::cout << "pre_blockID = " << pre_blockID << std::endl;
                            std::cout << "pre_char = " << pre_char << std::endl;
                            std::cout << "cur_char = " << cur_char << std::endl;
                            Logger::output_error(__FILE__,__LINE__);
                        }


                        if(lstar_sub_blockId_seq->get_reverse() & MASK_U16[1]) --g_rank, --t_rel_rank;

                        lstar_sub_blockId_seq->next_remove_reverse();

                        --counter;
                    }

                    lstar_sub_char_seq->next_remove_reverse();

                }

                /// next r_bucket
                if(cur_rel_bkt != 0)--cur_rel_bkt;
                else break;
            }

            delete m_pq;
        }
        else if(alpha_number == 1)
        {

            alphabet_type cur_char(cur_alpha_block_beg), pre_char(0);

            uint16 pre_blockID(0);

            InduceSortSameBkt< uint16 > * t_sortSameBkt = new InduceSortSameBkt< uint16 >(K_1024 * 1024 * 3);
            uint64 t_rank = 0;


            if(!bks_seqs[cur_alpha_block_id]->is_empty())
            {

                bks_seqs[cur_alpha_block_id]->start_read();

                while( !bks_seqs[cur_alpha_block_id]->is_eof())
                {

                    if(bks_seqs[cur_alpha_block_id]->get().second & MASK_U16[1] ) --g_rank;

                    is_SStar = bks_seqs[cur_alpha_block_id]->get().second & MASK_U16[2];

                    block_id = bks_seqs[cur_alpha_block_id]->get().second & MASK_U16[3];

                    ///if true, put this substring to 'lstar_sub_seq'
                    if(is_SStar)
                    {
                        if(sstar_rank != g_rank) ++name;

                        /// the first bit item must be true
                        if( !m_LMS_rPos_repeat_bit_seqs[block_id]->get() )
                        {
                            Logger::output_error(__FILE__,__LINE__);
                        }

                        m_LMS_rPos_repeat_bit_seqs[block_id]->next_remove();

                        /// push the first LMS position
                        s1_blocks[block_id]->push_back(s1_pair_type(m_LMS_rPos_seqs[block_id]->get(), name));
                        m_LMS_rPos_seqs[block_id]->next_remove();

                        compute_alphabet_blocks(
                            alpha_block_vec,
                            total_alpha_block_size,
                            alpha_block_max_size,
                            cur_alpha_size,
                            pre_name,
                            name,
                            ((sstar_rank != g_rank)? true : false)
                        );

                        /// push the position of the same LMS substrings
                        while(!m_LMS_rPos_repeat_bit_seqs[block_id]->is_eof())
                        {

                            if(m_LMS_rPos_repeat_bit_seqs[block_id]->get())break;
                            else
                            {
                                s1_blocks[block_id]->push_back(s1_pair_type(m_LMS_rPos_seqs[block_id]->get(), name));
                                m_LMS_rPos_seqs[block_id]->next_remove();

                                compute_alphabet_blocks(
                                    alpha_block_vec,
                                    total_alpha_block_size,
                                    alpha_block_max_size,
                                    cur_alpha_size,
                                    pre_name,
                                    name,
                                    false
                                );

                            }
                            m_LMS_rPos_repeat_bit_seqs[block_id]->next_remove();
                        }

                        sstar_rank = g_rank;

                    }
                    else   /// compute BWT
                    {

                        /// check whether to get next pair;
                        if( !cPair_seq[block_id].second )
                        {
                            cPair_seq[block_id] = m_sub_s_cbwt_seqs[block_id]->get();
                            m_sub_s_cbwt_seqs[block_id]->next_remove();
                        }

                        pre_char = cPair_seq[block_id].first;
                        --cPair_seq[block_id].second;

                        is_SStar = m_sub_s_bit_seqs[block_id]->get();
                        m_sub_s_bit_seqs[block_id]->next_remove();

                        if(is_SStar) block_id = ( block_id | MASK_U16[2] );

                        pre_blockID = get_induce_alpha_blockID(cur_alpha_block_id,pre_char);

                        if ( pre_blockID < cur_alpha_block_id )
                        {
                            if (g_rank != bks_ranks[pre_blockID]) block_id = (block_id | MASK_U16[1]);
                            bks_seqs[pre_blockID]->push_back(bks_pair_type(pre_char - m_alpha_blocks_info[pre_blockID].m_beg_alpha,block_id));
                            bks_ranks[pre_blockID] = g_rank;
                        }
                        else if (pre_blockID == cur_alpha_block_id)
                        {
                            t_sortSameBkt->push( ( (g_rank != t_rank) ? (block_id | MASK_U16[1]) : block_id )   );
                            t_rank = g_rank;
                        }
                        else
                        {
                            std::cout << "block_id = " << (block_id & MASK_U16[3]) << std::endl;
                            std::cout << "pre_blockID = " << pre_blockID << std::endl;
                            std::cout << "pre_char = " << pre_char << std::endl;
                            std::cout << "cur_char = " << cur_char << std::endl;
                            Logger::output_error(__FILE__,__LINE__);
                        }

                    }


                    bks_seqs[cur_alpha_block_id]->next_remove();
                }

                delete bks_seqs[cur_alpha_block_id];

                bks_seqs[cur_alpha_block_id] = nullptr;



            }
            else
            {
                delete bks_seqs[cur_alpha_block_id];

                bks_seqs[cur_alpha_block_id] = nullptr;
            }

            if (t_sortSameBkt->write_size())
            {

                t_sortSameBkt->swapBuf();

                while (true)
                {

                    t_rank = uint64_MAX;

                    while (t_sortSameBkt->read_size())
                    {
                        block_id = t_sortSameBkt->pop();

                        if (block_id & MASK_U16[1]) --g_rank;

                        is_SStar = block_id & MASK_U16[2];

                        block_id = block_id & MASK_U16[3];

                        if(is_SStar)
                        {
                            if(sstar_rank != g_rank) ++name;

                            /// the first bit item must be true
                            if( !m_LMS_rPos_repeat_bit_seqs[block_id]->get() )
                            {
                                Logger::output_error(__FILE__,__LINE__);
                            }

                            m_LMS_rPos_repeat_bit_seqs[block_id]->next_remove();

                            /// push the first LMS position
                            s1_blocks[block_id]->push_back(s1_pair_type(m_LMS_rPos_seqs[block_id]->get(), name));
                            m_LMS_rPos_seqs[block_id]->next_remove();

                            compute_alphabet_blocks(
                                alpha_block_vec,
                                total_alpha_block_size,
                                alpha_block_max_size,
                                cur_alpha_size,
                                pre_name,
                                name,
                                ((sstar_rank != g_rank)? true : false)
                            );

                            /// push the position of the same LMS substrings
                            while(!m_LMS_rPos_repeat_bit_seqs[block_id]->is_eof())
                            {

                                if(m_LMS_rPos_repeat_bit_seqs[block_id]->get())break;
                                else
                                {
                                    s1_blocks[block_id]->push_back(s1_pair_type(m_LMS_rPos_seqs[block_id]->get(), name));
                                    m_LMS_rPos_seqs[block_id]->next_remove();

                                    compute_alphabet_blocks(
                                        alpha_block_vec,
                                        total_alpha_block_size,
                                        alpha_block_max_size,
                                        cur_alpha_size,
                                        pre_name,
                                        name,
                                        false
                                    );

                                }
                                m_LMS_rPos_repeat_bit_seqs[block_id]->next_remove();
                            }

                            sstar_rank = g_rank;

                        }
                        else   /// compute BWT
                        {

                            /// check whether to get next pair;
                            if( !cPair_seq[block_id].second )
                            {
                                cPair_seq[block_id] = m_sub_s_cbwt_seqs[block_id]->get();
                                m_sub_s_cbwt_seqs[block_id]->next_remove();
                            }

                            pre_char = cPair_seq[block_id].first;
                            --cPair_seq[block_id].second;

                            is_SStar = m_sub_s_bit_seqs[block_id]->get();
                            m_sub_s_bit_seqs[block_id]->next_remove();

                            if(is_SStar) block_id = ( block_id | MASK_U16[2] );

                            pre_blockID = get_induce_alpha_blockID(cur_alpha_block_id,pre_char);

                            if ( pre_blockID < cur_alpha_block_id )
                            {
                                if (g_rank != bks_ranks[pre_blockID]) block_id = (block_id | MASK_U16[1]);
                                bks_seqs[pre_blockID]->push_back(bks_pair_type(pre_char - m_alpha_blocks_info[pre_blockID].m_beg_alpha,block_id));
                                bks_ranks[pre_blockID] = g_rank;
                            }
                            else if (pre_blockID == cur_alpha_block_id)
                            {
                                t_sortSameBkt->push( ( (g_rank != t_rank) ? (block_id | MASK_U16[1]) : block_id )   );
                                t_rank = g_rank;
                            }
                            else
                            {
                                std::cout << "block_id = " << (block_id & MASK_U16[3]) << std::endl;
                                std::cout << "pre_blockID = " << pre_blockID << std::endl;
                                std::cout << "pre_char = " << pre_char << std::endl;
                                std::cout << "cur_char = " << cur_char << std::endl;
                                Logger::output_error(__FILE__,__LINE__);
                            }

                        }




                    }

                    if (t_sortSameBkt->write_size())
                        t_sortSameBkt->swapBuf();
                    else
                    {
                        std::cout << "break.\n";
                        break;
                    }

                }
            }

            delete t_sortSameBkt;

            while(!lstar_sub_char_seq->is_eof() && lstar_sub_char_seq->get_reverse().first == cur_char)
            {

                //compressed_pair_type t_cPair(lstar_sub_char_seq->get_reverse());
                compress_type counter(lstar_sub_char_seq->get_reverse().second);

                while(counter)
                {

                    block_id = lstar_sub_blockId_seq->get_reverse() & MASK_U16[3];


                    /// check whether to get next pair;
                    if( !cPair_seq[block_id].second )
                    {
                        cPair_seq[block_id] = m_sub_s_cbwt_seqs[block_id]->get();
                        m_sub_s_cbwt_seqs[block_id]->next_remove();
                    }

                    pre_char = cPair_seq[block_id].first;
                    --cPair_seq[block_id].second;

                    is_SStar = m_sub_s_bit_seqs[block_id]->get();
                    m_sub_s_bit_seqs[block_id]->next_remove();

                    if(is_SStar) block_id = ( block_id | MASK_U16[2] );

                    pre_blockID = get_induce_alpha_blockID(cur_alpha_block_id,pre_char);

                    if ( pre_blockID < cur_alpha_block_id )
                    {
                        if (g_rank != bks_ranks[pre_blockID]) block_id = (block_id | MASK_U16[1]);
                        bks_seqs[pre_blockID]->push_back(bks_pair_type(pre_char - m_alpha_blocks_info[pre_blockID].m_beg_alpha, block_id));
                        bks_ranks[pre_blockID] = g_rank;
                    }
                    else
                    {
                        std::cout << "block_id = " << (block_id & MASK_U16[3]) << std::endl;
                        std::cout << "pre_blockID = " << pre_blockID << std::endl;
                        std::cout << "pre_char = " << pre_char << std::endl;
                        std::cout << "cur_char = " << cur_char << std::endl;
                        Logger::output_error(__FILE__,__LINE__);
                    }


                    if(lstar_sub_blockId_seq->get_reverse() & MASK_U16[1]) --g_rank;

                    lstar_sub_blockId_seq->next_remove_reverse();

                    --counter;
                }

                lstar_sub_char_seq->next_remove_reverse();

            }



        }
        else
        {

            Logger::output_error(__FILE__,__LINE__);
        }


        if(cur_alpha_block_id == 0) break;
        else --cur_alpha_block_id;
    }

    std::cout << "------------------------------------------------------------------------------------------------------------------------------.\n";

    std::cout << "S-type substrings sorting on level_" << m_level << " is over.\n";

    std::cout << "The current PDU = " << Logger::max_pdu << ".\n";

    std::cout << "------------------------------------------------------------------------------------------------------------------------------.\n";

    /// destruct
    {
        while(bks_seqs.empty())
        {
            if(bks_seqs.back())
            {
                delete bks_seqs.back();
            }
            bks_seqs.pop_back();
        }

        if(!lstar_sub_char_seq->is_eof()) /// the leftmost item, i.e., the first item in lstar_sub_seq must be not accessed
        {
            std::cout << "Leftmost pair of lstar = < " << lstar_sub_char_seq->get_reverse().first << " , " << lstar_sub_char_seq->get_reverse().second - 0 << std::endl;

            lstar_sub_char_seq->next_remove_reverse();

            //  Logger::output_error(__FILE__,__LINE__);

            if(lstar_sub_char_seq->is_eof())
            {
                delete lstar_sub_char_seq;

                std::cout << "deletet the lstar_sub_char_seq.\n";
            }
            else
            {
                Logger::output_error(__FILE__, __LINE__);
            }
        }



        if(lstar_sub_blockId_seq->is_eof())
        {
            delete lstar_sub_blockId_seq;
        }
        else
        {
            Logger::output_error(__FILE__,__LINE__);

        }



        while(!m_LMS_rPos_repeat_bit_seqs.empty())
        {
            if(m_LMS_rPos_repeat_bit_seqs.back()->is_eof())delete m_LMS_rPos_repeat_bit_seqs.back();
            else
            {
                Logger::output_error(__FILE__,__LINE__);
            }
            m_LMS_rPos_repeat_bit_seqs.pop_back();
        }

        while(!m_LMS_rPos_seqs.empty())
        {
            if(m_LMS_rPos_seqs.back()->is_eof())delete m_LMS_rPos_seqs.back();
            else
            {
                Logger::output_error(__FILE__,__LINE__);
            }
            m_LMS_rPos_seqs.pop_back();
        }


        while(!m_sub_s_cbwt_seqs.empty())
        {
            if(m_sub_s_cbwt_seqs.back()->is_eof())delete m_sub_s_cbwt_seqs.back();
            else
            {
                Logger::output_error(__FILE__,__LINE__);
            }
            m_sub_s_cbwt_seqs.pop_back();
        }

        while(!m_sub_s_bit_seqs.empty())
        {
            if(m_sub_s_bit_seqs.back()->is_eof())delete m_sub_s_bit_seqs.back();
            else
            {
                Logger::output_error(__FILE__,__LINE__);
            }
            m_sub_s_bit_seqs.pop_back();
        }

    }

    /**< compute the tail elements */
    /// Notice that the size of last alphabet plus the total_alpha_block_size exceeds the alpha_block_max_size,
    cur_name = name;

    if(pre_name == cur_name)
    {
        alpha_block_vec.push_back(Alpha_Block<offset_type>(pre_name, cur_name, cur_alpha_size));
    }
    else
    {
        if(cur_alpha_size >= alpha_block_max_size)
        {
            alpha_block_vec.push_back(Alpha_Block<offset_type>(pre_name, cur_name - 1, total_alpha_block_size));
            alpha_block_vec.push_back(Alpha_Block<offset_type>(cur_name, cur_name, cur_alpha_size));
        }
        else
        {
            alpha_block_vec.push_back(Alpha_Block<offset_type>(pre_name, cur_name, cur_alpha_size + total_alpha_block_size));
        }
    }

    /*	check whether the beginning alphabet is less or equal to the ending alphabet in each bucket-section	*/
    while(true)
    {
        bool is_error(false);
        for(uint64 i = 0; i < alpha_block_vec.size(); ++i)
        {
            if(alpha_block_vec[i].m_beg_alpha > alpha_block_vec[i].m_end_alpha)
            {
                std::cout << "alpha_block_vec[" << i <<"] is error.\n";
                is_error = true;
            }
        }
        if(is_error)
        {
            std::cout << "Alphabet blocks information array is error.\n";

            Logger::output_error(__FILE__,__LINE__);
        }
        else
        {
            break;
        }
    }


    /**< write the LMS alphabet block infor to disk, filename is alpha_blocks_info_Level_$level.dat */
    FILE * fot = fopen(("alpha_blocks_info_" + std::to_string(m_level) + ".dat").c_str(), "wb");

    fwrite(alpha_block_vec.data(), sizeof(alpha_block_vec[0]), alpha_block_vec.size(), fot);

    fclose(fot);

    Timer::add_reduce_sort_time(Timer::get_wall_time() - beg_time);

    Logger::add_red_iv( Logger::cur_iv - red_iv);

    Logger::add_red_ov( Logger::cur_ov - red_ov);

    std::cout << "Current Function: sortSub_4L() over.\n ";

    Timer::show_time();


    return name;
}


/// \brief merge block-wise sorted results
///
template<typename alphabet_type, typename offset_type, typename compress_type, typename relative_offset_type>
void DSAComputation<alphabet_type, offset_type, compress_type, relative_offset_type>::mergeSortedSuffixGlobal()
{


    if(m_alpha <= alphabet_type(255))
    {
        /// sort suffixes For small alphabet ( <= 1B)
        sortSuf_4S();
    }
    else
    {
        /// sort suffixes For Large alphabet ( > 1B)
        sortSuf_4L();
    }


    /// append the sentinel
    if(m_level)sa1_reverse->push_back(offset_type(m_s_len - 1));


}


template<typename alphabet_type, typename offset_type, typename compress_type, typename relative_offset_type>
void DSAComputation<alphabet_type, offset_type, compress_type, relative_offset_type>::sortSuf_4S()
{
    std::cout << "Current Function: sortSuf_4S().\n";

    Timer::show_time();

    uint64 ind_iv = Logger::cur_iv;

    uint64 ind_ov = Logger::cur_ov;

    double beg_time = Timer::get_wall_time();

    alphabet_type cur_char(0),pre_char(0);

    offset_type cur_pos(0);

    uint16 block_id = 0;

    bool is_LStar = false;

    /// count the number of L-suffixes in each bucket
    std::vector<uint64> l_suf_char_vector(m_alpha + 1, 0);

    /// store the sorted L-suffixes positions
    MyVector<offset_type> * l_suf_pos_seq = new MyVector<offset_type>();
    l_suf_pos_seq->start_write();

    /// indicate the L*-suffixes
    bit_vector_type * l_suf_bit_seq = new bit_vector_type();
    l_suf_bit_seq->start_write();


    /// the item of bucket-section: position
    std::vector< MyVector<offset_type> * > bks_seqs(m_alpha + 1, nullptr); /// bks_seqs(m_alpha + 1, nullptr);

    /// the type of L-suffixes,  use 1(true) to denote the L*-suffixes
    std::vector< bit_vector_type * > bit_type_seqs(m_alpha + 1, nullptr); /// indicate L*/S* type

    std::cout << "Initialize BKS_seq and bit_type_seq to start_write() "<< std::endl;

    for(uint64 i = 0; i <= m_alpha; ++i)
    {
        bks_seqs[i] = new MyVector<offset_type>();
        bks_seqs[i]->start_write();

        bit_type_seqs[i] = new bit_vector_type();
        bit_type_seqs[i]->start_write();

    }

    std::cout << "Initialize t_sortSameBkst()" << std::endl;

    /// induce the item belonging to the same bucket
    typedef Pair<offset_type, bool> same_bkt_pair_type;

    InduceSortSameBkt< same_bkt_pair_type > * t_sortSameBkt = new InduceSortSameBkt< same_bkt_pair_type >(K_1024 * 1024 * 1);

    /// initialize cPair_seq
    typedef Pair<alphabet_type, compress_type> cPair_type;

    /// allocate a temporary value for each compressed BWT pair
    std::vector<cPair_type> bwt_cPair;

    std::cout << "Make L-type BWT and bit sequences to start_read." << std::endl;

    for(uint64 i = 0; i < m_blocks_info.size(); i++) /// include the leftmost block
    {

        m_suf_l_cbwt_seqs[i]->start_read();

        /// assign the value to each cPair
        if(!m_suf_l_cbwt_seqs[i]->is_empty())
        {
            bwt_cPair.push_back( m_suf_l_cbwt_seqs[i]->get() );
            m_suf_l_cbwt_seqs[i]->next_remove();
        }
        else
        {
            std::cout << "the current level_" << m_level << "m_suf_l_cbwt_seqs["<< i << "]->is_empty() is true.\n ";

            std::cout << "The file : " << __FILE__ << ", the line : " << __LINE__ << std::endl;
        }

        m_suf_l_bit_seqs[i]->start_read_remove();

    }


    std::cout << "Make S* relative position sequence to start_read()." << std::endl;

    for(uint64 i = 0; i < m_blocks_info.size(); i++) /// include the leftmost block
    {
        m_LMS_rPos_seqs[i]->start_read();
    }


    std::cout << "start to read s1_blockID_seq." << std::endl;

    s1_blockID_seq->start_read();



    std::cout << "Start to sort L-type suffixes on level " << m_level << ".\n";

    std::cout << "The current PDU = " << Logger::max_pdu << ".\n";

    /**************************************************************/

    /// \note get BWT

    MyVector<uint8> * eLBWT = nullptr;
    MyVector<uint8> * eSBWT = nullptr;

    if(!m_level && is_getBWT)
    {
        std::cout << "Line: " << __LINE__ << ", Initialize the L- and S-type BWTs." << std::endl;


        for(uint64 i = 0; i < m_blocks_info.size(); i++)
        {

            ((MyVector<uint8> *)LStarPreChar_seqs[i])->start_read();
            ((MyVector<uint8> *)SStarPreChar_seqs[i])->start_read();
        }

        std::cout << "m_blocks_info.size() =  " << m_blocks_info.size() << " at line " << __LINE__ << std::endl;

        eLBWT = new MyVector<uint8>(), eLBWT->start_write();
        eSBWT = new MyVector<uint8>(), eSBWT->start_write();

        std::cout << "init eLBWT and eSBWT to write at line " << __LINE__ << std::endl;
    }

    //记录第一个元素所属的快和在块中的位置

//    int  firstCharIndex(0);
//    alphabet_type firstChar(0);
//    bool firstCharType;

    /**************************************************************/

    /// sort from 0 to m_alpha
    while( cur_char <= m_alpha)
    {

        std::cout << "cur_char = " << cur_char - 0 << std::endl;
        std::cout.flush();

        t_sortSameBkt->reset();

        if (!bks_seqs[cur_char]->is_empty())
        {
            /// step 1. induce sort the L-bucket
            bks_seqs[cur_char]->start_read();
            bit_type_seqs[cur_char]->start_read_remove();

            while (!bks_seqs[cur_char]->is_eof())
            {
                /// 1.1 get cur_position
                cur_pos = bks_seqs[cur_char]->get();
                bks_seqs[cur_char]->next_remove();

                /// 1.2: determine is_LStar
                is_LStar = bit_type_seqs[cur_char]->get();
                bit_type_seqs[cur_char]->next_remove();

                /// 1.3 append the l-type suffix
                ++l_suf_char_vector[cur_char];
                l_suf_pos_seq->push_back(cur_pos);
                l_suf_bit_seq->push_back(is_LStar);


                if(!m_level && is_getBWT) /// \note getBWT
                {
                    if(cur_pos == 0)
                    {

                        std::cout << "Line: " << __LINE__ << ", the current position is zero." << std::endl;

                        uint8 * buf = new uint8[1];
                        UtilityFunctions::fileRead<uint8>(buf,m_s,m_s_len-1,1);
                        eLBWT->push_back(buf[0]);

                        std::cout << "The last char is " << char(buf[0]) << std::endl;

                        delete [] buf;

                    }

                    if(is_LStar)
                    {

                        block_id = getBlockId(cur_pos);

                        eLBWT->push_back( ((MyVector<uint8> *)LStarPreChar_seqs[block_id])->get() );
                        ((MyVector<uint8> *)LStarPreChar_seqs[block_id])->next_remove();
                    }
                }


                /// 1.4: if cur_pos != 0 and is_LStar == false, compute preceding suffix
                if(cur_pos && !is_LStar) /// \note the leftmost block
                {
                    ///  get block_id
                    block_id = getBlockId(cur_pos);

                    /// check whether to get next pair;
                    if( !bwt_cPair[block_id].second )
                    {
                        bwt_cPair[block_id] = m_suf_l_cbwt_seqs[block_id]->get();
                        m_suf_l_cbwt_seqs[block_id]->next_remove();
                    }

                    pre_char = bwt_cPair[block_id].first;
                    --bwt_cPair[block_id].second;

                    if(!m_level && is_getBWT)eLBWT->push_back(pre_char);/// \note getBWT

                    is_LStar = m_suf_l_bit_seqs[block_id]->get();
                    m_suf_l_bit_seqs[block_id]->next_remove();

                    if (pre_char > cur_char)
                    {
                        bks_seqs[pre_char]->push_back(cur_pos - 1);
                        bit_type_seqs[pre_char]->push_back( is_LStar );
                    }
                    else if (pre_char == cur_char)
                    {
                        t_sortSameBkt->push( same_bkt_pair_type(cur_pos - 1, is_LStar));
                    }
                    else
                    {
                        Logger::output_error(__FILE__,__LINE__);
                    }

                }

            }

            /// free disk space of bks_seqs
            if (bks_seqs[cur_char]->is_eof())
                delete bks_seqs[cur_char], bks_seqs[cur_char] = nullptr;
            else
                Logger::output_error(__FILE__,__LINE__);

            /// free disk space of bit_type_seqs
            if (bit_type_seqs[cur_char]->is_eof())
                delete bit_type_seqs[cur_char], bit_type_seqs[cur_char] = nullptr;
            else
                Logger::output_error(__FILE__,__LINE__);

            /// step2: sort the same bkt
            if (t_sortSameBkt->write_size())
            {
                t_sortSameBkt->swapBuf();

                same_bkt_pair_type t_pair;

                while (true)
                {
                    while (t_sortSameBkt->read_size())
                    {

                        t_pair = t_sortSameBkt->pop();

                        ++l_suf_char_vector[cur_char];
                        l_suf_pos_seq->push_back(t_pair.first);
                        l_suf_bit_seq->push_back(t_pair.second);

                        if(!m_level && is_getBWT) /// \note getBWT
                        {
                            if(t_pair.first == 0)
                            {
                                std::cout << "Line: " << __LINE__ << ", the current position is zero." << std::endl;

                                uint8 * buf = new uint8[1];
                                UtilityFunctions::fileRead<uint8>(buf,m_s,m_s_len-1,1);
                                eLBWT->push_back(buf[0]);

                                std::cout << "The last char is " << char(buf[0]) << std::endl;

                                std::cout << __LINE__ << std::endl;

                                delete [] buf;

                            }

                            if(t_pair.second)
                            {
                                block_id = getBlockId(t_pair.first);

                                eLBWT->push_back( ((MyVector<uint8> *)LStarPreChar_seqs[block_id])->get() );
                                ((MyVector<uint8> *)LStarPreChar_seqs[block_id])->next_remove();
                            }
                        }


                        if(t_pair.first && !t_pair.second)
                        {

                            ///  get block_id
                            block_id = getBlockId(t_pair.first);

                            /// check whether to get next pair;
                            if( !bwt_cPair[block_id].second )
                            {
                                bwt_cPair[block_id] = m_suf_l_cbwt_seqs[block_id]->get();
                                m_suf_l_cbwt_seqs[block_id]->next_remove();
                            }

                            pre_char = bwt_cPair[block_id].first;
                            --bwt_cPair[block_id].second;

                            if(!m_level && is_getBWT)eLBWT->push_back(pre_char);/// \note getBWT

                            is_LStar = m_suf_l_bit_seqs[block_id]->get();
                            m_suf_l_bit_seqs[block_id]->next_remove();

                            if (pre_char > cur_char)
                            {
                                bks_seqs[pre_char]->push_back(t_pair.first - 1);
                                bit_type_seqs[pre_char]->push_back( is_LStar );
                            }
                            else if (pre_char == cur_char)
                            {
                                t_sortSameBkt->push( same_bkt_pair_type(t_pair.first - 1, is_LStar));
                            }
                            else
                            {
                                Logger::output_error(__FILE__,__LINE__);
                            }

                        }

                    }

                    if (t_sortSameBkt->write_size())
                        t_sortSameBkt->swapBuf();
                    else
                        break;
                }
            }

        }
        else   /// if the current bks_seqs is empty, free it
        {
            delete bks_seqs[cur_char];
            bks_seqs[cur_char] = nullptr;

            delete bit_type_seqs[cur_char];
            bit_type_seqs[cur_char] = nullptr;
        }

        /// 3.induce sort the LMS-bucket

        std::vector<uint64>::iterator it_beg = m_LMS_counter[cur_char]->begin(), it_end = m_LMS_counter[cur_char]->end();

        //block_id = 0;

        while( it_beg != it_end)
        {

            uint64 cur_counter = *it_beg;

            while(cur_counter)
            {

                block_id = s1_blockID_seq->get();
                s1_blockID_seq->next_remove();

                cur_pos = m_LMS_rPos_seqs[block_id]->get() + m_blocks_info[block_id]->m_beg_pos;
                m_LMS_rPos_seqs[block_id]->next_remove();


                /// check whether to get next pair;
                if( !bwt_cPair[block_id].second )
                {
                    bwt_cPair[block_id] = m_suf_l_cbwt_seqs[block_id]->get();
                    m_suf_l_cbwt_seqs[block_id]->next_remove();
                }

                pre_char = bwt_cPair[block_id].first;
                bwt_cPair[block_id].second--;

                is_LStar = m_suf_l_bit_seqs[block_id]->get();
                m_suf_l_bit_seqs[block_id]->next_remove();

                if (pre_char > cur_char)
                {
                    bks_seqs[pre_char]->push_back(cur_pos - 1);
                    bit_type_seqs[pre_char]->push_back(is_LStar);
                }
                else
                {
                    Logger::output_error(__FILE__,__LINE__);
                }

                cur_counter--;
            }

            it_beg++;

        }

        /// next bucket
        if(cur_char == m_alpha)break;
        else ++cur_char;
    }


    /// destructor
    {
        while (!m_suf_l_cbwt_seqs.empty())
        {
            if (!m_suf_l_cbwt_seqs.back()->is_eof())
            {
                Logger::output_error(__FILE__,__LINE__);
            }
            else
            {
                delete m_suf_l_cbwt_seqs.back();
            }

            m_suf_l_cbwt_seqs.pop_back();
        }

        while (!m_suf_l_bit_seqs.empty())
        {
            if (!m_suf_l_bit_seqs.back()->is_eof())
            {
                Logger::output_error(__FILE__,__LINE__);
            }
            else
            {
                delete m_suf_l_bit_seqs.back();
            }

            m_suf_l_bit_seqs.pop_back();
        }

        if(s1_blockID_seq->is_eof()) delete s1_blockID_seq;
        else Logger::output_error(__FILE__,__LINE__);

        while(!m_LMS_rPos_seqs.empty())
        {
            if(m_LMS_rPos_seqs.back()->is_eof()) delete m_LMS_rPos_seqs.back(), m_LMS_rPos_seqs.pop_back();
            else Logger::output_error(__FILE__,__LINE__);
        }

    }

    std::cout << "------------------------------------------------------------------------------------------------------------------------------.\n";

    std::cout << "L-type suffixes sorting on level_" << m_level << " is over.\n";

    std::cout << "The current PDU = " << Logger::max_pdu << ".\n";

    std::cout << "------------------------------------------------------------------------------------------------------------------------------.\n";

    std::cout << "Start to sort S-type suffixes on level " << m_level << ".\n";

    std::cout.flush();

    if(m_level)
    {
        bool is_SStar = false;

        /// redefine bucket-section type
        bks_seqs.clear();
        bit_type_seqs.clear();

        std::cout << "reset the bks_seqs and bit_type_seqs" << std::endl;

        for(uint64 i = 0; i <= m_alpha; ++i)
        {
            bks_seqs.push_back(new MyVector<offset_type>());
            bks_seqs.back()->start_write();

            bit_type_seqs.push_back(new bit_vector_type());
            bit_type_seqs.back()->start_write();

        }

        /// allocate a temporary value for each compressed BWT pair
        std::cout << "Initialize S-type BWT and bit sequences." << std::endl;

        bwt_cPair.clear();;

        for(uint64 i = 0; i < m_blocks_info.size(); i++) /// exclude the leftmost block
        {

            m_suf_s_cbwt_seqs[i]->start_read();

            /// assign the value to each cPair
            if(!m_suf_s_cbwt_seqs[i]->is_empty())
            {
                bwt_cPair.push_back( m_suf_s_cbwt_seqs[i]->get() );
                m_suf_s_cbwt_seqs[i]->next_remove();
            }
            else
            {
                std::cout << "the current level_" << m_level << "m_suf_s_cbwt_seqs["<< i << "]->is_empty() is true.\n ";

                std::cout << "The file : " << __FILE__ << ", the line : " << __LINE__ << std::endl;

            }

            m_suf_s_bit_seqs[i]->start_read_remove();

        }

        std::cout << "Initialize the sa1_reverse to start_write()." << std::endl;

        sa1_reverse = new MyVector<offset_type>();
        sa1_reverse->start_write();

        std::cout << "Initialize l_suf_bit_seq and l_suf_pos_seq to start_read_reverse()." << std::endl;
        l_suf_bit_seq->start_read_remove_reverse();
        l_suf_pos_seq->start_read_reverse();


        cur_char = m_alpha;

        while(cur_char >= alphabet_type(0))
        {
            t_sortSameBkt->reset();

            if (!bks_seqs[cur_char]->is_empty())
            {
                bks_seqs[cur_char]->start_read();

                bit_type_seqs[cur_char]->start_read_remove();

                while (!bks_seqs[cur_char]->is_eof())
                {
                    /// 1.1 get cur_position
                    cur_pos = bks_seqs[cur_char]->get();
                    bks_seqs[cur_char]->next_remove();

                    /// 1.2: determine is_LStar
                    is_SStar = bit_type_seqs[cur_char]->get();
                    bit_type_seqs[cur_char]->next_remove();

                    /// 1.3: write the suffix
                    sa1_reverse->push_back(cur_pos);


                    if(!is_SStar)
                    {
                        ///  get block_id
                        block_id = getBlockId(cur_pos);

                        /// check whether to get next pair;
                        if( !bwt_cPair[block_id].second )
                        {
                            bwt_cPair[block_id] = m_suf_s_cbwt_seqs[block_id]->get();
                            m_suf_s_cbwt_seqs[block_id]->next_remove();
                        }

                        pre_char = bwt_cPair[block_id].first;
                        bwt_cPair[block_id].second--;

                        is_SStar = m_suf_s_bit_seqs[block_id]->get();
                        m_suf_s_bit_seqs[block_id]->next_remove();

                        if (pre_char < cur_char)
                        {
                            bks_seqs[pre_char]->push_back(cur_pos - 1);
                            bit_type_seqs[pre_char]->push_back( is_SStar );
                        }
                        else if (pre_char == cur_char)
                        {
                            t_sortSameBkt->push( same_bkt_pair_type(cur_pos - 1, is_SStar));
                        }
                        else
                        {
                            Logger::output_error(__FILE__,__LINE__);
                        }

                    }

                }

                /// free disk space of bks_seqs
                if (bks_seqs[cur_char]->is_eof())
                    delete bks_seqs[cur_char], bks_seqs[cur_char] = nullptr;
                else
                    Logger::output_error(__FILE__,__LINE__);

                /// free disk space of bit_type_seqs
                if (bit_type_seqs[cur_char]->is_eof())
                    delete bit_type_seqs[cur_char], bit_type_seqs[cur_char] = nullptr;
                else
                    Logger::output_error(__FILE__,__LINE__);



                /// step2: sort the same bkt
                if (t_sortSameBkt->write_size())
                {
                    t_sortSameBkt->swapBuf();

                    same_bkt_pair_type t_pair;

                    while (true)
                    {
                        while (t_sortSameBkt->read_size())
                        {

                            t_pair = t_sortSameBkt->pop();

                            sa1_reverse->push_back(t_pair.first);

                            if(!t_pair.second)
                            {

                                ///  get block_id
                                block_id = getBlockId(t_pair.first);

                                /// check whether to get next pair;
                                if( !bwt_cPair[block_id].second )
                                {
                                    bwt_cPair[block_id] = m_suf_s_cbwt_seqs[block_id]->get();
                                    m_suf_s_cbwt_seqs[block_id]->next_remove();
                                }

                                pre_char = bwt_cPair[block_id].first;
                                bwt_cPair[block_id].second--;

                                is_SStar = m_suf_s_bit_seqs[block_id]->get();
                                m_suf_s_bit_seqs[block_id]->next_remove();

                                if (pre_char < cur_char)
                                {
                                    bks_seqs[pre_char]->push_back(t_pair.first - 1);
                                    bit_type_seqs[pre_char]->push_back( is_SStar );
                                }
                                else if (pre_char == cur_char)
                                {
                                    t_sortSameBkt->push( same_bkt_pair_type(t_pair.first - 1, is_SStar));
                                }
                                else
                                {
                                    Logger::output_error(__FILE__,__LINE__);
                                }

                            }

                        }

                        if (t_sortSameBkt->write_size())
                            t_sortSameBkt->swapBuf();
                        else
                            break;
                    }
                }

            }
            else
            {
                delete bks_seqs[cur_char];
                bks_seqs[cur_char] = nullptr;

                delete bit_type_seqs[cur_char];
                bit_type_seqs[cur_char] = nullptr;
            }


            /// step3: Induce sort the LStar-type
            while (l_suf_char_vector[cur_char])
            {
                sa1_reverse->push_back(l_suf_pos_seq->get_reverse());

                if(l_suf_bit_seq->get_reverse())
                {
                    cur_pos = l_suf_pos_seq->get_reverse();

                    block_id = getBlockId(cur_pos);

                    /// check whether to get next pair;
                    if( !bwt_cPair[block_id].second )
                    {
                        bwt_cPair[block_id] = m_suf_s_cbwt_seqs[block_id]->get();
                        m_suf_s_cbwt_seqs[block_id]->next_remove();
                    }

                    pre_char = bwt_cPair[block_id].first;
                    bwt_cPair[block_id].second--;

                    is_SStar = m_suf_s_bit_seqs[block_id]->get();
                    m_suf_s_bit_seqs[block_id]->next_remove();

                    if (pre_char < cur_char)
                    {
                        bks_seqs[pre_char]->push_back(cur_pos - 1);
                        bit_type_seqs[pre_char]->push_back(is_SStar);
                    }
                    else
                    {
                        Logger::output_error(__FILE__,__LINE__);
                    }


                }

                l_suf_char_vector[cur_char]--;
                l_suf_pos_seq->next_remove_reverse();
                l_suf_bit_seq->next_remove_reverse();

            }

            if(!cur_char)
                break;
            else --cur_char;
        }

        std::cout << "------------------------------------------------------------------------------------------------------------------------------.\n";

        std::cout << "S-type suffixes sorting on level_" << m_level << " is over.\n";

        std::cout << "The current PDU = " << Logger::max_pdu << ".\n";

        std::cout << "------------------------------------------------------------------------------------------------------------------------------.\n";

        /// destruct
        {

            delete t_sortSameBkt;

            if(l_suf_bit_seq->is_eof())
                delete l_suf_bit_seq;
            else
                Logger::output_error(__FILE__,__LINE__);

            if(l_suf_pos_seq->is_eof())
                delete l_suf_pos_seq;
            else
                Logger::output_error(__FILE__,__LINE__);


            while(!m_suf_s_cbwt_seqs.empty())
            {
                if(m_suf_s_cbwt_seqs.back()->is_eof())
                    delete m_suf_s_cbwt_seqs.back();
                else
                {
                    Logger::output_error(__FILE__,__LINE__);
                }
                m_suf_s_cbwt_seqs.pop_back();
            }

            while(!m_suf_s_bit_seqs.empty())
            {
                if(m_suf_s_bit_seqs.back()->is_eof())
                    delete m_suf_s_bit_seqs.back();
                else
                {
                    Logger::output_error(__FILE__,__LINE__);
                }
                m_suf_s_bit_seqs.pop_back();
            }

        }

        Timer::show_time();

        Timer::add_induce_sort_time(Timer::get_wall_time() - beg_time);

        Logger::add_ind_iv( Logger::cur_iv - ind_iv);

        Logger::add_ind_ov(Logger::cur_ov - ind_ov);

        std::cerr << "finish inducing S-type suffixes.\n";

        std::cout << " the cur_level = " << m_level << std::endl;


    }
    else
    {
        /// count the number of L-suffixes in each bucket
        std::vector<uint64> s_suf_char_vector(m_alpha + 1, 0);

        std::vector<uint64> l_suf_char_vector_1(m_alpha + 1, 0);

        for(uint64 i = 0; i < l_suf_char_vector.size(); ++i) l_suf_char_vector_1[i] = l_suf_char_vector[i];

        /// store the sorted L-suffixes positions
        MyVector<offset_type> * s_suf_pos_seq = new MyVector<offset_type>();
        s_suf_pos_seq->start_write();

        bool is_SStar = false;

        /// redefine bucket-section type
        bks_seqs.clear();
        bit_type_seqs.clear();

        std::cout << "reset the bks_seqs and bit_type_seqs" << std::endl;

        for(uint64 i = 0; i <= m_alpha; ++i)
        {
            bks_seqs.push_back(new MyVector<offset_type>());
            bks_seqs.back()->start_write();

            bit_type_seqs.push_back(new bit_vector_type());
            bit_type_seqs.back()->start_write();

        }

        /// allocate a temporary value for each compressed BWT pair
        std::cout << "Initialize S-type BWT and bit sequences." << std::endl;

        bwt_cPair.clear();;

        for(uint64 i = 0; i < m_blocks_info.size(); i++) /// exclude the leftmost block
        {

            m_suf_s_cbwt_seqs[i]->start_read();

            /// assign the value to each cPair
            if(!m_suf_s_cbwt_seqs[i]->is_empty())
            {
                bwt_cPair.push_back( m_suf_s_cbwt_seqs[i]->get() );
                m_suf_s_cbwt_seqs[i]->next_remove();
            }
            else
            {
                std::cout << "the current level_" << m_level << "m_suf_s_cbwt_seqs["<< i << "]->is_empty() is true.\n ";

                std::cout << "The file : " << __FILE__ << ", the line : " << __LINE__ << std::endl;

            }

            m_suf_s_bit_seqs[i]->start_read_remove();

        }

        l_suf_bit_seq->start_read_remove_reverse();
        l_suf_pos_seq->start_read_reverse();

        cur_char = m_alpha;

        while(cur_char >= alphabet_type(0))
        {
            t_sortSameBkt->reset();

            if (!bks_seqs[cur_char]->is_empty())
            {
                bks_seqs[cur_char]->start_read();

                bit_type_seqs[cur_char]->start_read_remove();

                while (!bks_seqs[cur_char]->is_eof())
                {
                    /// 1.1 get cur_position
                    cur_pos = bks_seqs[cur_char]->get();
                    bks_seqs[cur_char]->next_remove();

                    /// 1.2: determine is_LStar
                    is_SStar = bit_type_seqs[cur_char]->get();
                    bit_type_seqs[cur_char]->next_remove();

                    /// 1.3: write the suffix
                    /*sa1_reverse->push_back(cur_pos);*/
                    s_suf_char_vector[cur_char]++;
                    s_suf_pos_seq->push_back(cur_pos);


                    if(!m_level && is_getBWT) /// \note getBWT
                    {
                        if(cur_pos == 0)
                        {
                            std::cout << "Line: " << __LINE__ << ", the current position is zero." << std::endl;

                            uint8 * buf = new uint8[1];
                            UtilityFunctions::fileRead<uint8>(buf,m_s,m_s_len-1,1);
                            eSBWT->push_back(buf[0]);

                            std::cout << "The last char is " << buf[0]-0 << std::endl;

                            /*if(is_SStar)
                            {
                                std::cout << __LINE__ << ", error occurs" << std::endl;
                                std::cin.get();
                            }*/

                            delete [] buf;

                        }
                        else
                        {

                            if(is_SStar)
                            {
                                block_id = getBlockId(cur_pos);

                                if(cur_pos == m_blocks_info[block_id]->m_end_pos) block_id--;

                                eSBWT->push_back( ((MyVector<uint8> *)SStarPreChar_seqs[ block_id ])->get() );
                                ((MyVector<uint8> *)SStarPreChar_seqs[block_id])->next_remove();
                            }

                        }
                    }

                    if(!is_SStar)
                    {
                        ///  get block_id
                        block_id = getBlockId(cur_pos);

                        /// check whether to get next pair;
                        if( !bwt_cPair[block_id].second )
                        {
                            bwt_cPair[block_id] = m_suf_s_cbwt_seqs[block_id]->get();
                            m_suf_s_cbwt_seqs[block_id]->next_remove();
                        }

                        pre_char = bwt_cPair[block_id].first;
                        bwt_cPair[block_id].second--;

                        if(!m_level && is_getBWT)eSBWT->push_back(pre_char);/// \note getBWT

                        is_SStar = m_suf_s_bit_seqs[block_id]->get();
                        m_suf_s_bit_seqs[block_id]->next_remove();

                        if (pre_char < cur_char)
                        {
                            bks_seqs[pre_char]->push_back(cur_pos - 1);
                            bit_type_seqs[pre_char]->push_back( is_SStar );
                        }
                        else if (pre_char == cur_char)
                        {
                            t_sortSameBkt->push( same_bkt_pair_type(cur_pos - 1, is_SStar));
                        }
                        else
                        {
                            Logger::output_error(__FILE__,__LINE__);
                        }

                    }

                }

                /// free disk space of bks_seqs
                if (bks_seqs[cur_char]->is_eof())
                    delete bks_seqs[cur_char], bks_seqs[cur_char] = nullptr;
                else
                    Logger::output_error(__FILE__,__LINE__);

                /// free disk space of bit_type_seqs
                if (bit_type_seqs[cur_char]->is_eof())
                    delete bit_type_seqs[cur_char], bit_type_seqs[cur_char] = nullptr;
                else
                    Logger::output_error(__FILE__,__LINE__);



                /// step2: sort the same bkt
                if (t_sortSameBkt->write_size())
                {
                    t_sortSameBkt->swapBuf();

                    same_bkt_pair_type t_pair;

                    while (true)
                    {
                        while (t_sortSameBkt->read_size())
                        {

                            t_pair = t_sortSameBkt->pop();

                            /*sa1_reverse->push_back(t_pair.first);*/
                            s_suf_char_vector[cur_char]++;
                            s_suf_pos_seq->push_back(t_pair.first);

                            if(!m_level && is_getBWT) /// \note getBWT
                            {
                                if(t_pair.first == 0)
                                {
                                    std::cout << "Line: " << __LINE__ << ", the current position is zero." << std::endl;

                                    uint8 * buf = new uint8[1];
                                    UtilityFunctions::fileRead<uint8>(buf,m_s,m_s_len-1,1);
                                    eSBWT->push_back(buf[0]);

                                    std::cout << "The last char is " << buf[0]-0 << std::endl;

                                    /*if(t_pair.second)
                                    {
                                        std::cout << __LINE__ << ", error occurs" << std::endl;
                                        std::cin.get();
                                    }*/

                                    delete [] buf;

                                }
                                else
                                {

                                    if(t_pair.second) // is_SStar == true
                                    {
                                        block_id = getBlockId(t_pair.first);

                                        if(t_pair.first == m_blocks_info[block_id]->m_end_pos) block_id--;

                                        eSBWT->push_back( ((MyVector<uint8> *)SStarPreChar_seqs[block_id])->get() );
                                        ((MyVector<uint8> *)SStarPreChar_seqs[block_id])->next_remove();
                                    }

                                }
                            }

                            if(!t_pair.second)
                            {

                                ///  get block_id
                                block_id = getBlockId(t_pair.first);

                                /// check whether to get next pair;
                                if( !bwt_cPair[block_id].second )
                                {
                                    bwt_cPair[block_id] = m_suf_s_cbwt_seqs[block_id]->get();
                                    m_suf_s_cbwt_seqs[block_id]->next_remove();
                                }

                                pre_char = bwt_cPair[block_id].first;
                                bwt_cPair[block_id].second--;

                                if(!m_level && is_getBWT)eSBWT->push_back(pre_char);/// \note getBWT

                                is_SStar = m_suf_s_bit_seqs[block_id]->get();
                                m_suf_s_bit_seqs[block_id]->next_remove();

                                if (pre_char < cur_char)
                                {
                                    bks_seqs[pre_char]->push_back(t_pair.first - 1);
                                    bit_type_seqs[pre_char]->push_back( is_SStar );
                                }
                                else if (pre_char == cur_char)
                                {
                                    t_sortSameBkt->push( same_bkt_pair_type(t_pair.first - 1, is_SStar));
                                }
                                else
                                {
                                    Logger::output_error(__FILE__,__LINE__);
                                }

                            }

                        }

                        if (t_sortSameBkt->write_size())
                            t_sortSameBkt->swapBuf();
                        else
                            break;
                    }
                }

            }
            else
            {
                delete bks_seqs[cur_char];
                bks_seqs[cur_char] = nullptr;

                delete bit_type_seqs[cur_char];
                bit_type_seqs[cur_char] = nullptr;
            }


            /// step3: Induce sort the LStar-type
            while (l_suf_char_vector_1[cur_char])
            {
                //sa1_reverse->push_back(l_suf_pos_seq->get_reverse());

                if(l_suf_bit_seq->get_reverse())
                {
                    cur_pos = l_suf_pos_seq->get_reverse();

                    block_id = getBlockId(cur_pos);

                    /// check whether to get next pair;
                    if( !bwt_cPair[block_id].second )
                    {
                        bwt_cPair[block_id] = m_suf_s_cbwt_seqs[block_id]->get();
                        m_suf_s_cbwt_seqs[block_id]->next_remove();
                    }

                    pre_char = bwt_cPair[block_id].first;
                    bwt_cPair[block_id].second--;

                    is_SStar = m_suf_s_bit_seqs[block_id]->get();
                    m_suf_s_bit_seqs[block_id]->next_remove();

                    if (pre_char < cur_char)
                    {
                        bks_seqs[pre_char]->push_back(cur_pos - 1);
                        bit_type_seqs[pre_char]->push_back(is_SStar);
                    }
                    else
                    {
                        Logger::output_error(__FILE__,__LINE__);
                    }


                }

                l_suf_char_vector_1[cur_char]--;
                /*l_suf_pos_seq->next_remove_reverse();*/
                /*l_suf_bit_seq->next_remove_reverse();*/

                l_suf_pos_seq->next_reverse();
                l_suf_bit_seq->next_remove_reverse();

            }

            if(!cur_char)
                break;
            else --cur_char;
        }

        std::cout << "------------------------------------------------------------------------------------------------------------------------------.\n";

        std::cout << "S-type suffixes sorting on level_" << m_level << " is over.\n";

        std::cout << "The current PDU = " << Logger::max_pdu << ".\n";

        std::cout << "------------------------------------------------------------------------------------------------------------------------------.\n";

        /// destruct
        {

            delete t_sortSameBkt;

            if(l_suf_bit_seq->is_eof())
                delete l_suf_bit_seq;
            else
                Logger::output_error(__FILE__,__LINE__);

            /*if(l_suf_pos_seq->is_eof())
                delete l_suf_pos_seq;
            else
                Logger::output_error(__FILE__,__LINE__);*/


            while(!m_suf_s_cbwt_seqs.empty())
            {
                if(m_suf_s_cbwt_seqs.back()->is_eof())
                    delete m_suf_s_cbwt_seqs.back();
                else
                {
                    Logger::output_error(__FILE__,__LINE__);
                }
                m_suf_s_cbwt_seqs.pop_back();
            }

            while(!m_suf_s_bit_seqs.empty())
            {
                if(m_suf_s_bit_seqs.back()->is_eof())
                    delete m_suf_s_bit_seqs.back();
                else
                {
                    Logger::output_error(__FILE__,__LINE__);
                }
                m_suf_s_bit_seqs.pop_back();
            }

        }

        std::cerr << "finish inducing S-type suffixes.\n";

        std::cout << " the cur_level = " << m_level << std::endl;

        Timer::show_time();

        Timer::add_induce_sort_time(Timer::get_wall_time() - beg_time);

        Logger::add_ind_iv( Logger::cur_iv - ind_iv);

        Logger::add_ind_ov(Logger::cur_ov - ind_ov);


        double beg_time = Timer::get_wall_time();

        uint64 wrt_sa_iv = Logger::cur_iv;

        uint64 wrt_sa_ov = Logger::cur_ov;

        l_suf_pos_seq->start_read();
        s_suf_pos_seq->start_read_reverse();


        std::string sa_name = m_s + ".sa5";

        FILE * fp;

        uint64 buf_size;

        uint40 * sa_buf;

        uint64 index = 0;


        if(is_getSA)
        {
            std::cout << "write the final SA in ascending order.\n";

            fp = fopen(sa_name.c_str(),"wb");

            buf_size = 20 * K_1024 * sizeof(uint40);

            sa_buf = new uint40[ buf_size ];

        }
        else
        {

            delete l_suf_pos_seq;
            delete s_suf_pos_seq;

            std::cout << "*********************Don't write the final SA to disk.*********************\n";
        }

        if(!is_getBWT)
        {
            std::cout << "*********************Don't write the final BWT to disk.*********************\n";
        }



//////////////////////////////////////////////////////////////////////////////
        std::string bwt_name = m_s + ".bwt";
        uint64 bwt_buf_size = 20 * K_1024;
        uint8 * bwt_buf;
        FILE * bwt_fp;
        uint32 bwt_index = 0;

        std::vector<uint64> l_bwt_number_vector(l_suf_char_vector);
        std::vector<uint64> s_bwt_number_vector(s_suf_char_vector);

        if(!m_level && is_getBWT)
        {
            bwt_buf = new uint8[bwt_buf_size];
            bwt_fp = fopen(bwt_name.c_str(),"wb");

            eLBWT->start_read();
            eSBWT->start_read_reverse();

            std::cout << "eLBWT->size()  = " << eLBWT->size() << std::endl;
            std::cout << "eSBWT->size()  = " << eSBWT->size() << std::endl;

            std::cout << "eLBWT + eSBWT  = " <<  eLBWT->size() + eSBWT->size() << std::endl;
        }
///////////////////////////////////////////////////////////////////////////////

        for(uint64 i = 0; i <= m_alpha; i++)
        {
            if(is_getSA)
            {

                while(l_suf_char_vector[i])
                {
                    sa_buf[index++] = l_suf_pos_seq->get();

                    l_suf_pos_seq->next_remove();

                    l_suf_char_vector[i]--;

                    if (index == buf_size)
                    {

                        double bg = Timer::get_wall_time();

                        fwrite(sa_buf,sizeof(uint40),index,fp);

                        Logger::addPDU(index * sizeof(uint40));

                        Logger::addOV(index * sizeof(uint40));

                        index = 0;
                    }

                }


                while(s_suf_char_vector[i])
                {
                    sa_buf[index++] = s_suf_pos_seq->get_reverse();

                    s_suf_pos_seq->next_remove_reverse();

                    s_suf_char_vector[i]--;

                    if (index == buf_size)
                    {

                        double bg = Timer::get_wall_time();

                        fwrite(sa_buf,sizeof(uint40),index,fp);

                        Logger::addPDU(index * sizeof(uint40));

                        Logger::addOV(index * sizeof(uint40));

                        index = 0;
                    }

                }

            }

///////////////////////////////////////////////////////////////////////////////

            if(!m_level && is_getBWT)
            {

                while(l_bwt_number_vector[i])
                {
                    bwt_buf[bwt_index++] = eLBWT->get();

                    eLBWT->next_remove();

                    l_bwt_number_vector[i]--;

                    if (bwt_index == bwt_buf_size)
                    {

                        fwrite(bwt_buf,sizeof(uint8),bwt_index,bwt_fp);

                        Logger::addPDU(bwt_index * sizeof(uint8));

                        Logger::addOV(bwt_index * sizeof(uint8));

                        bwt_index = 0;
                    }

                }

                while(s_bwt_number_vector[i])
                {
                    bwt_buf[bwt_index++] = eSBWT->get_reverse();

                    eSBWT->next_remove_reverse();

                    s_bwt_number_vector[i]--;

                    if (bwt_index == bwt_buf_size)
                    {

                        fwrite(bwt_buf,sizeof(uint8),bwt_index,bwt_fp);

                        Logger::addPDU(bwt_index * sizeof(uint8));

                        Logger::addOV(bwt_index * sizeof(uint8));

                        bwt_index = 0;
                    }

                }

            }

///////////////////////////////////////////////////////////////////////////////

        }

        if(is_getSA)
        {

            fwrite(sa_buf, sizeof(uint40), index, fp);

            Logger::addPDU(index * sizeof(uint40));

            Logger::addOV(index * sizeof(uint40));

            fclose(fp);

            delete [] sa_buf;

            if(l_suf_pos_seq->is_eof())delete l_suf_pos_seq;
            else Logger::output_error(__FILE__,__LINE__);

            if(s_suf_pos_seq->is_eof())delete s_suf_pos_seq;
            else Logger::output_error(__FILE__,__LINE__);

            std::cout << "The task of writing the final SA to disk is over.";

        }

        if(!m_level && is_getBWT)
        {
            fwrite(bwt_buf, sizeof(uint8), bwt_index, bwt_fp);

            Logger::addPDU(bwt_index * sizeof(uint8));

            Logger::addOV(bwt_index * sizeof(uint8));

            fclose(bwt_fp);

            delete [] bwt_buf;

            std::cout << "The task of writing the final BWT into disk is over.\n";

            std::cout << "Destruct eLBWT and eSBWT.\n";

            if(!eLBWT->is_eof())
            {
                std::cout << "eLBWT seqs is error."  << std::endl;

                std::cin.get();
            }

            if(!eSBWT->is_eof())
            {
                std::cout << "eSBWT seqs is error."  << std::endl;

                std::cout << "eSBWT->read_size() = " << eSBWT->read_size() << std::endl;
                std::cout << "eSBWT->size() = " << eSBWT->size() << std::endl;

                std::cin.get();
            }

            delete eLBWT, eSBWT;

            std::cout << "Destruct LStarPreChar_seqs and SStarPreChar_seqs." << std::endl;

            for(int i = 0; i < m_blocks_info.size(); ++i)
            {
                if(  ((MyVector<uint8> *)LStarPreChar_seqs[i])->is_empty() )
                    std::cout << "LStarPreChar_seqs[" << i << "] is empty." << std::endl;

                if(  !((MyVector<uint8> *)LStarPreChar_seqs[i])->is_eof() )
                    std::cout << "Error occurs in LStarPreChar_seqs[" << i << "]." << std::endl;

                delete ((MyVector<uint8> *)LStarPreChar_seqs[i]);

                if(  ((MyVector<uint8> *)SStarPreChar_seqs[i])->is_empty() )
                    std::cout << "SStarPreChar_seqs[" << i << "] is empty." << std::endl;

                if(  !((MyVector<uint8> *)SStarPreChar_seqs[i])->is_eof() )
                    std::cout << "Error occurs in SStarPreChar_seqs[" << i << "]." << std::endl;

                delete ((MyVector<uint8> *)SStarPreChar_seqs[i]);
            }
        }




        Timer::add_write_sa_time(Timer::get_wall_time() - beg_time);

        Logger::add_wrt_sa_iv(Logger::cur_iv - wrt_sa_iv);

        Logger::add_wrt_sa_ov(Logger::cur_ov - wrt_sa_ov);

        sa1_reverse = nullptr;

    }



    return ;
}

/// \brief sort S-type suffix 4 Large alphabet ( >= 1B )
template<typename alphabet_type, typename offset_type, typename compress_type, typename relative_offset_type>
void DSAComputation<alphabet_type, offset_type, compress_type, relative_offset_type>::sortSuf_4L()
{
    std::cout << "Current Function: sortSuf_4L().\n";

    Timer::show_time();

    uint64 ind_iv = Logger::cur_iv;

    uint64 ind_ov = Logger::cur_ov;

    double beg_time = Timer::get_wall_time();

    /// < r_alphabet, position, bool >
    typedef Triple<relative_offset_type, offset_type, bool> pa_triple_type;

    /// \brief < r_alphabet, position >
    typedef Pair<relative_offset_type, offset_type> bks_pair_type;

    /// < r_alphabet, count,  pos, bool>
    typedef quadruple< relative_offset_type, relative_offset_type, offset_type, bool > pq_quadruple_type;


    uint16 block_id;


    bool is_LStar = false;


    /// compress max count
    compress_type c_max_count(std::numeric_limits<compress_type>::max()), cChar_cnt(0); /// the number of first compressed L*-character

    alphabet_type pre_cChar(std::numeric_limits<alphabet_type>::max()); /// first compressed L*-character

    std::cout << "Initialize the L-type suffix sequence" << std::endl;

    /// compress L-type suffix sequence
    MyVector< compressed_pair_type > * l_suf_char_seq = new MyVector< compressed_pair_type >();
    l_suf_char_seq->start_write();

    MyVector< offset_type > * l_suf_pos_seq = new MyVector< offset_type >();
    l_suf_pos_seq->start_write();

    bit_vector_type * l_suf_bit_seq = new bit_vector_type();
    l_suf_bit_seq->start_write();


    std::cout << "Initialize the bkt_seqs and bit_type_seq." << std::endl;

    /// define bucket-section type
    std::vector< MyVector< bks_pair_type > * > bks_seqs;

    std::vector< bit_vector_type * > bit_type_seqs;

    std::cout << "m_alpha_blocks_number = " << m_alpha_blocks_number << std::endl;

    for(uint64 i = 0; i < m_alpha_blocks_number; ++i)
    {
        bks_seqs.push_back(new MyVector< bks_pair_type >());
        bks_seqs.back()->start_write();

        bit_type_seqs.push_back(new bit_vector_type());
        bit_type_seqs.back()->start_write();

    }

    std::cout << "Initialize the LMS char and position sequence." << std::endl;

    std::vector<compressed_pair_type> lms_cPair(m_blocks_info.size());

    for(uint64 i = 0; i < m_blocks_info.size(); i++) /// exclude the leftmost block
    {
        m_cLMS_seqs[i]->start_read();

        /// assign the value to each cPair
        if(!m_cLMS_seqs[i]->is_empty())
        {
            lms_cPair[i] = m_cLMS_seqs[i]->get(); /// LMS char
            m_cLMS_seqs[i]->next_remove();
        }
        else
        {
            std::cout << "the current level_" << m_level << "m_cLMS_seqs["<< i << "]->is_empty() is true.\n ";

            std::cout << "The file : " << __FILE__ << ", the line : " << __LINE__ << std::endl;
        }

        m_LMS_rPos_seqs[i]->start_read(); /// LMS position

    }



    std::cout << "Initialize L-type BWT and type sequence." << std::endl;

    /// allocate a temporary value for each compressed BWT pair
    std::vector<compressed_pair_type> bwt_cPair(m_blocks_info.size());

    typedef Pair<offset_type,bool> bwt_pair_type;

    std::vector < MyVector<bwt_pair_type> * > bwt_pair_seq;

    for(uint64 i = 0; i < m_blocks_info.size(); i++) /// exclude the leftmost block
    {

        m_suf_l_cbwt_seqs[i]->start_read();

        /// assign the value to each cPair
        if(!m_suf_l_cbwt_seqs[i]->is_empty())
        {
            bwt_cPair[i] = m_suf_l_cbwt_seqs[i]->get();
            m_suf_l_cbwt_seqs[i]->next_remove();
        }
        else
        {
            std::cout << "the current level_" << m_level << "m_suf_l_cbwt_seqs["<< i << "]->is_empty() is true.\n ";

            std::cout << "The file : " << __FILE__ << ", the line : " << __LINE__ << std::endl;
        }

        m_suf_l_bit_seqs[i]->start_read_remove();

    }



    std::cout << "start to read s1_blockID_seq." << std::endl;

    s1_blockID_seq->start_read();


    std::cout << "Start to sort L-type suffixes on level " << m_level << ".\n";

    std::cout << "The current PDU = " << Logger::max_pdu << ".\n";

    uint16 cur_alpha_block_id(0);

    typedef std::priority_queue< pq_quadruple_type, std::vector<pq_quadruple_type>, TupleDscCmp2< pq_quadruple_type > > min_pq_type;

    std::vector<uint16> block_id_vec;

    /// sort from 0 to m_alpha
    while( cur_alpha_block_id < m_alpha_blocks_number)
    {

        std::cout << "cur_alpha_block_id = " << cur_alpha_block_id << std::endl;

        offset_type cur_alpha_block_beg(m_alpha_blocks_info[cur_alpha_block_id].m_beg_alpha);

        offset_type cur_alpha_block_end(m_alpha_blocks_info[cur_alpha_block_id].m_end_alpha);

        relative_offset_type alpha_number = cur_alpha_block_end - cur_alpha_block_beg + 1;

        std::cout << "cur_alpha_block_beg = " << cur_alpha_block_beg << std::endl;

        std::cout << "cur_alpha_block_end = " << cur_alpha_block_end << std::endl;

        std::cout << "alpha_number = " << alpha_number << std::endl;


        if(alpha_number > 1)
        {

            std::cout << "the current alpha block is not single block.\n";

            std::vector< pa_triple_type > pa_item_ary;

            typename std::vector< pa_triple_type >::iterator it_beg, it_end;

            if(!bks_seqs[cur_alpha_block_id]->is_empty())
            {
                std::cout << "bks_seqs size is " << bks_seqs[cur_alpha_block_id]->size() << std::endl;

                std::cout << "Load BKS items in to PA." << std::endl;


                std::cout.flush();

                //pa_item_ary.resize( bks_seqs[cur_alpha_block_id]->size() );
                pa_item_ary.reserve( bks_seqs[cur_alpha_block_id]->size() );
                pa_item_ary.clear();

                std::cout << "PA size is " << pa_item_ary.size() * sizeof(pa_triple_type) / K_1024 << "MB. \n";
                std::cout << "pa_item_ary.resize is done.\n";

                std::cout.flush();

                bks_seqs[cur_alpha_block_id]->start_read();

                bit_type_seqs[cur_alpha_block_id]->start_read_remove();

                //relative_offset_type index(0);

                while( !bks_seqs[cur_alpha_block_id]->is_eof() )
                {

                    pa_item_ary.push_back( pa_triple_type(bks_seqs[cur_alpha_block_id]->get().first
                                                          , bks_seqs[cur_alpha_block_id]->get().second
                                                          , bit_type_seqs[cur_alpha_block_id]->get())
                                         );

                    bks_seqs[cur_alpha_block_id]->next_remove();
                    bit_type_seqs[cur_alpha_block_id]->next_remove();

                    //++index;
                }

                delete bks_seqs[cur_alpha_block_id];
                bks_seqs[cur_alpha_block_id] = nullptr;

                delete bit_type_seqs[cur_alpha_block_id];
                bit_type_seqs[cur_alpha_block_id] = nullptr;


                std::cout << "start to sort PA.\n";
                std::cout.flush();

                std::stable_sort(pa_item_ary.begin(), pa_item_ary.end(), TupleAscCmp1<pa_triple_type>()); /// bug, using double working space

                std::cout << "Sorting done.\n";
                std::cout.flush();

            }
            else
            {
                std::cout << "The BKS is null." << std::endl;

                delete bks_seqs[cur_alpha_block_id];
                bks_seqs[cur_alpha_block_id] = nullptr;

                delete bit_type_seqs[cur_alpha_block_id];
                bit_type_seqs[cur_alpha_block_id] = nullptr;
            }

            it_beg = pa_item_ary.begin(), it_end = pa_item_ary.end();

            relative_offset_type pq_input_cnt(0);

            min_pq_type * m_pq = new min_pq_type();

            alphabet_type cur_char(0), pre_char(0);

            offset_type cur_pos, MIN_VALUE(0);

            uint16 pre_blockID(0);

            relative_offset_type cur_rel_bkt(0);

            std::cout << "start to sort the current block.\n";
            std::cout.flush();

            while(cur_rel_bkt < alpha_number)
            {
                cur_char = cur_alpha_block_beg + cur_rel_bkt;

                while ( it_beg != it_end && (*it_beg).first == cur_rel_bkt )
                {

                    cur_pos = (*it_beg).second;

                    is_LStar = (*it_beg).third;

                    {
                        /// L-type suffixes

                        l_suf_bit_seq->push_back(is_LStar);
                        l_suf_pos_seq->push_back(cur_pos);

                        /// compute the compressed character of L-type suffix
                        if(cur_char != pre_cChar)
                        {

                            l_suf_char_seq->push_back(compressed_pair_type(pre_cChar,cChar_cnt));

                            pre_cChar = cur_char, cChar_cnt = 1;

                        }
                        else
                        {

                            if(cChar_cnt != c_max_count)
                            {
                                cChar_cnt++;
                            }
                            else
                            {
                                l_suf_char_seq->push_back(compressed_pair_type(cur_char,cChar_cnt));
                                cChar_cnt = 1;
                            }

                        }

                    }


                    /// 1.4: is_LStar == false
                    if(cur_pos != MIN_VALUE && !is_LStar)
                    {
                        block_id = getBlockId(cur_pos);

                        //block_id_vec.push_back(block_id);

                        /// check whether to get next pair;
                        if( !bwt_cPair[block_id].second )
                        {
                            bwt_cPair[block_id] = m_suf_l_cbwt_seqs[block_id]->get();
                            m_suf_l_cbwt_seqs[block_id]->next_remove();

                            if(bwt_cPair[block_id].second == 0)
                            {
                                std::cout << "error, at line " << __LINE__ << std::endl;
                                std::cin.get();
                            }
                        }

                        pre_char = bwt_cPair[block_id].first;
                        --bwt_cPair[block_id].second;

                        is_LStar = m_suf_l_bit_seqs[block_id]->get();
                        m_suf_l_bit_seqs[block_id]->next_remove();


                        pre_blockID = get_alpha_blockID(cur_alpha_block_id,pre_char);


                        if ( pre_blockID > cur_alpha_block_id )
                        {
                            bks_seqs[pre_blockID]->push_back(bks_pair_type(pre_char - m_alpha_blocks_info[pre_blockID].m_beg_alpha,cur_pos - 1));
                            bit_type_seqs[pre_blockID]->push_back(is_LStar);
                        }
                        else if (pre_blockID == cur_alpha_block_id)
                        {
                            m_pq->push( pq_quadruple_type(relative_offset_type(pre_char - cur_alpha_block_beg), ++pq_input_cnt, cur_pos - 1, is_LStar) );
                            if( (pq_input_cnt % (100 * K_1024)) == 0 )
                            {
                                std::cout<< "PQ push into " << pq_input_cnt / K_1024   << " M items. at line " << __LINE__ << std::endl;
                                std::cout.flush();
                            }

                        }
                        else
                        {
                            Logger::output_error(__FILE__,__LINE__);
                        }

                    }

                    it_beg++;


                }


                /// step2: sort PQ
                while(!m_pq->empty() && (m_pq->top().first == cur_rel_bkt) )
                {

                    cur_pos = m_pq->top().third;

                    is_LStar = m_pq->top().forth;

                    {
                        /// L-type suffixes

                        l_suf_bit_seq->push_back(is_LStar);
                        l_suf_pos_seq->push_back(cur_pos);

                        /// compute the compressed character of L-type suffix
                        if(cur_char != pre_cChar)
                        {

                            l_suf_char_seq->push_back(compressed_pair_type(pre_cChar,cChar_cnt));

                            pre_cChar = cur_char, cChar_cnt = 1;

                        }
                        else
                        {

                            if(cChar_cnt != c_max_count)
                            {
                                cChar_cnt++;
                            }
                            else
                            {
                                l_suf_char_seq->push_back(compressed_pair_type(cur_char,cChar_cnt));
                                cChar_cnt = 1;
                            }

                        }

                    }


                    /// is_LStar == false
                    if(cur_pos != MIN_VALUE  && !is_LStar)
                    {
                        block_id = getBlockId(cur_pos);

                        //if(cur_alpha_block_id == 1)block_id_vec.push_back(block_id);

                        /// check whether to get next pair;
                        if( !bwt_cPair[block_id].second )
                        {
                            bwt_cPair[block_id] = m_suf_l_cbwt_seqs[block_id]->get();
                            m_suf_l_cbwt_seqs[block_id]->next_remove();

                            if( bwt_cPair[block_id].second == 0 )
                            {
                                std::cout << "error, at line " << __LINE__ << std::endl;
                                std::cin.get();
                            }
                        }

                        pre_char = bwt_cPair[block_id].first;
                        --bwt_cPair[block_id].second;

                        is_LStar = m_suf_l_bit_seqs[block_id]->get();
                        m_suf_l_bit_seqs[block_id]->next_remove();

                        pre_blockID = get_alpha_blockID(cur_alpha_block_id,pre_char);

                        if ( pre_blockID > cur_alpha_block_id )
                        {
                            bks_seqs[pre_blockID]->push_back(bks_pair_type(pre_char - m_alpha_blocks_info[pre_blockID].m_beg_alpha,cur_pos - 1));
                            bit_type_seqs[pre_blockID]->push_back(is_LStar);
                        }
                        else if (pre_blockID == cur_alpha_block_id)
                        {
                            m_pq->push( pq_quadruple_type(relative_offset_type(pre_char - cur_alpha_block_beg), ++pq_input_cnt, cur_pos - 1, is_LStar) );
                            if( (pq_input_cnt % (100 * K_1024)) == 0 )
                            {
                                std::cout<< "PQ push into " << pq_input_cnt / K_1024   << "M items. at line " << __LINE__ << std::endl;
                                std::cout.flush();
                            }
                        }
                        else
                        {
                            Logger::output_error(__FILE__,__LINE__);
                        }

                    }

                    m_pq->pop();
                }


                /// step3: induce sort the LMS-bucket
                while(!s1_blockID_seq->is_eof())
                {

                    block_id = s1_blockID_seq->get();

                    /// get lms_cPair from block_id
                    if( !lms_cPair[block_id].second )
                    {
                        lms_cPair[block_id] = m_cLMS_seqs[block_id]->get();
                        m_cLMS_seqs[block_id]->next_remove();

                        if( lms_cPair[block_id].second == 0 )
                        {
                            std::cout << "error, at line " << __LINE__ << std::endl;
                            std::cin.get();
                        }
                    }

                    if(lms_cPair[block_id].first == cur_char)
                    {

                        --lms_cPair[block_id].second;

                        cur_pos = m_LMS_rPos_seqs[block_id]->get() + m_blocks_info[block_id]->m_beg_pos; /// \note cur_pos is global
                        m_LMS_rPos_seqs[block_id]->next_remove();


                        if( !bwt_cPair[block_id].second )
                        {
                            bwt_cPair[block_id] = m_suf_l_cbwt_seqs[block_id]->get();
                            m_suf_l_cbwt_seqs[block_id]->next_remove();

                            if( bwt_cPair[block_id].second == 0 )
                            {
                                std::cout << "error, at line " << __LINE__ << std::endl;
                                std::cin.get();
                            }
                        }

                        pre_char = bwt_cPair[block_id].first;
                        --bwt_cPair[block_id].second;

                        is_LStar = m_suf_l_bit_seqs[block_id]->get();
                        m_suf_l_bit_seqs[block_id]->next_remove();

                        pre_blockID = get_alpha_blockID(cur_alpha_block_id,pre_char);

                        if ( pre_blockID > cur_alpha_block_id )
                        {
                            bks_seqs[pre_blockID]->push_back(bks_pair_type(pre_char - m_alpha_blocks_info[pre_blockID].m_beg_alpha,cur_pos - 1));
                            bit_type_seqs[pre_blockID]->push_back(is_LStar);
                        }
                        else if (pre_blockID == cur_alpha_block_id)
                        {
                            m_pq->push( pq_quadruple_type(relative_offset_type(pre_char - cur_alpha_block_beg), ++pq_input_cnt, cur_pos - 1, is_LStar) );
                            if( (pq_input_cnt % (100 * K_1024)) == 0 )
                            {
                                std::cout<< "PQ push into " << pq_input_cnt / K_1024   << "M items. at line " << __LINE__ << std::endl;
                                std::cout<< "m_pq.size() = " << m_pq->size() << std::endl;

                                std::cout.flush();
                            }
                        }
                        else
                        {
                            std::cout << "cur_char = " << cur_char << ", pre_char = " << pre_char << std::endl;
                            std::cout << "pre_blockID = " << pre_blockID  << std::endl;
                            std::cout << "block_id = " << block_id << std::endl;
                            Logger::output_error(__FILE__,__LINE__);
                        }

                        s1_blockID_seq->next_remove(); /// \note error-prone

                    }
                    else
                    {
                        break; ///  jump out the loop
                    }


                }

                /// next r_bucket
                ++cur_rel_bkt;
            }

            if(!m_pq->empty())
            {
                std::cout << "error, at line " << __LINE__ << std::endl;
                std::cin.get();
            }
            delete m_pq;
            m_pq = nullptr;

            std::cout << "Free PQ."<< std::endl;
            std::cout.flush();
        }
        else if( alpha_number == 1)
        {
            std::cout << "start to sort the singleton-block, the singleton-block ID = " << cur_alpha_block_id << std::endl;

            alphabet_type cur_char(cur_alpha_block_beg), pre_char(0);

            offset_type cur_pos(0);

            uint16 pre_blockID(0);

            typedef Pair< offset_type, bool> same_bkt_pair_type;

            InduceSortSameBkt< same_bkt_pair_type > * t_sortSameBkt = new InduceSortSameBkt< same_bkt_pair_type >(K_1024 * 1024 * 3);


            /// step 1: sort the bks_seqs item
            if(!bks_seqs[cur_alpha_block_id]->is_empty())
            {

                bks_seqs[cur_alpha_block_id]->start_read();
                bit_type_seqs[cur_alpha_block_id]->start_read_remove();

                while(!bks_seqs[cur_alpha_block_id]->is_eof())
                {
                    cur_pos = bks_seqs[cur_alpha_block_id]->get().second;
                    bks_seqs[cur_alpha_block_id]->next_remove();

                    is_LStar = bit_type_seqs[cur_alpha_block_id]->get();
                    bit_type_seqs[cur_alpha_block_id]->next_remove();


                    {
                        /// L-type suffixes

                        l_suf_bit_seq->push_back(is_LStar);
                        l_suf_pos_seq->push_back(cur_pos);

                        /// compute the compressed character of L-type suffix
                        if(cur_char != pre_cChar)
                        {

                            l_suf_char_seq->push_back(compressed_pair_type(pre_cChar,cChar_cnt));

                            pre_cChar = cur_char, cChar_cnt = 1;

                        }
                        else
                        {

                            if(cChar_cnt != c_max_count)
                            {
                                cChar_cnt++;
                            }
                            else
                            {
                                l_suf_char_seq->push_back(compressed_pair_type(cur_char,cChar_cnt));
                                cChar_cnt = 1;
                            }

                        }

                    }


                    /// 1.4: is_LStar == false
                    if( cur_pos && !is_LStar)
                    {
                        block_id = getBlockId(cur_pos);

                        /// check whether to get next pair;
                        if( !bwt_cPair[block_id].second )
                        {
                            bwt_cPair[block_id] = m_suf_l_cbwt_seqs[block_id]->get();
                            m_suf_l_cbwt_seqs[block_id]->next_remove();

                            if(bwt_cPair[block_id].second == 0)
                            {
                                std::cout << "error, at line " << __LINE__ << std::endl;
                                std::cin.get();
                            }
                        }

                        pre_char = bwt_cPair[block_id].first;
                        --bwt_cPair[block_id].second;

                        is_LStar = m_suf_l_bit_seqs[block_id]->get();
                        m_suf_l_bit_seqs[block_id]->next_remove();

                        pre_blockID = get_alpha_blockID(cur_alpha_block_id,pre_char);

                        if ( pre_blockID > cur_alpha_block_id )
                        {
                            bks_seqs[pre_blockID]->push_back(bks_pair_type(pre_char - m_alpha_blocks_info[pre_blockID].m_beg_alpha,cur_pos - 1));
                            bit_type_seqs[pre_blockID]->push_back(is_LStar);
                        }
                        else if (pre_blockID == cur_alpha_block_id)
                        {
                            t_sortSameBkt->push( same_bkt_pair_type(cur_pos - 1, is_LStar) );
                        }
                        else
                        {
                            Logger::output_error(__FILE__,__LINE__);
                        }

                    }


                }

                if(bks_seqs[cur_alpha_block_id]->is_eof())
                {
                    delete bks_seqs[cur_alpha_block_id];
                    bks_seqs[cur_alpha_block_id] = nullptr;

                    delete bit_type_seqs[cur_alpha_block_id];
                    bit_type_seqs[cur_alpha_block_id] = nullptr;

                }
                else
                {
                    Logger::output_error(__FILE__,__LINE__);
                }

            }
            else
            {
                delete bks_seqs[cur_alpha_block_id];
                bks_seqs[cur_alpha_block_id] = nullptr;

                delete bit_type_seqs[cur_alpha_block_id];
                bit_type_seqs[cur_alpha_block_id] = nullptr;
            }

            /// step 2: sort the same bucket item

            if (t_sortSameBkt->write_size())
            {

                t_sortSameBkt->swapBuf();

                while (true)
                {
                    same_bkt_pair_type t_pair;

                    while (t_sortSameBkt->read_size())
                    {
                        t_pair = t_sortSameBkt->pop();

                        cur_pos = t_pair.first;

                        is_LStar = t_pair.second;

                        {
                            /// L-type suffixes

                            l_suf_bit_seq->push_back(is_LStar);
                            l_suf_pos_seq->push_back(cur_pos);

                            /// compute the compressed character of L-type suffix
                            if(cur_char != pre_cChar)
                            {

                                l_suf_char_seq->push_back(compressed_pair_type(pre_cChar,cChar_cnt));

                                pre_cChar = cur_char, cChar_cnt = 1;

                            }
                            else
                            {

                                if(cChar_cnt != c_max_count)
                                {
                                    cChar_cnt++;
                                }
                                else
                                {
                                    l_suf_char_seq->push_back(compressed_pair_type(cur_char,cChar_cnt));
                                    cChar_cnt = 1;
                                }

                            }

                        }


                        /// 1.4: is_LStar == false
                        if( cur_pos && !is_LStar)
                        {
                            block_id = getBlockId(cur_pos);

                            /// check whether to get next pair;
                            if( !bwt_cPair[block_id].second )
                            {
                                bwt_cPair[block_id] = m_suf_l_cbwt_seqs[block_id]->get();
                                m_suf_l_cbwt_seqs[block_id]->next_remove();
                            }

                            pre_char = bwt_cPair[block_id].first;
                            --bwt_cPair[block_id].second;

                            is_LStar = m_suf_l_bit_seqs[block_id]->get();
                            m_suf_l_bit_seqs[block_id]->next_remove();

                            pre_blockID = get_alpha_blockID(cur_alpha_block_id,pre_char);

                            if ( pre_blockID > cur_alpha_block_id )
                            {
                                bks_seqs[pre_blockID]->push_back(bks_pair_type(pre_char - m_alpha_blocks_info[pre_blockID].m_beg_alpha,cur_pos - 1));
                                bit_type_seqs[pre_blockID]->push_back(is_LStar);
                            }
                            else if (pre_blockID == cur_alpha_block_id)
                            {
                                t_sortSameBkt->push( same_bkt_pair_type(cur_pos - 1, is_LStar) );
                            }
                            else
                            {
                                Logger::output_error(__FILE__,__LINE__);
                            }

                        }

                    }

                    if (t_sortSameBkt->write_size())
                        t_sortSameBkt->swapBuf();
                    else
                    {
                        std::cout << "break.\n";
                        break;
                    }

                }

            }

            delete t_sortSameBkt;

            /// step 3: sort the LMS item
            while(!s1_blockID_seq->is_eof())
            {

                block_id = s1_blockID_seq->get();

                /// get lms_cPair from block_id
                if( !lms_cPair[block_id].second )
                {
                    lms_cPair[block_id] = m_cLMS_seqs[block_id]->get();
                    m_cLMS_seqs[block_id]->next_remove();
                }

                if(lms_cPair[block_id].first == cur_char)
                {

                    --lms_cPair[block_id].second;

                    cur_pos = m_LMS_rPos_seqs[block_id]->get() + m_blocks_info[block_id]->m_beg_pos; /// \note cur_pos is global
                    m_LMS_rPos_seqs[block_id]->next_remove();

                    if( !bwt_cPair[block_id].second )
                    {
                        bwt_cPair[block_id] = m_suf_l_cbwt_seqs[block_id]->get();
                        m_suf_l_cbwt_seqs[block_id]->next_remove();
                    }

                    pre_char = bwt_cPair[block_id].first;
                    --bwt_cPair[block_id].second;

                    is_LStar = m_suf_l_bit_seqs[block_id]->get();
                    m_suf_l_bit_seqs[block_id]->next_remove();

                    pre_blockID = get_alpha_blockID(cur_alpha_block_id,pre_char);

                    if ( pre_blockID > cur_alpha_block_id )
                    {
                        bks_seqs[pre_blockID]->push_back(bks_pair_type(pre_char - m_alpha_blocks_info[pre_blockID].m_beg_alpha,cur_pos - 1));
                        bit_type_seqs[pre_blockID]->push_back(is_LStar);
                    }
                    else
                    {
                        Logger::output_error(__FILE__,__LINE__);
                    }

                    s1_blockID_seq->next_remove(); /// \note error-prone

                }
                else
                {

                    break; ///  jump out the loop
                }


            }


        }
        else
        {
            Logger::output_error(__FILE__,__LINE__);
        }
        /// next bucket
        ++cur_alpha_block_id;
    }

    /// push the last item of compressed L-suffix pair
    l_suf_char_seq->push_back(compressed_pair_type(pre_cChar,cChar_cnt)); /// \note error-prone


    /// destructor
    {
        std::cout << "destruct m_suf_l_cbwt_seqs.\n";
        while (!m_suf_l_cbwt_seqs.empty())
        {
            if (!m_suf_l_cbwt_seqs.back()->is_eof())
            {
                Logger::output_error(__FILE__,__LINE__);
            }
            else
            {
                delete m_suf_l_cbwt_seqs.back();
            }

            m_suf_l_cbwt_seqs.pop_back();
        }

        std::cout << "destruct m_suf_l_bit_seqs.\n";
        while (!m_suf_l_bit_seqs.empty())
        {
            if (!m_suf_l_bit_seqs.back()->is_eof())
            {
                Logger::output_error(__FILE__,__LINE__);
            }
            else
            {
                delete m_suf_l_bit_seqs.back();
            }

            m_suf_l_bit_seqs.pop_back();
        }

        std::cout << "destruct s1_blockID_seqs.\n";
        delete s1_blockID_seq;

        std::cout << "destruct m_LMS_rPos_seqs.\n";
        while(!m_LMS_rPos_seqs.empty())
        {

            if(m_LMS_rPos_seqs.back()->is_eof())
            {

                delete m_LMS_rPos_seqs.back();

                m_LMS_rPos_seqs.pop_back();
            }
            else
            {

                Logger::output_error(__FILE__,__LINE__);
            }

        }

        std::cout << "destruct m_cLMS_seqs.\n";
        while(!m_cLMS_seqs.empty())
        {

            if(m_cLMS_seqs.back()->is_eof())
            {

                delete m_cLMS_seqs.back();

                m_cLMS_seqs.pop_back();
            }
            else
            {

                Logger::output_error(__FILE__,__LINE__);
            }

        }

    }

    std::cout << "------------------------------------------------------------------------------------------------------------------------------.\n";

    std::cout << "L-type suffix sorting on level_" << m_level << " is over.\n";

    std::cout << "The current PDU = " << Logger::max_pdu << ".\n";

    std::cout << "------------------------------------------------------------------------------------------------------------------------------.\n";

    std::cout << "Start to sort S-type suffix on level " << m_level << ".\n";



    Logger::output_ram_use();


    bool is_SStar = false;

    /// redefine bucket-section type
    bks_seqs.clear();
    bit_type_seqs.clear();

    for(uint64 i = 0; i < m_alpha_blocks_number; ++i)
    {
        bks_seqs.push_back(new MyVector< bks_pair_type >() );
        bks_seqs.back()->start_write();

        bit_type_seqs.push_back(new bit_vector_type());
        bit_type_seqs.back()->start_write();

    }

    /// allocate a temporary value for each compressed BWT pair
    bwt_cPair.clear();

    for(uint64 i = 0; i < m_blocks_info.size(); i++) /// include the leftmost block
    {

        m_suf_s_cbwt_seqs[i]->start_read();

        /// assign the value to each cPair
        if(!m_suf_s_cbwt_seqs[i]->is_empty())
        {
            bwt_cPair.push_back( m_suf_s_cbwt_seqs[i]->get() );
            m_suf_s_cbwt_seqs[i]->next_remove();
        }
        else
        {
            std::cout << "the current level_" << m_level << "m_suf_s_cbwt_seqs["<< i << "]->is_empty() is true.\n ";

            std::cout << "The file : " << __FILE__ << ", the line : " << __LINE__ << std::endl;
        }

        m_suf_s_bit_seqs[i]->start_read_remove();

    }


    /// start to read the L-type suffix reversely
    std::cout << "Initialize the L-type suffixes reversely.\n";

    l_suf_bit_seq->start_read_remove_reverse();
    l_suf_pos_seq->start_read_reverse();
    l_suf_char_seq->start_read_reverse();

    /// initialize the compressed L-suffix character
    compressed_pair_type t_LSuf_cPair = l_suf_char_seq->get_reverse();
    l_suf_char_seq->next_remove_reverse();

    std::cout << "Initialize the sa1_reverse to start_write()." << std::endl;

    sa1_reverse = new MyVector<offset_type>();
    sa1_reverse->start_write();

    typedef std::priority_queue< pq_quadruple_type, std::vector<pq_quadruple_type>, TupleAscCmp2< pq_quadruple_type > > max_pq_type;

    cur_alpha_block_id = m_alpha_blocks_number - 1;

    std::cout << "start to sort the S-type suffixes.\n";

    while(cur_alpha_block_id >= 0)
    {

        std::cout << "cur_alpha_block_id = " << cur_alpha_block_id << std::endl;

        offset_type cur_alpha_block_beg(m_alpha_blocks_info[cur_alpha_block_id].m_beg_alpha);

        offset_type cur_alpha_block_end(m_alpha_blocks_info[cur_alpha_block_id].m_end_alpha);

        relative_offset_type alpha_number = cur_alpha_block_end - cur_alpha_block_beg + 1;

        std::cout << "cur_alpha_block_beg = " << cur_alpha_block_beg << std::endl;

        std::cout << "cur_alpha_block_end = " << cur_alpha_block_end << std::endl;

        std::cout << "alpha_number = " << alpha_number << std::endl;


        if(alpha_number > 1)
        {

            std::vector< pa_triple_type > pa_item_ary;

            typename std::vector< pa_triple_type >::iterator it_beg, it_end;

            if(!bks_seqs[cur_alpha_block_id]->is_empty())
            {

                pa_item_ary.resize( bks_seqs[cur_alpha_block_id]->size() );

                bks_seqs[cur_alpha_block_id]->start_read();

                bit_type_seqs[cur_alpha_block_id]->start_read_remove();

                uint64 index(0);

                while( !bks_seqs[cur_alpha_block_id]->is_eof() )
                {

                    pa_item_ary[index++] = pa_triple_type(bks_seqs[cur_alpha_block_id]->get().first
                                                          , bks_seqs[cur_alpha_block_id]->get().second
                                                          , bit_type_seqs[cur_alpha_block_id]->get()
                                                         );

                    bks_seqs[cur_alpha_block_id]->next_remove();
                    bit_type_seqs[cur_alpha_block_id]->next_remove();

                }

                delete bks_seqs[cur_alpha_block_id];
                bks_seqs[cur_alpha_block_id] = nullptr;

                delete bit_type_seqs[cur_alpha_block_id];
                bit_type_seqs[cur_alpha_block_id] = nullptr;

                std::stable_sort(pa_item_ary.begin(), pa_item_ary.end(), TupleDscCmp1<pa_triple_type>());

            }
            else
            {

                delete bks_seqs[cur_alpha_block_id];
                bks_seqs[cur_alpha_block_id] = nullptr;

                delete bit_type_seqs[cur_alpha_block_id];
                bit_type_seqs[cur_alpha_block_id] = nullptr;
            }

            it_beg = pa_item_ary.begin(), it_end = pa_item_ary.end();

            relative_offset_type cur_rel_bkt(alpha_number - 1); /// relative alpha

            relative_offset_type pq_input_cnt(std::numeric_limits<relative_offset_type>::max());

            alphabet_type cur_char(0), pre_char(0);

            offset_type cur_pos;

            uint16 pre_blockID(0);

            max_pq_type * m_pq = new max_pq_type();

            std::cout << "cur_rel_bkt = " << cur_rel_bkt << std::endl;


            while(cur_rel_bkt >= 0)
            {
                cur_char = cur_alpha_block_beg + cur_rel_bkt;

                while ( it_beg != it_end && (*it_beg).first == cur_rel_bkt )
                {
                    cur_pos = (*it_beg).second;

                    is_SStar = (*it_beg).third;

                    sa1_reverse->push_back(cur_pos);


                    if(!is_SStar)
                    {
                        block_id = getBlockId(cur_pos);

                        /// check whether to get next pair;
                        if( !bwt_cPair[block_id].second )
                        {
                            bwt_cPair[block_id] = m_suf_s_cbwt_seqs[block_id]->get();
                            m_suf_s_cbwt_seqs[block_id]->next_remove();
                        }

                        pre_char = bwt_cPair[block_id].first;
                        --bwt_cPair[block_id].second;

                        is_SStar = m_suf_s_bit_seqs[block_id]->get();
                        m_suf_s_bit_seqs[block_id]->next_remove();

                        pre_blockID = get_induce_alpha_blockID(cur_alpha_block_id,pre_char);

                        if ( pre_blockID < cur_alpha_block_id )
                        {
                            bks_seqs[pre_blockID]->push_back( bks_pair_type( pre_char - m_alpha_blocks_info[pre_blockID].m_beg_alpha, cur_pos - 1));
                            bit_type_seqs[pre_blockID]->push_back(is_SStar);
                        }
                        else if (pre_blockID == cur_alpha_block_id)
                        {
                            m_pq->push( pq_quadruple_type(relative_offset_type(pre_char - cur_alpha_block_beg), --pq_input_cnt, cur_pos - 1, is_SStar) );
                            if( (pq_input_cnt % (100* K_1024)) == 0 )
                            {
                                std::cout<< "PQ push into " << pq_input_cnt / K_1024   << "M items. at line " << __LINE__ << std::endl;
                                std::cout.flush();
                            }
                        }
                        else
                        {
                            Logger::output_error(__FILE__,__LINE__);
                        }

                    }

                    it_beg++;


                }


                /// step2: sort PQ
                while(!m_pq->empty() && m_pq->top().first == cur_rel_bkt)
                {

                    cur_pos = m_pq->top().third;

                    is_SStar = m_pq->top().forth;

                    sa1_reverse->push_back(cur_pos);

                    if(!is_SStar)
                    {
                        block_id = getBlockId(cur_pos);

                        /// check whether to get next pair;
                        if( !bwt_cPair[block_id].second )
                        {
                            bwt_cPair[block_id] = m_suf_s_cbwt_seqs[block_id]->get();
                            m_suf_s_cbwt_seqs[block_id]->next_remove();
                        }

                        pre_char = bwt_cPair[block_id].first;
                        --bwt_cPair[block_id].second;

                        is_SStar = m_suf_s_bit_seqs[block_id]->get();
                        m_suf_s_bit_seqs[block_id]->next_remove();

                        pre_blockID = get_induce_alpha_blockID(cur_alpha_block_id,pre_char);

                        if ( pre_blockID < cur_alpha_block_id )
                        {
                            bks_seqs[pre_blockID]->push_back( bks_pair_type( pre_char - m_alpha_blocks_info[pre_blockID].m_beg_alpha, cur_pos - 1));
                            bit_type_seqs[pre_blockID]->push_back(is_SStar);
                        }
                        else if (pre_blockID == cur_alpha_block_id)
                        {
                            m_pq->push( pq_quadruple_type(relative_offset_type(pre_char - cur_alpha_block_beg), --pq_input_cnt, cur_pos - 1, is_SStar) );
                            if( (pq_input_cnt % (100* K_1024)) == 0 )
                            {
                                std::cout<< "PQ push into " << pq_input_cnt / K_1024   << "M items. at line " << __LINE__ << std::endl;
                                std::cout.flush();
                            }
                        }
                        else
                        {
                            Logger::output_error(__FILE__,__LINE__);
                        }

                    }


                    m_pq->pop();
                }


                /// step3: induce sort the LMS-bucket
                while(!l_suf_pos_seq->is_eof())
                {
                    if(!t_LSuf_cPair.second)
                    {
                        t_LSuf_cPair = l_suf_char_seq->get_reverse();
                        l_suf_char_seq->next_remove_reverse();
                    }

                    if(t_LSuf_cPair.first == cur_char)
                    {

                        --t_LSuf_cPair.second;

                        cur_pos = l_suf_pos_seq->get_reverse();
                        l_suf_pos_seq->next_remove_reverse(); /// next position

                        is_LStar = l_suf_bit_seq->get_reverse();
                        l_suf_bit_seq->next_remove_reverse();

                        sa1_reverse->push_back(cur_pos);

                        if(is_LStar)
                        {

                            block_id = getBlockId(cur_pos);

                            if( !bwt_cPair[block_id].second )
                            {
                                bwt_cPair[block_id] = m_suf_s_cbwt_seqs[block_id]->get();
                                m_suf_s_cbwt_seqs[block_id]->next_remove();
                            }

                            pre_char = bwt_cPair[block_id].first;
                            --bwt_cPair[block_id].second;

                            is_SStar = m_suf_s_bit_seqs[block_id]->get();
                            m_suf_s_bit_seqs[block_id]->next_remove();

                            pre_blockID = get_induce_alpha_blockID(cur_alpha_block_id,pre_char);

                            if ( pre_blockID < cur_alpha_block_id )
                            {
                                bks_seqs[pre_blockID]->push_back( bks_pair_type( pre_char - m_alpha_blocks_info[pre_blockID].m_beg_alpha, cur_pos - 1));
                                bit_type_seqs[pre_blockID]->push_back(is_SStar);
                            }
                            else if (pre_blockID == cur_alpha_block_id)
                            {
                                m_pq->push( pq_quadruple_type(relative_offset_type(pre_char - cur_alpha_block_beg), --pq_input_cnt, cur_pos - 1, is_SStar) );
                                if( (pq_input_cnt % (100* K_1024)) == 0 )
                                {
                                    std::cout<< "PQ push into " << pq_input_cnt / K_1024   << "M items. at line " << __LINE__ << std::endl;
                                    std::cout.flush();
                                }
                            }
                            else
                            {
                                Logger::output_error(__FILE__,__LINE__);
                            }

                        }


                    }
                    else
                    {

                        break;
                    }

                }

                /// next r_bucket
                if(cur_rel_bkt != 0)--cur_rel_bkt;
                else break;
            }

            if(!m_pq->empty())
            {
                std::cout << "error, at line " << __LINE__ << std::endl;
                std::cin.get();
            }

            delete m_pq;

        }
        else if(alpha_number == 1)
        {
            alphabet_type cur_char(cur_alpha_block_beg), pre_char;

            offset_type cur_pos;

            uint16 block_id,pre_blockID;

            typedef Pair<offset_type, bool> same_bkt_pair_type;

            InduceSortSameBkt<same_bkt_pair_type> * t_sortSameBkt = new InduceSortSameBkt<same_bkt_pair_type>(K_1024 * 1024 * 3);

            /// step 1: sort the bks_seqs item
            if(!bks_seqs[cur_alpha_block_id]->is_empty())
            {

                bks_seqs[cur_alpha_block_id]->start_read();
                bit_type_seqs[cur_alpha_block_id]->start_read_remove();

                while(!bks_seqs[cur_alpha_block_id]->is_eof())
                {
                    cur_pos = bks_seqs[cur_alpha_block_id]->get().second;
                    bks_seqs[cur_alpha_block_id]->next_remove();

                    is_SStar = bit_type_seqs[cur_alpha_block_id]->get();
                    bit_type_seqs[cur_alpha_block_id]->next_remove();

                    sa1_reverse->push_back(cur_pos);

                    /// 1.4: is_LStar == false
                    if(!is_SStar)
                    {
                        block_id = getBlockId(cur_pos);

                        /// check whether to get next pair;
                        if( !bwt_cPair[block_id].second )
                        {
                            bwt_cPair[block_id] = m_suf_s_cbwt_seqs[block_id]->get();
                            m_suf_s_cbwt_seqs[block_id]->next_remove();
                        }

                        pre_char = bwt_cPair[block_id].first;
                        --bwt_cPair[block_id].second;

                        is_SStar = m_suf_s_bit_seqs[block_id]->get();
                        m_suf_s_bit_seqs[block_id]->next_remove();

                        pre_blockID = get_induce_alpha_blockID(cur_alpha_block_id,pre_char);

                        if ( pre_blockID < cur_alpha_block_id )
                        {
                            bks_seqs[pre_blockID]->push_back(bks_pair_type(pre_char - m_alpha_blocks_info[pre_blockID].m_beg_alpha,cur_pos - 1));
                            bit_type_seqs[pre_blockID]->push_back(is_SStar);
                        }
                        else if (pre_blockID == cur_alpha_block_id)
                        {
                            t_sortSameBkt->push( same_bkt_pair_type(cur_pos - 1, is_SStar) );
                        }
                        else
                        {
                            Logger::output_error(__FILE__,__LINE__);
                        }

                    }


                }

                if(bks_seqs[cur_alpha_block_id]->is_eof())
                {
                    delete bks_seqs[cur_alpha_block_id];
                    bks_seqs[cur_alpha_block_id] = nullptr;

                    delete bit_type_seqs[cur_alpha_block_id];
                    bit_type_seqs[cur_alpha_block_id] = nullptr;

                }
                else
                {
                    Logger::output_error(__FILE__,__LINE__);
                }

            }
            else
            {
                delete bks_seqs[cur_alpha_block_id];
                bks_seqs[cur_alpha_block_id] = nullptr;

                delete bit_type_seqs[cur_alpha_block_id];
                bit_type_seqs[cur_alpha_block_id] = nullptr;
            }

            /// step 2: sort the same bucket item

            if (t_sortSameBkt->write_size())
            {

                t_sortSameBkt->swapBuf();

                while (true)
                {
                    same_bkt_pair_type t_pair;

                    while (t_sortSameBkt->read_size())
                    {
                        t_pair = t_sortSameBkt->pop();

                        cur_pos = t_pair.first;

                        is_SStar = t_pair.second;

                        sa1_reverse->push_back(cur_pos);


                        if(!is_SStar)
                        {
                            block_id = getBlockId(cur_pos);

                            /// check whether to get next pair;
                            if( !bwt_cPair[block_id].second )
                            {
                                bwt_cPair[block_id] = m_suf_s_cbwt_seqs[block_id]->get();
                                m_suf_s_cbwt_seqs[block_id]->next_remove();
                            }

                            pre_char = bwt_cPair[block_id].first;
                            --bwt_cPair[block_id].second;

                            is_SStar = m_suf_s_bit_seqs[block_id]->get();
                            m_suf_s_bit_seqs[block_id]->next_remove();

                            pre_blockID = get_induce_alpha_blockID(cur_alpha_block_id,pre_char);

                            if ( pre_blockID < cur_alpha_block_id )
                            {
                                bks_seqs[pre_blockID]->push_back(bks_pair_type(pre_char - m_alpha_blocks_info[pre_blockID].m_beg_alpha,cur_pos - 1));
                                bit_type_seqs[pre_blockID]->push_back(is_SStar);
                            }
                            else if (pre_blockID == cur_alpha_block_id)
                            {
                                t_sortSameBkt->push( same_bkt_pair_type(cur_pos - 1, is_SStar) );
                            }
                            else
                            {
                                Logger::output_error(__FILE__,__LINE__);
                            }

                        }

                    }

                    if (t_sortSameBkt->write_size())
                        t_sortSameBkt->swapBuf();
                    else
                    {
                        std::cout << "break.\n";
                        break;
                    }

                }

            }

            delete t_sortSameBkt;


            /// step3: induce sort the LMS-bucket
            while(!l_suf_pos_seq->is_eof())
            {
                if(!t_LSuf_cPair.second)
                {

                    t_LSuf_cPair = l_suf_char_seq->get_reverse();
                    l_suf_char_seq->next_remove_reverse();
                }

                if(t_LSuf_cPair.first == cur_char)
                {

                    --t_LSuf_cPair.second;

                    cur_pos = l_suf_pos_seq->get_reverse();
                    l_suf_pos_seq->next_remove_reverse();

                    is_LStar = l_suf_bit_seq->get_reverse();
                    l_suf_bit_seq->next_remove_reverse();

                    sa1_reverse->push_back(cur_pos);

                    if(is_LStar)
                    {

                        block_id = getBlockId(cur_pos);

                        if( !bwt_cPair[block_id].second )
                        {
                            bwt_cPair[block_id] = m_suf_s_cbwt_seqs[block_id]->get();
                            m_suf_s_cbwt_seqs[block_id]->next_remove();
                        }

                        pre_char = bwt_cPair[block_id].first;
                        --bwt_cPair[block_id].second;

                        is_SStar = m_suf_s_bit_seqs[block_id]->get();
                        m_suf_s_bit_seqs[block_id]->next_remove();

                        pre_blockID = get_induce_alpha_blockID(cur_alpha_block_id,pre_char);

                        if ( pre_blockID < cur_alpha_block_id )
                        {
                            bks_seqs[pre_blockID]->push_back( bks_pair_type( pre_char - m_alpha_blocks_info[pre_blockID].m_beg_alpha, cur_pos - 1));
                            bit_type_seqs[pre_blockID]->push_back(is_SStar);
                        }
                        else
                        {
                            Logger::output_error(__FILE__,__LINE__);
                        }

                    }


                }
                else
                {

                    break;
                }

            }



        }
        else
        {
            Logger::output_error(__FILE__,__LINE__);
        }


        if(cur_alpha_block_id == 0) break;
        else --cur_alpha_block_id;
    }

    std::cout << "------------------------------------------------------------------------------------------------------------------------------.\n";

    std::cout << "S-type suffixes sorting on level_" << m_level << " is over.\n";

    std::cout << "The current PDU = " << Logger::max_pdu << ".\n";

    std::cout << "------------------------------------------------------------------------------------------------------------------------------.\n";

    /// destruct
    {
        if(l_suf_bit_seq->is_eof()) delete l_suf_bit_seq;
        else
        {
            Logger::output_error(__FILE__,__LINE__);
        }

        if(l_suf_pos_seq->is_eof()) delete l_suf_pos_seq;
        else
        {
            Logger::output_error(__FILE__,__LINE__);
        }

        if(l_suf_char_seq->is_eof()) delete l_suf_char_seq;
        else
        {
            /// the last compressed pair is useless
            std::cout << " compressed pair of l_suf_char_seq is " << l_suf_char_seq->get_reverse().first << ", " << l_suf_char_seq->get_reverse().second << std::endl;

            Logger::output_separator_line(" l_suf_char_seq->is_eof() is false. This situation is  normal.");

            delete l_suf_char_seq;
            //Logger::output_error(__FILE__,__LINE__);
        }


        while(!m_suf_s_cbwt_seqs.empty())
        {
            if(m_suf_s_cbwt_seqs.back()->is_eof())delete m_suf_s_cbwt_seqs.back();
            else
            {
                Logger::output_error(__FILE__,__LINE__);
            }
            m_suf_s_cbwt_seqs.pop_back();
        }

        while(!m_suf_s_bit_seqs.empty())
        {
            if(m_suf_s_bit_seqs.back()->is_eof())delete m_suf_s_bit_seqs.back();
            else
            {
                Logger::output_error(__FILE__,__LINE__);
            }
            m_suf_s_bit_seqs.pop_back();
        }

    }

    std::cerr << "finish inducing S-type suffixes.\n";

    Timer::show_time();

    Timer::add_induce_sort_time(Timer::get_wall_time() - beg_time);

    Logger::add_ind_iv( Logger::cur_iv - ind_iv);

    Logger::add_ind_ov(Logger::cur_ov - ind_ov);

    std::cout << " the cur_level = " << m_level << std::endl;



    return ;
}


/// \brief run
///
/// execute the body
template<typename alphabet_type, typename offset_type, typename compress_type, typename relative_offset_type>
void DSAComputation<alphabet_type, offset_type, compress_type, relative_offset_type>::run()
{
    double beg_time_level_0 = 0;

    if (m_level == 0) beg_time_level_0 = Timer::get_wall_time();

    std::cout << "Before computing blocks information." << std::endl;

    compute_blocks_info(m_s_len - 1);

    std::cout << "After computing blocks information." << std::endl;

    if(!m_blocks_info.empty())
    {

//       std::cout << "Output the block information.\n";

        for (uint64 i = 0; i < m_blocks_info.size(); i++)
        {
            m_blocks_info[i]->m_size = m_blocks_info[i]->m_end_pos - m_blocks_info[i]->m_beg_pos + 1;

            std::cout << i << " blocks : m_beg_pos = " << m_blocks_info[i]->m_beg_pos << " , m_end_pos = " << m_blocks_info[i]->m_end_pos;

            std::cout << " \n m_lms_num = " << m_blocks_info[i]->m_lms_num << "\n m_size = " << m_blocks_info[i]->m_size << " is_single = " << m_blocks_info[i]->is_single() <<std::endl;

            std::cout << "is_left_block = " << m_blocks_info[i]->is_left_block() << std::endl;


        }


    }

    std::cout << "Before computing the bwt" << std::endl;

    double t_time = Timer::get_wall_time();

    uint64 cpt_bwt_iv = Logger::cur_iv;

    uint64 cpt_bwt_ov = Logger::cur_ov;

    compute_bwt();

    std::cout << "After computing the bwt" << std::endl;

    Logger::add_cpt_bwt_iv(Logger::cur_iv - cpt_bwt_iv);

    Logger::add_cpt_bwt_ov(Logger::cur_ov - cpt_bwt_ov);

    Timer::add_bwt_time(Timer::get_wall_time() - t_time);


    std::cout << "sampling position block ID.\n";

    compute_block_id_of_samplings();

//    std::cout << "test getBlockID().\n";
//
//    for(int i = 0; i < m_blocks_info.size(); ++i){
//
//        std::cout << i << "th block's end position = " << getBlockId( m_blocks_info[i]->m_end_pos ) - 0 << std::endl;
//        std::cout << i << "th block's beg position = " << getBlockId( m_blocks_info[i]->m_beg_pos ) - 0 << std::endl;
//
//    }

    //std::cin.get();
    /// sampling the alphabet-block ID
    if(m_level)
    {

        compute_alphabet_blocks_id_of_samplings();

        std::cout << "sampling the alphabet-block ID." << std::endl;

    }

    std::cout << "start to conduct mergeSortedSStarGlobal(). \n";

    double time_beg = Timer::get_wall_time();

    std::cout << "Before sorting substrings" << std::endl;

    Logger::output_ram_use();

    bool is_unique = mergeSortedSStarGlobal();

    std::cout << "After sorting substrings" << std::endl;

    Logger::output_ram_use();

    if (!is_unique)
    {
        /**< UtilityFunctions::getFileLength() return value is the byte number */
        uint64 s1_size = UtilityFunctions::getFileLength(m_s1);

        std::cout << "alphabet on next level = " << alphabetSet_vec.back() << std::endl;

        //std::cin.get();

        uint64 name = alphabetSet_vec.back();

        /// the alphabet size on next level
        uint64 alphabet_size(0);

        if(name >= (uint32_MAX - 1))
        {
            alphabet_size = 5;
            s1_size /= 5;
        }
        else if(name >= (uint24_MAX - 1))
        {
            alphabet_size = 4;
            s1_size /= 4;
        }
        else if(name >= (uint16_MAX - 1))
        {
            alphabet_size = 4;
            s1_size /= 4;
        }
        else if(name >= (uint8_MAX - 1))
        {
            alphabet_size = 2;
            s1_size /= 2;
        }
        else alphabet_size = 1;

        std::cout << "s1_size = " << s1_size << std::endl;

        offset_type sa_size = s1_size * sizeof(offset_type);
        offset_type input_size = s1_size * alphabet_size;
        uint64 bucket_size = (name) * sizeof(offset_type);/// may be lager than 2^40-1
        offset_type type_size = s1_size / 8;

        std::cout << "MAX_MEM = " << MAX_MEM << std::endl;

        std::cout << "(sa_size + input_size + bucket_size + type_size) = " << (sa_size + input_size + bucket_size + type_size) << std::endl;

        if(MAX_MEM > (sa_size + input_size + bucket_size + type_size)) std::cout << "MAX_MEM > (sa_size + input_size + bucket_size + type_size) is true.\n";
        else std::cout << "MAX_MEM > (sa_size + input_size + bucket_size + type_size) is false.\n";

        if(MAX_MEM > (sa_size + input_size + bucket_size + type_size))
            // if (MAX_MEM >= (s1_size * (alphabet_size + alphabet_size + sizeof(offset_type)) + s1_size / 8 ) )
        {

            std::string alpha_blocks_file_name = "alpha_blocks_info_" + std::to_string(m_level) + ".dat";

            std::remove(alpha_blocks_file_name.c_str());

            std::cout << " SAIS in Memory.\n";

            double beg_time = Timer::get_wall_time();

            uint64 cpt_sa_iv = Logger::cur_iv;

            uint64 cpt_sa_ov = Logger::cur_ov;

            if(alphabet_size == 1)SAIS<uint8,offset_type>(m_s1, sa1_reverse);
            else if(alphabet_size == 2)SAIS<uint16,offset_type>(m_s1, sa1_reverse);
            else if(alphabet_size == 3)SAIS<uint32,offset_type>(m_s1, sa1_reverse);
            else if(alphabet_size == 4)SAIS<uint32,offset_type>(m_s1, sa1_reverse);
            else if(alphabet_size == 5)SAIS<uint40,offset_type>(m_s1, sa1_reverse);
            else {
                std::cout << "error. at line " << __LINE__ << std::endl;
                std::cin.get();
            }
            std::cout << "SAIS in RAM is over.\n";

            Timer::add_sais_time(Timer::get_wall_time() - beg_time);

            Timer::prog_mid_time = Timer::get_wall_time();

            Logger::add_cpt_sa_iv(Logger::cur_iv - cpt_sa_iv);

            Logger::add_cpt_sa_ov(Logger::cur_ov - cpt_sa_ov);

            Logger::delPDU(UtilityFunctions::getFileLength(m_s1));

            std::remove(m_s1.c_str());


        }
        else
        {
            if(alphabet_size == 1)
            {
                std::cout << "recursion begin, each alphabet occupys 4 bytes.\n";

                std::cout << " alphabetSet_vec[" << m_level + 1 << "] = " << alphabetSet_vec[m_level + 1] << std::endl;

                DSAComputation<uint8, offset_type, compress_type, relative_offset_type> dsac_recurse(m_s1, m_level + 1, sa1_reverse, alphabetSet_vec[m_level + 1]);
                dsac_recurse.run();

            }
            else if(alphabet_size == 2)
            {
                std::cout << "recursion begin, each alphabet occupys 2 bytes.\n";

                std::cout << " alphabetSet_vec[" << m_level + 1 << "] = " << alphabetSet_vec[m_level + 1] << std::endl;

                DSAComputation<uint16, offset_type, compress_type, relative_offset_type> dsac_recurse(m_s1, m_level + 1, sa1_reverse, alphabetSet_vec[m_level + 1]);
                dsac_recurse.run();

            }
            else if(alphabet_size == 3)
            {
                std::cout << "recursion begin, each alphabet occupys 4 bytes.\n";

                std::cout << " alphabetSet_vec[" << m_level + 1 << "] = " << alphabetSet_vec[m_level + 1] << std::endl;

                DSAComputation<uint32, offset_type, compress_type, relative_offset_type> dsac_recurse(m_s1, m_level + 1, sa1_reverse, alphabetSet_vec[m_level + 1]);
                dsac_recurse.run();

            }
            else if(alphabet_size == 4)
            {

                std::cout << "recursion begin, each alphabet occupys 4 bytes.\n";

                std::cout << " alphabetSet_vec[" << m_level + 1 << "] = " << alphabetSet_vec[m_level + 1] << std::endl;

                DSAComputation<uint32, offset_type, compress_type, relative_offset_type> dsac_recurse(m_s1, m_level + 1, sa1_reverse, alphabetSet_vec[m_level + 1]);
                dsac_recurse.run();
            }
            else
            {

                std::cout << "recursion begin, each alphabet occupys 5 bytes.\n";

                DSAComputation<uint40, offset_type, compress_type, relative_offset_type> dsac_recurse(m_s1, m_level + 1, sa1_reverse, alphabetSet_vec[m_level + 1]);
                dsac_recurse.run();
            }


        }

    }
    else
    {
        double beg_time = Timer::get_wall_time();

        uint64 get_s1_iv = Logger::cur_iv;

        uint64 get_s1_ov = Logger::cur_ov;

        uint64 name = alphabetSet_vec.back(), counter(0);

        std::cout << "the name for next level " << name << std::endl;

        std::cout << "is_unique is true, directly compute the m_s1 as sa.\n";

        typedef Pair<alphabet_type, offset_type> pair_type; /// <name, r_pos>

        typedef TupleDscCmp1<pair_type> pair_comparator_dsc;

        typedef MySorter<pair_type, pair_comparator_dsc> sorter_type;

        sorter_type * sa_sorter = new sorter_type(MAX_MEM / 3 * 2);


        /// exclude the leftmost block
        for (uint64 i = m_blocks_info.size() - 2; i >= 0; --i)
        {

            ///allocate the write buffer in MEM, size equal to the elements number in block
            uint64 t_buf_size = m_blocks_info[i]->m_end_pos - m_blocks_info[i]->m_beg_pos + 1;

            std::cout << "t_buf_size = " << t_buf_size << std::endl;

            alphabet_type * name_vector = new alphabet_type[t_buf_size];

            for(uint64 id = 0; id < t_buf_size; ++id) name_vector[id] = 0;

            s1_blocks[i]->start_read();

            while (!s1_blocks[i]->is_eof())
            {
                name_vector[s1_blocks[i]->get().first] = name - s1_blocks[i]->get().second + 1;

                s1_blocks[i]->next_remove();
            }

            delete s1_blocks[i];

            for(uint64 id = 0; id < t_buf_size; id++)
                if( name_vector[id] ) sa_sorter->push(pair_type( name_vector[id], counter++) ); /// if not equal to zero

            delete []name_vector;

            if (i == 0) break;
        }


        /// \note at the end, append the sentinel to sorter.
        sa_sorter->push(pair_type(0,counter));

        std::cout << "sa_sorter.size() = " << sa_sorter->size();

        std::cout << "start to sort.\n";

        sa_sorter->sort();

        std::cout << "sort end.\n";

        sa1_reverse = new offset_vector_type();

        sa1_reverse->start_write();

        while (!sa_sorter->empty())
        {

            sa1_reverse->push_back((*(*sa_sorter)).second);

            ++(*sa_sorter);

        }


        std::cout << "sa1_reverse is ready when is_unique is true.\n";

        delete sa_sorter;

        sa_sorter = nullptr;

        Timer::add_generate_s1_time(Timer::get_wall_time() - beg_time);

        Logger::add_get_s1_iv(Logger::cur_iv - get_s1_iv);

        Logger::add_get_s1_ov(Logger::cur_ov - get_s1_ov);


    }

    std::string alpha_blocks_file_name = "alpha_blocks_info_" + std::to_string(m_level - 1) + ".dat";
    if(fopen(alpha_blocks_file_name.c_str(),"rb"))
        std::remove(alpha_blocks_file_name.c_str());


    std::cout << "======================= Begin to Solution induction phase ====================== " << std::endl;

    std::cout << "The current level = " << m_level << std::endl;

    std::cout << "begin to execute compute_lms_blocks_id_of_samplings()." << std::endl;

    /// step 1 : sampling LMS position number
    compute_lms_blocks_id_of_samplings();

    std::cout << " m_lms_size = " << m_lms_size << std::endl;


    Logger::output_separator_line("start to compute the BWT in induction stage.");

    /// \note compute BWT
    {

        /// reset the m_LMS_counter
        if(m_alpha <= alphabet_type(255) )
        {

            if(!m_LMS_counter.empty()) m_LMS_counter.clear();

            /// reset the m_LMS_counter
            for(uint64 i(0); i <= m_alpha; ++i)
                m_LMS_counter.push_back(new std::vector<uint64>());

        }

        /// reset m_LMS_rPos_seqs
        m_LMS_rPos_seqs.clear();

        /// reset m_cLMS_Seqs
        m_cLMS_seqs.clear();



        /// the LMS index number in each input block
        std::vector<MyVector<relative_offset_type> * > lms_order_seqs;

        /// 1. scan sa from left to right to classify the position into different seqs
        for (uint64 i = 0; i < m_blocks_info.size(); i++)
        {
            lms_order_seqs.push_back(new MyVector<relative_offset_type>());
            lms_order_seqs.back()->start_write();
        }

        std::cout << " start to sampling LMS position number.\n";


        double beg_time = Timer::get_wall_time();

        uint64 cpt_lms_pos_iv = Logger::cur_iv;

        uint64 cpt_lms_pos_ov = Logger::cur_ov;


        uint16 blockId = 0;

        sa1_reverse->start_read_reverse();

        std::cout << " sa1_reverse size = " << sa1_reverse->size() << std::endl;

        Logger::output_separator_line("scan sa1_reverse from small to large to compute s1_blockID_seq.");

        s1_blockID_seq = new MyVector<uint16>();

        s1_blockID_seq->start_write();

        /// step1 : sampling LMS position number
        while (!sa1_reverse->is_eof())
        {
            blockId = getLmsBlockId(sa1_reverse->get_reverse());

            s1_blockID_seq->push_back(blockId);

            if (m_blocks_info[blockId]->m_lms_num == 1)
                lms_order_seqs[blockId]->push_back(0);
            else
                lms_order_seqs[blockId]->push_back(sa1_reverse->get_reverse() - m_lms_blocks_index[m_blocks_info.size() - 1 - blockId - 1]);

            sa1_reverse->next_remove_reverse();
        }

        delete sa1_reverse;

        sa1_reverse = nullptr;

        std::cout << " sampling LMS position number close.\n";

        Timer::add_compute_lms_pos_time(Timer::get_wall_time() - beg_time);

        Logger::add_cpt_lms_pos_iv(Logger::cur_iv - cpt_lms_pos_iv);

        Logger::add_cpt_lms_pos_ov(Logger::cur_ov - cpt_lms_pos_ov);


        /// compute BWTs
        /*for(uint32 i = 0; i < m_blocks_info.size(); i++)*/
        {
            /// compute BWTs


            /// \note if (getBWT $$ !m_level) is true, then perform the task as follows:

            {
                if(!m_level && is_getBWT)
                {

                    int block_num = m_blocks_info.size();

                    for(int i = 0; i < block_num; i++)
                    {

                        LStarPreChar_seqs.push_back(new MyVector<uint8>());
                        ((MyVector<uint8> *)LStarPreChar_seqs.back())->start_write();

                        SStarPreChar_seqs.push_back(new MyVector<uint8>());
                        ((MyVector<uint8> *)SStarPreChar_seqs.back())->start_write();

                    }

                }

            }



            cpt_bwt_iv = Logger::cur_iv;

            cpt_bwt_ov = Logger::cur_ov;

            beg_time = Timer::get_wall_time();


            /// include the leftmost block
            uint64 loop = m_blocks_info.size() / 4; /// \note four thread

            int is_remainder = (m_blocks_info.size() % 4); /// \note four thread

            for(uint64 i = 0; i < loop; ++i)
            {

                std::cout << "start to compute the BWTs of loop number = " << i << std::endl;

                if(sizeof(alphabet_type) == 1)
                {

                    /// thread 1;
                    bool is_sentinel_1 = (( (m_level == 0) && (i == 0) ) ? true : false);

                    for(uint64 j = 0; j <= m_alpha; ++j)
                        m_LMS_counter[j]->push_back(0);

                    m_suf_l_cbwt_seqs.push_back( new compressed_pair_vector_type() );
                    m_suf_l_cbwt_seqs.back()->start_write();

                    m_suf_l_bit_seqs.push_back( new bit_vector_type() );
                    m_suf_l_bit_seqs.back()->start_write();


                    m_suf_s_cbwt_seqs.push_back( new compressed_pair_vector_type() );
                    m_suf_s_cbwt_seqs.back()->start_write();

                    m_suf_s_bit_seqs.push_back( new bit_vector_type() );
                    m_suf_s_bit_seqs.back()->start_write();

                    m_LMS_rPos_seqs.push_back( new relative_offset_vector_type() );
                    m_LMS_rPos_seqs.back()->start_write();

                    uint64 id_1 = i * 4;

                    std::cout << "id_1 = " << id_1 << std::endl;

                    uint64 beg_pos_1 = m_blocks_info[id_1]->m_beg_pos;
                    uint64 end_pos_1 = m_blocks_info[id_1]->m_end_pos;

                    bool is_single_1 = m_blocks_info[id_1]->is_single();



                    /// thread 2;

                    bool is_sentinel_2 = false;

                    for(uint64 j = 0; j <= m_alpha; ++j)
                        m_LMS_counter[j]->push_back(0);

                    m_suf_l_cbwt_seqs.push_back( new compressed_pair_vector_type() );
                    m_suf_l_cbwt_seqs.back()->start_write();

                    m_suf_l_bit_seqs.push_back( new bit_vector_type() );
                    m_suf_l_bit_seqs.back()->start_write();

                    m_suf_s_cbwt_seqs.push_back( new compressed_pair_vector_type() );
                    m_suf_s_cbwt_seqs.back()->start_write();

                    m_suf_s_bit_seqs.push_back( new bit_vector_type() );
                    m_suf_s_bit_seqs.back()->start_write();

                    m_LMS_rPos_seqs.push_back( new relative_offset_vector_type() );
                    m_LMS_rPos_seqs.back()->start_write();

                    uint64 id_2 = i * 4 + 1;

                    std::cout << "id_2 = " << id_2 << std::endl;

                    uint64 beg_pos_2 = m_blocks_info[id_2]->m_beg_pos;
                    uint64 end_pos_2 = m_blocks_info[id_2]->m_end_pos;

                    bool is_single_2 = m_blocks_info[id_2]->is_single();



                    /// thread 3;

                    bool is_sentinel_3 = false;

                    for(uint64 j = 0; j <= m_alpha; ++j)
                        m_LMS_counter[j]->push_back(0);

                    m_suf_l_cbwt_seqs.push_back( new compressed_pair_vector_type() );
                    m_suf_l_cbwt_seqs.back()->start_write();

                    m_suf_l_bit_seqs.push_back( new bit_vector_type() );
                    m_suf_l_bit_seqs.back()->start_write();

                    m_suf_s_cbwt_seqs.push_back( new compressed_pair_vector_type() );
                    m_suf_s_cbwt_seqs.back()->start_write();

                    m_suf_s_bit_seqs.push_back( new bit_vector_type() );
                    m_suf_s_bit_seqs.back()->start_write();

                    m_LMS_rPos_seqs.push_back( new relative_offset_vector_type() );
                    m_LMS_rPos_seqs.back()->start_write();


                    uint64 id_3 = i * 4 + 2;

                    std::cout << "id_3 = " << id_3 << std::endl;

                    uint64 beg_pos_3 = m_blocks_info[id_3]->m_beg_pos;
                    uint64 end_pos_3 = m_blocks_info[id_3]->m_end_pos;

                    bool is_single_3 = m_blocks_info[id_3]->is_single();


                    /// thread 4;

                    bool is_sentinel_4 = false;

                    for(uint64 j = 0; j <= m_alpha; ++j)
                        m_LMS_counter[j]->push_back(0);

                    m_suf_l_cbwt_seqs.push_back( new compressed_pair_vector_type() );
                    m_suf_l_cbwt_seqs.back()->start_write();

                    m_suf_l_bit_seqs.push_back( new bit_vector_type() );
                    m_suf_l_bit_seqs.back()->start_write();

                    m_suf_s_cbwt_seqs.push_back( new compressed_pair_vector_type() );
                    m_suf_s_cbwt_seqs.back()->start_write();

                    m_suf_s_bit_seqs.push_back( new bit_vector_type() );
                    m_suf_s_bit_seqs.back()->start_write();

                    m_LMS_rPos_seqs.push_back( new relative_offset_vector_type() );
                    m_LMS_rPos_seqs.back()->start_write();


                    uint64 id_4 = i * 4 + 3;

                    std::cout << "id_4 = " << id_4 << std::endl;

                    uint64 beg_pos_4 = m_blocks_info[id_4]->m_beg_pos;
                    uint64 end_pos_4 = m_blocks_info[id_4]->m_end_pos;

                    bool is_single_4 = m_blocks_info[id_4]->is_single();


                    /// thread 1
                    std::thread th1(UtilityFunctions::th_induceSort_bwt_1B<alphabet_type,compress_type,relative_offset_type>
                                    , m_s
                                    , beg_pos_1
                                    , end_pos_1
                                    , m_alpha
                                    , m_level
                                    , id_1
                                    , std::ref(m_LMS_counter)
                                    , std::ref(m_suf_l_cbwt_seqs[id_1])
                                    , std::ref(m_suf_l_bit_seqs[id_1])
                                    , std::ref(m_suf_s_cbwt_seqs[id_1])
                                    , std::ref(m_suf_s_bit_seqs[id_1])
                                    , std::ref(m_LMS_rPos_seqs[id_1])
                                    , is_single_1
                                    , is_sentinel_1
                                    , std::ref(lms_order_seqs[id_1])
                                   );

                    //th1.join(); //2024-03-04

                    /// thread 2
                    std::thread th2(UtilityFunctions::th_induceSort_bwt_1B<alphabet_type,compress_type,relative_offset_type>
                                    , m_s
                                    , beg_pos_2
                                    , end_pos_2
                                    , m_alpha
                                    , m_level
                                    , id_2
                                    , std::ref(m_LMS_counter)
                                    , std::ref(m_suf_l_cbwt_seqs[id_2])
                                    , std::ref(m_suf_l_bit_seqs[id_2])
                                    , std::ref(m_suf_s_cbwt_seqs[id_2])
                                    , std::ref(m_suf_s_bit_seqs[id_2])
                                    , std::ref(m_LMS_rPos_seqs[id_2])
                                    , is_single_2
                                    , is_sentinel_2
                                    , std::ref(lms_order_seqs[id_2])
                                   );

                    /// thread 3
                    std::thread th3(UtilityFunctions::th_induceSort_bwt_1B<alphabet_type,compress_type,relative_offset_type>
                                    , m_s
                                    , beg_pos_3
                                    , end_pos_3
                                    , m_alpha
                                    , m_level
                                    , id_3
                                    , std::ref(m_LMS_counter)
                                    , std::ref(m_suf_l_cbwt_seqs[id_3])
                                    , std::ref(m_suf_l_bit_seqs[id_3])
                                    , std::ref(m_suf_s_cbwt_seqs[id_3])
                                    , std::ref(m_suf_s_bit_seqs[id_3])
                                    , std::ref(m_LMS_rPos_seqs[id_3])
                                    , is_single_3
                                    , is_sentinel_3
                                    , std::ref(lms_order_seqs[id_3])
                                   );
                    /// thread 4
                    std::thread th4(UtilityFunctions::th_induceSort_bwt_1B<alphabet_type,compress_type,relative_offset_type>
                                    , m_s
                                    , beg_pos_4
                                    , end_pos_4
                                    , m_alpha
                                    , m_level
                                    , id_4
                                    , std::ref(m_LMS_counter)
                                    , std::ref(m_suf_l_cbwt_seqs[id_4])
                                    , std::ref(m_suf_l_bit_seqs[id_4])
                                    , std::ref(m_suf_s_cbwt_seqs[id_4])
                                    , std::ref(m_suf_s_bit_seqs[id_4])
                                    , std::ref(m_LMS_rPos_seqs[id_4])
                                    , is_single_4
                                    , is_sentinel_4
                                    , std::ref(lms_order_seqs[id_4])
                                   );

                    th1.join();
                    th2.join();
                    th3.join();
                    th4.join();

                }
                else
                {

                    /// thread 1;
                    bool is_sentinel_1 = false;

                    m_cLMS_seqs.push_back(new compressed_pair_vector_type());
                    m_cLMS_seqs.back()->start_write();

                    m_suf_l_cbwt_seqs.push_back( new compressed_pair_vector_type() );
                    m_suf_l_cbwt_seqs.back()->start_write();

                    m_suf_l_bit_seqs.push_back( new bit_vector_type() );
                    m_suf_l_bit_seqs.back()->start_write();


                    m_suf_s_cbwt_seqs.push_back( new compressed_pair_vector_type() );
                    m_suf_s_cbwt_seqs.back()->start_write();

                    m_suf_s_bit_seqs.push_back( new bit_vector_type() );
                    m_suf_s_bit_seqs.back()->start_write();

                    m_LMS_rPos_seqs.push_back( new relative_offset_vector_type() );
                    m_LMS_rPos_seqs.back()->start_write();


                    uint64 id_1 = i * 4;

                    std::cout << "id_1 = " << id_1 << std::endl;

                    uint64 beg_pos_1 = m_blocks_info[id_1]->m_beg_pos;
                    uint64 end_pos_1 = m_blocks_info[id_1]->m_end_pos;

                    bool is_single_1 = m_blocks_info[id_1]->is_single();


                    /// thread 2;

                    bool is_sentinel_2 = false;

                    m_cLMS_seqs.push_back(new compressed_pair_vector_type());
                    m_cLMS_seqs.back()->start_write();

                    m_suf_l_cbwt_seqs.push_back( new compressed_pair_vector_type() );
                    m_suf_l_cbwt_seqs.back()->start_write();

                    m_suf_l_bit_seqs.push_back( new bit_vector_type() );
                    m_suf_l_bit_seqs.back()->start_write();

                    m_suf_s_cbwt_seqs.push_back( new compressed_pair_vector_type() );
                    m_suf_s_cbwt_seqs.back()->start_write();

                    m_suf_s_bit_seqs.push_back( new bit_vector_type() );
                    m_suf_s_bit_seqs.back()->start_write();

                    m_LMS_rPos_seqs.push_back( new relative_offset_vector_type() );
                    m_LMS_rPos_seqs.back()->start_write();

                    uint64 id_2 = i * 4 + 1;

                    std::cout << "id_2 = " << id_2 << std::endl;

                    uint64 beg_pos_2 = m_blocks_info[id_2]->m_beg_pos;
                    uint64 end_pos_2 = m_blocks_info[id_2]->m_end_pos;

                    bool is_single_2 = m_blocks_info[id_2]->is_single();


                    /// thread 3;

                    bool is_sentinel_3 = false;

                    m_cLMS_seqs.push_back(new compressed_pair_vector_type());
                    m_cLMS_seqs.back()->start_write();

                    m_suf_l_cbwt_seqs.push_back( new compressed_pair_vector_type() );
                    m_suf_l_cbwt_seqs.back()->start_write();

                    m_suf_l_bit_seqs.push_back( new bit_vector_type() );
                    m_suf_l_bit_seqs.back()->start_write();

                    m_suf_s_cbwt_seqs.push_back( new compressed_pair_vector_type() );
                    m_suf_s_cbwt_seqs.back()->start_write();

                    m_suf_s_bit_seqs.push_back( new bit_vector_type() );
                    m_suf_s_bit_seqs.back()->start_write();

                    m_LMS_rPos_seqs.push_back( new relative_offset_vector_type() );
                    m_LMS_rPos_seqs.back()->start_write();

                    uint64 id_3 = i * 4 + 2;

                    std::cout << "id_3 = " << id_3 << std::endl;

                    uint64 beg_pos_3 = m_blocks_info[id_3]->m_beg_pos;
                    uint64 end_pos_3 = m_blocks_info[id_3]->m_end_pos;

                    bool is_single_3 = m_blocks_info[id_3]->is_single();


                    /// thread 4;

                    bool is_sentinel_4 = false;

                    m_cLMS_seqs.push_back(new compressed_pair_vector_type());
                    m_cLMS_seqs.back()->start_write();

                    m_suf_l_cbwt_seqs.push_back( new compressed_pair_vector_type() );
                    m_suf_l_cbwt_seqs.back()->start_write();

                    m_suf_l_bit_seqs.push_back( new bit_vector_type() );
                    m_suf_l_bit_seqs.back()->start_write();

                    m_suf_s_cbwt_seqs.push_back( new compressed_pair_vector_type() );
                    m_suf_s_cbwt_seqs.back()->start_write();

                    m_suf_s_bit_seqs.push_back( new bit_vector_type() );
                    m_suf_s_bit_seqs.back()->start_write();

                    m_LMS_rPos_seqs.push_back( new relative_offset_vector_type() );
                    m_LMS_rPos_seqs.back()->start_write();

                    uint64 id_4 = i * 4 + 3;

                    std::cout << "id_4 = " << id_4 << std::endl;

                    uint64 beg_pos_4 = m_blocks_info[id_4]->m_beg_pos;
                    uint64 end_pos_4 = m_blocks_info[id_4]->m_end_pos;

                    bool is_single_4 = m_blocks_info[id_4]->is_single();


                    /// thread 1
                    std::thread th1(UtilityFunctions::th_induceSort_bwt_4B<alphabet_type,compress_type,relative_offset_type>
                                    , m_s
                                    , beg_pos_1
                                    , end_pos_1
                                    , m_alpha
                                    , m_level
                                    , id_1
                                    , std::ref(m_cLMS_seqs[id_1])
                                    , std::ref(m_suf_l_cbwt_seqs[id_1])
                                    , std::ref(m_suf_l_bit_seqs[id_1])
                                    , std::ref(m_suf_s_cbwt_seqs[id_1])
                                    , std::ref(m_suf_s_bit_seqs[id_1])
                                    , std::ref(m_LMS_rPos_seqs[id_1])
                                    , is_single_1
                                    , is_sentinel_1
                                    , std::ref(lms_order_seqs[id_1])
                                   );

                    //std::cin.get();
                    //th1.join(); //2024-03-04

                    /// thread 2
                    std::thread th2(UtilityFunctions::th_induceSort_bwt_4B<alphabet_type,compress_type,relative_offset_type>
                                    , m_s
                                    , beg_pos_2
                                    , end_pos_2
                                    , m_alpha
                                    , m_level
                                    , id_2
                                    , std::ref(m_cLMS_seqs[id_2])
                                    , std::ref(m_suf_l_cbwt_seqs[id_2])
                                    , std::ref(m_suf_l_bit_seqs[id_2])
                                    , std::ref(m_suf_s_cbwt_seqs[id_2])
                                    , std::ref(m_suf_s_bit_seqs[id_2])
                                    , std::ref(m_LMS_rPos_seqs[id_2])
                                    , is_single_2
                                    , is_sentinel_2
                                    , std::ref(lms_order_seqs[id_2])
                                   );
                    /// thread 3
                    std::thread th3(UtilityFunctions::th_induceSort_bwt_4B<alphabet_type,compress_type,relative_offset_type>
                                    , m_s
                                    , beg_pos_3
                                    , end_pos_3
                                    , m_alpha
                                    , m_level
                                    , id_3
                                    , std::ref(m_cLMS_seqs[id_3])
                                    , std::ref(m_suf_l_cbwt_seqs[id_3])
                                    , std::ref(m_suf_l_bit_seqs[id_3])
                                    , std::ref(m_suf_s_cbwt_seqs[id_3])
                                    , std::ref(m_suf_s_bit_seqs[id_3])
                                    , std::ref(m_LMS_rPos_seqs[id_3])
                                    , is_single_3
                                    , is_sentinel_3
                                    , std::ref(lms_order_seqs[id_3])
                                   );
                    /// thread 4
                    std::thread th4(UtilityFunctions::th_induceSort_bwt_4B<alphabet_type,compress_type,relative_offset_type>
                                    , m_s
                                    , beg_pos_4
                                    , end_pos_4
                                    , m_alpha
                                    , m_level
                                    , id_4
                                    , std::ref(m_cLMS_seqs[id_4])
                                    , std::ref(m_suf_l_cbwt_seqs[id_4])
                                    , std::ref(m_suf_l_bit_seqs[id_4])
                                    , std::ref(m_suf_s_cbwt_seqs[id_4])
                                    , std::ref(m_suf_s_bit_seqs[id_4])
                                    , std::ref(m_LMS_rPos_seqs[id_4])
                                    , is_single_4
                                    , is_sentinel_4
                                    , std::ref(lms_order_seqs[id_4])
                                   );
                    th1.join();
                    th2.join();
                    th3.join();
                    th4.join();

//                    std::cout << "th1 and th2 is closed.\n";

                }

            }

            if(is_remainder == 1)
            {
                std::cout << "********************************* is_remainder = 1 *********************************" << std::endl;
                if(sizeof(alphabet_type) == 1)
                {
                    bool is_sentinel_1 = (( (m_level == 0) && (m_LMS_rPos_seqs.size() == 0) ) ? true : false);

                    for(uint64 j = 0; j <= m_alpha; ++j)
                        m_LMS_counter[j]->push_back(0);

                    m_suf_l_cbwt_seqs.push_back( new compressed_pair_vector_type() );
                    m_suf_l_cbwt_seqs.back()->start_write();

                    m_suf_l_bit_seqs.push_back( new bit_vector_type() );
                    m_suf_l_bit_seqs.back()->start_write();


                    m_suf_s_cbwt_seqs.push_back( new compressed_pair_vector_type() );
                    m_suf_s_cbwt_seqs.back()->start_write();

                    m_suf_s_bit_seqs.push_back( new bit_vector_type() );
                    m_suf_s_bit_seqs.back()->start_write();

                    m_LMS_rPos_seqs.push_back( new relative_offset_vector_type() );
                    m_LMS_rPos_seqs.back()->start_write();

                    uint64 id = m_blocks_info.size()-1;

                    std::cout << "id = " << id << std::endl;

                    uint64 beg_pos_1 = m_blocks_info[id]->m_beg_pos;
                    uint64 end_pos_1 = m_blocks_info[id]->m_end_pos;

                    bool is_single_1 = m_blocks_info[id]->is_single();

                    UtilityFunctions::th_induceSort_bwt_1B<alphabet_type,compress_type,relative_offset_type>(
                        m_s
                        , beg_pos_1
                        , end_pos_1
                        , m_alpha
                        , m_level
                        , id
                        , std::ref(m_LMS_counter)
                        , std::ref(m_suf_l_cbwt_seqs[id])
                        , std::ref(m_suf_l_bit_seqs[id])
                        , std::ref(m_suf_s_cbwt_seqs[id])
                        , std::ref(m_suf_s_bit_seqs[id])
                        , std::ref(m_LMS_rPos_seqs[id])
                        , is_single_1
                        , is_sentinel_1
                        , std::ref(lms_order_seqs[id])
                    );

                }
                else
                {
                    bool is_sentinel_1 = false;

                    m_cLMS_seqs.push_back(new compressed_pair_vector_type());
                    m_cLMS_seqs.back()->start_write();

                    m_suf_l_cbwt_seqs.push_back( new compressed_pair_vector_type() );
                    m_suf_l_cbwt_seqs.back()->start_write();

                    m_suf_l_bit_seqs.push_back( new bit_vector_type() );
                    m_suf_l_bit_seqs.back()->start_write();

                    m_suf_s_cbwt_seqs.push_back( new compressed_pair_vector_type() );
                    m_suf_s_cbwt_seqs.back()->start_write();

                    m_suf_s_bit_seqs.push_back( new bit_vector_type() );
                    m_suf_s_bit_seqs.back()->start_write();

                    m_LMS_rPos_seqs.push_back( new relative_offset_vector_type() );
                    m_LMS_rPos_seqs.back()->start_write();

                    uint64 id = m_blocks_info.size()-1;

                    std::cout << "id = " << id << std::endl;

                    uint64 beg_pos_1 = m_blocks_info[id]->m_beg_pos;
                    uint64 end_pos_1 = m_blocks_info[id]->m_end_pos;

                    bool is_single_1 = m_blocks_info[id]->is_single();

                    UtilityFunctions::th_induceSort_bwt_4B<alphabet_type,compress_type,relative_offset_type>(
                        m_s
                        , beg_pos_1
                        , end_pos_1
                        , m_alpha
                        , m_level
                        , id
                        , std::ref(m_cLMS_seqs[id])
                        , std::ref(m_suf_l_cbwt_seqs[id])
                        , std::ref(m_suf_l_bit_seqs[id])
                        , std::ref(m_suf_s_cbwt_seqs[id])
                        , std::ref(m_suf_s_bit_seqs[id])
                        , std::ref(m_LMS_rPos_seqs[id])
                        , is_single_1
                        , is_sentinel_1
                        , std::ref(lms_order_seqs[id])
                    );
                }


            }

            if(is_remainder == 2)
            {
                std::cout << "********************************* is_remainder = 2 *********************************" << std::endl;
                if(sizeof(alphabet_type) == 1)
                {
                    /// thread 1
                    bool is_sentinel_1 = (( (m_level == 0) && (m_LMS_rPos_seqs.size() == 0) ) ? true : false);

                    for(uint64 j = 0; j <= m_alpha; ++j)
                        m_LMS_counter[j]->push_back(0);

                    m_suf_l_cbwt_seqs.push_back( new compressed_pair_vector_type() );
                    m_suf_l_cbwt_seqs.back()->start_write();

                    m_suf_l_bit_seqs.push_back( new bit_vector_type() );
                    m_suf_l_bit_seqs.back()->start_write();


                    m_suf_s_cbwt_seqs.push_back( new compressed_pair_vector_type() );
                    m_suf_s_cbwt_seqs.back()->start_write();

                    m_suf_s_bit_seqs.push_back( new bit_vector_type() );
                    m_suf_s_bit_seqs.back()->start_write();

                    m_LMS_rPos_seqs.push_back( new relative_offset_vector_type() );
                    m_LMS_rPos_seqs.back()->start_write();

                    uint64 id_1 = m_blocks_info.size()-2;
                    std::cout << "id = " << id_1 << std::endl;

                    uint64 beg_pos_1 = m_blocks_info[id_1]->m_beg_pos;
                    uint64 end_pos_1 = m_blocks_info[id_1]->m_end_pos;

                    bool is_single_1 = m_blocks_info[id_1]->is_single();

                    /// thread 2
                    bool is_sentinel_2 = false;

                    for(uint64 j = 0; j <= m_alpha; ++j)
                        m_LMS_counter[j]->push_back(0);

                    m_suf_l_cbwt_seqs.push_back( new compressed_pair_vector_type() );
                    m_suf_l_cbwt_seqs.back()->start_write();

                    m_suf_l_bit_seqs.push_back( new bit_vector_type() );
                    m_suf_l_bit_seqs.back()->start_write();


                    m_suf_s_cbwt_seqs.push_back( new compressed_pair_vector_type() );
                    m_suf_s_cbwt_seqs.back()->start_write();

                    m_suf_s_bit_seqs.push_back( new bit_vector_type() );
                    m_suf_s_bit_seqs.back()->start_write();

                    m_LMS_rPos_seqs.push_back( new relative_offset_vector_type() );
                    m_LMS_rPos_seqs.back()->start_write();

                    uint64 id_2 = m_blocks_info.size()-1;
                    std::cout << "id = " << id_2 << std::endl;

                    uint64 beg_pos_2 = m_blocks_info[id_2]->m_beg_pos;
                    uint64 end_pos_2 = m_blocks_info[id_2]->m_end_pos;

                    bool is_single_2 = m_blocks_info[id_2]->is_single();


                    /// induce thread 1
                    std::thread th1(UtilityFunctions::th_induceSort_bwt_1B<alphabet_type,compress_type,relative_offset_type>
                                    , m_s
                                    , beg_pos_1
                                    , end_pos_1
                                    , m_alpha
                                    , m_level
                                    , id_1
                                    , std::ref(m_LMS_counter)
                                    , std::ref(m_suf_l_cbwt_seqs[id_1])
                                    , std::ref(m_suf_l_bit_seqs[id_1])
                                    , std::ref(m_suf_s_cbwt_seqs[id_1])
                                    , std::ref(m_suf_s_bit_seqs[id_1])
                                    , std::ref(m_LMS_rPos_seqs[id_1])
                                    , is_single_1
                                    , is_sentinel_1
                                    , std::ref(lms_order_seqs[id_1])
                                   );

                    /// induce thread 2
                    std::thread th2(UtilityFunctions::th_induceSort_bwt_1B<alphabet_type,compress_type,relative_offset_type>
                                    , m_s
                                    , beg_pos_2
                                    , end_pos_2
                                    , m_alpha
                                    , m_level
                                    , id_2
                                    , std::ref(m_LMS_counter)
                                    , std::ref(m_suf_l_cbwt_seqs[id_2])
                                    , std::ref(m_suf_l_bit_seqs[id_2])
                                    , std::ref(m_suf_s_cbwt_seqs[id_2])
                                    , std::ref(m_suf_s_bit_seqs[id_2])
                                    , std::ref(m_LMS_rPos_seqs[id_2])
                                    , is_single_2
                                    , is_sentinel_2
                                    , std::ref(lms_order_seqs[id_2])
                                   );

                    th1.join();
                    th2.join();

                }
                else
                {
                    /// induce  thread 1
                    bool is_sentinel_1 = false;

                    m_cLMS_seqs.push_back(new compressed_pair_vector_type());
                    m_cLMS_seqs.back()->start_write();

                    m_suf_l_cbwt_seqs.push_back( new compressed_pair_vector_type() );
                    m_suf_l_cbwt_seqs.back()->start_write();

                    m_suf_l_bit_seqs.push_back( new bit_vector_type() );
                    m_suf_l_bit_seqs.back()->start_write();

                    m_suf_s_cbwt_seqs.push_back( new compressed_pair_vector_type() );
                    m_suf_s_cbwt_seqs.back()->start_write();

                    m_suf_s_bit_seqs.push_back( new bit_vector_type() );
                    m_suf_s_bit_seqs.back()->start_write();

                    m_LMS_rPos_seqs.push_back( new relative_offset_vector_type() );
                    m_LMS_rPos_seqs.back()->start_write();

                    uint64 id_1 = m_blocks_info.size()-2;

                    uint64 beg_pos_1 = m_blocks_info[id_1]->m_beg_pos;
                    uint64 end_pos_1 = m_blocks_info[id_1]->m_end_pos;

                    bool is_single_1 = m_blocks_info[id_1]->is_single();

                    /// induce thread 2
                    bool is_sentinel_2 = false;

                    m_cLMS_seqs.push_back(new compressed_pair_vector_type());
                    m_cLMS_seqs.back()->start_write();

                    m_suf_l_cbwt_seqs.push_back( new compressed_pair_vector_type() );
                    m_suf_l_cbwt_seqs.back()->start_write();

                    m_suf_l_bit_seqs.push_back( new bit_vector_type() );
                    m_suf_l_bit_seqs.back()->start_write();

                    m_suf_s_cbwt_seqs.push_back( new compressed_pair_vector_type() );
                    m_suf_s_cbwt_seqs.back()->start_write();

                    m_suf_s_bit_seqs.push_back( new bit_vector_type() );
                    m_suf_s_bit_seqs.back()->start_write();

                    m_LMS_rPos_seqs.push_back( new relative_offset_vector_type() );
                    m_LMS_rPos_seqs.back()->start_write();

                    uint64 id_2 = m_blocks_info.size()-1;

                    uint64 beg_pos_2 = m_blocks_info[id_2]->m_beg_pos;
                    uint64 end_pos_2 = m_blocks_info[id_2]->m_end_pos;

                    bool is_single_2 = m_blocks_info[id_2]->is_single();


                    /// induce thread 1
                    std::thread th1(UtilityFunctions::th_induceSort_bwt_4B<alphabet_type,compress_type,relative_offset_type>
                                    , m_s
                                    , beg_pos_1
                                    , end_pos_1
                                    , m_alpha
                                    , m_level
                                    , id_1
                                    , std::ref(m_cLMS_seqs[id_1])
                                    , std::ref(m_suf_l_cbwt_seqs[id_1])
                                    , std::ref(m_suf_l_bit_seqs[id_1])
                                    , std::ref(m_suf_s_cbwt_seqs[id_1])
                                    , std::ref(m_suf_s_bit_seqs[id_1])
                                    , std::ref(m_LMS_rPos_seqs[id_1])
                                    , is_single_1
                                    , is_sentinel_1
                                    , std::ref(lms_order_seqs[id_1])
                                   );

                    /// induce thread 2
                    std::thread th2(UtilityFunctions::th_induceSort_bwt_4B<alphabet_type,compress_type,relative_offset_type>
                                    , m_s
                                    , beg_pos_2
                                    , end_pos_2
                                    , m_alpha
                                    , m_level
                                    , id_2
                                    , std::ref(m_cLMS_seqs[id_2])
                                    , std::ref(m_suf_l_cbwt_seqs[id_2])
                                    , std::ref(m_suf_l_bit_seqs[id_2])
                                    , std::ref(m_suf_s_cbwt_seqs[id_2])
                                    , std::ref(m_suf_s_bit_seqs[id_2])
                                    , std::ref(m_LMS_rPos_seqs[id_2])
                                    , is_single_2
                                    , is_sentinel_2
                                    , std::ref(lms_order_seqs[id_2])
                                   );

                    th1.join();
                    th2.join();

                }


            }

            if(is_remainder == 3)
            {
                std::cout << "********************************* is_remainder = 3 *********************************" << std::endl;
                if(sizeof(alphabet_type) == 1)
                {
                    /// thread 1
                    bool is_sentinel_1 = (( (m_level == 0) && (m_LMS_rPos_seqs.size() == 0) ) ? true : false);

                    for(uint64 j = 0; j <= m_alpha; ++j)
                        m_LMS_counter[j]->push_back(0);

                    m_suf_l_cbwt_seqs.push_back( new compressed_pair_vector_type() );
                    m_suf_l_cbwt_seqs.back()->start_write();

                    m_suf_l_bit_seqs.push_back( new bit_vector_type() );
                    m_suf_l_bit_seqs.back()->start_write();


                    m_suf_s_cbwt_seqs.push_back( new compressed_pair_vector_type() );
                    m_suf_s_cbwt_seqs.back()->start_write();

                    m_suf_s_bit_seqs.push_back( new bit_vector_type() );
                    m_suf_s_bit_seqs.back()->start_write();

                    m_LMS_rPos_seqs.push_back( new relative_offset_vector_type() );
                    m_LMS_rPos_seqs.back()->start_write();

                    uint64 id_1 = m_blocks_info.size()-3;
                    std::cout << "id = " << id_1 << std::endl;

                    uint64 beg_pos_1 = m_blocks_info[id_1]->m_beg_pos;
                    uint64 end_pos_1 = m_blocks_info[id_1]->m_end_pos;

                    bool is_single_1 = m_blocks_info[id_1]->is_single();

                    /// thread 2
                    bool is_sentinel_2 = false;

                    for(uint64 j = 0; j <= m_alpha; ++j)
                        m_LMS_counter[j]->push_back(0);

                    m_suf_l_cbwt_seqs.push_back( new compressed_pair_vector_type() );
                    m_suf_l_cbwt_seqs.back()->start_write();

                    m_suf_l_bit_seqs.push_back( new bit_vector_type() );
                    m_suf_l_bit_seqs.back()->start_write();


                    m_suf_s_cbwt_seqs.push_back( new compressed_pair_vector_type() );
                    m_suf_s_cbwt_seqs.back()->start_write();

                    m_suf_s_bit_seqs.push_back( new bit_vector_type() );
                    m_suf_s_bit_seqs.back()->start_write();

                    m_LMS_rPos_seqs.push_back( new relative_offset_vector_type() );
                    m_LMS_rPos_seqs.back()->start_write();

                    uint64 id_2 = m_blocks_info.size()-2;
                    std::cout << "id = " << id_2 << std::endl;

                    uint64 beg_pos_2 = m_blocks_info[id_2]->m_beg_pos;
                    uint64 end_pos_2 = m_blocks_info[id_2]->m_end_pos;

                    bool is_single_2 = m_blocks_info[id_2]->is_single();


                    /// thread 3
                    bool is_sentinel_3 = false;

                    for(uint64 j = 0; j <= m_alpha; ++j)
                        m_LMS_counter[j]->push_back(0);

                    m_suf_l_cbwt_seqs.push_back( new compressed_pair_vector_type() );
                    m_suf_l_cbwt_seqs.back()->start_write();

                    m_suf_l_bit_seqs.push_back( new bit_vector_type() );
                    m_suf_l_bit_seqs.back()->start_write();


                    m_suf_s_cbwt_seqs.push_back( new compressed_pair_vector_type() );
                    m_suf_s_cbwt_seqs.back()->start_write();

                    m_suf_s_bit_seqs.push_back( new bit_vector_type() );
                    m_suf_s_bit_seqs.back()->start_write();

                    m_LMS_rPos_seqs.push_back( new relative_offset_vector_type() );
                    m_LMS_rPos_seqs.back()->start_write();

                    uint64 id_3 = m_blocks_info.size()-1;
                    std::cout << "id = " << id_3 << std::endl;

                    uint64 beg_pos_3 = m_blocks_info[id_3]->m_beg_pos;
                    uint64 end_pos_3 = m_blocks_info[id_3]->m_end_pos;

                    bool is_single_3 = m_blocks_info[id_3]->is_single();


                    /// induce thread 1
                    std::thread th1(UtilityFunctions::th_induceSort_bwt_1B<alphabet_type,compress_type,relative_offset_type>
                                    , m_s
                                    , beg_pos_1
                                    , end_pos_1
                                    , m_alpha
                                    , m_level
                                    , id_1
                                    , std::ref(m_LMS_counter)
                                    , std::ref(m_suf_l_cbwt_seqs[id_1])
                                    , std::ref(m_suf_l_bit_seqs[id_1])
                                    , std::ref(m_suf_s_cbwt_seqs[id_1])
                                    , std::ref(m_suf_s_bit_seqs[id_1])
                                    , std::ref(m_LMS_rPos_seqs[id_1])
                                    , is_single_1
                                    , is_sentinel_1
                                    , std::ref(lms_order_seqs[id_1])
                                   );

                    /// induce thread 2
                    std::thread th2(UtilityFunctions::th_induceSort_bwt_1B<alphabet_type,compress_type,relative_offset_type>
                                    , m_s
                                    , beg_pos_2
                                    , end_pos_2
                                    , m_alpha
                                    , m_level
                                    , id_2
                                    , std::ref(m_LMS_counter)
                                    , std::ref(m_suf_l_cbwt_seqs[id_2])
                                    , std::ref(m_suf_l_bit_seqs[id_2])
                                    , std::ref(m_suf_s_cbwt_seqs[id_2])
                                    , std::ref(m_suf_s_bit_seqs[id_2])
                                    , std::ref(m_LMS_rPos_seqs[id_2])
                                    , is_single_2
                                    , is_sentinel_2
                                    , std::ref(lms_order_seqs[id_2])
                                   );
                    /// induce thread 3
                    std::thread th3(UtilityFunctions::th_induceSort_bwt_1B<alphabet_type,compress_type,relative_offset_type>
                                    , m_s
                                    , beg_pos_3
                                    , end_pos_3
                                    , m_alpha
                                    , m_level
                                    , id_3
                                    , std::ref(m_LMS_counter)
                                    , std::ref(m_suf_l_cbwt_seqs[id_3])
                                    , std::ref(m_suf_l_bit_seqs[id_3])
                                    , std::ref(m_suf_s_cbwt_seqs[id_3])
                                    , std::ref(m_suf_s_bit_seqs[id_3])
                                    , std::ref(m_LMS_rPos_seqs[id_3])
                                    , is_single_3
                                    , is_sentinel_3
                                    , std::ref(lms_order_seqs[id_3])
                                   );

                    th1.join();
                    th2.join();
                    th3.join();
                }
                else
                {
                    /// induce  thread 1
                    bool is_sentinel_1 = false;

                    m_cLMS_seqs.push_back(new compressed_pair_vector_type());
                    m_cLMS_seqs.back()->start_write();

                    m_suf_l_cbwt_seqs.push_back( new compressed_pair_vector_type() );
                    m_suf_l_cbwt_seqs.back()->start_write();

                    m_suf_l_bit_seqs.push_back( new bit_vector_type() );
                    m_suf_l_bit_seqs.back()->start_write();

                    m_suf_s_cbwt_seqs.push_back( new compressed_pair_vector_type() );
                    m_suf_s_cbwt_seqs.back()->start_write();

                    m_suf_s_bit_seqs.push_back( new bit_vector_type() );
                    m_suf_s_bit_seqs.back()->start_write();

                    m_LMS_rPos_seqs.push_back( new relative_offset_vector_type() );
                    m_LMS_rPos_seqs.back()->start_write();

                    uint64 id_1 = m_blocks_info.size()-3;

                    uint64 beg_pos_1 = m_blocks_info[id_1]->m_beg_pos;
                    uint64 end_pos_1 = m_blocks_info[id_1]->m_end_pos;

                    bool is_single_1 = m_blocks_info[id_1]->is_single();

                    /// induce thread 2
                    bool is_sentinel_2 = false;

                    m_cLMS_seqs.push_back(new compressed_pair_vector_type());
                    m_cLMS_seqs.back()->start_write();

                    m_suf_l_cbwt_seqs.push_back( new compressed_pair_vector_type() );
                    m_suf_l_cbwt_seqs.back()->start_write();

                    m_suf_l_bit_seqs.push_back( new bit_vector_type() );
                    m_suf_l_bit_seqs.back()->start_write();

                    m_suf_s_cbwt_seqs.push_back( new compressed_pair_vector_type() );
                    m_suf_s_cbwt_seqs.back()->start_write();

                    m_suf_s_bit_seqs.push_back( new bit_vector_type() );
                    m_suf_s_bit_seqs.back()->start_write();

                    m_LMS_rPos_seqs.push_back( new relative_offset_vector_type() );
                    m_LMS_rPos_seqs.back()->start_write();

                    uint64 id_2 = m_blocks_info.size()-2;

                    uint64 beg_pos_2 = m_blocks_info[id_2]->m_beg_pos;
                    uint64 end_pos_2 = m_blocks_info[id_2]->m_end_pos;

                    bool is_single_2 = m_blocks_info[id_2]->is_single();


                    /// induce thread 3
                    bool is_sentinel_3 = false;

                    m_cLMS_seqs.push_back(new compressed_pair_vector_type());
                    m_cLMS_seqs.back()->start_write();

                    m_suf_l_cbwt_seqs.push_back( new compressed_pair_vector_type() );
                    m_suf_l_cbwt_seqs.back()->start_write();

                    m_suf_l_bit_seqs.push_back( new bit_vector_type() );
                    m_suf_l_bit_seqs.back()->start_write();

                    m_suf_s_cbwt_seqs.push_back( new compressed_pair_vector_type() );
                    m_suf_s_cbwt_seqs.back()->start_write();

                    m_suf_s_bit_seqs.push_back( new bit_vector_type() );
                    m_suf_s_bit_seqs.back()->start_write();

                    m_LMS_rPos_seqs.push_back( new relative_offset_vector_type() );
                    m_LMS_rPos_seqs.back()->start_write();

                    uint64 id_3 = m_blocks_info.size()-1;

                    uint64 beg_pos_3 = m_blocks_info[id_3]->m_beg_pos;
                    uint64 end_pos_3 = m_blocks_info[id_3]->m_end_pos;

                    bool is_single_3 = m_blocks_info[id_3]->is_single();


                    /// induce thread 1
                    std::thread th1(UtilityFunctions::th_induceSort_bwt_4B<alphabet_type,compress_type,relative_offset_type>
                                    , m_s
                                    , beg_pos_1
                                    , end_pos_1
                                    , m_alpha
                                    , m_level
                                    , id_1
                                    , std::ref(m_cLMS_seqs[id_1])
                                    , std::ref(m_suf_l_cbwt_seqs[id_1])
                                    , std::ref(m_suf_l_bit_seqs[id_1])
                                    , std::ref(m_suf_s_cbwt_seqs[id_1])
                                    , std::ref(m_suf_s_bit_seqs[id_1])
                                    , std::ref(m_LMS_rPos_seqs[id_1])
                                    , is_single_1
                                    , is_sentinel_1
                                    , std::ref(lms_order_seqs[id_1])
                                   );

                    /// induce thread 2
                    std::thread th2(UtilityFunctions::th_induceSort_bwt_4B<alphabet_type,compress_type,relative_offset_type>
                                    , m_s
                                    , beg_pos_2
                                    , end_pos_2
                                    , m_alpha
                                    , m_level
                                    , id_2
                                    , std::ref(m_cLMS_seqs[id_2])
                                    , std::ref(m_suf_l_cbwt_seqs[id_2])
                                    , std::ref(m_suf_l_bit_seqs[id_2])
                                    , std::ref(m_suf_s_cbwt_seqs[id_2])
                                    , std::ref(m_suf_s_bit_seqs[id_2])
                                    , std::ref(m_LMS_rPos_seqs[id_2])
                                    , is_single_2
                                    , is_sentinel_2
                                    , std::ref(lms_order_seqs[id_2])
                                   );
                    /// induce thread 2
                    std::thread th3(UtilityFunctions::th_induceSort_bwt_4B<alphabet_type,compress_type,relative_offset_type>
                                    , m_s
                                    , beg_pos_3
                                    , end_pos_3
                                    , m_alpha
                                    , m_level
                                    , id_3
                                    , std::ref(m_cLMS_seqs[id_3])
                                    , std::ref(m_suf_l_cbwt_seqs[id_3])
                                    , std::ref(m_suf_l_bit_seqs[id_3])
                                    , std::ref(m_suf_s_cbwt_seqs[id_3])
                                    , std::ref(m_suf_s_bit_seqs[id_3])
                                    , std::ref(m_LMS_rPos_seqs[id_3])
                                    , is_single_3
                                    , is_sentinel_3
                                    , std::ref(lms_order_seqs[id_3])
                                   );

                    th1.join();
                    th2.join();
                    th3.join();

                }


            }



            Timer::add_bwt_time(Timer::get_wall_time() - beg_time);

            Logger::add_cpt_bwt_iv(Logger::cur_iv - cpt_bwt_iv);

            Logger::add_cpt_bwt_ov(Logger::cur_ov - cpt_bwt_ov);
        }


        while(!lms_order_seqs.empty())
        {
            if(lms_order_seqs.back() != nullptr)
            {
                if(lms_order_seqs.back()->is_eof())
                {
                    delete lms_order_seqs.back();
                }
                else
                {
                    std::cout << "It is the leftmost block." << std::endl;
                    std::cout << "File: " << __FILE__ << "Line number:" << __LINE__ << std::endl;
                }
            }
            lms_order_seqs.pop_back();
        }
    }

    if(m_level)
    {
        Logger::delPDU(UtilityFunctions::getFileLength(m_s));
        std::remove(m_s.c_str());
    }
    std::cout << "start to induce the suffixes.\n";

    std::cout << "Before sorting suffixes" << std::endl;

    Logger::output_ram_use();

    mergeSortedSuffixGlobal();

    std::cout << "After sorting suffixes" << std::endl;

    Logger::output_ram_use();

    return;
}




/** \brief 此处采样的前提是，至少有2块，如果只有1快的话直接返回
 *
 * \param
 * \param
 * \return
 *
 */
template<typename alphabet_type, typename offset_type, typename compress_type, typename relative_offset_type>
void DSAComputation<alphabet_type, offset_type, compress_type, relative_offset_type>::compute_lms_blocks_id_of_samplings()
{
    if (m_blocks_info.size() == 1)
    {

        std::cout << " the block number is only one, it is impossible, return.\n";

        return;
    }

    /// the left most block
    m_lms_blocks_index.push_back(0);

    m_lms_size = 1;

    if (m_blocks_info.size() > 1)
    {

        for (uint64 i = m_blocks_info.size() - 2; i >= 0; --i)
        {

            m_lms_size = m_lms_size + m_blocks_info[i]->m_lms_num;

            m_lms_blocks_index.push_back(m_lms_size - 1);

            if (i == 0) break;
        }
    }

    /// bug is ocurren, change uint32 to uint64
    //uint64 average = m_lms_size / m_blocks_info.size();
    uint64 average = m_lms_size / m_blocks_info.size();

    uint64 tp = average;

    uint64 block_id = 0;

    m_lms_blocks_of_samplings.push_back(block_id);

    while (block_id < m_blocks_info.size() - 1)
    {
        if (average > m_lms_blocks_index[block_id] && average <= m_lms_blocks_index[block_id + 1])
        {
            m_lms_blocks_of_samplings.push_back(block_id + 1);

            average += tp;
        }
        else
        {

            ++block_id;
        }
    }


    return;
}


template<typename alphabet_type, typename offset_type, typename compress_type, typename relative_offset_type>
uint16 DSAComputation<alphabet_type, offset_type, compress_type, relative_offset_type>::getLmsBlockId(const offset_type & _pos)
{

    if ((uint64)_pos == 0) return m_blocks_info.size() - 1;

    uint64 average = m_lms_size / m_blocks_info.size();

    uint16 block_id = m_lms_blocks_of_samplings[_pos / average];

    if (block_id > 0) --block_id;

    while (block_id < m_blocks_info.size() - 1)
    {
        if (_pos > m_lms_blocks_index[block_id] && _pos <= m_lms_blocks_index[block_id + 1])
        {
            return (m_blocks_info.size() - (block_id + 1) - 1);
        }
        else
        {
            ++block_id;
        }
    }

    std::cout << __LINE__ << " getLmsBlockID() is error.\n";
    std::cout << " _pos = " << _pos << ". \n";
    std::cout << " average = " << average << std::endl;
    std::cout << "block_id = " << m_lms_blocks_of_samplings[_pos / average] << std::endl;
    std::cout << "_pos / average = " << _pos / average << std::endl;
    std::cout << "list m_lms_blocks_of_samplings :";
    for (uint64 i = 0; i < m_lms_blocks_of_samplings.size(); ++i)
    {
        std::cout << i << " : " << m_lms_blocks_of_samplings[i] << std::endl;
    }

    return 0;

}



/// \brief compute samplings' block_id
///
/// \note samplings are positioned at 0, 1 * m_capacity, 2 * m_capacity, ... m_s_len - 1
template<typename alphabet_type, typename offset_type, typename compress_type, typename relative_offset_type>
void DSAComputation<alphabet_type, offset_type, compress_type, relative_offset_type>::compute_block_id_of_samplings()
{

    uint64 sampling_pos = 0;

    ///  size() - 1
    uint64 block_id = m_blocks_info.size() - 1;

    while (sampling_pos < m_s_len - 1)
    {

        while (true)
        {

            if (sampling_pos <= m_blocks_info[block_id]->m_end_pos)
            {

                m_block_id_of_samplings.push_back(block_id);

                break;
            }

            --block_id;
        }

        //sampling_pos += m_blocks_info[0]->m_capacity;
        sampling_pos += m_block_size;
    }

    sampling_pos = m_s_len - 1; // the sentinel is also a sampling

    while (true)
    {

        if (sampling_pos <= m_blocks_info[block_id]->m_end_pos)
        {

            m_block_id_of_samplings.push_back(block_id);

            break;
        }

        --block_id;
    }

    return;
}

/// \brief return the block idx of the character at the given position
///
template<typename alphabet_type, typename offset_type, typename compress_type, typename relative_offset_type>
uint32 DSAComputation<alphabet_type, offset_type, compress_type, relative_offset_type>::getBlockId(const offset_type & _pos)
{

    //uint64 sampling_id = static_cast<uint64>(_pos / m_blocks_info[0]->m_capacity); // _pos >= sampling_id * m_capacity
    uint64 sampling_id = static_cast<uint64>(_pos / m_block_size);

    uint64 block_id = m_block_id_of_samplings[sampling_id];

    while (true)
    {

        if (m_blocks_info[block_id]->m_end_pos >= _pos)
        {

            return block_id;
        }

        --block_id;
    }

    std::cerr << "getBlockId is wrong\n";

    return 0;
}

/// \brief return the block idx of the character at the given position
/// this function can return the right block ID, when the _pos = 0.
template<typename alphabet_type, typename offset_type, typename compress_type, typename relative_offset_type>
uint16 DSAComputation<alphabet_type, offset_type, compress_type, relative_offset_type>::getAlphaBlockId(const alphabet_type & _pos)
{

    uint64 alpha_average = m_alpha / m_alpha_blocks_number;

    uint64 sampling_id = static_cast<uint32>(_pos / alpha_average);

    uint16 block_id = m_alpha_blocks_id_of_samplings[sampling_id];

    while (true)
    {

        if (m_alpha_blocks_info[block_id].m_beg_alpha <= _pos && m_alpha_blocks_info[block_id].m_end_alpha >= _pos)
        {

            return block_id;
        }

        ++block_id;

    }

    std::cerr << "getAlphaBlockId is wrong\n";
    std::cin.get();

    return 0;
}

/// \brief return the alpha block no empty id in each text block
template<typename alphabet_type, typename offset_type, typename compress_type, typename relative_offset_type>
void DSAComputation<alphabet_type, offset_type, compress_type, relative_offset_type>::getAlphaBktID(uint32 block_id, uint32 * _alpha_size)
{

    if( _alpha_size[alpha_bkt_noEmpty_id[block_id]] != 0 ) return ;
    else
    {
        //if(alpha_bkt_noEmpty_id[block_id] == m_alpha_blocks_number) return m_alpha_blocks_number;

        uint64 i = alpha_bkt_noEmpty_id[block_id];

        for(; i < m_alpha_blocks_number; i++)
        {

            if( _alpha_size[i] != 0 )
            {

                alpha_bkt_noEmpty_id[block_id] = i;

                return ;
            }
        }
        return ;
        /*if(i == m_alpha_blocks_number){
            std::cout << "Line : " << __LINE__ << std::endl;
            std::cin.get();
        }*/

    }


}

template<typename alphabet_type, typename offset_type, typename compress_type, typename relative_offset_type>
uint32 DSAComputation<alphabet_type, offset_type, compress_type, relative_offset_type>::getLStarAlphaBktID(uint32 alpha_block_id, uint64 * _alpha_size)
{

    if( _alpha_size[alpha_block_id] != 0 ) return alpha_block_id;
    else
    {

        int32 i = alpha_block_id ;

        for(; i >= 0; i--)
        {

            if( _alpha_size[i] != 0 )
            {

                return i;
            }
        }

        return i;

    }


}

template<typename alphabet_type, typename offset_type, typename compress_type, typename relative_offset_type>
bool DSAComputation<alphabet_type, offset_type, compress_type, relative_offset_type>::get_build_way()
{
    return is_ram_build;
}


template<typename alphabet_type, typename offset_type, typename compress_type, typename relative_offset_type>
void DSAComputation<alphabet_type, offset_type, compress_type, relative_offset_type>::compute_alphabet_blocks(
    std::vector< Alpha_Block<offset_type> > & _alpha_block_vec
    , offset_type &_total_alpha_block_size
    , offset_type &_alpha_block_max_size
    , offset_type &_cur_alpha_size
    , offset_type &_pre_name
    , offset_type _cur_name
    , bool is_change)
{
    if(is_change)
    {

        if (_total_alpha_block_size + _cur_alpha_size < _alpha_block_max_size)
        {

            _total_alpha_block_size += _cur_alpha_size;

            _cur_alpha_size = 0;

        }
        else if (_total_alpha_block_size + _cur_alpha_size == _alpha_block_max_size)
        {

            _alpha_block_vec.push_back(Alpha_Block<offset_type>(_pre_name, _cur_name - 1, _total_alpha_block_size + _cur_alpha_size));

            _cur_alpha_size = 0;

            _total_alpha_block_size = 0;

            _pre_name = _cur_name;

        }
        else
        {

            if(_cur_alpha_size >= _alpha_block_max_size)
            {

                if(_pre_name != offset_type(_cur_name - 1)) /// error-prone : '--_cur_name'
                {
                    _alpha_block_vec.push_back(Alpha_Block<offset_type>(_pre_name, _cur_name - 2, _total_alpha_block_size));
                    _alpha_block_vec.push_back(Alpha_Block<offset_type>(_cur_name - 1, _cur_name - 1, _cur_alpha_size));
                }
                else
                {
                    _alpha_block_vec.push_back(Alpha_Block<offset_type>(_cur_name - 1, _cur_name - 1, _cur_alpha_size));
                }
                _pre_name = _cur_name;
                _total_alpha_block_size = 0;
                _cur_alpha_size = 0;

            }
            else
            {
                _alpha_block_vec.push_back(Alpha_Block<offset_type>(_pre_name, _cur_name - 2, _total_alpha_block_size));
                _pre_name = _cur_name-1;
                _total_alpha_block_size = _cur_alpha_size;
                _cur_alpha_size = 0;
            }

        }
    }

    ++_cur_alpha_size;
}//end


template<typename alphabet_type, typename offset_type, typename compress_type, typename relative_offset_type>
void DSAComputation<alphabet_type, offset_type, compress_type, relative_offset_type>::compute_alphabet_blocks_id_of_samplings()
{

    for(uint64 block_id = 0; block_id < m_alpha_blocks_number; block_id++)
    {
        m_alpha_blocks_id_of_samplings.push_back(m_alpha_blocks_info[block_id].m_beg_alpha);

        std::cout << "m_alpha_blocks_id_of_samplings[" << block_id << "] = " << m_alpha_blocks_id_of_samplings.back() << std::endl;

        m_alpha_blocks_id_of_samplings.push_back(m_alpha_blocks_info[block_id].m_end_alpha);

        std::cout << "m_alpha_blocks_id_of_samplings[" << block_id << "] = " << m_alpha_blocks_id_of_samplings.back() << std::endl;
    }


    return;
}

template<typename alphabet_type, typename offset_type, typename compress_type, typename relative_offset_type>
uint16 DSAComputation<alphabet_type, offset_type, compress_type, relative_offset_type>::get_alpha_blockID(const alphabet_type & _beg_block_id, const alphabet_type & _key)
{

    //uint16 mid(0), _beg(_beg_block_id << 1), _end((m_alpha_blocks_number << 1) - 1);
    uint16 mid(0), _beg(0), _end((m_alpha_blocks_number << 1) - 1);

    while( _beg < _end )
    {

        mid = (_beg + _end ) >> 1;

        if(m_alpha_blocks_id_of_samplings[mid] < _key)
        {

            _beg = mid + 1;

        }
        else if(m_alpha_blocks_id_of_samplings[mid] > _key)
        {

            _end = mid - 1;

        }
        else return (mid >> 1);

    }

    return (_beg >> 1);

}


template<typename alphabet_type, typename offset_type, typename compress_type, typename relative_offset_type>
uint16 DSAComputation<alphabet_type, offset_type, compress_type, relative_offset_type>::get_induce_alpha_blockID(const alphabet_type & _end_block_id, const alphabet_type & _key)
{

    //uint16 mid(0), _beg(0), _end(_end_block_id << 1);

    uint16 mid(0), _beg(0), _end((m_alpha_blocks_number << 1) - 1);

    while( _beg < _end )
    {

        mid = (_beg + _end ) >> 1;

        if(m_alpha_blocks_id_of_samplings[mid] < _key)
        {

            _beg = mid + 1;

        }
        else if(m_alpha_blocks_id_of_samplings[mid] > _key)
        {

            _end = mid - 1;

        }
        else return (mid >> 1);

    }

    return (_beg >> 1);
}

/// s1 alphabet of 2 bytes
template<typename alphabet_type, typename offset_type, typename compress_type, typename relative_offset_type>
void DSAComputation<alphabet_type, offset_type, compress_type, relative_offset_type>::generate_s1_1B( offset_type & _name)
{

    std::cout << "Generating s1 of 1B.\n";

    uint8 name(_name);

    FILE * fot = fopen(m_s1.c_str(), "wb");

    /// exclude the leftmost block
    for (uint64 i = m_blocks_info.size() - 2; i >= 0; --i)
    {

        ///allocate the write buffer in MEM, size equal to the elements number in block
        uint64 t_buf_size = m_blocks_info[i]->m_end_pos - m_blocks_info[i]->m_beg_pos + 1;

        uint8 * name_buf_ary = new uint8[t_buf_size];

        for (uint64 k(0); k < t_buf_size; ++k) name_buf_ary[k] = 0;

        s1_blocks[i]->start_read();

        while (!s1_blocks[i]->is_eof())
        {
            name_buf_ary[s1_blocks[i]->get().first] = name - s1_blocks[i]->get().second + 1;;

            //  bkt_set_counter[ (name - s1_blocks[i]->get().second + 1) >> 20 ];

            s1_blocks[i]->next_remove();
        }

        delete s1_blocks[i];

        /// compact the name_buf_ary
        {
            uint64 index = 0, j = 0;

            while (j < t_buf_size)
            {
                if ( name_buf_ary[j] ) name_buf_ary[index++] = name_buf_ary[j];

                ++j;
            }

            fwrite(name_buf_ary, sizeof(uint8), index, fot);

            Logger::addPDU(index * sizeof(uint8));

            Logger::addOV(index * sizeof(uint8));


            fflush(fot);
        }

        delete[] name_buf_ary;

        if (i == 0) break;
    }


    //append the sentinel
    uint8 sentinel(0);

    fwrite(&sentinel, sizeof(uint8), 1, fot);

    fclose(fot);

    s1_blocks.clear();


}


/// s1 alphabet of 2 bytes
template<typename alphabet_type, typename offset_type, typename compress_type, typename relative_offset_type>
void DSAComputation<alphabet_type, offset_type, compress_type, relative_offset_type>::generate_s1_2B( offset_type & _name)
{

    std::cout << "Generating s1 of 2B.\n";

    uint16 name(_name);

    FILE * fot = fopen(m_s1.c_str(), "wb");

    /// exclude the leftmost block
    for (uint64 i = m_blocks_info.size() - 2; i >= 0; --i)
    {

        ///allocate the write buffer in MEM, size equal to the elements number in block
        uint64 t_buf_size = m_blocks_info[i]->m_end_pos - m_blocks_info[i]->m_beg_pos + 1;

        uint16 * name_buf_ary = new uint16[t_buf_size];

        for (uint64 k(0); k < t_buf_size; ++k) name_buf_ary[k] = 0;

        s1_blocks[i]->start_read();

        while (!s1_blocks[i]->is_eof())
        {
            name_buf_ary[s1_blocks[i]->get().first] = name - s1_blocks[i]->get().second + 1;;

            //  bkt_set_counter[ (name - s1_blocks[i]->get().second + 1) >> 20 ];

            s1_blocks[i]->next_remove();
        }

        delete s1_blocks[i];

        /// compact the name_buf_ary
        {
            uint64 index = 0, j = 0;

            while (j < t_buf_size)
            {
                if ( name_buf_ary[j] ) name_buf_ary[index++] = name_buf_ary[j];

                ++j;
            }

            fwrite(name_buf_ary, sizeof(uint16), index, fot);

            Logger::addPDU(index * sizeof(uint16));

            Logger::addOV(index * sizeof(uint16));


            fflush(fot);
        }

        delete[] name_buf_ary;

        if (i == 0) break;
    }


    //append the sentinel
    uint16 sentinel(0);

    fwrite(&sentinel, sizeof(uint16), 1, fot);

    fclose(fot);

    s1_blocks.clear();


}


/// s1 alphabet of 3 bytes
template<typename alphabet_type, typename offset_type, typename compress_type, typename relative_offset_type>
void DSAComputation<alphabet_type, offset_type, compress_type, relative_offset_type>::generate_s1_3B( offset_type & _name)
{

    std::cout << "Generating s1 of 4B.\n";

    uint32 name(_name);

    FILE * fot = fopen(m_s1.c_str(), "wb");

    /// exclude the leftmost block
    for (uint32 i = m_blocks_info.size() - 2; i >= 0; --i)
    {

        ///allocate the write buffer in MEM, size equal to the elements number in block
        uint64 t_buf_size = m_blocks_info[i]->m_end_pos - m_blocks_info[i]->m_beg_pos + 1;

        uint32 * name_buf_ary = new uint32[t_buf_size];

        for (uint64 k(0); k < t_buf_size; ++k) name_buf_ary[k] = 0;

        s1_blocks[i]->start_read();

        while (!s1_blocks[i]->is_eof())
        {
            name_buf_ary[s1_blocks[i]->get().first] = name - s1_blocks[i]->get().second + 1;

            //  bkt_set_counter[ (name - s1_blocks[i]->get().second + 1) >> 20 ];

            s1_blocks[i]->next_remove();
        }

        delete s1_blocks[i];

        /// compact the name_buf_ary
        {
            uint64 index = 0, j = 0;

            while (j < t_buf_size)
            {
                if ( name_buf_ary[j] ) name_buf_ary[index++] = name_buf_ary[j];

                ++j;
            }

            fwrite(name_buf_ary, sizeof(uint32), index, fot);

            Logger::addPDU(index * sizeof(uint32));

            Logger::addOV(index * sizeof(uint32));


            fflush(fot);
        }

        delete[] name_buf_ary;

        if (i == 0) break;
    }


    //append the sentinel
    uint32 sentinel(0);

    fwrite(&sentinel, sizeof(uint32), 1, fot);

    fclose(fot);

    s1_blocks.clear();


}


/// s1 alphabet of 4 bytes
template<typename alphabet_type, typename offset_type, typename compress_type, typename relative_offset_type>
void DSAComputation<alphabet_type, offset_type, compress_type, relative_offset_type>::generate_s1_4B( offset_type & _name)
{

    std::cout << "Generating s1 of 4B.\n";

    uint32 name(_name);

    FILE * fot = fopen(m_s1.c_str(), "wb");

    /// exclude the leftmost block
    for (uint32 i = m_blocks_info.size() - 2; i >= 0; --i)
    {

        ///allocate the write buffer in MEM, size equal to the elements number in block
        uint64 t_buf_size = m_blocks_info[i]->m_end_pos - m_blocks_info[i]->m_beg_pos + 1;

        uint32 * name_buf_ary = new uint32[t_buf_size];

        for (uint64 k(0); k < t_buf_size; ++k) name_buf_ary[k] = 0;

        s1_blocks[i]->start_read();

        while (!s1_blocks[i]->is_eof())
        {
            name_buf_ary[s1_blocks[i]->get().first] = name - s1_blocks[i]->get().second + 1;;

            //  bkt_set_counter[ (name - s1_blocks[i]->get().second + 1) >> 20 ];

            s1_blocks[i]->next_remove();
        }

        delete s1_blocks[i];

        /// compact the name_buf_ary
        {
            uint64 index = 0, j = 0;

            while (j < t_buf_size)
            {
                if ( name_buf_ary[j] ) name_buf_ary[index++] = name_buf_ary[j];

                ++j;
            }

            fwrite(name_buf_ary, sizeof(uint32), index, fot);

            Logger::addPDU(index * sizeof(uint32));

            Logger::addOV(index * sizeof(uint32));


            fflush(fot);
        }

        delete[] name_buf_ary;

        if (i == 0) break;
    }


    //append the sentinel
    uint32 sentinel(0);

    fwrite(&sentinel, sizeof(uint32), 1, fot);

    fclose(fot);

    s1_blocks.clear();


}

/// s1 alphabet of 5 bytes
template<typename alphabet_type, typename offset_type, typename compress_type, typename relative_offset_type>
void DSAComputation<alphabet_type, offset_type, compress_type, relative_offset_type>::generate_s1_5B( offset_type & _name )
{


    std::cout << "Generating s1 of 5B.\n";

    FILE * fot = fopen(m_s1.c_str(), "wb");

    uint16 t_blockId = uint16_MAX;

    /// exclude the leftmost block
    for (uint32 i = m_blocks_info.size() - 2; i >= 0; --i)
    {

        ///allocate the write buffer in MEM, size equal to the elements number in block
        uint64 t_buf_size = m_blocks_info[i]->m_end_pos - m_blocks_info[i]->m_beg_pos + 1;

        offset_type * name_buf_ary = new offset_type[t_buf_size];

        for (relative_offset_type k(0); k < t_buf_size; ++k) name_buf_ary[k] = 0;

        s1_blocks[i]->start_read();

        while (!s1_blocks[i]->is_eof())
        {
            name_buf_ary[s1_blocks[i]->get().first] = _name - s1_blocks[i]->get().second + 1;;

            //  bkt_set_counter[ (name - s1_blocks[i]->get().second + 1) >> 20 ];

            s1_blocks[i]->next_remove();
        }

        delete s1_blocks[i];

        /// compact the name_buf_ary
        {
            uint64 index = 0, j = 0;

            while (j < t_buf_size)
            {
                if ( name_buf_ary[j] ) name_buf_ary[index++] = name_buf_ary[j];

                ++j;
            }

            fwrite(name_buf_ary, sizeof(offset_type), index, fot);

            Logger::addPDU(index * sizeof(offset_type));

            Logger::addOV(index * sizeof(offset_type));


            fflush(fot);
        }

        delete[] name_buf_ary;

        if (i == 0) break;
    }

    //append the sentinel
    offset_type sentinel(0);

    fwrite(&sentinel, sizeof(offset_type), 1, fot);

    fclose(fot);


    s1_blocks.clear();

}




#endif
