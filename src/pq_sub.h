#ifndef PQ_SUB
#define PQ_SUB

#include<deque>

#include"tuple_sorter.h"

///priority queue in MEM
template<typename alphabet_type, typename offset_type>
class PQL_SUB{

private:


    /**< Triple < char, rank, pos > */
    typedef quadruple<alphabet_type,alphabet_type,alphabet_type,offset_type> quadruple_type;

    typedef TupleDscCmp3< quadruple_type > tuple_Asc_Com3_type;

    std::vector< std::vector<quadruple_type> * > alpha_vec;

    /*uint64 cur_alpha,pre_alpha;

    uint64 cur_rank,pre_rank;*/

    uint32 alpha_size;

    uint32 alpha_block_max_size;

    uint32 alpha_block_number;

    uint32 m_size;

    uint32 cur_block_number;

    //uint32 m_read;

public:

/** \brief constructor
 *
 * \param _alpha : Alphabet_Block.m_end_alpha - Alphabet_Block.m_beg_alpha + 1
 * \param
 * \return
 *
 */
    PQL_SUB(const uint32 & _alpha):
    //cur_alpha(0)
    alpha_size(_alpha)
    ,alpha_block_max_size(0)
    ,alpha_block_number(0)
    ,m_size(0)
  //  ,cur_rank(0)
   // ,pre_rank(uint64_MAX)
   // ,pre_alpha(uint64_MAX)
    ,cur_block_number(0)
    {
        alpha_block_number = 4;

        alpha_block_max_size = alpha_size / alpha_block_number;

        uint32 difference = alpha_size - alpha_block_max_size * alpha_block_number ;

        alpha_block_number = alpha_block_number +  (difference % alpha_block_max_size  == 0  ? difference / alpha_block_max_size : difference / alpha_block_max_size + 1);

        for(uint32 i = 0; i < alpha_block_number; ++i){

            alpha_vec.push_back( new std::vector<quadruple_type>() );
        }

    }

    ~PQL_SUB(){

        for(uint32 i = 0; i < alpha_block_number; ++i){

            if(alpha_vec[i] != nullptr ){

                delete alpha_vec[i];
            }
        }
    }

    quadruple_type top(){

        return alpha_vec[cur_block_number]->front();
    }

    void push( const quadruple_type & _quadruple_type){

/*        pre_alpha = cur_alpha;

        pre_rank = cur_rank;

        cur_alpha = _triple.first;

        cur_rank = _triple.second;*/

        alpha_vec[ _quadruple_type.first / alpha_block_max_size ]->push_back(_quadruple_type);

        if(_quadruple_type.first / alpha_block_max_size == cur_block_number)
            std::push_heap(alpha_vec[ cur_block_number ]->begin(), alpha_vec[ cur_block_number ]->end(), tuple_Asc_Com3_type() );

        ++m_size;
    }

    void pop(){

        std::pop_heap(alpha_vec[cur_block_number]->begin(), alpha_vec[cur_block_number]->end(), tuple_Asc_Com3_type() );

        alpha_vec[cur_block_number]->pop_back();

        if(alpha_vec[cur_block_number]->empty()){

            delete alpha_vec[cur_block_number];

            alpha_vec[cur_block_number] = nullptr;

            ++cur_block_number;

            if( cur_block_number < alpha_block_number)
                std::make_heap(alpha_vec[cur_block_number]->begin(), alpha_vec[cur_block_number]->end(), tuple_Asc_Com3_type());

        }

        --m_size;
    }

    /*bool is_diff(){

        return (cur_alpha != pre_alpha) || (cur_alpha == pre_alpha && cur_rank != pre_rank);
    }*/

    bool is_empty(){

        return m_size == 0;
    }


    uint32 size(){

        return m_size;
    }



};


#endif // PQ_SUB
