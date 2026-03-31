/******************************************************************************
 * sortLMS.h
 *
 ******************************************************************************
 * Ling Bo Han <hanlb (at) mail2.sysu.edu.cn>
 * Yi Wu <wu.yi.christian (at) gmail.com>
 * Ge Nong <issng (at) mail.sysu.edu.cn>
 ******************************************************************************
*/
#ifndef _SORTLMS_H
#define _SORTLMS_H

#include<vector>
#include<algorithm>

#include"common.h"
#include"vector.h"
#include"tuple.h"

template<typename alphabet_type, typename offset_type = uint32>
class lmsSorter {

private:

    //	typedef Pair<offset_type , offset_type> pair_pos_id; // pair(seqs_id,pos_relative)

    typedef Pair<alphabet_type, offset_type> pair_lms_type ;

    typedef MyVector<pair_lms_type> lms_vector_type; // LMS priority queue

    std::vector<lms_vector_type * > m_lms_seqs; // pointer

    alphabet_type m_cur_bkt;

    uint8 m_seqs_id;

    uint32 m_read_size;

    uint32 m_seqs_num;

    uint32 m_size;


public:

    lmsSorter( const std::vector<lms_vector_type * >  & _lms_seqs)
        :    m_lms_seqs(_lms_seqs)
        ,m_cur_bkt(std::numeric_limits<alphabet_type>::max())
        ,m_seqs_id(_lms_seqs.size()-1)
        ,m_read_size(0)
        ,m_seqs_num(_lms_seqs.size())
        ,m_size(0) {

        for(uint32 i = 0; i < m_seqs_num; i++) {

            m_lms_seqs[i]->start_read();

            m_size += m_lms_seqs[i]->size();

        }

        // std::cout << "m_size = " << m_size << std::endl;

        // std::cout << "m_seqs_num = " << m_seqs_num << std::endl;

    }

    /** \brief compute the m_cur_bkt and m_seqs_id
     *
     * \param void
     * \param
     * \return void
     *
     */
    void initBkt() {

        m_cur_bkt =  std::numeric_limits<alphabet_type>::max();

        for(int32 i = m_seqs_num - 1; i >= 0; --i) {

            if( !m_lms_seqs[i]->is_eof() && m_lms_seqs[i]->get().first < m_cur_bkt ) {

                m_seqs_id = i;

                m_cur_bkt = m_lms_seqs[i]->get().first;

            }
        }

        return;
    }


    /** \brief compute the next LMS sequence ID in BlockInfo vector
     *
     * \param void
     * \param
     * \return void
     *
     */

    void computeBkt() {

        for(int32 i = m_seqs_id; i >= 0; --i) {

            if( !m_lms_seqs[i]->is_eof() && m_lms_seqs[i]->get().first == m_cur_bkt ) {

                m_seqs_id = i;

                return ;

            }
        }

        initBkt();

        return;
    }


    alphabet_type getLmsBkt() const {

        return m_cur_bkt;
    }

    uint8 getSeqsId() const {

        return m_seqs_id;
    }

    uint32 getLmsPos() {

        return m_lms_seqs[m_seqs_id]->get().second;
    }

    void next() {

        //  std::cout << " m_seqs_id = " << m_seqs_id - 0 << std::endl;

        ++m_read_size;

        m_lms_seqs[m_seqs_id]->next_remove();

        if(!is_empty()) {

            if(m_lms_seqs[m_seqs_id]->is_eof()) {

                if( !m_seqs_id ) {

                    initBkt();

                } else {

                    --m_seqs_id;

                    computeBkt();
                }

                //     std::cout << "1. m_seqs_id = " << m_seqs_id - 0 << std::endl;

                return ;

            }

            if( !m_lms_seqs[m_seqs_id]->is_eof() && m_lms_seqs[m_seqs_id]->get().first != m_cur_bkt ) {

                if( !m_seqs_id ) {

                    initBkt();

                } else {

                    --m_seqs_id;

                    computeBkt();

                    //      std::cout << "2. m_seqs_id = " << m_seqs_id - 0 << std::endl;
                }

                return ;

            }
        }

        return ;

    }

    bool is_empty() const {

        return m_size == m_read_size;
    }

};


template < typename alphabet_type, typename offset_type = uint32>
class LMSHeapSorter {

private:

    typedef Pair<alphabet_type, offset_type> pair_lms_type ;

    typedef MyVector< pair_lms_type > lms_vector_type;

    typedef Pair< alphabet_type, uint16> pair_type;

    std::vector< lms_vector_type * >  m_lms_seqs; // pointer

    std::vector< pair_type > pair_vector;

    uint64 m_size;

    uint64 m_read_size;

public:
    LMSHeapSorter( std::vector<lms_vector_type * > & _lms_seqs )
        :m_lms_seqs(_lms_seqs)
        ,pair_vector(_lms_seqs.size())
        ,m_size(0)
        ,m_read_size(0) {

        std::cout << " m_lms_seqs size in SortLMS.h :" << pair_vector.size() << std::endl;


        for(uint32 i = 0; i < _lms_seqs.size(); i++) {

            m_lms_seqs[i]->start_read();

            m_size += m_lms_seqs[i]->size();

            if(!m_lms_seqs[i]->is_eof()) {

                pair_vector[i].first = m_lms_seqs[i]->get().first;

                pair_vector[i].second = i;

            }

        }

      //  __DEBUG_INFO__ ;
    }

    ~LMSHeapSorter(){

//        while(!m_lms_seqs.empty()){
//
//            if(m_lms_seqs.back())delete m_lms_seqs.back();
//
//            m_lms_seqs.pop_back();
//        }

    }

    void initBkt(){

        //__DEBUG_INFO__ ;

        std::make_heap(pair_vector.begin(),pair_vector.end(),TupleAscCmp2_LMS< pair_type >() );

    }


    alphabet_type getLmsBkt() const {

        return pair_vector[0].first;
    }

    uint16 getSeqsId() const {

        return pair_vector[0].second;
    }

    uint32 getLmsPos() {

        return m_lms_seqs[pair_vector[0].second]->get().second;
    }

    void next(){

        ++m_read_size;

        m_lms_seqs[pair_vector[0].second]->next_remove();

        if(!is_empty()){

            if( !m_lms_seqs[pair_vector[0].second]->is_eof()){

                if( m_lms_seqs[pair_vector[0].second]->get().first != pair_vector[0].first ){

                    pair_vector[0].first = m_lms_seqs[pair_vector[0].second]->get().first;

                    initBkt();
                }
            }else{

                std::pop_heap(pair_vector.begin(),pair_vector.end(),TupleAscCmp2_LMS< pair_type >() );

                //delete m_lms_seqs[pair_vector.back().second];

                //m_lms_seqs[pair_vector.back().second] = nullptr;

//                std::cout << " LMS seqs is nullpt or not : " << m_lms_seqs[pair_vector.back().second] << std::endl;
//
//                std::cin.get();

                pair_vector.pop_back();
            }
        }else {

            if(!pair_vector.empty()){

                //delete m_lms_seqs[pair_vector.back().second];

              //  m_lms_seqs[pair_vector.back().second] = nullptr;

//                std::cout << " LMS seqs is nullpt or not : " << m_lms_seqs[pair_vector.back().second] << std::endl;
//
//                std::cin.get();

                pair_vector.pop_back();

            }

        }

    }

    bool is_empty() const {

        return m_size == m_read_size;
    }

};


//template < typename alphabet_type, typename offset_type = uint32>
/// for level_1
class LMSHeapSorter_1 {

private:

    typedef Pair<uint24, uint32> pair_lms_type ;

    typedef MyVector< pair_lms_type > lms_vector_type;

    typedef Pair< uint40, uint16> pair_type;

    std::vector< lms_vector_type * >  m_lms_seqs; // pointer

    std::vector< pair_type > pair_vector;

    //std::vector< std::vector<uint32> * > alpha_size_in_each_block; /// 各个文本块中，各个相对字符的个数

    ///为各个文本块分配一个计数器，记录各个字符块的元素个数
    std::vector< uint32 * > alpha_bkt_size_in_each_block;
    ///记录各个文本块中，非空字符块的ID
    std::vector< uint32 > alpha_bkt_noEmpty_id;

    Alpha_Block<uint64> * m_alpha_blocks_info;

    uint32 m_alpha_blocks_number;/// 字符块个数

    uint64 m_size;

    uint64 m_read_size;

public:
    LMSHeapSorter_1( std::vector<lms_vector_type * > & _lms_seqs
                    , std::vector< uint32 * > & _alpha_bkt_size_in_each_block
                    , Alpha_Block<uint64> * _alpha_blocks_info
                    , uint32 _alpha_blocks_number )
        :m_lms_seqs(_lms_seqs)
        ,alpha_bkt_size_in_each_block(_alpha_bkt_size_in_each_block)
        ,m_alpha_blocks_info(_alpha_blocks_info)
        ,m_alpha_blocks_number(_alpha_blocks_number)
        ,pair_vector(_lms_seqs.size())
        ,m_size(0)
        ,m_read_size(0) {

        std::cout << "m_lms_seqs.size() = " << m_lms_seqs.size() << std::endl;

        std::cout << "m_alpha_blocks_number= " << m_alpha_blocks_number << std::endl;


        /*for(int i = 0; i < m_lms_seqs.size(); ++i){
                std::cout << "the " << i << "th text block : ";
            for(int j = 0; j < m_alpha_blocks_number; j++){

                std::cout << j << " : " << alpha_bkt_size_in_each_block[i][j] << " ";
            }
            std::cout << std::endl;
        }*/

        /*for(int j = 0; j < m_alpha_blocks_number; j++){

            std::cout << j << " : " << m_alpha_blocks_info[j].m_beg_alpha << "," << m_alpha_blocks_info[j].m_end_alpha << std::endl;
        }*/

        std::cout << " m_lms_seqs size in LMSHeapSorter_1 in sortLMS.h :" << pair_vector.size() << std::endl;

        alpha_bkt_noEmpty_id.clear();

        for(uint32 i = 0; i < m_lms_seqs.size(); i++) alpha_bkt_noEmpty_id.push_back(0);

        for(uint32 i = 0; i < m_lms_seqs.size(); i++) {

            m_lms_seqs[i]->start_read();

            m_size += m_lms_seqs[i]->size();

            getAlphaBktID(i,alpha_bkt_size_in_each_block[i]);

            if(!m_lms_seqs[i]->is_eof()) {

                pair_vector[i].first = uint64(m_lms_seqs[i]->get().first) + m_alpha_blocks_info[alpha_bkt_noEmpty_id[i]].m_beg_alpha;

                pair_vector[i].second = i;

               // std::cout << "m_lms_seqs[" << i << "] = " << m_lms_seqs[i]->get().first << "  " << std::endl;

            }

        }

      //  for(uint32 i = 0; i < _lms_seqs.size(); i++) std::cout << i << alpha_bkt_noEmpty_id[i] << std::endl;



      //  __DEBUG_INFO__ ;
    }

    ~LMSHeapSorter_1(){

//        while(!m_lms_seqs.empty()){
//
//            if(m_lms_seqs.back())delete m_lms_seqs.back();
//
//            m_lms_seqs.pop_back();
//        }

    }

    void initBkt(){


        std::make_heap(pair_vector.begin(),pair_vector.end(),TupleAscCmp2_LMS< pair_type >() );

    }


    uint40 getLmsBkt() const {

        return pair_vector[0].first;
    }

    uint16 getSeqsId() const {

        return pair_vector[0].second;
    }

    uint32 getLmsPos() {

        return m_lms_seqs[pair_vector[0].second]->get().second;
    }

    void next(){

        ++m_read_size;

        m_lms_seqs[pair_vector[0].second]->next_remove();

        if(alpha_bkt_size_in_each_block[ pair_vector[0].second ][alpha_bkt_noEmpty_id[ pair_vector[0].second  ]])
            alpha_bkt_size_in_each_block[pair_vector[0].second][alpha_bkt_noEmpty_id[pair_vector[0].second]]--;

        getAlphaBktID(pair_vector[0].second,alpha_bkt_size_in_each_block[pair_vector[0].second]);

        if(!is_empty()){

            if( !m_lms_seqs[pair_vector[0].second]->is_eof()){

                if( uint64(m_lms_seqs[pair_vector[0].second]->get().first) + m_alpha_blocks_info[alpha_bkt_noEmpty_id[ pair_vector[0].second ]].m_beg_alpha
                   != pair_vector[0].first ){

                    pair_vector[0].first = uint64(m_lms_seqs[pair_vector[0].second]->get().first) + m_alpha_blocks_info[alpha_bkt_noEmpty_id[ pair_vector[0].second ]].m_beg_alpha;

                    initBkt();
                   // std::cout << "-------------" << pair_vector[0].first << std::endl;
                }
            }else{

                std::pop_heap(pair_vector.begin(),pair_vector.end(),TupleAscCmp2_LMS< pair_type >() );

                pair_vector.pop_back();
            }
        }else {

            if(!pair_vector.empty()){

                pair_vector.pop_back();

            }

        }

    }

    bool is_empty() const {

        return m_size == m_read_size;
    }

    void getAlphaBktID(uint16 block_id, uint32 * _alpha_size){

        if( _alpha_size[alpha_bkt_noEmpty_id[block_id]] != 0 ) return;
        else{
            //if(alpha_bkt_noEmpty_id[block_id] == m_alpha_blocks_number) return m_alpha_blocks_number;

            uint32 i = alpha_bkt_noEmpty_id[block_id];

            for(; i < m_alpha_blocks_number; i++){

                if( _alpha_size[i] != 0 ) {

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

};


///20190219
class SStarSorter {

private:

   typedef Pair<uint24,uint40> pair_type;

   uint32 bkt_id;

   uint64 total;

   uint64 read_size;

  // bool is_reverse;

   std::vector< MyVector<Pair<uint24,uint40> > *> bkt_vector;

public:
    SStarSorter(std::vector< MyVector<Pair<uint24,uint40> > *> _bkt_vector):bkt_vector(_bkt_vector) {

        bkt_id = 0;

        total = 0;

        for(uint32 i = 0; i < bkt_vector.size(); i++) {

            std::cout << i << " block size = " << bkt_vector[i]->size();
            total += bkt_vector[i]->size();
        }
        read_size = 0;

        std::cout << "total = " << total << std::endl;

        std::cout << "bkt_vector.size() = " << bkt_vector.size() << std::endl;


        if(read_size != total)initialize();
        else {

            for(uint32 i = 0; i < bkt_vector.size(); i++) delete bkt_vector[bkt_id];

            std::cout << "the SStar_sorter is empty.\n";
        }

    }

    pair_type get(){
        //read_size++;
        return bkt_vector[bkt_id]->get();
    }

    void next_remove(){
        read_size++;
        bkt_vector[bkt_id]->next_remove();

        if(bkt_vector[bkt_id]->is_eof()){
            delete bkt_vector[bkt_id];
            ++bkt_id;
            if( bkt_id < bkt_vector.size() )initialize();
            else return;
        }
    }

    bool is_empty(){

        return (total == 0) || (total == read_size);
    }

    void initialize(){

        while(bkt_vector[bkt_id]->is_empty()){

            delete bkt_vector[bkt_id];

            bkt_id++;

            if(bkt_id == bkt_vector.size())break;
        }

        if(bkt_id == bkt_vector.size()) return;

        std::cout << "bkt_id = " << bkt_id << std::endl;

        bkt_vector[bkt_id]->start_read();

        uint32 p_size = bkt_vector[bkt_id]->size();

        pair_type * p_buf = new pair_type[p_size];

        uint32 index = 0;

        while(!bkt_vector[bkt_id]->is_eof()){

            p_buf[index++] = bkt_vector[bkt_id]->get();

            bkt_vector[bkt_id]->next_remove();
        }

        delete bkt_vector[bkt_id];

        std::cout << "p_size = " << p_size << std::endl;

        std::sort(p_buf, p_buf + p_size, TupleAscCmp2<pair_type>());

        bkt_vector[bkt_id] = new MyVector<pair_type>();

        bkt_vector[bkt_id]->start_write();

        for(uint32 i = 0; i < p_size; i++){

            bkt_vector[bkt_id]->push_back(p_buf[i]);

            if(i < 100)std::cout << "p_buf[" << i << "] = < " << p_buf[i].first << ", " << p_buf[i].second << " >\n";

        }

        delete p_buf;

        bkt_vector[bkt_id]->start_read();

    }

    uint24 getLmsBkt(){

        return bkt_vector[bkt_id]->get().first;
    }

    uint32 get_bkt_id(){

        return bkt_id;
    }

};

///20190220
class sorted_LStar_reader_reverse {

private:

   typedef Pair<uint24,uint40> pair_type;

   int32 bkt_id; // signed integer

   uint64 total;

   uint64 read_size;

   std::vector< MyVector<Pair<uint24,uint40> > *> bkt_vector;

public:
    sorted_LStar_reader_reverse(std::vector< MyVector<Pair<uint24,uint40> > *> _bkt_vector):bkt_vector(_bkt_vector) {

        bkt_id = bkt_vector.size() - 1;

        total = 0;

        for(int32 i = bkt_vector.size() - 1; i >= 0; i--) {
            if(bkt_vector[bkt_id]->size()) {
                bkt_id = i;
                total += bkt_vector[i]->size();
                for( i = i - 1; i >= 0; i--) {
                    total += bkt_vector[i]->size();
                }
            }
        }
        read_size = 0;

        std::cout << "sorted_LStar_reader_reverse.total = " << total << std::endl;

        std::cout << "bkt_vector.size() = " << bkt_vector.size() << std::endl;

        if(read_size != total)initialize();
        else {

            for(uint32 i = 0; i < bkt_vector.size(); i++) delete bkt_vector[bkt_id];

            std::cout << "the sorted_LStar_reader_reverse is empty.\n";
        }

    }

    pair_type get_reverse(){

        return bkt_vector[bkt_id]->get_reverse();
    }

    void next_remove_reverse(){

        read_size++;

        bkt_vector[bkt_id]->next_remove_reverse();

        if(bkt_vector[bkt_id]->is_eof()){

            delete bkt_vector[bkt_id];

            --bkt_id;

            if( bkt_id >= 0 )initialize();
            else return;
        }
    }

    bool is_eof(){

        return (total == read_size);
    }

    void initialize(){

        while(bkt_vector[bkt_id]->is_empty()){

            delete bkt_vector[bkt_id];

            bkt_id--;

            if( bkt_id < 0 )break;
        }

        if(bkt_id < 0) return;

        bkt_vector[bkt_id]->start_read_reverse();

    }

    uint32 get_bkt_id(){

        return bkt_id;
    }

};

///20190220
class sorted_LStar_reader {

private:

   typedef Pair<uint24,uint40> pair_type;

   int32 bkt_id; // signed integer

   uint64 total;

   uint64 read_size;

   std::vector< MyVector<Pair<uint24,uint40> > *> bkt_vector;

public:
    sorted_LStar_reader(std::vector< MyVector<Pair<uint24,uint40> > *> _bkt_vector):bkt_vector(_bkt_vector) {

        bkt_id = 0;

        total = 0;

        for(uint32 i = 0; i < bkt_vector.size(); i++) total += bkt_vector[i]->size();

        std::cout << "sorted_LStar_reader.total = " << total << std::endl;

        std::cout << "bkt_vector.size() = " << bkt_vector.size() << std::endl;

        read_size = 0;

        if(read_size != total)initialize();
        else {

           // for(uint32 i = 0; i < bkt_vector.size(); i++) delete bkt_vector[bkt_id];

            std::cout << "the LStarReader is empty.\n";
        }

    }

    pair_type get(){

        return bkt_vector[bkt_id]->get();
    }

    void next(){

        read_size++;

        bkt_vector[bkt_id]->next();

        if(bkt_vector[bkt_id]->is_eof()){

            //delete bkt_vector[bkt_id];

            ++bkt_id;

            if( bkt_id < bkt_vector.size() )initialize();
            else return;
        }
    }

    bool is_eof(){

        return (total == read_size);
    }

    void initialize(){

        while(bkt_vector[bkt_id]->is_empty()){

            //delete bkt_vector[bkt_id];

            bkt_id++;

            if( bkt_id == bkt_vector.size())break;
        }

        if(bkt_id == bkt_vector.size()) return;

        bkt_vector[bkt_id]->start_read();

    }

};



/// 20190506
/**< High level, used for get the least compressed LMS character */
template < typename alphabet_type, typename compress_type = uint8>
class c_LMSHeapSorter {

private:

    typedef Pair<alphabet_type, compress_type> pair_lms_type ;
    typedef MyVector< pair_lms_type > lms_vector_type;
    std::vector< lms_vector_type * >  m_cLMS_seqs; /// pointer

    typedef Triple< alphabet_type, uint16, compress_type> triple_type;
    std::vector< triple_type > triple_vector;

    uint64 m_size;

    uint64 m_read_size;

public:
    c_LMSHeapSorter( std::vector<lms_vector_type * > & _cLMS_seqs )
        :m_cLMS_seqs(_cLMS_seqs)
        ,triple_vector(_cLMS_seqs.size())
        ,m_size(0)
        ,m_read_size(0) {

        std::cout << "m_lms_seqs size in SortLMS.h :" << triple_vector.size() << std::endl;


        for(uint32 i = 0; i < _cLMS_seqs.size(); i++) {

            m_cLMS_seqs[i]->start_read();

            m_size += m_cLMS_seqs[i]->size();

            /// if the rightmost block is singleton and the sentinel is the unique LMS char, the current cLMS_seqs will be omitted
            if(!m_cLMS_seqs[i]->is_eof()) {

                triple_vector[i].first = m_cLMS_seqs[i]->get().first;

                triple_vector[i].second = i;

                triple_vector[i].third = m_cLMS_seqs[i]->get().second;

            }
        }

/*         /// if the rightmost block is singleton block, the unique LMS char of the sentinel should be pop_back
 *         if( m_cLMS_seqs[0]->size() == 1 ){
 *
 *             /// the block 0 or leftmost block is singleton block
 *             if(triple_vector[0].first == 0 && triple_vector[0].third == 0){ /// \note reference to utility.h -> compute_singleton_bwt()
 *
 *                 triple_vector[0].first = triple_vector.back().first;
 *
 *                 triple_vector[0].second = triple_vector.back().second;
 *
 *                 triple_vector[0].third = triple_vector.back().third;
 *
 *                 triple_vector.pop_back();
 *
 *                 Logger::output_separator_line("The rightmost block is singleton block.");
 *
 *             }
 *
 *         }
 */

    }

    ~c_LMSHeapSorter(){

        for(uint32 i = 0; i < m_cLMS_seqs.size(); i++) {

            delete m_cLMS_seqs[i];
        }

        m_cLMS_seqs.clear();

        std::cout << "destructor of LMS_sorter.\n";
    }

    void initBkt(){

        std::make_heap(triple_vector.begin(), triple_vector.end(), TripleAscCmp2< triple_type >() );
    }


    alphabet_type get_bkt() const {

        return triple_vector[0].first;
    }

    uint16 get_block_id() const {

        return triple_vector[0].second;
    }

    uint64 get_bkt_num() const {

        return triple_vector[0].third;
    }

    void next(){

        assert(triple_vector[0].third);

        if(--triple_vector[0].third)return ; /// return ;

        ++m_read_size;

        m_cLMS_seqs[triple_vector[0].second]->next_remove();

        if(!is_empty()){

            if( !m_cLMS_seqs[triple_vector[0].second]->is_eof() ){


                if( m_cLMS_seqs[triple_vector[0].second]->get().first != triple_vector[0].first ){

                    triple_vector[0].first = m_cLMS_seqs[triple_vector[0].second]->get().first;

                    //triple_vector[0].second = triple_vector[0].second;

                    triple_vector[0].third = m_cLMS_seqs[triple_vector[0].second]->get().second;

                    initBkt();

                }else{

                    triple_vector[0].first = m_cLMS_seqs[triple_vector[0].second]->get().first;

                    //triple_vector[0].second = triple_vector[0].second;

                    triple_vector[0].third = m_cLMS_seqs[triple_vector[0].second]->get().second;
                }

            }else{

                std::pop_heap(triple_vector.begin(),triple_vector.end(),TripleAscCmp2< triple_type >() );

                triple_vector.pop_back();
            }
        }else {

            if(!triple_vector.empty()){

                triple_vector.pop_back();

            }

        }

    }

    bool is_empty() const {

        return m_size == m_read_size;
    }

    uint64 get_size(){

        return m_size ;
    }

};



#endif // _SORTLMS_H
