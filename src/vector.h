//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// Copyright (c) 2017, Sun Yat-sen University.
/// All rights reserved.
/// \file vector.h
/// \brief A self-defined external-memory vector designed for read/write operations I/O-efficiently.
///
/// The vector is implemented by primitive I/O functions
/// The vector provides interfaces for scanning elements rightward and leftward, but it doesn't support random access operations.
/// The vector supports two read modes: read-only and read-remove.
///
/// \author Yi Wu
/// \date 2017.7
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


#ifndef _VECTOR_H
#define _VECTOR_H


#include <fstream>
#include <cstdio>
#include <string>
#include <bitset>

#include <thread>
#include <chrono>
#include <functional>
#include <atomic>
#include <future>
#include <mutex>
#include <fstream>


#include "common.h"
#include "logger.h"
#include "tuple.h"
//#include "bitArray.h"

#define STATISTICS_COLLECTION


//#define THREAD

#ifndef THREAD

/// \brief Definition of a virtual vector consisting of one or multiple physical vectors.
///
template<typename element_type>
class MyVector
{

private:

    uint64 m_size; ///< number of elements in the virtual vector

    uint64 m_read; ///< number of elements already read from the virtual vector

private:

    /// \brief self-defined RAM buffer
    ///
    struct MyBuf
    {

    private:

        const uint32 m_capacity; ///< the capacity of the buffer, specified by VEC_BUF_RAM in common.h

        element_type *m_data; ///< handler to the payload of the buffer

        uint32 m_size;  ///< number of elements in the buffer

        uint32 m_read; ///< number of elements read from the buffer


    public:

        /// \brief ctor
        ///
        MyBuf() : m_capacity(VEC_BUF_RAM / sizeof(element_type))
        {

            m_data = new element_type[m_capacity]; // allocate RAM space for buffer
        }

        /// \brief dtor
        ///
        ~MyBuf()
        {

            delete[] m_data;

            m_data = nullptr;
        }

        /// \brief read a data block from the file
        ///
        /// \param _file file handler
        /// \param _offset offset
        /// \param _num number of elements to be read
        void read_block(FILE*& _file, const long _offset, const uint32 _num)
        {

            fseek(_file, _offset * sizeof(element_type), SEEK_SET); // seek read position
            double bg = Timer::get_wall_time();
            fread(m_data, sizeof(element_type), _num, _file); // sequentially read
            Timer::add_read_time(Timer::get_wall_time() - bg);
            Logger::addIV(_num * sizeof(element_type));

            m_size = _num;
        }

        /// \brief start read elements from the buffer
        ///
        void start_read()
        {

            m_read = 0;
        }

        /// \brief check if all the elements are read already
        ///
        bool is_eof() const
        {

            return m_size == m_read;
        }

        /// \brief get an element from the buffer (forwardly)
        ///
        /// \note check if eof before calling this function
        const element_type& get() const
        {

            return m_data[m_read];
        }

        /// \brief move forwardly
        ///
        void next()
        {

            ++m_read; // move to next element
        }

        /// \brief get an element from the buffer (reversely)
        ///
        const element_type& get_reverse() const
        {

            return m_data[m_size - 1 - m_read]; // e.g., m_read = 0, then read m_data[m_size - 1] from the buffer
        }

        /// \brief move reversely
        ///
        void next_reverse()
        {

            ++m_read; // call next()
        }


        /// \brief put elements into the buffer
        void start_write()
        {

            m_size = 0;
        }

        /// \brief check if the buffer is already full
        ///
        bool full() const
        {

            return m_size == m_capacity;
        }

        /// \brief put an element to the buffer
        ///
        /// \note check if full() before calling this funciton
        void put(const element_type& _elem)
        {

            m_data[m_size++] = _elem;
        }

        /// \brief write a data block to the file
        ///
        /// \note sequentially write into the block
        void write_block(FILE*& _file)
        {

            if(!_file)
            {

                std::cout << " _file is nullptr." << std::endl;

                Logger::output_error(__FILE__,__LINE__);

            }
            double bg = Timer::get_wall_time();

            fwrite(m_data, sizeof(element_type), m_size, _file);

            Timer::add_write_time(Timer::get_wall_time()-bg);
            Logger::addPDU(m_size * sizeof(element_type));

            Logger::addOV(m_size * sizeof(element_type));

        }


        /// \brief check if empty
        ///
        bool empty() const
        {

            return m_size == 0;
        }

        /// \brief return capacity
        ///
        const uint32& capacity() const
        {

            return m_capacity;
        }

        const uint32& size() const
        {

            return m_size;
        }
    };

    /// \brief self-define phyiscal vector
    ///
    struct MyPhiVector
    {

    private:

        std::string m_fname; ///< name of the file associated with the vector

        FILE* m_file; ///< handler

        const uint32 m_capacity; ///< capacity of the vector

        uint32 m_size; ///< number of elements in the vector

        uint32 m_read; ///< number of elements already read from the vector

        MyBuf*& m_buf; ///< a RAM buffer for facilitating I/O operations on the vectorcal


    public:

        /// \brief ctor
        ///
        MyPhiVector(MyBuf*&_buf) : m_capacity(PHI_VEC_EM / sizeof(element_type)), m_buf(_buf), m_file(nullptr)
        {

            //m_fname = "tmp_dsais1n_" + std::to_string(global_file_idx) + ".dat";
            m_fname = "tmp_file_" + random_string_hash() + ".dat";

            ++global_file_idx; // plus one each time to keep unique
        }


        ~MyPhiVector()
        {

            if(m_file){

                //std::cout << "close I/O stream.\n";

                fclose(m_file);
            }

            FILE * fp = fopen(m_fname.c_str(), "rb");

            if(fp){

                Logger::delPDU(m_size * sizeof(element_type)); /// error-prone

                fclose(fp);
            }


            std::remove(m_fname.c_str());

        }

        /// \brief return a random string for the file name
        std::string random_string_hash()
        {

//            const uint64_t hash = (uint64_t)rand() * RAND_MAX + rand();
//
//            return std::to_string(hash);

            // 使用现代随机数生成器，避免rand()的缺陷
            static std::mt19937_64 rng(std::chrono::steady_clock::now().time_since_epoch().count());
            std::uniform_int_distribution<uint64_t> dist;
            uint64_t hash = dist(rng);

            // 转十六进制，降低碰撞概率
            std::stringstream ss;
            ss << std::hex << hash;
            return ss.str();


        }

        /// \brief prepare for writing
        ///
        void start_write()
        {

            if(m_file){

                std::cout << "close I/O stream before start wirte.\n";

                fclose(m_file);
            }

            double bg = Timer::get_wall_time();
            m_file = fopen(m_fname.c_str(), "wb");
            Timer::add_write_time(Timer::get_wall_time()-bg);

            m_size = 0;

            m_buf->start_write();
        }

        /// \brief check if the vector is full
        ///
        bool full() const
        {

            return m_size == m_capacity;
        }

        /// \brief put an element into the vector
        ///
        /// \note check if full before calling this function
        void push_back(const element_type& _elem)
        {

            if (m_buf->full())   // buffer is full
            {

                m_buf->write_block(m_file);

                m_buf->start_write(); // clear up the buffer
            }

            m_buf->put(_elem), ++m_size; /// put an element into the buffer
        }

        /// \brief ending writing process
        ///
        /// \note call the function after the vector is full
        void end_write()
        {

            if (!m_buf->empty())   // flush the remaining elements in the buffer
            {

                m_buf->write_block(m_file);

            }

            fclose(m_file); // close file

            m_file = nullptr;
        }

        /// \brief prepare for reading forwardly
        ///
        void start_read()
        {

            if(m_file){

                std::cout << "close I/O stream before start read.\n";

                fclose(m_file);
            }

            double bg = Timer::get_wall_time();
            m_file = fopen(m_fname.c_str(), "rb");
            Timer::add_read_time(Timer::get_wall_time() - bg);
            m_read = 0;

            //clock_t beg = clock();
            m_buf->read_block(m_file, m_read, std::min(m_size - m_read, m_buf->capacity()));

            m_buf->start_read();
        }

        /// \brief check if no more to read
        ///
        bool is_eof()
        {

            return m_read == m_size;
        }

        /// \brief prepare for reading reversely
        ///
        void start_read_reverse()
        {

            if(m_file){

                std::cout << "close I/O stream before start read reverse.\n";

                fclose(m_file);
            }

            m_file = fopen(m_fname.c_str(), "rb");

            m_read = 0;

            m_buf->read_block(m_file, m_size - m_read - std::min(m_size - m_read, m_buf->capacity()), std::min(m_size - m_read, m_buf->capacity()));

            m_buf->start_read();
        }

        /// \brief get an element from the vector forwardly
        ///
        /// \note check if eof before calling the function
        const element_type& get() const
        {

            return m_buf->get();
        }

        /// \brief get an element from the vector reversely
        ///
        /// \note check if eof before calling the function
        const element_type& get_reverse() const
        {

            return m_buf->get_reverse();
        }


        /// \brief move to next element forwardly
        ///
        void next()
        {

            ++m_read, m_buf->next();

            if (m_buf->is_eof())   // the elements in the buffer are already processed
            {

                m_buf->read_block(m_file, m_read, std::min(m_size - m_read, m_buf->capacity()));

                m_buf->start_read();
            }
        }

        /// \brief move to next element reversely
        ///
        void next_reverse()
        {

            ++m_read, m_buf->next_reverse();

            if (m_buf->is_eof())   // the elements in the buffer are already processed
            {

                m_buf->read_block(m_file, m_size - m_read - std::min(m_size - m_read, m_buf->capacity()), std::min(m_size - m_read, m_buf->capacity()));

                m_buf->start_read();
            }
        }

        /// \brief finish reading
        ///
        /// \note call the function after is_eof() == true
        void end_read()
        {

            fclose(m_file);

            m_file = nullptr;
        }

        /// \brief remove the physical file
        ///
        void remove_file()
        {

            std::remove(m_fname.c_str());

            Logger::delPDU(m_size * sizeof(element_type));

        }

        /// \brief get m_size
        ///
        const uint32& size() const
        {

            return m_size;
        }

        const uint32& capacity() const
        {

            return m_capacity;
        }

        std::string get_name()
        {

            return m_fname;
        }

        void set_size( uint32 size)
        {
            m_size = size;
        }

    };

private:

    MyBuf* m_buf; ///< handler to RAM buffer

    std::vector<MyPhiVector*> m_phi_vectors; ///< handlers to physical vectors

    uint32 m_phi_vector_write_idx; ///< index of the vector being written

    uint32 m_phi_vector_read_idx; ///< index of the vector being read

    bool m_flag; ///< set false if changing from writting to reading

public:

    /// \brief ctor
    ///
    MyVector()
    {

        m_buf = new MyBuf(); // create the RAM buffer

        m_phi_vectors.clear();

        // start_write();

        m_flag = false;

        m_size = 0;

        m_read = 0;
    }


    // write by block manner
    MyVector(element_type * T, uint32 size)
    {

        m_buf = new MyBuf();

        m_phi_vectors.clear();

        // compute the block number
        uint32 block_number = PHI_VEC_EM / sizeof(element_type);

        uint32 loop = (size % block_number) == 0 ? (size / block_number) : (size / block_number) + 1;

        FILE * fp;

        for (uint32 i = 0; i < loop; i++)
        {

            MyPhiVector * phy_vec = new MyPhiVector(m_buf);

            fp = fopen((phy_vec->get_name()).c_str(),"wb");

            //uint32 write_size = (i + 1) * block_number < size ? block_number : size % block_number;

            uint32 write_size;

            if((i+1)*block_number < size)
            {
                write_size = block_number;
            }
            else{

                write_size = size - i * block_number;
            }

            //std::cout << "write_size = " << write_size << std::endl;

            if(write_size == 0)Logger::output_error(__FILE__,__LINE__);

            fwrite(T+i*block_number,sizeof(element_type),write_size,fp);

            //std::cout << " T = " << (uint64)(T + i*block_number) << std::endl;

            fclose(fp);

            phy_vec->set_size(write_size);

            m_phi_vectors.push_back(phy_vec);

            m_size += write_size;
        }


        m_flag = true;
    }

    /// \brief dtor
    ///
    ~MyVector()
    {

        delete m_buf;

        m_buf = nullptr;

        for (uint32 i = 0; i < m_phi_vectors.size(); i++)
        {

            // std::cout << " delete file name : " << m_phi_vectors[i]->get_name() << " ." << std::endl;

            if(m_phi_vectors[i])
            {

                delete m_phi_vectors[i];

                m_phi_vectors[i]=nullptr;
            }

        }
    }


    /// \brief start writting
    ///
    void start_write()
    {

        m_size = 0;

        m_phi_vectors.push_back(new MyPhiVector(m_buf)); // create a physical vector at the beginning

        m_phi_vector_write_idx = 0;

        m_phi_vectors[m_phi_vector_write_idx]->start_write();
    }

    /// \brief check if the vector is empty
    ///
    bool is_empty()
    {

        return m_size == 0;
    }

    /// \brief append an element to the end of the virtual vector
    ///
    /// \note
    void push_back(const element_type & _value)
    {

        if (m_phi_vectors[m_phi_vector_write_idx]->full())
        {

            m_phi_vectors[m_phi_vector_write_idx]->end_write();

            m_phi_vectors.push_back(new MyPhiVector(m_buf));

            ++m_phi_vector_write_idx;

            m_phi_vectors[m_phi_vector_write_idx]->start_write();
        }

        m_phi_vectors[m_phi_vector_write_idx]->push_back(_value), ++m_size;

        //	std::cerr << "phi idx: " << m_phi_vector_write_idx << " m_size: " << m_size << std::endl;
    }

    /// \brief finish writing the virtual vector
    ///
    void end_write()
    {

        //m_flag = true;

        m_phi_vectors[m_phi_vector_write_idx]->end_write(); // finish writing the last physical vector
    }

    /// \brief start reading the virtual vector forwardly
    ///
    void start_read()
    {

        if (m_flag == false)
        {

            end_write(); // finish writing the final vector

            m_flag = true;
        }

        m_read = 0;

        m_phi_vector_read_idx = 0;

        m_phi_vectors[m_phi_vector_read_idx]->start_read();
    }

    /// \brief start reading the virtual vector reversely
    ///
    void start_read_reverse()
    {

        if (m_flag == false)
        {

            end_write();

            m_flag = true;
        }

        m_read = 0;

        m_phi_vector_read_idx = m_phi_vectors.size() - 1;

        m_phi_vectors[m_phi_vector_read_idx]->start_read_reverse();
    }

    /// \brief check if all the elements are already read
    ///
    bool is_eof()
    {

        return m_size == m_read;
    }

    /// \brief get an element from the virtual vector forwardly
    ///
    const element_type& get() const
    {

        return m_phi_vectors[m_phi_vector_read_idx]->get();
    }

    /// \brief get an element from the virtual vector reversely
    ///
    const element_type& get_reverse() const
    {

        return m_phi_vectors[m_phi_vector_read_idx]->get_reverse();
    }

    /// \brief move to next + remove
    ///
    void next_remove()
    {

        ++m_read, m_phi_vectors[m_phi_vector_read_idx]->next();

        if (m_phi_vectors[m_phi_vector_read_idx]->is_eof())
        {

            m_phi_vectors[m_phi_vector_read_idx]->end_read();

            m_phi_vectors[m_phi_vector_read_idx]->remove_file();

            delete m_phi_vectors[m_phi_vector_read_idx];

            m_phi_vectors[m_phi_vector_read_idx] = nullptr;

            if (!is_eof())
            {

                ++m_phi_vector_read_idx;

                m_phi_vectors[m_phi_vector_read_idx]->start_read();
            }
        }
    }

    /// \brief move to next element
    ///
    void next()
    {

        ++m_read, m_phi_vectors[m_phi_vector_read_idx]->next();

        if (m_phi_vectors[m_phi_vector_read_idx]->is_eof())
        {

            m_phi_vectors[m_phi_vector_read_idx]->end_read();

            if (!is_eof())
            {

                ++m_phi_vector_read_idx;

                m_phi_vectors[m_phi_vector_read_idx]->start_read();
            }
        }
    }

    /// \brief move to next + remove +reverse
    ///
    void next_remove_reverse()
    {

        ++m_read, m_phi_vectors[m_phi_vector_read_idx]->next_reverse();

        if (m_phi_vectors[m_phi_vector_read_idx]->is_eof())
        {

            m_phi_vectors[m_phi_vector_read_idx]->end_read();

            m_phi_vectors[m_phi_vector_read_idx]->remove_file();

            delete m_phi_vectors[m_phi_vector_read_idx];

            m_phi_vectors[m_phi_vector_read_idx] = nullptr;

            if (!is_eof())
            {

                --m_phi_vector_read_idx;

                m_phi_vectors[m_phi_vector_read_idx]->start_read_reverse();
            }
        }
    }

    /// \brief next + reverse
    ///
    void next_reverse()
    {

        ++m_read, m_phi_vectors[m_phi_vector_read_idx]->next_reverse();

        if (m_phi_vectors[m_phi_vector_read_idx]->is_eof())
        {

            m_phi_vectors[m_phi_vector_read_idx]->end_read();

            if (!is_eof())
            {

                --m_phi_vector_read_idx;

                m_phi_vectors[m_phi_vector_read_idx]->start_read_reverse();
            }
        }
    }

    /// \brief report status
    ///
    void report()
    {

        std::cerr << "number of phi vecs: " << m_phi_vectors.size() << std::endl;
    }

    /// \brief return size
    ///
    uint64 size() const
    {

        return m_size;
    }

    uint64 read_size() const
    {

        return m_read;
    }
};


#else

template<typename element_type>
void write_buf(FILE * _fp, element_type * _buf, uint64 _size)//, std::promise<uint32> *& _prom)
{
    //std::lock_guard<std::mutex> lock(io_mu);
    //io_mu.lock();
    if(!_fp)
    {
        std::cout << "File : " << __FILE__ << " , Line number : " << __LINE__ << " , fp is null.\n";
        std::cin.get();
    }

    if(_size)
    {

        fwrite(_buf,sizeof(element_type),_size,_fp);

        fflush(_fp);
    }

    /*std::cout << " thread run, buf address = " << (uint64)_buf << std::endl;
    std::cout << " _fp = " << _fp << std::endl;
    for(int i = 0; i < _size; ++i){
        std::cout << "buf[" << i << "] = " << _buf[i]  - 0 << std::endl;
    }*/

    //_prom->set_value(100);

    //io_mu.unlock();

    //std::this_thread::sleep_for(std::chrono::milliseconds(1000));

    //  std::cout << "write thread ID : " << std::this_thread::get_id() << " sleep 1 second.\n";

    return;
}
template<typename element_type>
void read_buf(FILE * _fp, element_type * _buf, uint32 _size, uint64 _offset)//, std::promise<uint32> *& _prom)
{
    // std::lock_guard<std::mutex> lock(io_mu);
    //io_mu.lock();
    if(!_fp)
    {
        std::cout << "File : " << __FILE__ << " , Line number : " << __LINE__ << " , fp is null.\n";
        std::cin.get();
    }

    if(_size)
    {

        fseek(_fp, _offset * sizeof(element_type), SEEK_SET); // seek read position

        fread(_buf, sizeof(element_type), _size, _fp); // sequentially read
    }

    //_prom->set_value(200);

    //io_mu.unlock();

    //std::this_thread::sleep_for(std::chrono::milliseconds(1000));

    // std::cout << "read thread ID : " << std::this_thread::get_id() << " sleep 1 second.\n";


    return;
}


/// \brief Definition of a virtual vector consisting of one or multiple physical vectors.
///
template<typename element_type>
class MyVector
{

private:

    uint64 m_size; ///< number of elements in the virtual vector

    uint64 m_read; ///< number of elements already read from the virtual vector

private:

    /// \brief self-defined RAM buffer
    ///
    struct MyBuf
    {

    private:

        const uint32 m_capacity; ///< the capacity of the buffer, specified by VEC_BUF_RAM in common.h

        //element_type *m_data; ///< handler to the payload of the buffer

        element_type * m_data[2]; /// handler to the payload of the buffer

        //uint8 is_ready[2]; /// 标记buf的状态，1 = true, 0 = false

        uint32 m_size;  ///< number of elements in the buffer

        uint32 m_read; ///< number of elements read from the buffer

        uint8 cur_buf_id; ///< only 0 or 1, means m_data[0] or m_data[1], in memory

        //std::future<uint32> *fut;

        // std::promise<uint32> *prom;

        std::thread * pt_thread;

        // bool write_mark;///< the first write mark;

    public:

        /// \brief constructor
        ///
        MyBuf() : m_capacity(VEC_BUF_RAM / sizeof(element_type))
        {

            cur_buf_id = 0;
//            fp_switch = true;

            m_data[0] = new element_type[m_capacity]; // allocate RAM space for buffer
            m_data[1] = new element_type[m_capacity]; // allocate RAM space for buffer

            //  fut = new std::future<uint32>();
            //  prom = new std::promise<uint32>();
            // *fut = prom->get_future();
            //std::async(std::launch::async, func, std::ref(prom));
            // prom->set_value(100);

            pt_thread = new std::thread([]()
            {
                return;
            });
        }

        /// \brief dtor
        ///
        ~MyBuf()
        {

            delete m_data[0];
            delete m_data[1];

            m_data[0] = nullptr;
            m_data[1] = nullptr;

            /*delete fut;
            delete prom;*/

            pt_thread->join();
            delete pt_thread;
        }

        void pre_read(FILE *& _file, uint64 _size, bool _is_reverse)
        {



            uint32 read_size = 0;

            uint64 offset = 0;

            if(_is_reverse)
            {
                read_size = (_size > m_capacity) ? m_capacity : _size;
                offset = _size - read_size;
            }
            else
            {
                read_size = (_size > m_capacity) ? m_capacity : _size;
                offset = 0;
            }

            cur_buf_id = 0;

            pt_thread->join();

            delete pt_thread;

            pt_thread = new std::thread(read_buf<element_type>, _file, m_data[1], read_size, offset);

            //std::async(std::launch::async, read_buf<element_type>, _fp, m_data[1], read_size, offset,std::ref(prom));
        }


        /// \brief read a data block from the file
        /// \param _file file handler
        /// \param _size the total size of the PHI_VEC_EM
        /// \param _offset offset
        /// \param _num number of elements to be read
        void read_block(FILE*& _file, const uint32 _size, const uint64 _offset, const uint32 _num, bool _is_reverse)
        {

            uint32 read_size = 0;

            uint64 offset = 0;

            if(_is_reverse)
            {
                offset = (_offset >= m_capacity) ? (_offset - m_capacity) : 0;
                read_size = (offset != 0) ? m_capacity : _offset;
            }
            else
            {
                offset = (_offset + m_capacity) <= _size ? (_offset + m_capacity) : _size;
                read_size = (offset + m_capacity) <= _size ? m_capacity : (_size - offset);
            }


            pt_thread->join();

            delete pt_thread;

            pt_thread = new std::thread(read_buf<element_type>, _file, m_data[cur_buf_id], read_size, offset);

            switch_buf_id();

            Logger::addIV(_num * sizeof(element_type));

            m_size = _num;
        }

        /// \brief write a data block to the file
        ///
        /// \note sequentially write into the block
        void write_block(FILE*& _file)
        {

            if(!_file)
            {
                std::cout << " _file is nullptr." << std::endl;

                std::cin.get();
            }
            ///step 1: get the cur_buf_id,
            ///step 2: get the other id, test thread is closed,
            /// if closed, create the new thread to write.
            /*uint32 rtv = fut->get();

            if(rtv == 100)
            {
                delete fut;
                delete prom;
                fut = new std::future<uint32>();
                prom = new std::promise<uint32>();
                *fut = prom->get_future();
                std::async(std::launch::async, write_buf<element_type>, _file, m_data[cur_buf_id], m_size, std::ref(prom));
            }
            else
            {
                std::cout << "File : " << __FILE__ << " , Line number : " << __LINE__ << " , thread exception.\n";
                std::cin.get();
            }*/
            pt_thread->join();

            delete pt_thread;

            pt_thread = new std::thread(write_buf<element_type>, _file, m_data[cur_buf_id], m_size);

            switch_buf_id();

            //fwrite(m_data[cur_buf_id], sizeof(element_type), m_size, _file);

            Logger::addPDU(m_size * sizeof(element_type));

            Logger::addOV(m_size * sizeof(element_type));

        }

        void end_write()
        {

            pt_thread->join();

            delete pt_thread;

            pt_thread = new std::thread([]()
            {
                return;
            });

            cur_buf_id = 0;
        }

        void end_read()
        {

            pt_thread->join();

            delete pt_thread;

            pt_thread = new std::thread([]()
            {
                return;
            });

        }

        void switch_buf_id()
        {
            if(cur_buf_id == 1)cur_buf_id = 0;
            else cur_buf_id = 1;
        }

        /// \brief start read elements from the buffer
        ///
        void start_read()
        {

            m_read = 0;
        }

        /// \brief check if all the elements are read already
        ///
        bool is_eof() const
        {

            return m_size == m_read;
        }

        /// \brief get an element from the buffer (forwardly)
        ///
        /// \note check if eof before calling this function
        const element_type& get() const
        {

            return m_data[cur_buf_id][m_read];
        }

        /// \brief move forwardly
        ///
        void next()
        {

            ++m_read; // move to next element
        }

        /// \brief get an element from the buffer (reversely)
        ///
        const element_type& get_reverse() const
        {

            return m_data[cur_buf_id][m_size - 1 - m_read]; // e.g., m_read = 0, then read m_data[m_size - 1] from the buffer
        }

        /// \brief move reversely
        ///
        void next_reverse()
        {

            ++m_read; // call next()
        }


        /// \brief put elements into the buffer
        void start_write()
        {

            m_size = 0;
        }

        /// \brief check if the buffer is already full
        ///
        bool full() const
        {

            return m_size == m_capacity;
        }

        /// \brief put an element to the buffer
        ///
        /// \note check if full() before calling this funciton
        void put(const element_type& _elem)
        {

            m_data[cur_buf_id][m_size++] = _elem;
        }


        /// \brief check if empty
        ///
        bool empty() const
        {

            return m_size == 0;
        }

        /// \brief return capacity
        ///
        const uint32& capacity() const
        {

            return m_capacity;
        }

        const uint32& size() const
        {

            return m_size;
        }


    };

    /// \brief self-define phyiscal vector
    ///
    struct MyPhiVector
    {

    private:

        std::string m_fname; ///< name of the file associated with the vector

        FILE* m_file; ///< handler

        const uint32 m_capacity; ///< capacity of the vector

        uint32 m_size; ///< number of elements in the vector

        uint32 m_read; ///< number of elements already read from the vector

        MyBuf*& m_buf; ///< a RAM buffer for facilitating I/O operations on the vector


    public:

        /// \brief ctor
        ///
        MyPhiVector(MyBuf*&_buf) : m_capacity(PHI_VEC_EM / sizeof(element_type)), m_buf(_buf)
        {

            m_fname = "tmp_dsais1n_" + std::to_string(global_file_idx) + ".dat";

//            if(global_file_idx == 287){
//
//                std::cout << "global_file_idx == 287.\n";
//                std::cin.get();
//            }

            ++global_file_idx; // plus one each time to keep unique
        }


        ~MyPhiVector()
        {

            if(fopen(m_fname.c_str(), "rb"))Logger::delPDU(m_size * sizeof(element_type));

            std::remove(m_fname.c_str());
        }
        /// \brief prepare for writing
        ///
        void start_write()
        {

            m_file = fopen(m_fname.c_str(), "wb");

            m_size = 0;

            m_buf->start_write();
        }

        /// \brief check if the vector is full
        ///
        bool full() const
        {

            return m_size == m_capacity;
        }

        /// \brief put an element into the vector
        ///
        /// \note check if full before calling this function
        void push_back(const element_type& _elem)
        {

            if (m_buf->full())   // buffer is full
            {

                m_buf->write_block(m_file);

                m_buf->start_write(); // clear up the buffer
            }

            m_buf->put(_elem), ++m_size; /// put an element into the buffer
        }

        /// \brief ending writing process
        ///
        /// \note call the function after the vector is full
        void end_write()
        {
            // std::cout << " end_write() in physical block : m_file = " << (uint64)m_file << std::endl;

            if(!m_file)
            {
                std::cout << __FILE__ << " , line number : " << __LINE__ << " , m_file is empty()";
                std::cin.get();
            }


            if (!m_buf->empty())   // flush the remaining elements in the buffer
            {

                m_buf->write_block(m_file);

            }

            m_buf->end_write();

            fclose(m_file); // close file

        }

        /// \brief prepare for reading forwardly
        ///
        void start_read()
        {


            m_file = fopen(m_fname.c_str(), "rb");

            m_read = 0;

            /// for thread
            m_buf->pre_read(m_file, m_size,false);

            m_buf->read_block(m_file, m_size, m_read, std::min(m_size - m_read, m_buf->capacity()), false);

            m_buf->start_read();
        }

        /// \brief check if no more to read
        ///
        bool is_eof()
        {

            return m_read == m_size;
        }

        /// \brief prepare for reading reversely
        ///
        void start_read_reverse()
        {

            m_file = fopen(m_fname.c_str(), "rb");

            m_read = 0;

            /// for thread
            m_buf->pre_read(m_file, m_size,true);

            m_buf->read_block(m_file, m_size, m_size - m_read - std::min(m_size - m_read, m_buf->capacity()), std::min(m_size - m_read, m_buf->capacity()),true);

            m_buf->start_read();
        }

        /// \brief get an element from the vector forwardly
        ///
        /// \note check if eof before calling the function
        const element_type& get() const
        {

            return m_buf->get();
        }

        /// \brief get an element from the vector reversely
        ///
        /// \note check if eof before calling the function
        const element_type& get_reverse() const
        {

            return m_buf->get_reverse();
        }


        /// \brief move to next element forwardly
        ///
        void next()
        {

            ++m_read, m_buf->next();

            if (m_buf->is_eof())   // the elements in the buffer are already processed
            {

                m_buf->read_block(m_file, m_size, m_read, std::min(m_size - m_read, m_buf->capacity()), false);

                m_buf->start_read();
            }
        }

        /// \brief move to next element reversely
        ///
        void next_reverse()
        {

            ++m_read, m_buf->next_reverse();

            if (m_buf->is_eof())   // the elements in the buffer are already processed
            {

                m_buf->read_block(m_file, m_size, m_size - m_read - std::min(m_size - m_read, m_buf->capacity()), std::min(m_size - m_read, m_buf->capacity()), true);

                m_buf->start_read();
            }
        }

        /// \brief finish reading
        ///
        /// \note call the function after is_eof() == true
        void end_read()
        {
            m_buf->end_read();
            fclose(m_file);
        }

        /// \brief remove the physical file
        ///
        void remove_file()
        {

            std::remove(m_fname.c_str());

            Logger::delPDU(m_size * sizeof(element_type));

        }

        /// \brief get m_size
        ///
        const uint32& size() const
        {

            return m_size;
        }

        const uint32& capacity() const
        {

            return m_capacity;
        }

        std::string get_name()
        {

            return m_fname;
        }

        void set_size( uint32 size)
        {
            m_size = size;
        }

    };

private:

    MyBuf* m_buf; ///< handler to RAM buffer

    std::vector<MyPhiVector*> m_phi_vectors; ///< handlers to physical vectors

    uint32 m_phi_vector_write_idx; ///< index of the vector being written

    uint32 m_phi_vector_read_idx; ///< index of the vector being read

    bool m_flag; ///< set false if changing from writting to reading

public:

    /// \brief ctor
    ///
    MyVector()
    {

        m_buf = new MyBuf(); // create the RAM buffer

        m_phi_vectors.clear();

        // start_write();

        m_flag = false;

        m_size = 0;

        m_read = 0;
    }


    /// \brief dtor
    ///
    ~MyVector()
    {

        delete m_buf;

        m_buf = nullptr;

        for (uint32 i = 0; i < m_phi_vectors.size(); i++)
        {

            // std::cout << " delete file name : " << m_phi_vectors[i]->get_name() << " ." << std::endl;

            if(m_phi_vectors[i])
            {

                delete m_phi_vectors[i];

                m_phi_vectors[i]=nullptr;
            }

        }
    }


    /// \brief start writting
    ///
    void start_write()
    {

        m_size = 0;

        m_phi_vectors.push_back(new MyPhiVector(m_buf)); // create a physical vector at the beginning

        m_phi_vector_write_idx = 0;

        m_phi_vectors[m_phi_vector_write_idx]->start_write();
    }

    /// \brief check if the vector is empty
    ///
    bool is_empty()
    {

        return m_size == 0;
    }

    /// \brief append an element to the end of the virtual vector
    ///
    /// \note
    void push_back(const element_type & _value)
    {

        if (m_phi_vectors[m_phi_vector_write_idx]->full())
        {

            m_phi_vectors[m_phi_vector_write_idx]->end_write();

            m_phi_vectors.push_back(new MyPhiVector(m_buf));

            ++m_phi_vector_write_idx;

            m_phi_vectors[m_phi_vector_write_idx]->start_write();
        }

        m_phi_vectors[m_phi_vector_write_idx]->push_back(_value), ++m_size;

        //	std::cerr << "phi idx: " << m_phi_vector_write_idx << " m_size: " << m_size << std::endl;
    }

    /// \brief finish writing the virtual vector
    ///
    void end_write()
    {

        //m_flag = true;

        m_phi_vectors[m_phi_vector_write_idx]->end_write(); // finish writing the last physical vector
    }

    /// \brief start reading the virtual vector forwardly
    ///
    void start_read()
    {


        if (m_flag == false)
        {

            end_write(); // finish writing the final vector

            m_flag = true;

            //std::cout << " end_write closed.\n";
        }

        m_read = 0;

        m_phi_vector_read_idx = 0;

        m_phi_vectors[m_phi_vector_read_idx]->start_read();
    }

    /// \brief start reading the virtual vector reversely
    ///
    void start_read_reverse()
    {

        if (m_flag == false)
        {

            end_write();

            m_flag = true;
        }

        m_read = 0;

        m_phi_vector_read_idx = m_phi_vectors.size() - 1;

        m_phi_vectors[m_phi_vector_read_idx]->start_read_reverse();
    }

    /// \brief check if all the elements are already read
    ///
    bool is_eof()
    {

        return m_size == m_read;
    }

    /// \brief get an element from the virtual vector forwardly
    ///
    const element_type& get() const
    {

        return m_phi_vectors[m_phi_vector_read_idx]->get();
    }

    /// \brief get an element from the virtual vector reversely
    ///
    const element_type& get_reverse() const
    {

        return m_phi_vectors[m_phi_vector_read_idx]->get_reverse();
    }

    /// \brief move to next + remove
    ///
    void next_remove()
    {

        ++m_read, m_phi_vectors[m_phi_vector_read_idx]->next();

        if (m_phi_vectors[m_phi_vector_read_idx]->is_eof())
        {

            m_phi_vectors[m_phi_vector_read_idx]->end_read();

            m_phi_vectors[m_phi_vector_read_idx]->remove_file();

            delete m_phi_vectors[m_phi_vector_read_idx];

            m_phi_vectors[m_phi_vector_read_idx] = nullptr;

            if (!is_eof())
            {

                ++m_phi_vector_read_idx;

                m_phi_vectors[m_phi_vector_read_idx]->start_read();
            }
        }
    }

    /// \brief move to next element
    ///
    void next()
    {

        ++m_read, m_phi_vectors[m_phi_vector_read_idx]->next();

        if (m_phi_vectors[m_phi_vector_read_idx]->is_eof())
        {

            m_phi_vectors[m_phi_vector_read_idx]->end_read();

            if (!is_eof())
            {

                ++m_phi_vector_read_idx;

                m_phi_vectors[m_phi_vector_read_idx]->start_read();
            }
        }
    }

    /// \brief move to next + remove +reverse
    ///
    void next_remove_reverse()
    {

        ++m_read, m_phi_vectors[m_phi_vector_read_idx]->next_reverse();

        if (m_phi_vectors[m_phi_vector_read_idx]->is_eof())
        {

            m_phi_vectors[m_phi_vector_read_idx]->end_read();

            m_phi_vectors[m_phi_vector_read_idx]->remove_file();

            delete m_phi_vectors[m_phi_vector_read_idx];

            m_phi_vectors[m_phi_vector_read_idx] = nullptr;

            if (!is_eof())
            {

                --m_phi_vector_read_idx;

                m_phi_vectors[m_phi_vector_read_idx]->start_read_reverse();
            }
        }
    }

    /// \brief next + reverse
    ///
    void next_reverse()
    {

        ++m_read, m_phi_vectors[m_phi_vector_read_idx]->next_reverse();

        if (m_phi_vectors[m_phi_vector_read_idx]->is_eof())
        {

            m_phi_vectors[m_phi_vector_read_idx]->end_read();

            if (!is_eof())
            {

                --m_phi_vector_read_idx;

                m_phi_vectors[m_phi_vector_read_idx]->start_read_reverse();
            }
        }
    }

    /// \brief report status
    ///
    void report()
    {

        std::cerr << "number of phi vecs: " << m_phi_vectors.size() << std::endl;
    }

    /// \brief return size
    ///
    uint64 size() const
    {

        return m_size;
    }

    uint64 read_size() const
    {

        return m_read;
    }


};

#endif // THREAD



/// template < uint64, N = 64> , <uint32 ,N = 32> ,
template<typename element_type, int N = 64>
class BitVectorWapper
{
private:

    /**< 外存Vector对象 */
    MyVector<element_type> * m_vector;

    /**< N bit的二进制数 */
    std::bitset<N> m_buf;

    /**< 读取,写入位置 */
    uint8 m_pos;

    /**< bit 总数 */
    uint64 m_size;

    /**< 读取的bit 数 */
    uint64 m_read_size;

    bool is_write;
public:

    /**< 构造函数，新生产一个外存vector成员，其余数据成员职位零 */
    BitVectorWapper():m_vector(nullptr),m_buf(0),m_pos(0),m_size(0),m_read_size(0),is_write(false)
    {

        m_vector = new MyVector<element_type> ();
    }

    ///destructor
    ~BitVectorWapper()
    {

        if(m_vector)
        {

            delete m_vector;

            m_vector = nullptr;
        }
    }

    /**< 在写之前调用，将缓存N bit置为0 */
    void start_write()
    {

        m_vector->start_write();

        m_size = 0;

        m_read_size = 0;

        m_buf.reset();

        m_pos = 0;
    }

    /** \brief push_back( bool _t)
     *
     * \param bool
     * \param
     * \return void
     *
     */

    void push_back( bool _t)
    {

        if(_t)
            m_buf.set( N - m_pos - 1);
        ++m_size;
        ++m_pos;
        if( m_pos%N == 0 )
        {
            m_vector->push_back(m_buf.to_ullong());
            m_pos = 0;
            m_buf.reset();
        }
    }

    /** \brief 写入结束调用
     *
     * \param void
     * \param
     * \return void
     *
     */
    void write_end()   ///revised 20181225
    {

        if(!is_write)
        {

            if(m_pos)
            {

                m_vector->push_back(m_buf.to_ullong());
            }

            is_write = true;
        }

        m_buf.reset();

        m_pos = 0;
    }

    /** \brief 顺序读之前调用,在执行此函数之前需判断是否为空
     *
     * \param void
     * \param
     * \return void
     *
     */

    void start_read()
    {

        write_end();

        m_vector->start_read();

        m_buf.reset();

        m_buf = m_vector->get();

        m_vector->next();

        m_pos = 0;

        m_read_size = 0;
    }

    /** \brief 顺序一次性读之前调用,在执行此函数之前需判断是否为空
     *
     * \param void
     * \param
     * \return void
     *
     */

    void start_read_remove()
    {

        write_end();

        m_vector->start_read();

        m_buf.reset();

        m_buf = m_vector->get();

        m_vector->next_remove();

        m_pos = 0;

        m_read_size = 0;
    }


    /** \brief 逆序读之前调用，在执行此函数之前需判断是否为空
     *
     * \param void
     * \param
     * \return void
     *
     */

    void start_read_reverse()
    {

        write_end();

        m_vector->start_read_reverse();

        m_buf.reset();

        m_buf = m_vector->get_reverse();

        m_vector->next_reverse();

        m_pos = ( (m_size % N == 0) ? 0 : (N - (m_size % N)) );

        m_read_size = 0;
    }

    /** \brief 逆序一次性读之前调用，在执行此函数之前需判断是否为空
     *
     * \param void
     * \param
     * \return void
     *
     */

    void start_read_remove_reverse()
    {

        write_end();

        m_vector->start_read_reverse();

        m_buf.reset();

        m_buf = m_vector->get_reverse();

        m_vector->next_remove_reverse();

        m_pos = ( (m_size % N == 0) ? 0 : (N - (m_size % N)) )  ;

        m_read_size = 0;
    }


    /** \brief 取得当前类型(顺序),该方法不能用于自增, 否则就不能取数据多次
     *
     * \param void
     * \param
     * \return bool
     *
     */

    bool get()
    {

        //++m_read_size;

        return m_buf[ N - 1 - m_pos];
    }


    /** \brief 取得当前类型(逆序),get_reverse 方法不能用于自增, 否则就不能get()数据多次
     *
     * \param void
     * \param
     * \return bool
     *
     */

    bool get_reverse()
    {

        //++m_read_size;

        return m_buf[m_pos];
    }

    /** \brief 指向下一个读取位置(顺序)
     *
     * \param void
     * \param
     * \return void
     *
     */

    void next()
    {

        ++m_read_size;

        ++m_pos;

        if( (m_pos % N == 0) && !m_vector->is_eof())
        {

            m_buf = m_vector->get();

            m_vector->next();

            m_pos = 0;
        }
    }

    /** \brief 指向下一个读取位置，同时删除当前值
     *
     * \param void
     * \param
     * \return void
     *
     */
    void next_remove()
    {

        ++m_read_size;

        ++m_pos;

        if( (m_pos % N == 0) && !m_vector->is_eof())
        {

            m_buf = m_vector->get();

            m_vector->next_remove();

            m_pos = 0;
        }
    }


    /** \brief 指向下一个逆序读取位置(逆序)
     *
     * \param void
     * \param
     * \return void
     *
     */
    void next_reverse()
    {

        ++m_read_size;

        ++m_pos;

        if(m_pos % N == 0 && !m_vector->is_eof())
        {

            m_buf = m_vector->get_reverse();

            m_vector->next_reverse();

            m_pos = 0;
        }
    }


    /** \brief 指向下一个逆序读取位置，同时删除当前值
     *
     * \param void
     * \param
     * \return void
     *
     */
    void next_remove_reverse()
    {

        ++m_read_size;

        ++m_pos;

        if(m_pos % N == 0 && !m_vector->is_eof())
        {

            m_buf = m_vector->get_reverse();

            m_vector->next_remove_reverse();

            m_pos = 0;
        }
    }

    /** \brief 判断是否读完毕
     *
     * \param void
     * \param
     * \return bool 读完返回true, 否则返回false
     *
     */
    bool is_eof()
    {

        return m_size == m_read_size;
    }


    /** \brief 返回bit总数
     *
     * \param void
     * \param
     * \return uint64
     *
     */
    uint64 size()
    {

        return m_size ;
    }

    /** \brief 判断bit总数是否为0
     *
     * \param void
     * \param
     * \return bool
     *
     */
    bool is_empty()
    {

        return m_size == 0 ;
    }

    /** \brief 返回读取的bit总数
     *
     * \param void
     * \param
     * \return uint64
     *
     */
    uint64 read_size()
    {

        return m_read_size;
    }

};

// <T1, T2, T3> = <position, block_id, diff>, just for Level 0.
template < typename T >
class InduceSortSameBkt
{

private:

    /**< 定义三元组，T1是块内位置，T2是块号，T3是标志位（标记相邻的两个元组是否相同） */
    //typedef Triple< T1, T2, T3> triple_type;

    /**< 在内存申请的2个缓存区 */
    std::vector< T *> m_buf_array;

    /**< 如果内存缓存区满，将内存缓存区数据写入外存 */
    std::vector< MyVector< T > * > m_disk_vector;

    /**< 内存缓存区的两个数组的位置索引 */
    std::vector< uint32 >  m_pos_array;

    /**< 分别记录读元素和写元素的个数 */
    std::vector< uint32 > m_size;

    /**< 读和写操作切换标记 */
    uint8 read_buf, write_buf;

    /**< 内存缓存数组元素初始化个数 */
    uint64 m_buf_length;

    /// return obj
    T tp;


public:
    /** \brief 构造函数
     *
     * \param m_buf_array(2,nullptr):2个内存数组指针，初始化为空指针
     * \param m_disk_vector(2,nullptr): 2个外存vector指针，初始化为空，如果内存数组满，则将数据交换到外存
     * \param m_pos_array(3,0): 前两个内存数组的位置标记，第3个用于pop()方法中，用于记录内存中元素个数（余数个）
     * \param m_size(2,0): 分别记录读写元素的个数
     * \param read_buf: 标记读数组
     * \param write_buf: 标记写数组
     * \param m_buf_length: 内存数组长度
     * \return null
     *
     */

    InduceSortSameBkt(uint64 _mem_size):m_buf_array(2,nullptr)
        ,m_disk_vector(2,nullptr)
        ,m_pos_array(3,0)
        ,m_size(2,0)
        ,read_buf(1)
        ,write_buf(0)
        ,m_buf_length(0)
    {

        m_buf_length = _mem_size/sizeof(T)/2;

        m_buf_array[0] = new T[m_buf_length];

        m_buf_array[1] = new T[m_buf_length];

    }

    /** \brief 析构函数，将2个内存数组释放，如果外存vector指针不为空，先delete然后赋值为nullptr
     *
     * \return void
     *
     */

    ~InduceSortSameBkt()
    {

        if(m_buf_array[0])
            delete [] m_buf_array[0];

        if(m_buf_array[1])
            delete [] m_buf_array[1];

        if(m_disk_vector[0])
        {

            if(m_disk_vector[0]->is_eof())
                delete m_disk_vector[0];
            else
            {
                std::cout<<" Error ocurrence in Function ~InduceSortSameBkt() in class InduceSortSameBkt in vector.h.\n";
            }

        }

        if(m_disk_vector[1])
        {

            if(m_disk_vector[1]->is_eof())
                delete m_disk_vector[1];
            else
            {
                std::cout<<" Error ocurrence in Function ~InduceSortSameBkt() in class InduceSortSameBkt in vector.h.\n";
            }
        }

        std::cout << "~InduceSortSameBkt() \n";

    }

    /** \brief 写元素，先写内存数组，满则写出至外存
     *
     * \param _tp，三元组引用Triple<T1,T2,T3>
     * \param
     * \return
     *
     */

    void push(T const & _tp)
    {

        if( m_pos_array[write_buf] != 0 && m_pos_array[write_buf] % m_buf_length == 0)
        {

            if( m_disk_vector[write_buf] == nullptr )
            {

                m_disk_vector[write_buf] = new MyVector< T >();

                m_disk_vector[write_buf]->start_write();

            }

            for(uint32 i=0; i < m_pos_array[write_buf]; i++)
            {

                m_disk_vector[write_buf]->push_back( *(m_buf_array[write_buf] + i) );
            }

            m_pos_array[write_buf] = 0;
        }

        *(m_buf_array[write_buf] + m_pos_array[write_buf]) = _tp;

        ++m_pos_array[write_buf];

        ++m_size[write_buf];

    }

    /** \brief 交换读写对象，将写标记更换为读，将读标记更换为写，在写操作执行完成后，读操作执行前调用
     *
     * \param 无
     * \return void
     *
     */
    void swapBuf()
    {

        if(m_disk_vector[write_buf])
            m_disk_vector[write_buf]->start_read();

        m_pos_array[2] = m_pos_array[write_buf];

        m_pos_array[write_buf] = 0;

        m_pos_array[read_buf] = 0;

        if(read_buf)
        {

            read_buf = 0;

            write_buf = 1;

        }
        else
        {

            read_buf = 1;

            write_buf = 0;

        }
    }

    /** \brief 出元素，如有外存元素，先出外存元素，再出内存数组中的元素（内存元素是从左向右读，从0位置开始向右）
     *
     * \param 无
     * \return Triple<T1,T2,T3>
     *
     */
    T pop()
    {

        //T tp;

        if( m_disk_vector[read_buf] != nullptr )
        {

            if(!m_disk_vector[read_buf]->is_eof())
            {

                tp = m_disk_vector[read_buf]->get();

                m_disk_vector[read_buf]->next_remove();

            }
            else
            {

                delete m_disk_vector[read_buf];

                m_disk_vector[read_buf] = nullptr;

                if( m_pos_array[read_buf] < m_pos_array[2] )
                {

                    tp = *(m_buf_array[read_buf] + m_pos_array[read_buf]);

                    ++m_pos_array[read_buf];
                }
            }

        }
        else if( m_pos_array[read_buf] < m_pos_array[2] )
        {

            tp = *(m_buf_array[read_buf] + m_pos_array[read_buf]);

            ++m_pos_array[read_buf];
        }

        if(m_size[read_buf])
            --m_size[read_buf];

        return tp;
    }


    /*void next(){

        return ;
    }*/

    /** \brief 返回读元素个数
     *
     * \param 无
     * \return uint32，元素个数
     *
     */

    uint32 read_size()const
    {

        return m_size[read_buf];
    }

    /** \brief 返回写元素个数
     *
     * \param 无
     * \return uint32，写元素个数
     *
     */
    uint32 write_size()const
    {

        return m_size[write_buf];
    }

    /** \brief 重置对象，在调用对象之前执行，同一个对象多次使用
     *
     * \param 无
     * \param
     * \return 无
     *
     */
    void reset()
    {

        if(m_disk_vector[0])
        {

            if(m_disk_vector[0]->is_eof())
                delete m_disk_vector[0];
            else
            {
                std::cout<<" Error ocurrence in Function ~InduceSortSameBkt() in class InduceSortSameBkt in vector.h.\n";
            }

        }

        if(m_disk_vector[1])
        {

            if(m_disk_vector[1]->is_eof())
                delete m_disk_vector[1];
            else
            {
                std::cout<<" Error ocurrence in Function ~InduceSortSameBkt() in class InduceSortSameBkt in vector.h.\n";
            }
        }

        m_pos_array[0] = 0;

        m_pos_array[1] = 0;

        m_pos_array[2] = 0;

        m_size[0] = 0;

        m_size[1] = 0;

        read_buf = 1;

        write_buf = 0;
    }

};



class InduceSortSameBit
{

private:

    typedef BitVectorWapper< uint64 > bit_vector_type;
    /**< 不需要内存buffer，BitVectorWapper自带buffer */
    std::vector< bit_vector_type * > m_disk_vector;

    /**< 分别记录读元素和写元素的个数 */
    std::vector< uint32 > m_size;

    /**< 读和写操作切换标记 */
    uint8 read_buf, write_buf;

    /// return obj
    bool is_true;


public:
    /** \brief 构造函数
     *
     * \param m_disk_vector(2,nullptr): 2个外存vector指针，初始化为空，如果内存数组满，则将数据交换到外存
     * \param m_size(2,0): 分别记录读写元素的个数
     * \param read_buf: 标记读数组
     * \param write_buf: 标记写数组
     * \param m_buf_length: 内存数组长度
     * \return null
     *
     */

    InduceSortSameBit():
        m_disk_vector(2,nullptr)
        ,m_size(2,0)
        ,read_buf(1)
        ,write_buf(0)
    {

    }

    /** \brief 析构函数，将2个内存数组释放，如果外存vector指针不为空，先delete然后赋值为nullptr
     *
     * \return void
     *
     */

    ~InduceSortSameBit()
    {

        if(m_disk_vector[0])
        {
            delete m_disk_vector[0];
        }

        if(m_disk_vector[1])
        {
            delete m_disk_vector[1];
        }

        std::cout << "~InduceSortSameBkt() \n";

    }

    /** \brief 写元素，先写内存数组，满则写出至外存
     *
     * \param _tp，三元组引用Triple<T1,T2,T3>
     * \param
     * \return
     *
     */

    void push(bool const & _tp)
    {
        m_disk_vector[write_buf]->push_back(_tp);
        m_size[write_buf]++;
    }

    /** \brief 交换读写对象，将写标记更换为读，将读标记更换为写，在写操作执行完成后，读操作执行前调用
     *
     * \param 无
     * \return void
     *
     */
    void swapBuf()
    {

        if(read_buf)
        {

            read_buf = 0;

            write_buf = 1;

        }
        else
        {

            read_buf = 1;

            write_buf = 0;

        }

        m_disk_vector[read_buf]->start_read();

        if(m_disk_vector[write_buf])delete m_disk_vector[write_buf];
        m_disk_vector[write_buf] = new bit_vector_type();
        m_disk_vector[write_buf]->start_write();
    }

    /** \brief 出元素，如有外存元素，先出外存元素，再出内存数组中的元素（内存元素是从左向右读，从0位置开始向右）
     *
     * \param 无
     * \return Triple<T1,T2,T3>
     *
     */
    bool pop()
    {
        bool is_ture = m_disk_vector[read_buf]->get();
        m_disk_vector[read_buf]->next_remove();
        m_size[read_buf]--;
        return is_ture;
    }


    /*void next(){

        return ;
    }*/

    /** \brief 返回读元素个数
     *
     * \param 无
     * \return uint32，元素个数
     *
     */

    uint32 read_size()const
    {

        return m_size[read_buf];
    }

    /** \brief 返回写元素个数
     *
     * \param 无
     * \return uint32，写元素个数
     *
     */
    uint32 write_size()const
    {

        return m_size[write_buf];
    }

    /** \brief 重置对象，在调用对象之前执行，同一个对象多次使用
     *
     * \param 无
     * \param
     * \return 无
     *
     */
    void reset()
    {

        if(m_disk_vector[0])
        {

            delete m_disk_vector[0];
            m_disk_vector[0] = nullptr;
        }

        if(m_disk_vector[1])
        {

            delete m_disk_vector[1];
            m_disk_vector[1] = nullptr;
        }


        m_size[0] = 0;

        m_size[1] = 0;

        read_buf = 1;

        write_buf = 0;

        m_disk_vector[write_buf] = new bit_vector_type();
        m_disk_vector[write_buf]->start_write();
    }

};



#endif

