///////////////////////////////////////////////////////////////////////////////////////////////////////////
/// Copyright (c) 2017, Sun Yat-sen University.
/// All rights reserved.
/// \file sorter.h
/// \brief A self-defined EM sorter for implementing space-efficient I/O operations.
///
/// The sorter is based on stxxl::vector. It first splits elements into blocks and sorts each block in RAM.
/// Afterward, it merges the block-wise results by a loser tree.
///
/// \author Yi Wu
/// \date 2017.7
/////////////////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef _SORTER_H
#define _SORTER_H

#include "vector.h"
#include "losertree.h"
///STL parallel sort Library
#include <parallel/algorithm>

/// \brief definition of I/O operations on self-defined disk-based sorter
///
/// All elements are naturally splitted into multiple blocks, the elements in each block are sorted in RAM and forwarded to EM for later use.
/// Assuming that O(n) <= m^2, a loser tree is applied to merging the blockwise results into a whole.
/// \note any two elements are assumed to be different from each other
template<typename element_type, typename comparator_type>
class MySorter{

// alias
private:

	typedef MyVector<element_type> element_vector_type; // vector type

public:

	/// \brief merge elements by a loser tree
	///
	struct ElementCompare{

	public:

		std::vector<element_vector_type*>& m_block; ///< handlers to blocks in EM

		uint32 m_block_num; ///< number of blocks

		std::vector<element_type> m_head; ///< smallest elements among thoes unvisited in each block

		std::vector<bool> m_exists; ///< indicator for existence

	public:

		/// \brief ctor
		///
		ElementCompare(std::vector<element_vector_type*>& _block, const uint32 _block_num) : m_block(_block), m_block_num(_block_num) {

			m_exists.resize(m_block_num, true);

			m_head.resize(m_block_num);

			for (uint32 i = 0; i < m_block_num; ++i) {

				m_block[i]->start_read();

				m_head[i] = m_block[i]->get();

				m_block[i]->next_remove();
			}
		}

		/// \brief fetch current smallest of the specified block
		///
		void fetch(const uint32 _block_idx) {

			if (!m_block[_block_idx]->is_eof()) {

				m_head[_block_idx] = m_block[_block_idx]->get();

				m_exists[_block_idx] = true;

				m_block[_block_idx]->next_remove();
			}
			else {

				m_exists[_block_idx] = false;
			}
		}

		/// \brief check if there remains elements to be read from the specified block
		///
		bool exists(const uint32 _block_idx) const {

			return m_exists[_block_idx];
		}

		/// \brief compare two elements
		///
		/// \note any two elements must be different from each other (assumption)
		int operator() (const uint32 _block_idx_a, const uint32 _block_idx_b) const {

			static comparator_type ct = comparator_type();

			return ct(m_head[_block_idx_a], m_head[_block_idx_b]) ? -1 : 1;
		}
	};

	const uint64 m_block_capacity; ///< capacity for each block

	std::vector<element_type> *m_ram_block; ///< a block container in RAM

	std::vector<element_vector_type*> m_em_blocks; ///< blocks stored in EM

	uint64 m_size; ///< number of elements in all the blocks

	ElementCompare *m_cmp;  ///< instance of compare

	LoserTree3Way<ElementCompare> *m_ltree; ///< loser tree

	uint32 m_block_num; ///< number of blocks

	uint32 m_pre_loser; ///< block idx of previously popped element

	uint32 m_cur_loser; ///< block idx of currently popped element

	clock_t m_timer;

public:

	/// \brief ctor
	///
	MySorter(const uint64 _avail_mem) : m_block_capacity(_avail_mem / sizeof(element_type)) {

		m_ram_block = new std::vector<element_type>();

		m_ram_block->clear();

		m_em_blocks.clear();

		m_size = 0;
	}

	/// \brief push an element into the sorter
	///
	/// The elements are orgainzed into multiple blocks, where each block is sorted in RAM before it is redirected to EM.
	void push(const element_type & _value) {

		if (m_ram_block->size() == m_block_capacity) {

#ifdef PARALELL_SORT
			__gnu_parallel::sort(m_ram_block->begin(), m_ram_block->end(), comparator_type());
#else
			std::sort(m_ram_block->begin(), m_ram_block->end(), comparator_type());
#endif // PARALELL_SORT



			m_em_blocks.push_back(new element_vector_type()), m_em_blocks.back()->start_write();

			for (uint64 i = 0; i < m_ram_block->size(); ++i) {

				m_em_blocks[m_em_blocks.size() - 1]->push_back((*m_ram_block)[i]);
			}

			m_size += m_ram_block->size();

			m_ram_block->clear();
		}

		m_ram_block->push_back(_value);
	}

	/// \brief employ the loser-tree to merge the block-wise results
	///
	void sort() {

		// process the remaining elements in the ram block
		if (m_ram_block->size() != 0) {

#ifdef PARALELL_SORT
			__gnu_parallel::sort(m_ram_block->begin(), m_ram_block->end(), comparator_type());
#else
			std::sort(m_ram_block->begin(), m_ram_block->end(), comparator_type());
#endif // PARALELL_SORT




			m_em_blocks.push_back(new element_vector_type()), m_em_blocks.back()->start_write();

			for (uint64 i = 0; i < m_ram_block->size(); ++i) {

				m_em_blocks[m_em_blocks.size() - 1]->push_back((*m_ram_block)[i]);
			}

			m_size += m_ram_block->size();

			m_ram_block->clear();
		}

		// free m_ram_block
		delete m_ram_block; m_ram_block = nullptr;

		// initialize the loser tree
		m_block_num = m_em_blocks.size();

		m_cmp = new ElementCompare(m_em_blocks, m_block_num);

		m_pre_loser = std::numeric_limits<uint32>::max();

		m_ltree = new LoserTree3Way<ElementCompare>(*m_cmp);

		m_ltree->play_initial(m_block_num);
	}

	/// \brief get the top element
	///
	const element_type& operator*() {

		m_cur_loser = m_ltree->top();

		return m_cmp->m_head[m_cur_loser];
	}

	/// \brief forward
	///
	void operator++() {

		m_pre_loser = m_cur_loser;

		m_cmp->fetch(m_pre_loser);

		m_ltree->replay();
	}

	/// \brief check if empty
	///
	bool empty() const {

		return m_ltree->done();
	}

	/// \brief return the number of elements in the sorter
	///
	uint64 size() {

		return m_size;
	}

	~MySorter() {

		delete m_ltree; m_ltree = nullptr;

		delete m_cmp; m_cmp = nullptr;
	}


	void report() {

		std::cerr << "block num: " << m_em_blocks.size() << std::endl;

		std::cerr << "size: " << m_size << std::endl;
	}
};

#endif // _SORTER_H

