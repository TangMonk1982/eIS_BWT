/******************************************************************************
 * common.h
 *
 ******************************************************************************
 * Ling Bo Han <hanlb (at) mail2.sysu.edu.cn>
 * Yi Wu <wu.yi.christian (at) gmail.com>
 * Ge Nong <issng (at) mail.sysu.edu.cn>
 ******************************************************************************
*/
#ifndef __COMMON_H
#define __COMMON_H

#include<limits>
#include<vector>
#include<utility>
#include<time.h>
#include<stdlib.h>
#include<stdio.h>
#include<iostream>
//#include<cassert>
#include<algorithm>
#include<queue>
#include<chrono>
#include <random>
#include <assert.h>
#include <mutex>
#include <memory>

#include <cilk/cilk.h>
#include <cilk/cilk_api.h>
#include <cilk/holder.h>

#include"type.h"
#include"uint40.h"
//#include"uint24.h"
#include"logger.h"
#include"timer.h"
#include "tuple.h"



// LMS- L- and S-type
constexpr uint8 LML_TYPE = 3, LMS_TYPE = 2, S_TYPE = 1, L_TYPE = 0;

const uint8 uint8_MAX = std::numeric_limits<uint8>::max();

const uint16 uint16_MAX = std::numeric_limits<uint16>::max();

const uint32 uint24_MAX = 256 * 256 * 256;

const uint32 uint32_MAX = std::numeric_limits<uint32>::max();

const uint40 uint40_MAX = std::numeric_limits<uint40>::max();

const uint64 uint64_MAX = std::numeric_limits<uint64>::max();

// memory allocation
constexpr uint64 K_512 = 512 * 1024;

constexpr uint64 K_1024 = 1024 * 1024;

constexpr uint64 MAX_MEM = 24 * 1024 * K_1024;
//constexpr uint64 MAX_MEM = 3 * 1024 * K_1024;
//onstexpr uint64 MAX_MEM = 50 * K_1024;

static uint64 VEC_BUF_RAM = 2 * K_1024; // 2M

static uint64 PHI_VEC_EM = 200 * K_1024; // 20M

static uint32 global_file_idx = 0;

/**< alphabet set size on each level. */
static std::vector<uint40> alphabetSet_vec;

uint8 transTable[256] = {0};

/**<  MASK_U16[0]: is_diff = true and is_LStar = true; MASK_U16[1]: is_diff = true and is_LStar = false; MASK_U16[2]: is_diff = false and is_LStar = true; */
static const uint16 MASK_U16[] = { 0xc000, 0x8000, 0x4000, 0x3fff, 0x7fff };//MASK_U16[] = { 0xc000, 0x8000, 0x4000, 0x3fff, 0x7fff }

static bool global_unique = false;

// compute BWT
/// \note computing BWT

bool is_getBWT = true; // false or true

bool is_getSA = false; // false or true

// record the preceding char of the LStar suffix
std::vector< void * > LStarPreChar_seqs;
std::vector< void * > SStarPreChar_seqs;
//std::vector< void * > S_eBWT_seqs;


#endif // __COMMON_H






















