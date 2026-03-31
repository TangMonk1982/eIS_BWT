/******************************************************************************
 * sais.h
 *
 ******************************************************************************
 * Ling Bo Han <hanlb (at) mail2.sysu.edu.cn>
 * Yi Wu <wu.yi.christian (at) gmail.com>
 * Ge Nong <issng (at) mail.sysu.edu.cn>
 ******************************************************************************
*/

#ifndef MY_SAIS_H
#define MY_SAIS_H


#include "common.h"
#include "utility.h"
#include "vector.h"


//#define TEST_DEBUG2
template<typename alphabet_type, typename offset_type>
class SAComputation;

/// \brief portal to SA computation on RAM
///
template<typename alphabet_type, typename offset_type>
class SAIS{

private:

	typedef MyVector<offset_type> offset_vector_type;

public:

	SAIS( std::string & _s, offset_vector_type *& _sa_reverse);

};

/// \brief compute SA for the input string residing on EM
///
/// \note alphabet_type = offset_type > uint32
template<typename alphabet_type, typename offset_type>
SAIS<alphabet_type, offset_type>::SAIS(std::string & _s, offset_vector_type *& _sa_reverse) {

    std::cout << "Compute SA in RAM. \n";
	// load _s into ram
	// note that characters in _s are named in a condensed way
	//uint32 s_size = _s->size();
	offset_type s_size = UtilityFunctions::getFileLength(_s) / sizeof(alphabet_type);

	alphabet_type * s = new alphabet_type [s_size];

    UtilityFunctions::fileRead<alphabet_type>(s,_s,0,s_size);

    /// remove the _s string
    //std::remove(_s.c_str());

    //Logger::delPDU(s_size * sizeof(alphabet_type));

    std::cout << "read buf closed.\n";

	alphabet_type max_alpha = 0;

	for (offset_type i = 0; i < s_size; ++i) {

        max_alpha = std::max(max_alpha, s[i]);
	}

	std::cout << "max_alpha = " << max_alpha << std::endl;

	// compute sa
	offset_type *sa = new offset_type[s_size];

	std::cout << "start to compute the sa.\n";

	SAComputation<alphabet_type, offset_type>(s, s_size, max_alpha, sa);

	std::cout << "compute sa closed.\n";

	// output sa reversely
	_sa_reverse = new offset_vector_type();

	_sa_reverse->start_write();

	for (uint64 i = s_size - 1; i >= 0; --i) {

		_sa_reverse->push_back(sa[i]);

		if(i == 0) break;
	}

	delete [] s;

	delete [] sa;
}



/// \brief compute SA on RAM
///
template<typename alphabet_type, typename offset_type>
class SAComputation{

	const offset_type uint_MAX = std::numeric_limits<offset_type>::max();

public:

	SAComputation(alphabet_type * _s, const offset_type _sLen, const alphabet_type _alpha, offset_type * _sa);

	void getBuckets(alphabet_type * _s, const offset_type _sLen, offset_type * _bkt, const offset_type _bktNum, const bool _end);

	void induceL(alphabet_type * _s, BitWrapper & _t, offset_type * _sa, const offset_type _sLen, offset_type * _bkt, const alphabet_type _alpha);

	void induceS(alphabet_type * _s, BitWrapper & _t, offset_type * _sa, const offset_type _sLen, offset_type * _bkt, const alphabet_type _alpha);
};

/// \brief ctor
///
template<typename alphabet_type, typename offset_type>
SAComputation<alphabet_type, offset_type>::SAComputation(alphabet_type * _s, const offset_type _sLen, const alphabet_type _alpha, offset_type * _sa) {

	offset_type i, j, offset_zero(0), offset_one(1);

	char * t_buf = new char[_sLen / 8 + 1]; BitWrapper t(t_buf);

	//compute t
	t.set(_sLen - 1, S_TYPE); t.set(_sLen - 2, L_TYPE);
	for (i = _sLen - 3; ; --i) {
		t.set(i, (_s[i] < _s[i + 1] || (_s[i] == _s[i + 1] && t.get(i + 1) == S_TYPE)) ? S_TYPE : L_TYPE);
		if (i == offset_zero) break;
	}

	// sort all the S-substrings
	offset_type *bkt = new offset_type[_alpha + 1];
	getBuckets(_s, _sLen, bkt, _alpha + 1, true);


	for (i = 0; i < _sLen; ++i) _sa[i] = uint_MAX;//init sa
	for (i = _sLen - 3; i >= offset_one; --i) {	// find lms-chars in _s[1, _sLen - 3]
		if (t.get(i) && !t.get(i - 1)) {
			_sa[bkt[_s[i]]] = i;
			--bkt[_s[i]];
		}
	}
	_sa[0] = _sLen - 1;	//_s[_sLen - 1] is the sentinel smallest than any other characters

	induceL(_s, t, _sa, _sLen, bkt, _alpha); //induce l-type lms-prefix
	induceS(_s, t, _sa, _sLen, bkt, _alpha); //induce s-type lms-prefix

	delete[] bkt;

	// compact all the sorted substrings into the first n1 items of s
	offset_type s1_len = 0;
	for (i = 0; i < _sLen; ++i) {
		if (_sa[i] > offset_zero && t.get(_sa[i]) && !t.get(_sa[i] - 1)) {
			_sa[s1_len] = _sa[i];
			++s1_len;
		}
	}

	for (i = s1_len; i < _sLen; ++i) _sa[i] = uint_MAX; //init

	// find the lexicographic names of all substrings
	offset_type name = 0;
	offset_type prev = uint_MAX;
	for (i = 0; i< s1_len; ++i) {
		offset_type pos = _sa[i]; bool diff = false;
		for (offset_type d = 0; d < _sLen; ++d) {
			if (prev == uint_MAX || pos + d == _sLen - 1 || prev + d == _sLen - 1 || _s[pos + d] != _s[prev + d] || t.get(pos + d) != t.get(prev + d)) {
				diff = true;
				break;
			}
			else {
				if (d > offset_zero && ((t.get(pos + d) && !t.get(pos + d - 1)) || (t.get(prev + d) && !t.get(prev + d - 1)))) {
					break;
				}
			}
		}
		if (diff) {
			++name;
			prev = pos;
		}
		pos = pos / 2;
		_sa[s1_len + pos] = name - 1;
	}
	for (i = _sLen - 1, j = _sLen - 1; i >= s1_len; --i)
		if (_sa[i] != uint_MAX) _sa[j] = _sa[i],--j;

	// s1 is done now
	offset_type *sa1 = _sa;
	offset_type *s1 = _sa + _sLen - s1_len;
	// stage 2: solve the reduced problem

	// recurse if names are not yet unique
	if (name < s1_len) {
		SAComputation<offset_type, offset_type>(s1, s1_len, name - 1, sa1);
	}
	else { // generate the suffix array of s1 directly
		for (i = 0; i < s1_len; ++i) {
			sa1[s1[i]] = i;
		}
	}

	// stage 3: induce the result for the original problem
	// put all left-most S characters into their buckets
	bkt = new offset_type[_alpha + 1];
	getBuckets(_s, _sLen, bkt, _alpha + 1, true); // find ends of buckets

	j = 0;
	for (i = 1; i < _sLen; ++i) {
		if (t.get(i) && !t.get(i - 1)) s1[j] = i, ++j;// get p1
	}
	for (i = 0; i < s1_len; ++i) sa1[i] = s1[sa1[i]]; // get index in s1
	for (i = s1_len; i < _sLen; ++i) _sa[i] = uint_MAX; // init SA[n1..n-1]
	for (i = s1_len - 1; ; --i) {
		j = _sa[i]; _sa[i] = uint_MAX;
		_sa[bkt[_s[j]]] = j, --bkt[_s[j]];
		if (i == offset_zero) break;
	}
	induceL(_s, t, _sa, _sLen, bkt, _alpha);
	induceS(_s, t, _sa, _sLen, bkt, _alpha);

	delete[] bkt; bkt = nullptr;
	delete[] t_buf; t_buf = nullptr;
}


//@usage: compute bucket size.
template<typename alphabet_type, typename offset_type>
void SAComputation<alphabet_type, offset_type>::getBuckets(alphabet_type * _s, const offset_type _sLen, offset_type * _bkt, const offset_type _bktNum, const bool _end) {
	offset_type i;
	offset_type sum = 0;
	for (i = 0; i < _bktNum; ++i) _bkt[i] = 0; // clear all buckets
	for (i = 0; i < _sLen; ++i) ++_bkt[_s[i]]; // compute the size of each bucket
	for (i = 0; i < _bktNum; ++i) { sum += _bkt[i]; _bkt[i] = _end ? sum - 1 : sum - _bkt[i]; }
}

//@usage: induce L.
template<typename alphabet_type, typename offset_type>
void SAComputation<alphabet_type, offset_type>::induceL(alphabet_type * _s, BitWrapper & _t, offset_type * _sa, const offset_type _sLen, offset_type * _bkt, const alphabet_type _alpha) {
	offset_type i, j, offset_zero(0);
	getBuckets(_s, _sLen, _bkt, _alpha + 1, false); // find heads of buckets
	for (i = 0; i < _sLen; ++i) {
		if (_sa[i] != uint_MAX && _sa[i] != offset_zero) {
			j = _sa[i] - 1;
			if (!_t.get(j)) { // _sa[i] is unsigned, if _sa[i] == 0, then j does not exist.
				_sa[_bkt[_s[j]]] = j, ++_bkt[_s[j]];
			}
		}
	}
}

//@usage: induce S.
template<typename alphabet_type, typename offset_type>
void SAComputation<alphabet_type, offset_type>::induceS(alphabet_type * _s, BitWrapper & _t, offset_type * _sa, const offset_type _sLen, offset_type * _bkt, const alphabet_type _alpha) {
	offset_type i, j, offset_zero(0);
	getBuckets(_s, _sLen, _bkt, _alpha + 1, true); // find ends of buckets
	for (i = _sLen - 1; ; --i) {
		if (_sa[i] != uint_MAX && _sa[i] != offset_zero) {
			j = _sa[i] - 1;
			if (_t.get(j)) { // _sa[i] is unsigned, if _sa[i] == 0, then j == EMPTY, which does not exist.
				_sa[_bkt[_s[j]]] = j, --_bkt[_s[j]];
			}
		}
		if (i == offset_zero) break;
	}
}

#endif // MY_SAIS_H
