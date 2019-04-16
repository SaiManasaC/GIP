#ifndef UTILS_HPP_
#define UTILS_HPP_

#include <iostream>
#include <limits>
#include <string>
#include <vector>
#include <cmath>
#include <algorithm>
//#include <execution>
//#include <numeric>
#include <pthread.h>
#include <assert.h>
#include <unistd.h>
#include <sys/time.h>
#include <sys/resource.h>

#define KMER_LENGTH 14
#define REF_LEN_MAX 4294967294
//From "parallel_stable_sort.h"
#define SORT_CUT_OFF 500
#define MERGE_CUT_OFF 2000

#define __UTILS_CHAR_UINT8 	\
    static inline std::uint8_t charToUint8(const char c) {\
        switch (c) {\
            case 'A':\
            case 'a':\
                     return 0;\
            case 'C':\
            case 'c':\
                     return 1;\
            case 'G':\
            case 'g':\
                     return 2;\
            case 'T':\
            case 't':\
                     return 3;\
            default:\
                     return 4;\
        }\
    }

#define __UTILS_HASH_VALUE 	 \
static inline int hashValue(const std::uint8_t *seq, int windowSize) {\
	int i = 0;\
	int val = 0;\
	while (i < windowSize) {\
		if (seq[i] == (std::uint8_t) 4) {\
			return -1;\
		}\
		val = (val << 2) | ((int) seq[i]);\
		i++;\
	}\
	return val;\
}

struct InputArgs{
    std::string operation;
    std::string refFileName;
    std::string rd1FileName;
    std::string rd2FileName;
    std::string comFileName;

    //Need below only for compression
    std::uint32_t threadCount;
    std::uint32_t rd1Length;
    std::uint32_t rd2Length;
    //std::uint32_t kmerLength;

    InputArgs(){
        threadCount = 2;
        rd1Length   = 0;
        rd2Length   = 0;
        //kmerLength = 0;
    }
};

struct CompressionDataStructures{
    std::uint32_t ref_length; //Concatenated reference length excluding Ns
    std::uint8_t * ref_bases;
    std::vector<std::uint32_t> lookup_table;
    std::vector<std::uint32_t> occurrence_table;
    
    CompressionDataStructures(){
        ref_length = 0;
        ref_bases = NULL;
    }
};

struct DecompressionDataStructures{
    std::uint32_t ref_length; //Concatenated reference length excluding Ns
    std::uint8_t * ref_bases;
    
    DecompressionDataStructures(){
        ref_length = 0;
        ref_bases = NULL;
    }
};

double realtime();
double cputime();

#endif /* UTILS_HPP_ */






















































