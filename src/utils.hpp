#ifndef UTILS_HPP_
#define UTILS_HPP_

#include <iostream>
#include <limits>
#include <string>
#include <cstring>
#include <vector>
#include <cmath>
#include <algorithm>
#include <queue>
//#include <execution>
//#include <numeric>
#include <stdlib.h>
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
//Also, see strict_differences and candis_distance
#define EDIT_DISTANCE 7
#define READ_BLOCK_SIZE 2000
#define READ_LEN_MAX 160
#define READ_LEN_U32 15

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

struct Read{
	int length;
	std::uint8_t bases[READ_LEN_MAX];
	std::uint8_t rc_bases[READ_LEN_MAX];
	char charBases[READ_LEN_MAX];
};

// Can handle reads of length up to 158
// 158 bases x 3 bits/base + 6 spare bits
struct UnmappedRead{
    std::uint32_t fwd_read[READ_LEN_U32];
    std::uint32_t bwd_read[READ_LEN_U32];
    std::uint32_t location;
};

// Can handle reads of length up to 158
// 158 bases x 3 bits/base + 6 spare bits
// MS-Nibble holds number of differences
struct MappedRead{
    std::uint32_t read[READ_LEN_U32];
    std::uint32_t location;
    std::uint32_t base_location;
    int diff_location;
    std::uint32_t id;
    char tmp_cigar[16];
};

struct PairingInfoFwd {
    std::uint32_t id;
    std::uint32_t my_rel_posn;
};

struct PairingInfoBwd {
    std::uint32_t id;
    std::uint32_t my_location;
    std::uint32_t pe_location;
    std::uint32_t pe_rel_posn;
};

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
    std::uint32_t rdLength;
    //std::uint32_t kmerLength;

    InputArgs(){
        threadCount = 2;
        rd1Length   = 0;
        rd2Length   = 0;
        rdLength    = 0;
        //kmerLength = 0;
    }
};

struct CompressionDataStructures{
    std::uint32_t ref_length; //Concatenated reference length excluding Ns
    std::uint8_t * ref_bases;
    std::vector<std::uint32_t> lookup_table;
    std::vector<std::uint32_t> occurrence_table;
    std::vector<UnmappedRead> unmapped_reads;
    std::vector<MappedRead> fwd_mapped_reads;
    std::vector<MappedRead> bwd_mapped_reads;
    std::vector<PairingInfoFwd> fwd_pairing_info;
    std::vector<PairingInfoBwd> bwd_pairing_info;
    
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






















































