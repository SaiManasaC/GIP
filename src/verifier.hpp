#ifndef VERIFIER_HPP_
#define VERIFIER_HPP_

#include "utils.hpp"
#include <emmintrin.h>
#include <smmintrin.h>

//we consider alphabet set as {A, T, C, G, N}
#define ALPHABET_SIZE 5
//v_c and v_p are 8 and 16 respectively
#define V_CPU 8
#define V_PHI 16
#define CANDIDATES_COUNT (V_CPU*80)
#define PLACES_COUNT (V_CPU*80*3)

#define SAVE_LOCS_COUNT 16
#define SetBit(A,k)   (A[(k/32)] |=  (1 << (k%32)))
#define TestBit(A,k)  (A[(k/32)] &   (1 << (k%32)))
#define ClearBit(A,k) (A[(k/32)] &= ~(1 << (k%32)))

struct TwoTuple{
	uint32_t a;
	int b;
};

struct AlignmentStatistics{
    std::vector<std::size_t> as_unmapped_reads_n_count; //threadCount
    std::vector<std::size_t> as_unmapped_diffs_n_count; //threadCount
    std::vector<std::size_t> as_paired_read_count; //threadCount
    std::vector<std::size_t> as_paired_mapped_count; //threadCount
    std::vector<std::size_t> as_single_mapped_count; //threadCount
    std::vector<std::size_t> as_both_unmapped_count; //threadCount
    std::vector<std::size_t> as_can_unmapped_count; //threadCount
    std::vector<std::size_t> as_pla_unmapped_count; //threadCount
    std::vector<std::size_t> as_nil_unmapped_count; //threadCount
    
    //Includes two-mapped and one-mapped reads
    std::vector<std::size_t> as_differences_count; //threadCount * 16
    
    std::vector<std::size_t> as_less_than_length_count; //threadCount
    std::vector<std::size_t> as_equal_to_length_count; //threadCount
    std::vector<std::size_t> as_more_than_length_count; //threadCount
    std::vector<std::size_t> as_neg_diff_count_bef; //threadCount
    std::vector<std::size_t> as_neg_diff_count_aft; //threadCount
    std::vector<std::size_t> as_cigar_unequal_count; //threadCount
    std::vector<std::size_t> as_location_unequal_count; //threadCount
    std::vector<std::size_t> as_soft_clipped_count; //threadCount
    
    std::vector<std::size_t> as_paired_read_order; //threadCount * 8
};

struct VerifyArgsForThread{
    std::uint32_t vaft_thread_id;
    InputArgs * vaft_in_args;
    CompressionDataStructures * vaft_com_ds;
    AlignmentStatistics * asp;
    Read * read_1;
    Read * read_2;
    MappedRead * fwd_read;
    MappedRead * bwd_read;
    UnmappedRead * unm_read;
};

bool NJ_CPUVBM(VerifyArgsForThread *);
inline void CPUVMBMKernel(VerifyArgsForThread * vaftp, const uint32_t regIndex, const uint8_t* text, const int readLen, const uint32_t *candidates, const uint32_t candiNum, int16_t *errors, int16_t *locations);
int CPUAlignKernel(const uint8_t* pattern, const uint8_t* text, int match_site, char* cigar, const int readLen, int best_error); 
void CPUMDTag(const uint8_t* pattern, const uint8_t* text, uint32_t startLocation, char* cigar, int* my_diff_locs, int* my_diff_vals);
void uncompact_mapped_read (uint8_t *, uint32_t *, uint32_t);

#endif /* VERIFIER_HPP_ */

