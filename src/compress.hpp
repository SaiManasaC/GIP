#ifndef COMPRESS_HPP_
#define COMPRESS_HPP_

#include "utils.hpp"
#include "verifier.hpp"

struct CompressionStatistics{
    std::size_t cs_pe_mapped_count;
    
    std::vector<std::size_t> cs_cigar_resolved_count; //threadCount
    std::vector<std::size_t> cs_diff_more_count; //threadCount
    std::vector<std::size_t> cs_diff_more_value; //threadCount
    std::vector<std::size_t> cs_diff_same_count; //threadCount
    std::vector<std::size_t> cs_diff_less_count; //threadCount
    std::vector<std::size_t> cs_diff_less_value; //threadCount
    std::vector<std::size_t> cs_edit_count_befr; //threadCount
    std::vector<std::size_t> cs_subs_count_befr; //threadCount
    std::vector<std::size_t> cs_inss_count_befr; //threadCount
    std::vector<std::size_t> cs_dels_count_befr; //threadCount
    std::vector<std::size_t> cs_fwd_locns_count; //threadCount
    std::vector<std::size_t> cs_bwd_locns_count; //threadCount
    std::vector<std::size_t> cs_fwd_diff_counts_count; //threadCount
    std::vector<std::size_t> cs_bwd_diff_counts_count; //threadCount
    std::vector<std::size_t> cs_fwd_diff_posns_count; //threadCount
    std::vector<std::size_t> cs_bwd_diff_posns_count; //threadCount
    std::vector<std::size_t> cs_fwd_diff_values_count; //threadCount
    std::vector<std::size_t> cs_bwd_diff_values_count; //threadCount
    std::vector<std::size_t> cs_pe_rel_locns_count; //threadCount
    std::vector<std::size_t> cs_pe_rel_posns_count; //threadCount
    
    //Only two-mapped reads
    std::vector<std::size_t> cs_differences_count_befr; //threadCount * 16
    std::vector<std::size_t> cs_differences_count_aftr; //threadCount * 16

    std::vector<std::size_t> cs_edit_count_aftr; //threadCount
    std::vector<std::size_t> cs_subs_count_aftr; //threadCount
    std::vector<std::size_t> cs_inss_count_aftr; //threadCount
    std::vector<std::size_t> cs_dels_count_aftr; //threadCount
};

struct CompressArgsForThread{
    std::uint32_t caft_thread_id;
    InputArgs * caft_in_args;
    CompressionDataStructures * caft_com_ds;
    pthread_barrier_t * caft_barrier;
    CompressionStatistics * csp;
    //FILE ** caft_fp;
    std::uint8_t * caft_sur;
    bool caft_is_fwd_read;
    std::uint32_t caft_locn_diff;
    int64_t caft_pe_rel_locn_sum;
    int64_t caft_pe_rel_locn_cnt;
    int64_t caft_pe_rel_locn_mean;
    std::vector<std::uint8_t> caft_fwd_locns;
    std::vector<std::uint8_t> caft_bwd_locns;
    std::vector<std::uint8_t> caft_fwd_diff_counts;
    std::vector<std::uint8_t> caft_bwd_diff_counts;
    std::vector<std::uint8_t> caft_fwd_diff_posns;
    std::vector<std::uint8_t> caft_bwd_diff_posns;
    std::vector<std::uint8_t> caft_fwd_diff_values;
    std::vector<std::uint8_t> caft_bwd_diff_values;
    std::vector<std::uint8_t> caft_pe_rel_locns;
    std::vector<std::uint8_t> caft_pe_rel_posns;
};

int compress_reads(InputArgs&, CompressionDataStructures&);

#endif /* COMPRESS_HPP_ */

