#ifndef DECOMPRESS_HPP_
#define DECOMPRESS_HPP_

#include "utils.hpp"
#include "kseq.h"

struct DecompressionStatistics{
    std::size_t ds_pe_mapped_count;
    
    std::vector<std::size_t> ds_edit_count_befr; //threadCount
    std::vector<std::size_t> ds_subs_count_befr; //threadCount
    std::vector<std::size_t> ds_inss_count_befr; //threadCount
    std::vector<std::size_t> ds_dels_count_befr; //threadCount
    std::vector<std::size_t> ds_fwd_locns_count; //threadCount
    std::vector<std::size_t> ds_bwd_locns_count; //threadCount
    std::vector<std::size_t> ds_fwd_diff_counts_count; //threadCount
    std::vector<std::size_t> ds_bwd_diff_counts_count; //threadCount
    std::vector<std::size_t> ds_fwd_diff_posns_count; //threadCount
    std::vector<std::size_t> ds_bwd_diff_posns_count; //threadCount
    std::vector<std::size_t> ds_fwd_diff_values_count; //threadCount
    std::vector<std::size_t> ds_bwd_diff_values_count; //threadCount
    std::vector<std::size_t> ds_pe_rel_locns_count; //threadCount
    std::vector<std::size_t> ds_pe_rel_posns_count; //threadCount
    
    //Only two-mapped reads
    std::vector<std::size_t> ds_differences_count_befr; //threadCount * 16
};

struct DecompressArgsForThread{
    std::uint32_t daft_thread_id;
    InputArgs * daft_in_args;
    DecompressionDataStructures * daft_decom_ds;
    pthread_barrier_t * daft_barrier;
    DecompressionStatistics * dsp;
    std::vector<std::uint8_t> daft_fwd_locns;
    std::vector<std::uint8_t> daft_bwd_locns;
    std::vector<std::uint8_t> daft_fwd_diff_counts;
    std::vector<std::uint8_t> daft_bwd_diff_counts;
    std::vector<std::uint8_t> daft_fwd_diff_posns;
    std::vector<std::uint8_t> daft_bwd_diff_posns;
    std::vector<std::uint8_t> daft_fwd_diff_values;
    std::vector<std::uint8_t> daft_bwd_diff_values;
    std::vector<std::uint8_t> daft_pe_rel_locns;
    std::vector<std::uint8_t> daft_pe_rel_posns;
    std::vector<char> daft_fwd_reads;
};

int decompress_reads(InputArgs&, DecompressionDataStructures&);

#endif /* DECOMPRESS_HPP_ */

