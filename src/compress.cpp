#include "compress.hpp"
//software.intel.com/en-us/articles/a-parallel-stable-sort-using-c11-for-tbb-cilk-plus-and-openmp
#include "parallel_stable_sort.h"

__UTILS_CHAR_UINT8
__UTILS_HASH_VALUE

void initialize_compression_stats(InputArgs& in_args, CompressionStatistics * ics_csp) {
    (ics_csp->cs_cigar_resolved_count).resize(in_args.threadCount, 0);
    (ics_csp->cs_diff_more_count).resize(in_args.threadCount, 0);
    (ics_csp->cs_diff_more_value).resize(in_args.threadCount, 0);
    (ics_csp->cs_diff_same_count).resize(in_args.threadCount, 0);
    (ics_csp->cs_diff_less_count).resize(in_args.threadCount, 0);
    (ics_csp->cs_diff_less_value).resize(in_args.threadCount, 0);
    (ics_csp->cs_edit_count_befr).resize(in_args.threadCount, 0);
    (ics_csp->cs_subs_count_befr).resize(in_args.threadCount, 0);
    (ics_csp->cs_inss_count_befr).resize(in_args.threadCount, 0);
    (ics_csp->cs_dels_count_befr).resize(in_args.threadCount, 0);
    (ics_csp->cs_fwd_locns_count).resize(in_args.threadCount, 0);
    (ics_csp->cs_bwd_locns_count).resize(in_args.threadCount, 0);
    (ics_csp->cs_fwd_diff_counts_count).resize(in_args.threadCount, 0);
    (ics_csp->cs_bwd_diff_counts_count).resize(in_args.threadCount, 0);
    (ics_csp->cs_fwd_diff_posns_count).resize(in_args.threadCount, 0);
    (ics_csp->cs_bwd_diff_posns_count).resize(in_args.threadCount, 0);
    (ics_csp->cs_fwd_diff_values_count).resize(in_args.threadCount, 0);
    (ics_csp->cs_bwd_diff_values_count).resize(in_args.threadCount, 0);
    (ics_csp->cs_pe_rel_locns_count).resize(in_args.threadCount, 0);
    (ics_csp->cs_pe_rel_posns_count).resize(in_args.threadCount, 0);
    
    (ics_csp->cs_differences_count_befr).resize((in_args.threadCount * 16), 0);
    (ics_csp->cs_differences_count_aftr).resize((in_args.threadCount * 16), 0);

    (ics_csp->cs_edit_count_aftr).resize(in_args.threadCount, 0);
    (ics_csp->cs_subs_count_aftr).resize(in_args.threadCount, 0);
    (ics_csp->cs_inss_count_aftr).resize(in_args.threadCount, 0);
    (ics_csp->cs_dels_count_aftr).resize(in_args.threadCount, 0);
}

void display_compression_stats(InputArgs& in_args, CompressionStatistics * dcs_csp) {
    std::size_t total_cigar_resolved_count = 0;
    std::size_t total_diff_more_count = 0;
    std::size_t total_diff_more_value = 0;
    std::size_t total_diff_same_count = 0;
    std::size_t total_diff_less_count = 0;
    std::size_t total_diff_less_value = 0;
    std::size_t total_edit_count_befr = 0;
    std::size_t total_subs_count_befr = 0;
    std::size_t total_inss_count_befr = 0;
    std::size_t total_dels_count_befr = 0;
    std::size_t total_fwd_locns_count = 0;
    std::size_t total_bwd_locns_count = 0;
    std::size_t total_fwd_diff_counts_count = 0;
    std::size_t total_bwd_diff_counts_count = 0;
    std::size_t total_fwd_diff_posns_count = 0;
    std::size_t total_bwd_diff_posns_count = 0;
    std::size_t total_fwd_diff_values_count = 0;
    std::size_t total_bwd_diff_values_count = 0;
    std::size_t total_pe_rel_locns_count = 0;
    std::size_t total_pe_rel_posns_count = 0;
    std::size_t total_edit_count_aftr = 0;
    std::size_t total_subs_count_aftr = 0;
    std::size_t total_inss_count_aftr = 0;
    std::size_t total_dels_count_aftr = 0;
    for (std::uint32_t i = 0; i < in_args.threadCount; i++) {
        total_cigar_resolved_count += (dcs_csp->cs_cigar_resolved_count)[i];
        total_diff_more_count += (dcs_csp->cs_diff_more_count)[i];
        total_diff_more_value += (dcs_csp->cs_diff_more_value)[i];
        total_diff_same_count += (dcs_csp->cs_diff_same_count)[i];
        total_diff_less_count += (dcs_csp->cs_diff_less_count)[i];
        total_diff_less_value += (dcs_csp->cs_diff_less_value)[i];
        total_edit_count_befr += (dcs_csp->cs_edit_count_befr)[i];
        total_subs_count_befr += (dcs_csp->cs_subs_count_befr)[i];
        total_inss_count_befr += (dcs_csp->cs_inss_count_befr)[i];
        total_dels_count_befr += (dcs_csp->cs_dels_count_befr)[i];
        total_fwd_locns_count += (dcs_csp->cs_fwd_locns_count)[i];
        total_bwd_locns_count += (dcs_csp->cs_bwd_locns_count)[i];
        total_fwd_diff_counts_count += (dcs_csp->cs_fwd_diff_counts_count)[i];
        total_bwd_diff_counts_count += (dcs_csp->cs_bwd_diff_counts_count)[i];
        total_fwd_diff_posns_count += (dcs_csp->cs_fwd_diff_posns_count)[i];
        total_bwd_diff_posns_count += (dcs_csp->cs_bwd_diff_posns_count)[i];
        total_fwd_diff_values_count += (dcs_csp->cs_fwd_diff_values_count)[i];
        total_bwd_diff_values_count += (dcs_csp->cs_bwd_diff_values_count)[i];
        total_pe_rel_locns_count += (dcs_csp->cs_pe_rel_locns_count)[i];
        total_pe_rel_posns_count += (dcs_csp->cs_pe_rel_posns_count)[i];
        total_edit_count_aftr += (dcs_csp->cs_edit_count_aftr)[i];
        total_subs_count_aftr += (dcs_csp->cs_subs_count_aftr)[i];
        total_inss_count_aftr += (dcs_csp->cs_inss_count_aftr)[i];
        total_dels_count_aftr += (dcs_csp->cs_dels_count_aftr)[i];
    }
    std::size_t total_differences_count_befr[16] = {0};
    for (std::uint32_t i = 0; i < in_args.threadCount; i++) {
        for (std::uint32_t j = 0; j < 16; j++) {
            total_differences_count_befr[j] += (dcs_csp->cs_differences_count_befr)[i*16 + j];
        }
    }
    std::size_t total_differences_count_aftr[16] = {0};
    for (std::uint32_t i = 0; i < in_args.threadCount; i++) {
        for (std::uint32_t j = 0; j < 16; j++) {
            total_differences_count_aftr[j] += (dcs_csp->cs_differences_count_aftr)[i*16 + j];
        }
    }
    
    fprintf(stderr, "total_cigar_resolved_count   : %lu, %.2lf.\n", total_cigar_resolved_count,
    (((double) total_cigar_resolved_count)/((double) dcs_csp->cs_pe_mapped_count))*100);
    fprintf(stderr, "total_diff_more_count        : %lu, %.2lf.\n", total_diff_more_count,
    (((double) total_diff_more_count)/((double) dcs_csp->cs_pe_mapped_count))*100);
    fprintf(stderr, "total_diff_more_value        : %lu.\n", total_diff_more_value);
    fprintf(stderr, "total_diff_same_count        : %lu, %.2lf.\n", total_diff_same_count,
    (((double) total_diff_same_count)/((double) dcs_csp->cs_pe_mapped_count))*100);
    fprintf(stderr, "total_diff_less_count        : %lu, %.2lf.\n", total_diff_less_count,
    (((double) total_diff_less_count)/((double) dcs_csp->cs_pe_mapped_count))*100);
    fprintf(stderr, "total_diff_less_value        : %lu.\n", total_diff_less_value);
    fprintf(stderr, "total_edit_count_befr        : %lu, %.2lf.\n", total_edit_count_befr,
    (((double) total_edit_count_befr)/((double) dcs_csp->cs_pe_mapped_count))*100);
    fprintf(stderr, "total_subs_count_befr        : %lu, %.2lf.\n", total_subs_count_befr,
    (((double) total_subs_count_befr)/((double) total_edit_count_befr))*100);
    fprintf(stderr, "total_inss_count_befr        : %lu, %.2lf.\n", total_inss_count_befr,
    (((double) total_inss_count_befr)/((double) total_edit_count_befr))*100);
    fprintf(stderr, "total_dels_count_befr        : %lu, %.2lf.\n", total_dels_count_befr,
    (((double) total_dels_count_befr)/((double) total_edit_count_befr))*100);
    fprintf(stderr, "total_edit_count_aftr        : %lu, %.2lf.\n", total_edit_count_aftr,
    (((double) total_edit_count_aftr)/((double) dcs_csp->cs_pe_mapped_count))*100);
    fprintf(stderr, "total_subs_count_aftr        : %lu, %.2lf.\n", total_subs_count_aftr,
    (((double) total_subs_count_aftr)/((double) total_edit_count_aftr))*100);
    fprintf(stderr, "total_inss_count_aftr        : %lu, %.2lf.\n", total_inss_count_aftr,
    (((double) total_inss_count_aftr)/((double) total_edit_count_aftr))*100);
    fprintf(stderr, "total_dels_count_aftr        : %lu, %.2lf.\n", total_dels_count_aftr,
    (((double) total_dels_count_aftr)/((double) total_edit_count_aftr))*100);
    fprintf(stderr, "total_fwd_locns_count        : %lu, %.2lf.\n", total_fwd_locns_count,
    (((double) total_fwd_locns_count)/((double) dcs_csp->cs_pe_mapped_count))*200);
    fprintf(stderr, "total_bwd_locns_count        : %lu, %.2lf.\n", total_bwd_locns_count,
    (((double) total_bwd_locns_count)/((double) dcs_csp->cs_pe_mapped_count))*200);
    fprintf(stderr, "total_fwd_diff_counts_count  : %lu, %.2lf.\n", total_fwd_diff_counts_count,
    (((double) total_fwd_diff_counts_count)/((double) dcs_csp->cs_pe_mapped_count))*200);
    fprintf(stderr, "total_bwd_diff_counts_count  : %lu, %.2lf.\n", total_bwd_diff_counts_count,
    (((double) total_bwd_diff_counts_count)/((double) dcs_csp->cs_pe_mapped_count))*200);
    fprintf(stderr, "total_fwd_diff_posns_count   : %lu, %.2lf.\n", total_fwd_diff_posns_count,
    (((double) total_fwd_diff_posns_count)/((double) dcs_csp->cs_pe_mapped_count))*200);
    fprintf(stderr, "total_bwd_diff_posns_count   : %lu, %.2lf.\n", total_bwd_diff_posns_count,
    (((double) total_bwd_diff_posns_count)/((double) dcs_csp->cs_pe_mapped_count))*200);
    fprintf(stderr, "total_fwd_diff_values_count  : %lu, %.2lf.\n", total_fwd_diff_values_count,
    (((double) total_fwd_diff_values_count)/((double) dcs_csp->cs_pe_mapped_count))*200);
    fprintf(stderr, "total_bwd_diff_values_count  : %lu, %.2lf.\n", total_bwd_diff_values_count,
    (((double) total_bwd_diff_values_count)/((double) dcs_csp->cs_pe_mapped_count))*200);
    fprintf(stderr, "total_pe_rel_locns_count     : %lu, %.2lf.\n", total_pe_rel_locns_count,
    (((double) total_pe_rel_locns_count)/((double) dcs_csp->cs_pe_mapped_count))*200);
    fprintf(stderr, "total_pe_rel_posns_count     : %lu, %.2lf.\n", total_pe_rel_posns_count,
    (((double) total_pe_rel_posns_count)/((double) dcs_csp->cs_pe_mapped_count))*200);

    for (std::uint32_t j = 0; j < 16; j++) {
        fprintf(stdout, "total_differences_count_bef[%2u] : %lu, %.2lf.\n", j, total_differences_count_befr[j],
        (((double) total_differences_count_befr[j])/((double) dcs_csp->cs_pe_mapped_count))*100);
    }
    for (std::uint32_t j = 0; j < 16; j++) {
        fprintf(stdout, "total_differences_count_aft[%2u] : %lu, %.2lf.\n", j, total_differences_count_aftr[j],
        (((double) total_differences_count_aftr[j])/((double) dcs_csp->cs_pe_mapped_count))*100);
    }
}

void *compactReads1Thread(void *arg) {
    struct CompressArgsForThread * tap;
    tap = (struct CompressArgsForThread *) arg;
    std::uint32_t indices_per_thread = ceil(((double) tap->caft_com_ds->unmapped_reads.size())/((double) tap->caft_in_args->threadCount));
    std::uint32_t index_start = (tap->caft_thread_id) * indices_per_thread;
    std::uint32_t index_end = ((tap->caft_thread_id+1) * indices_per_thread) - 1;
    if ((tap->caft_com_ds->unmapped_reads.size() - 1) < index_end) {
        index_end = tap->caft_com_ds->unmapped_reads.size() - 1;
    }
    //std::cout << tap->caft_thread_id << " : " << index_start << " : " << index_end << std::endl;
    
    std::uint8_t * b1 = tap->caft_sur + ((std::size_t) index_start) * ((tap->caft_in_args->rdLength) * 2);
    std::uint32_t ri, k;
    int z, one_base;
    for (std::uint32_t i = index_start; i <= index_end; i++) {
        UnmappedRead & u1 = (tap->caft_com_ds->unmapped_reads)[i];
        for (ri = 0, k = 0; ri < (tap->caft_in_args->rdLength); ri++) {
            one_base = 0;
            z = 1;
            if (TestBit(u1.fwd_read, k)) {
                one_base |= z;
            }
            k++;
            z = (z << 1);
            if (TestBit(u1.fwd_read, k)) {
                one_base |= z;
            }
            k++;
            z = (z << 1);
            if (TestBit(u1.fwd_read, k)) {
                one_base |= z;
            }
            k++;
            //z = (z << 1); //Not necessary
            (*b1) = (uint8_t) one_base;
            b1++;
        }
        for (ri = 0, k = 0; ri < (tap->caft_in_args->rdLength); ri++) {
            one_base = 0;
            z = 1;
            if (TestBit(u1.bwd_read, k)) {
                one_base |= z;
            }
            k++;
            z = (z << 1);
            if (TestBit(u1.bwd_read, k)) {
                one_base |= z;
            }
            k++;
            z = (z << 1);
            if (TestBit(u1.bwd_read, k)) {
                one_base |= z;
            }
            k++;
            //z = (z << 1); //Not necessary
            (*b1) = (uint8_t) one_base;
            b1++;
        }
    }
    
    pthread_barrier_wait(tap->caft_barrier); //TODO[LATER] : Just to be safe
    
    indices_per_thread = ceil(((double) tap->caft_com_ds->fwd_mapped_reads.size())/((double) tap->caft_in_args->threadCount));
    index_start = (tap->caft_thread_id) * indices_per_thread;
    index_end = ((tap->caft_thread_id+1) * indices_per_thread) - 1;
    if ((tap->caft_com_ds->fwd_mapped_reads.size() - 1) < index_end) {
        index_end = tap->caft_com_ds->fwd_mapped_reads.size() - 1;
    }
    //std::cout << tap->caft_thread_id << " : " << index_start << " : " << index_end << std::endl;
    
    for (std::uint32_t i = index_start; i <= index_end; i++) {
#if !NDEBUG
        assert ((tap->caft_com_ds->fwd_mapped_reads)[i].id == (tap->caft_com_ds->bwd_mapped_reads)[i].id);
#endif
        (tap->caft_com_ds->bwd_pairing_info)[i].id = (tap->caft_com_ds->bwd_mapped_reads)[i].id;
        (tap->caft_com_ds->bwd_pairing_info)[i].my_location = (tap->caft_com_ds->bwd_mapped_reads)[i].location;
        (tap->caft_com_ds->bwd_pairing_info)[i].pe_location = (tap->caft_com_ds->fwd_mapped_reads)[i].location;
    }
    
    return NULL;
}

void populate_compaction_data(CompressArgsForThread * pcd_tap, int* pcd_diff_locs, int* pcd_diff_vals) {
    std::uint32_t threadID = pcd_tap->caft_thread_id;
    if (pcd_tap->caft_is_fwd_read) {
        if (pcd_tap->caft_locn_diff < 253) {
            (pcd_tap->caft_fwd_locns).push_back((std::uint8_t) pcd_tap->caft_locn_diff);
        } else if (pcd_tap->caft_locn_diff < 65536) {
            (pcd_tap->caft_fwd_locns).push_back((std::uint8_t) 253);
            (pcd_tap->caft_fwd_locns).push_back((std::uint8_t) ((pcd_tap->caft_locn_diff      ) & 0x000000ff));
            (pcd_tap->caft_fwd_locns).push_back((std::uint8_t) ((pcd_tap->caft_locn_diff >> 8 ) & 0x000000ff));
        } else if (pcd_tap->caft_locn_diff < 16777216) {
            (pcd_tap->caft_fwd_locns).push_back((std::uint8_t) 254);
            (pcd_tap->caft_fwd_locns).push_back((std::uint8_t) ((pcd_tap->caft_locn_diff      ) & 0x000000ff));
            (pcd_tap->caft_fwd_locns).push_back((std::uint8_t) ((pcd_tap->caft_locn_diff >> 8 ) & 0x000000ff));
            (pcd_tap->caft_fwd_locns).push_back((std::uint8_t) ((pcd_tap->caft_locn_diff >> 16) & 0x000000ff));
        } else {
            (pcd_tap->caft_fwd_locns).push_back((std::uint8_t) 255);
            (pcd_tap->caft_fwd_locns).push_back((std::uint8_t) ((pcd_tap->caft_locn_diff      ) & 0x000000ff));
            (pcd_tap->caft_fwd_locns).push_back((std::uint8_t) ((pcd_tap->caft_locn_diff >> 8 ) & 0x000000ff));
            (pcd_tap->caft_fwd_locns).push_back((std::uint8_t) ((pcd_tap->caft_locn_diff >> 16) & 0x000000ff));
            (pcd_tap->caft_fwd_locns).push_back((std::uint8_t) ((pcd_tap->caft_locn_diff >> 24) & 0x000000ff));
        }
        (pcd_tap->caft_fwd_diff_counts).push_back((std::uint8_t) pcd_diff_locs[0]);
        int pcd_diff_count = pcd_diff_locs[0];
        pcd_diff_locs[0] = 0;
        for (int i = 1; i <= pcd_diff_count; i++) {
            (pcd_tap->caft_fwd_diff_posns).push_back((std::uint8_t) (pcd_diff_locs[i] - pcd_diff_locs[i-1]));
            (pcd_tap->caft_fwd_diff_values).push_back((std::uint8_t) pcd_diff_vals[i]);
            (pcd_tap->csp->cs_edit_count_aftr)[threadID] += 1;
            if (pcd_diff_vals[i] < 4) {
                (pcd_tap->csp->cs_subs_count_aftr)[threadID] += 1;
            } else if (pcd_diff_vals[i] < 9) {
                (pcd_tap->csp->cs_inss_count_aftr)[threadID] += 1;
            } else {
                (pcd_tap->csp->cs_dels_count_aftr)[threadID] += 1;
            }
        }
        pcd_diff_locs[0] = pcd_diff_count;
    } else {
        if (pcd_tap->caft_locn_diff < 253) {
            (pcd_tap->caft_bwd_locns).push_back((std::uint8_t) pcd_tap->caft_locn_diff);
        } else if (pcd_tap->caft_locn_diff < 65536) {
            (pcd_tap->caft_bwd_locns).push_back((std::uint8_t) 253);
            (pcd_tap->caft_bwd_locns).push_back((std::uint8_t) ((pcd_tap->caft_locn_diff      ) & 0x000000ff));
            (pcd_tap->caft_bwd_locns).push_back((std::uint8_t) ((pcd_tap->caft_locn_diff >> 8 ) & 0x000000ff));
        } else if (pcd_tap->caft_locn_diff < 16777216) {
            (pcd_tap->caft_bwd_locns).push_back((std::uint8_t) 254);
            (pcd_tap->caft_bwd_locns).push_back((std::uint8_t) ((pcd_tap->caft_locn_diff      ) & 0x000000ff));
            (pcd_tap->caft_bwd_locns).push_back((std::uint8_t) ((pcd_tap->caft_locn_diff >> 8 ) & 0x000000ff));
            (pcd_tap->caft_bwd_locns).push_back((std::uint8_t) ((pcd_tap->caft_locn_diff >> 16) & 0x000000ff));
        } else {
            (pcd_tap->caft_bwd_locns).push_back((std::uint8_t) 255);
            (pcd_tap->caft_bwd_locns).push_back((std::uint8_t) ((pcd_tap->caft_locn_diff      ) & 0x000000ff));
            (pcd_tap->caft_bwd_locns).push_back((std::uint8_t) ((pcd_tap->caft_locn_diff >> 8 ) & 0x000000ff));
            (pcd_tap->caft_bwd_locns).push_back((std::uint8_t) ((pcd_tap->caft_locn_diff >> 16) & 0x000000ff));
            (pcd_tap->caft_bwd_locns).push_back((std::uint8_t) ((pcd_tap->caft_locn_diff >> 24) & 0x000000ff));
        }
        (pcd_tap->caft_bwd_diff_counts).push_back((std::uint8_t) pcd_diff_locs[0]);
        int pcd_diff_count = pcd_diff_locs[0];
        pcd_diff_locs[0] = 0;
        for (int i = 1; i <= pcd_diff_count; i++) {
            (pcd_tap->caft_bwd_diff_posns).push_back((std::uint8_t) (pcd_diff_locs[i] - pcd_diff_locs[i-1]));
            (pcd_tap->caft_bwd_diff_values).push_back((std::uint8_t) pcd_diff_vals[i]);
            (pcd_tap->csp->cs_edit_count_aftr)[threadID] += 1;
            if (pcd_diff_vals[i] < 4) {
                (pcd_tap->csp->cs_subs_count_aftr)[threadID] += 1;
            } else if (pcd_diff_vals[i] < 9) {
                (pcd_tap->csp->cs_inss_count_aftr)[threadID] += 1;
            } else {
                (pcd_tap->csp->cs_dels_count_aftr)[threadID] += 1;
            }
        }
        pcd_diff_locs[0] = pcd_diff_count;
    }
}
/* For every forward two aligned read and every reverse two aligned read append (to the lists) : 
(data stream 1)relative start location – location wrt. previous read 
(data stream 2)number of differences (between read and reference genome) – one byte per count
(data stream 3)relative difference positions – one byte per difference position 
(data stream 4)differing bases – 2 bits to capture substitution to N (1 possibility), deletion (4 possibilities), 
insertion of N (1 possibility), insertion of regular base (4 possibilities) */
void *compactReads2Thread(void *arg) {
    struct CompressArgsForThread * tap;
    tap = (struct CompressArgsForThread *) arg;
    std::uint32_t indices_per_thread = ceil(((double) tap->caft_com_ds->fwd_mapped_reads.size())/((double) tap->caft_in_args->threadCount));
    std::uint32_t index_start = (tap->caft_thread_id) * indices_per_thread;
    std::uint32_t index_end = ((tap->caft_thread_id+1) * indices_per_thread) - 1;
    if ((tap->caft_com_ds->fwd_mapped_reads.size() - 1) < index_end) {
        index_end = tap->caft_com_ds->fwd_mapped_reads.size() - 1;
    }
    //std::cout << tap->caft_thread_id << " : " << index_start << " : " << index_end << std::endl;
    
    (tap->caft_fwd_locns).reserve(indices_per_thread * 1.1);
    (tap->caft_bwd_locns).reserve(indices_per_thread * 1.1);
    (tap->caft_fwd_diff_counts).reserve(indices_per_thread);
    (tap->caft_bwd_diff_counts).reserve(indices_per_thread);
    (tap->caft_fwd_diff_posns).reserve(indices_per_thread * 1.1);
    (tap->caft_bwd_diff_posns).reserve(indices_per_thread * 1.1);
    (tap->caft_fwd_diff_values).reserve(indices_per_thread * 1.1);
    (tap->caft_bwd_diff_values).reserve(indices_per_thread * 1.1);
    
    std::uint32_t local_diffs_fwd[READ_LEN_MAX*3][EDIT_DISTANCE][10];
    memset(local_diffs_fwd, 0, READ_LEN_MAX*3*EDIT_DISTANCE*10*sizeof(std::uint32_t));
    std::uint8_t  local_bases_fwd[READ_LEN_MAX*3][EDIT_DISTANCE];
    memset(local_bases_fwd, 254, READ_LEN_MAX*3*EDIT_DISTANCE*sizeof(std::uint8_t));
    std::uint32_t local_count_fwd[READ_LEN_MAX*3][EDIT_DISTANCE];
    memset(local_count_fwd, 0, READ_LEN_MAX*3*EDIT_DISTANCE*sizeof(std::uint32_t));
    std::uint32_t local_usage_fwd[READ_LEN_MAX*3];
    memset(local_usage_fwd, 0, READ_LEN_MAX*3*sizeof(std::uint32_t));
    std::uint32_t curr_ref_buf_pos_fwd = 0;
    std::uint32_t temp_ref_buf_pos_fwd = 0;
    int temp_diff_cnt_fwd = 0;
    int prev_diff_loc_fwd = -1;
    
    std::uint32_t local_diffs_bwd[READ_LEN_MAX*3][EDIT_DISTANCE][10];
    memset(local_diffs_bwd, 0, READ_LEN_MAX*3*EDIT_DISTANCE*10*sizeof(std::uint32_t));
    std::uint8_t  local_bases_bwd[READ_LEN_MAX*3][EDIT_DISTANCE];
    memset(local_bases_bwd, 254, READ_LEN_MAX*3*EDIT_DISTANCE*sizeof(std::uint8_t));
    std::uint32_t local_count_bwd[READ_LEN_MAX*3][EDIT_DISTANCE];
    memset(local_count_bwd, 0, READ_LEN_MAX*3*EDIT_DISTANCE*sizeof(std::uint32_t));
    std::uint32_t local_usage_bwd[READ_LEN_MAX*3];
    memset(local_usage_bwd, 0, READ_LEN_MAX*3*sizeof(std::uint32_t));
    std::uint32_t curr_ref_buf_pos_bwd = 0;
    std::uint32_t temp_ref_buf_pos_bwd = 0;
    int temp_diff_cnt_bwd = 0;
    int prev_diff_loc_bwd = -1;
    
    int k, z, my_diffs;
    char cigar[READ_LEN_MAX*2];
    int diff_vals[32], aftr_diff_vals[32];
    int diff_locs[32], aftr_diff_locs[32];
    MappedRead * m1 = NULL, * m2 = NULL;
    Read r1, r2;
    std::uint32_t fwd_prev_locn = 0, bwd_prev_locn = 0;
    std::uint32_t threadID = tap->caft_thread_id;
    std::uint8_t * my_ref_bases = tap->caft_com_ds->ref_bases;
    for (std::uint32_t i = index_start; i <= index_end; i++) {
        m1 = &(tap->caft_com_ds->fwd_mapped_reads[i]);
        my_diffs = 0;
        k = 32*READ_LEN_U32 - 4;
        z = 1;
        if (TestBit(m1->read, k)) {
            my_diffs |= z;
        }
        k++;
        z = (z << 1);
        if (TestBit(m1->read, k)) {
            my_diffs |= z;
        }
        k++;
        z = (z << 1);
        if (TestBit(m1->read, k)) {
            my_diffs |= z;
        }
        k++;
        z = (z << 1);
        if (TestBit(m1->read, k)) {
            my_diffs |= z;
        }
        assert (my_diffs < 16);
        if (strcmp(m1->tmp_cigar, "Y")) {
            (tap->csp->cs_cigar_resolved_count)[threadID] += 1;
            if (my_diffs) {
                uncompact_mapped_read(r1.bases, m1->read, tap->caft_in_args->rdLength);
                //r1.length = (int) tap->caft_in_args->rdLength;
                CPUMDTag(my_ref_bases, r1.bases, m1->location, m1->tmp_cigar, diff_locs, diff_vals);
                if (diff_locs[0] > my_diffs) {
                    (tap->csp->cs_diff_more_count)[threadID] += 1;
                    (tap->csp->cs_diff_more_value)[threadID] += (diff_locs[0] - my_diffs);
                } else if (diff_locs[0] == my_diffs) {
                    (tap->csp->cs_diff_same_count)[threadID] += 1;
                } else {
                    (tap->csp->cs_diff_less_count)[threadID] += 1;
                    (tap->csp->cs_diff_less_value)[threadID] += (my_diffs - diff_locs[0]);
                }
                (tap->csp->cs_differences_count_befr)[threadID*16 + diff_locs[0]]++;
            } else {
                (tap->csp->cs_differences_count_befr)[threadID*16]++;
                (tap->csp->cs_diff_same_count)[threadID] += 1;
                diff_locs[0] = 0;
            }
        } else {
            uncompact_mapped_read(r1.bases, m1->read, tap->caft_in_args->rdLength);
            r1.length = (int) tap->caft_in_args->rdLength;
            CPUAlignKernel((my_ref_bases + m1->base_location), r1.bases, m1->diff_location, cigar, r1.length, my_diffs);
            CPUMDTag(my_ref_bases, r1.bases, m1->location, cigar, diff_locs, diff_vals);
                if (diff_locs[0] > my_diffs) {
                    (tap->csp->cs_diff_more_count)[threadID] += 1;
                    (tap->csp->cs_diff_more_value)[threadID] += (diff_locs[0] - my_diffs);
                } else if (diff_locs[0] == my_diffs) {
                    (tap->csp->cs_diff_same_count)[threadID] += 1;
                } else {
                    (tap->csp->cs_diff_less_count)[threadID] += 1;
                    (tap->csp->cs_diff_less_value)[threadID] += (my_diffs - diff_locs[0]);
                }
                (tap->csp->cs_differences_count_befr)[threadID*16 + diff_locs[0]]++;
        }
        
        aftr_diff_locs[0] = diff_locs[0];
        for (int i = 1; i <= diff_locs[0]; i++) {
            aftr_diff_locs[i] = diff_locs[i];
            aftr_diff_vals[i] = diff_vals[i];
            (tap->csp->cs_edit_count_befr)[threadID] += 1;
            if (diff_vals[i] < 4) {
                (tap->csp->cs_subs_count_befr)[threadID] += 1;
            } else if (diff_vals[i] < 9) {
                (tap->csp->cs_inss_count_befr)[threadID] += 1;
            } else {
                (tap->csp->cs_dels_count_befr)[threadID] += 1;
            }
        }
        if (tap->caft_in_args->updateReference) {
            temp_diff_cnt_fwd = 0;
            prev_diff_loc_fwd = -1;
            if (fwd_prev_locn) {
                for (std::uint32_t i = fwd_prev_locn; i < m1->location; i++) {
                    temp_ref_buf_pos_fwd = i % (READ_LEN_MAX*3);
//                     assert (temp_ref_buf_pos_fwd < (READ_LEN_MAX*3));
//                     assert (local_usage_fwd[temp_ref_buf_pos_fwd] < EDIT_DISTANCE);
                    for (std::uint32_t j = 0; j < local_usage_fwd[temp_ref_buf_pos_fwd]; j++) {
                        local_bases_fwd[temp_ref_buf_pos_fwd][j] = 254;
                        local_count_fwd[temp_ref_buf_pos_fwd][j] = 0;
                        for (std::uint32_t k = 0; k < 10; k++) {
                            local_diffs_fwd[temp_ref_buf_pos_fwd][j][k] = 0;
                        }
                    }
                    local_usage_fwd[temp_ref_buf_pos_fwd] = 0;
                }
            }
            
            curr_ref_buf_pos_fwd = m1->location % (READ_LEN_MAX*3);
            aftr_diff_locs[0] = 0;
            
            for (int i = 1; i <= diff_locs[0]; i++) {
                temp_ref_buf_pos_fwd = (curr_ref_buf_pos_fwd + diff_locs[i]) % (READ_LEN_MAX*3);
//                 assert (temp_ref_buf_pos_fwd < (READ_LEN_MAX*3));
                if (diff_locs[i] == prev_diff_loc_fwd) {
                    temp_diff_cnt_fwd++;
                } else {
                    prev_diff_loc_fwd = diff_locs[i];
                    temp_diff_cnt_fwd = 1;
                }
                local_usage_fwd[temp_ref_buf_pos_fwd] = (std::uint32_t) temp_diff_cnt_fwd;
                
                if (local_bases_fwd[temp_ref_buf_pos_fwd][temp_diff_cnt_fwd - 1] == ((std::uint8_t) diff_vals[i])) {
                    local_diffs_fwd[temp_ref_buf_pos_fwd][temp_diff_cnt_fwd - 1][diff_vals[i]] += 1;
                    local_count_fwd[temp_ref_buf_pos_fwd][temp_diff_cnt_fwd - 1] += 1;
                } else {
                    aftr_diff_locs[0] += 1;
                    aftr_diff_locs[aftr_diff_locs[0]] = diff_locs[i];
                    aftr_diff_vals[aftr_diff_locs[0]] = diff_vals[i];
                    local_diffs_fwd[temp_ref_buf_pos_fwd][temp_diff_cnt_fwd - 1][diff_vals[i]] += 1;
                    if (local_diffs_fwd[temp_ref_buf_pos_fwd][temp_diff_cnt_fwd - 1][diff_vals[i]] > (local_count_fwd[temp_ref_buf_pos_fwd][temp_diff_cnt_fwd - 1])) {
                        local_count_fwd[temp_ref_buf_pos_fwd][temp_diff_cnt_fwd - 1] = local_diffs_fwd[temp_ref_buf_pos_fwd][temp_diff_cnt_fwd - 1][diff_vals[i]];
                        local_bases_fwd[temp_ref_buf_pos_fwd][temp_diff_cnt_fwd - 1] = ((std::uint8_t) diff_vals[i]);
                    }
                }
            }
        }
        
        tap->caft_is_fwd_read = true;
        tap->caft_locn_diff = m1->location - fwd_prev_locn;
        (tap->csp->cs_differences_count_aftr)[threadID*16 + aftr_diff_locs[0]]++;
        populate_compaction_data(tap, aftr_diff_locs, aftr_diff_vals);
        fwd_prev_locn = m1->location;
        if (diff_locs[0] > EDIT_DISTANCE) {
            std::cerr << "Fwd ID : " << m1->id << ", " << diff_locs[0] << std::endl;
        }
        
        m2 = &(tap->caft_com_ds->bwd_mapped_reads[i]);
        my_diffs = 0;
        k = 32*READ_LEN_U32 - 4;
        z = 1;
        if (TestBit(m2->read, k)) {
            my_diffs |= z;
        }
        k++;
        z = (z << 1);
        if (TestBit(m2->read, k)) {
            my_diffs |= z;
        }
        k++;
        z = (z << 1);
        if (TestBit(m2->read, k)) {
            my_diffs |= z;
        }
        k++;
        z = (z << 1);
        if (TestBit(m2->read, k)) {
            my_diffs |= z;
        }
        assert (my_diffs < 16);
        if (strcmp(m2->tmp_cigar, "Y")) {
            (tap->csp->cs_cigar_resolved_count)[threadID] += 1;
            if (my_diffs) {
                uncompact_mapped_read(r2.bases, m2->read, tap->caft_in_args->rdLength);
                //r2.length = (int) tap->caft_in_args->rdLength;
                CPUMDTag(my_ref_bases, r2.bases, m2->location, m2->tmp_cigar, diff_locs, diff_vals);
                if (diff_locs[0] > my_diffs) {
                    (tap->csp->cs_diff_more_count)[threadID] += 1;
                    (tap->csp->cs_diff_more_value)[threadID] += (diff_locs[0] - my_diffs);
                } else if (diff_locs[0] == my_diffs) {
                    (tap->csp->cs_diff_same_count)[threadID] += 1;
                } else {
                    (tap->csp->cs_diff_less_count)[threadID] += 1;
                    (tap->csp->cs_diff_less_value)[threadID] += (my_diffs - diff_locs[0]);
                }
                (tap->csp->cs_differences_count_befr)[threadID*16 + diff_locs[0]]++;
            } else {
                (tap->csp->cs_differences_count_befr)[threadID*16]++;
                (tap->csp->cs_diff_same_count)[threadID] += 1;
                diff_locs[0] = 0;
            }
        } else {
            uncompact_mapped_read(r2.bases, m2->read, tap->caft_in_args->rdLength);
            r2.length = (int) tap->caft_in_args->rdLength;
            CPUAlignKernel((my_ref_bases + m2->base_location), r2.bases, m2->diff_location, cigar, r2.length, my_diffs);
            CPUMDTag(my_ref_bases, r2.bases, m2->location, cigar, diff_locs, diff_vals);
                if (diff_locs[0] > my_diffs) {
                    (tap->csp->cs_diff_more_count)[threadID] += 1;
                    (tap->csp->cs_diff_more_value)[threadID] += (diff_locs[0] - my_diffs);
                } else if (diff_locs[0] == my_diffs) {
                    (tap->csp->cs_diff_same_count)[threadID] += 1;
                } else {
                    (tap->csp->cs_diff_less_count)[threadID] += 1;
                    (tap->csp->cs_diff_less_value)[threadID] += (my_diffs - diff_locs[0]);
                }
                (tap->csp->cs_differences_count_befr)[threadID*16 + diff_locs[0]]++;
        }
        
        aftr_diff_locs[0] = diff_locs[0];
        for (int i = 1; i <= diff_locs[0]; i++) {
            aftr_diff_locs[i] = diff_locs[i];
            aftr_diff_vals[i] = diff_vals[i];
            (tap->csp->cs_edit_count_befr)[threadID] += 1;
            if (diff_vals[i] < 4) {
                (tap->csp->cs_subs_count_befr)[threadID] += 1;
            } else if (diff_vals[i] < 9) {
                (tap->csp->cs_inss_count_befr)[threadID] += 1;
            } else {
                (tap->csp->cs_dels_count_befr)[threadID] += 1;
            }
        }
        if (tap->caft_in_args->updateReference) {
            temp_diff_cnt_bwd = 0;
            prev_diff_loc_bwd = -1;
            if (bwd_prev_locn) {
                for (std::uint32_t i = bwd_prev_locn; i < m2->location; i++) {
                    temp_ref_buf_pos_bwd = i % (READ_LEN_MAX*3);
//                     assert (temp_ref_buf_pos_bwd < (READ_LEN_MAX*3));
//                     assert (local_usage_bwd[temp_ref_buf_pos_bwd] < EDIT_DISTANCE);
                    for (std::uint32_t j = 0; j < local_usage_bwd[temp_ref_buf_pos_bwd]; j++) {
                        local_bases_bwd[temp_ref_buf_pos_bwd][j] = 254;
                        local_count_bwd[temp_ref_buf_pos_bwd][j] = 0;
                        for (std::uint32_t k = 0; k < 10; k++) {
                            local_diffs_bwd[temp_ref_buf_pos_bwd][j][k] = 0;
                        }
                    }
                    local_usage_bwd[temp_ref_buf_pos_bwd] = 0;
                }
            }
            
            curr_ref_buf_pos_bwd = m2->location % (READ_LEN_MAX*3);
            aftr_diff_locs[0] = 0;
            
            for (int i = 1; i <= diff_locs[0]; i++) {
                temp_ref_buf_pos_bwd = (curr_ref_buf_pos_bwd + diff_locs[i]) % (READ_LEN_MAX*3);
//                 assert (temp_ref_buf_pos_bwd < (READ_LEN_MAX*3));
                if (diff_locs[i] == prev_diff_loc_bwd) {
                    temp_diff_cnt_bwd++;
                } else {
                    prev_diff_loc_bwd = diff_locs[i];
                    temp_diff_cnt_bwd = 1;
                }
                local_usage_bwd[temp_ref_buf_pos_bwd] = (std::uint32_t) temp_diff_cnt_bwd;
                
                if (local_bases_bwd[temp_ref_buf_pos_bwd][temp_diff_cnt_bwd - 1] == ((std::uint8_t) diff_vals[i])) {
                    local_diffs_bwd[temp_ref_buf_pos_bwd][temp_diff_cnt_bwd - 1][diff_vals[i]] += 1;
                    local_count_bwd[temp_ref_buf_pos_bwd][temp_diff_cnt_bwd - 1] += 1;
                } else {
                    aftr_diff_locs[0] += 1;
                    aftr_diff_locs[aftr_diff_locs[0]] = diff_locs[i];
                    aftr_diff_vals[aftr_diff_locs[0]] = diff_vals[i];
                    local_diffs_bwd[temp_ref_buf_pos_bwd][temp_diff_cnt_bwd - 1][diff_vals[i]] += 1;
                    if (local_diffs_bwd[temp_ref_buf_pos_bwd][temp_diff_cnt_bwd - 1][diff_vals[i]] > (local_count_bwd[temp_ref_buf_pos_bwd][temp_diff_cnt_bwd - 1])) {
                        local_count_bwd[temp_ref_buf_pos_bwd][temp_diff_cnt_bwd - 1] = local_diffs_bwd[temp_ref_buf_pos_bwd][temp_diff_cnt_bwd - 1][diff_vals[i]];
                        local_bases_bwd[temp_ref_buf_pos_bwd][temp_diff_cnt_bwd - 1] = ((std::uint8_t) diff_vals[i]);
                    }
                }
            }
        }
        
        tap->caft_is_fwd_read = false;
        tap->caft_locn_diff = m2->location - bwd_prev_locn;
        (tap->csp->cs_differences_count_aftr)[threadID*16 + aftr_diff_locs[0]]++;
        populate_compaction_data(tap, aftr_diff_locs, aftr_diff_vals);
        bwd_prev_locn = m2->location;
        if (diff_locs[0] > EDIT_DISTANCE) {
            std::cerr << "Bwd ID : " << m2->id << ", " << diff_locs[0] << std::endl;
        }
    }
    
    fwd_prev_locn = std::numeric_limits<std::uint32_t>::max();
    std::uint32_t locn_rel_posn = 0;
    for (std::uint32_t i = index_start; i <= index_end; i++) {
        if ((tap->caft_com_ds->fwd_mapped_reads)[i].location == fwd_prev_locn) {
            locn_rel_posn++;
        } else {
            fwd_prev_locn = (tap->caft_com_ds->fwd_mapped_reads)[i].location;
            locn_rel_posn = 0;
        }
        (tap->caft_com_ds->fwd_pairing_info)[i].id = (tap->caft_com_ds->fwd_mapped_reads)[i].id;
        (tap->caft_com_ds->fwd_pairing_info)[i].my_rel_posn = locn_rel_posn;
    }
    
    pthread_barrier_wait(tap->caft_barrier);
    
    //Process first location again
    if (threadID) {
        fwd_prev_locn = (tap->caft_com_ds->fwd_mapped_reads)[index_start - 1].location;
        locn_rel_posn = (tap->caft_com_ds->fwd_pairing_info)[index_start - 1].my_rel_posn;
        for (std::uint32_t i = index_start; i <= index_end; i++) {
            if ((tap->caft_com_ds->fwd_mapped_reads)[i].location == fwd_prev_locn) {
                locn_rel_posn++;
            } else {
                break;
            }
            (tap->caft_com_ds->fwd_pairing_info)[i].my_rel_posn = locn_rel_posn;
        }
    }
    
    (tap->csp->cs_fwd_locns_count)[threadID] += tap->caft_fwd_locns.size();
    (tap->csp->cs_bwd_locns_count)[threadID] += tap->caft_bwd_locns.size();
    (tap->csp->cs_fwd_diff_counts_count)[threadID] += tap->caft_fwd_diff_counts.size();
    (tap->csp->cs_bwd_diff_counts_count)[threadID] += tap->caft_bwd_diff_counts.size();
    (tap->csp->cs_fwd_diff_posns_count)[threadID] += tap->caft_fwd_diff_posns.size();
    (tap->csp->cs_bwd_diff_posns_count)[threadID] += tap->caft_bwd_diff_posns.size();
    (tap->csp->cs_fwd_diff_values_count)[threadID] += tap->caft_fwd_diff_values.size();
    (tap->csp->cs_bwd_diff_values_count)[threadID] += tap->caft_bwd_diff_values.size();
    
/*//Will free immediately
    std::vector<std::uint8_t>(tap->caft_fwd_locns).swap(tap->caft_fwd_locns);
    std::vector<std::uint8_t>(tap->caft_bwd_locns).swap(tap->caft_bwd_locns);
    std::vector<std::uint8_t>(tap->caft_fwd_diff_counts).swap(tap->caft_fwd_diff_counts);
    std::vector<std::uint8_t>(tap->caft_bwd_diff_counts).swap(tap->caft_bwd_diff_counts);
    std::vector<std::uint8_t>(tap->caft_fwd_diff_posns).swap(tap->caft_fwd_diff_posns);
    std::vector<std::uint8_t>(tap->caft_bwd_diff_posns).swap(tap->caft_bwd_diff_posns);
    std::vector<std::uint8_t>(tap->caft_fwd_diff_values).swap(tap->caft_fwd_diff_values);
    std::vector<std::uint8_t>(tap->caft_bwd_diff_values).swap(tap->caft_bwd_diff_values);
*/
    
    return NULL;
}
/* Create a list of tuples – (location of forward read, location of reverse read, position 
of forward read)
 */
void *compactReads3Thread(void *arg) {
    struct CompressArgsForThread * tap;
    tap = (struct CompressArgsForThread *) arg;
    std::uint32_t indices_per_thread = ceil(((double) tap->caft_com_ds->bwd_pairing_info.size())/((double) tap->caft_in_args->threadCount));
    std::uint32_t index_start = (tap->caft_thread_id) * indices_per_thread;
    std::uint32_t index_end = ((tap->caft_thread_id+1) * indices_per_thread) - 1;
    if ((tap->caft_com_ds->bwd_pairing_info.size() - 1) < index_end) {
        index_end = tap->caft_com_ds->bwd_pairing_info.size() - 1;
    }
    //std::cout << tap->caft_thread_id << " : " << index_start << " : " << index_end << std::endl;
    
    int64_t my_rel_locn_a, my_rel_locn_b;
    tap->caft_pe_rel_locn_sum = 0;
    tap->caft_pe_rel_locn_cnt = 0;
    for (std::uint32_t i = index_start; i <= index_end; i++) {
#if !NDEBUG
        assert ((tap->caft_com_ds->bwd_pairing_info)[i].id == (tap->caft_com_ds->fwd_pairing_info)[i].id);
#endif
        (tap->caft_com_ds->bwd_pairing_info)[i].pe_rel_posn = (tap->caft_com_ds->fwd_pairing_info)[i].my_rel_posn;
        my_rel_locn_a = (int64_t) (tap->caft_com_ds->bwd_pairing_info)[i].my_location;
        my_rel_locn_b = (int64_t) (tap->caft_com_ds->bwd_pairing_info)[i].pe_location;
        if ((my_rel_locn_a - my_rel_locn_b > 0) && (my_rel_locn_a - my_rel_locn_b < 512)) {
            tap->caft_pe_rel_locn_sum += (my_rel_locn_a - my_rel_locn_b);
            tap->caft_pe_rel_locn_cnt += 1;
        }
    }
    return NULL;
}
/* (data stream 5)locations of other ends – for every reverse read capture relative 
location of its forward read */
void *compactReads4Thread(void *arg) {
    struct CompressArgsForThread * tap;
    tap = (struct CompressArgsForThread *) arg;
    std::uint32_t indices_per_thread = ceil(((double) tap->caft_com_ds->bwd_pairing_info.size())/((double) tap->caft_in_args->threadCount));
    std::uint32_t index_start = (tap->caft_thread_id) * indices_per_thread;
    std::uint32_t index_end = ((tap->caft_thread_id+1) * indices_per_thread) - 1;
    if ((tap->caft_com_ds->bwd_pairing_info.size() - 1) < index_end) {
        index_end = tap->caft_com_ds->bwd_pairing_info.size() - 1;
    }
    //std::cout << tap->caft_thread_id << " : " << index_start << " : " << index_end << std::endl;
    
    (tap->caft_pe_rel_posns).reserve(indices_per_thread * 1.1);
    (tap->caft_pe_rel_locns).reserve(indices_per_thread * 1.8);
    
    std::uint32_t my_rel_posn;
    int64_t my_rel_locn_1, my_rel_locn_a, my_rel_locn_b;
    std::uint64_t my_rel_locn_2;
    for (std::uint32_t i = index_start; i <= index_end; i++) {
        my_rel_posn = (tap->caft_com_ds->bwd_pairing_info)[i].pe_rel_posn;
        if (my_rel_posn < 253) {
            (tap->caft_pe_rel_posns).push_back((std::uint8_t) my_rel_posn);
        } else if (my_rel_posn < 65536) {
            (tap->caft_pe_rel_posns).push_back((std::uint8_t) 253);
            (tap->caft_pe_rel_posns).push_back((std::uint8_t) ((my_rel_posn      ) & 0x000000ff));
            (tap->caft_pe_rel_posns).push_back((std::uint8_t) ((my_rel_posn >> 8 ) & 0x000000ff));
        } else if (my_rel_posn < 16777216) {
            (tap->caft_pe_rel_posns).push_back((std::uint8_t) 254);
            (tap->caft_pe_rel_posns).push_back((std::uint8_t) ((my_rel_posn      ) & 0x000000ff));
            (tap->caft_pe_rel_posns).push_back((std::uint8_t) ((my_rel_posn >> 8 ) & 0x000000ff));
            (tap->caft_pe_rel_posns).push_back((std::uint8_t) ((my_rel_posn >> 16) & 0x000000ff));
        } else {
            (tap->caft_pe_rel_posns).push_back((std::uint8_t) 255);
            (tap->caft_pe_rel_posns).push_back((std::uint8_t) ((my_rel_posn      ) & 0x000000ff));
            (tap->caft_pe_rel_posns).push_back((std::uint8_t) ((my_rel_posn >> 8 ) & 0x000000ff));
            (tap->caft_pe_rel_posns).push_back((std::uint8_t) ((my_rel_posn >> 16) & 0x000000ff));
            (tap->caft_pe_rel_posns).push_back((std::uint8_t) ((my_rel_posn >> 24) & 0x000000ff));
        }
        
        my_rel_locn_a = (int64_t) (tap->caft_com_ds->bwd_pairing_info)[i].my_location;
        my_rel_locn_b = (int64_t) (tap->caft_com_ds->bwd_pairing_info)[i].pe_location;
        my_rel_locn_1 = (my_rel_locn_a - my_rel_locn_b);
        my_rel_locn_1 = my_rel_locn_1 - tap->caft_pe_rel_locn_mean;
        if (my_rel_locn_1 < 0) {
            my_rel_locn_1 = my_rel_locn_1 * (-2);
        } else if (my_rel_locn_1 > 0) {
            my_rel_locn_1 = my_rel_locn_1 * 2 - 1;
        } else {
            //Do nothing
        }
#if !NDEBUG
        assert ((my_rel_locn_1 >= 0) && (my_rel_locn_1 < (4294967296*4)));
#endif
        my_rel_locn_2 = (std::uint64_t) my_rel_locn_1;
        if (my_rel_locn_2 < 252) {
            (tap->caft_pe_rel_locns).push_back((std::uint8_t) my_rel_locn_2);
        } else if (my_rel_locn_2 < 65536) {
            (tap->caft_pe_rel_locns).push_back((std::uint8_t) 252);
            (tap->caft_pe_rel_locns).push_back((std::uint8_t) ((my_rel_locn_2      ) & 0x00000000000000ff));
            (tap->caft_pe_rel_locns).push_back((std::uint8_t) ((my_rel_locn_2 >> 8 ) & 0x00000000000000ff));
        } else if (my_rel_locn_2 < 16777216) {
            (tap->caft_pe_rel_locns).push_back((std::uint8_t) 253);
            (tap->caft_pe_rel_locns).push_back((std::uint8_t) ((my_rel_locn_2      ) & 0x00000000000000ff));
            (tap->caft_pe_rel_locns).push_back((std::uint8_t) ((my_rel_locn_2 >> 8 ) & 0x00000000000000ff));
            (tap->caft_pe_rel_locns).push_back((std::uint8_t) ((my_rel_locn_2 >> 16) & 0x00000000000000ff));
        } else if (my_rel_locn_2 < 4294967296) {
            (tap->caft_pe_rel_locns).push_back((std::uint8_t) 254);
            (tap->caft_pe_rel_locns).push_back((std::uint8_t) ((my_rel_locn_2      ) & 0x00000000000000ff));
            (tap->caft_pe_rel_locns).push_back((std::uint8_t) ((my_rel_locn_2 >> 8 ) & 0x00000000000000ff));
            (tap->caft_pe_rel_locns).push_back((std::uint8_t) ((my_rel_locn_2 >> 16) & 0x00000000000000ff));
            (tap->caft_pe_rel_locns).push_back((std::uint8_t) ((my_rel_locn_2 >> 24) & 0x00000000000000ff));
        } else {
            (tap->caft_pe_rel_locns).push_back((std::uint8_t) 255);
            (tap->caft_pe_rel_locns).push_back((std::uint8_t) ((my_rel_locn_2      ) & 0x00000000000000ff));
            (tap->caft_pe_rel_locns).push_back((std::uint8_t) ((my_rel_locn_2 >> 8 ) & 0x00000000000000ff));
            (tap->caft_pe_rel_locns).push_back((std::uint8_t) ((my_rel_locn_2 >> 16) & 0x00000000000000ff));
            (tap->caft_pe_rel_locns).push_back((std::uint8_t) ((my_rel_locn_2 >> 24) & 0x00000000000000ff));
            (tap->caft_pe_rel_locns).push_back((std::uint8_t) ((my_rel_locn_2 >> 32) & 0x00000000000000ff));
        }
    }
    
    (tap->csp->cs_pe_rel_posns_count)[tap->caft_thread_id] += tap->caft_pe_rel_posns.size();
    (tap->csp->cs_pe_rel_locns_count)[tap->caft_thread_id] += tap->caft_pe_rel_locns.size();
    return NULL;
}
/* write single end reads*/
void write_se_data(FILE * wsd_fp, CompressArgsForThread * wsd_tap) {
    std::size_t num_write;
	std::size_t my_size;
    for (std::uint32_t i = 0; i < wsd_tap[0].caft_in_args->threadCount; i++) {
	    my_size = wsd_tap[i].caft_fwd_locns.size();
	    num_write = fwrite(&(my_size), 1, sizeof(std::size_t), wsd_fp);
	    if (num_write) {
	        //std::cout << "Bytes written : " << num_write << std::endl;
	    } else {
	        std::cerr << "Error writing to .cpct file" << std::endl;
	        assert (0);
	    }
	    wsd_tap[i].caft_com_ds->totalBytes += num_write;
	    num_write = fwrite(wsd_tap[i].caft_fwd_locns.data(), 1, wsd_tap[i].caft_fwd_locns.size(), wsd_fp);
	    if (num_write) {
	        //std::cout << "Wrote " << num_write << " bytes to .cpct file" << std::endl;
	    } else {
	        std::cerr << "Error writing to .cpct file" << std::endl;
	        assert (0);
	    }
	    wsd_tap[i].caft_com_ds->totalBytes += num_write;
	    std::vector<std::uint8_t>().swap(wsd_tap[i].caft_fwd_locns); //Free memory
    }
    for (std::uint32_t i = 0; i < wsd_tap[0].caft_in_args->threadCount; i++) {
	    my_size = wsd_tap[i].caft_bwd_locns.size();
	    num_write = fwrite(&(my_size), 1, sizeof(std::size_t), wsd_fp);
	    if (num_write) {
	        //std::cout << "Bytes written : " << num_write << std::endl;
	    } else {
	        std::cerr << "Error writing to .cpct file" << std::endl;
	        assert (0);
	    }
	    wsd_tap[i].caft_com_ds->totalBytes += num_write;
	    num_write = fwrite(wsd_tap[i].caft_bwd_locns.data(), 1, wsd_tap[i].caft_bwd_locns.size(), wsd_fp);
	    if (num_write) {
	        //std::cout << "Wrote " << num_write << " bytes to .cpct file" << std::endl;
	    } else {
	        std::cerr << "Error writing to .cpct file" << std::endl;
	        assert (0);
	    }
	    wsd_tap[i].caft_com_ds->totalBytes += num_write;
	    std::vector<std::uint8_t>().swap(wsd_tap[i].caft_bwd_locns); //Free memory
    }
    for (std::uint32_t i = 0; i < wsd_tap[0].caft_in_args->threadCount; i++) {
	    my_size = wsd_tap[i].caft_fwd_diff_counts.size();
	    num_write = fwrite(&(my_size), 1, sizeof(std::size_t), wsd_fp);
	    if (num_write) {
	        //std::cout << "Bytes written : " << num_write << std::endl;
	    } else {
	        std::cerr << "Error writing to .cpct file" << std::endl;
	        assert (0);
	    }
	    wsd_tap[i].caft_com_ds->totalBytes += num_write;
	    num_write = fwrite(wsd_tap[i].caft_fwd_diff_counts.data(), 1, wsd_tap[i].caft_fwd_diff_counts.size(), wsd_fp);
	    if (num_write) {
	        //std::cout << "Wrote " << num_write << " bytes to .cpct file" << std::endl;
	    } else {
	        std::cerr << "Error writing to .cpct file" << std::endl;
	        assert (0);
	    }
	    wsd_tap[i].caft_com_ds->totalBytes += num_write;
	    std::vector<std::uint8_t>().swap(wsd_tap[i].caft_fwd_diff_counts); //Free memory
    }
    for (std::uint32_t i = 0; i < wsd_tap[0].caft_in_args->threadCount; i++) {
	    my_size = wsd_tap[i].caft_bwd_diff_counts.size();
	    num_write = fwrite(&(my_size), 1, sizeof(std::size_t), wsd_fp);
	    if (num_write) {
	        //std::cout << "Bytes written : " << num_write << std::endl;
	    } else {
	        std::cerr << "Error writing to .cpct file" << std::endl;
	        assert (0);
	    }
	    wsd_tap[i].caft_com_ds->totalBytes += num_write;
	    num_write = fwrite(wsd_tap[i].caft_bwd_diff_counts.data(), 1, wsd_tap[i].caft_bwd_diff_counts.size(), wsd_fp);
	    if (num_write) {
	        //std::cout << "Wrote " << num_write << " bytes to .cpct file" << std::endl;
	    } else {
	        std::cerr << "Error writing to .cpct file" << std::endl;
	        assert (0);
	    }
	    wsd_tap[i].caft_com_ds->totalBytes += num_write;
	    std::vector<std::uint8_t>().swap(wsd_tap[i].caft_bwd_diff_counts); //Free memory
    }
    for (std::uint32_t i = 0; i < wsd_tap[0].caft_in_args->threadCount; i++) {
	    my_size = wsd_tap[i].caft_fwd_diff_posns.size();
	    num_write = fwrite(&(my_size), 1, sizeof(std::size_t), wsd_fp);
	    if (num_write) {
	        //std::cout << "Bytes written : " << num_write << std::endl;
	    } else {
	        std::cerr << "Error writing to .cpct file" << std::endl;
	        assert (0);
	    }
	    wsd_tap[i].caft_com_ds->totalBytes += num_write;
	    num_write = fwrite(wsd_tap[i].caft_fwd_diff_posns.data(), 1, wsd_tap[i].caft_fwd_diff_posns.size(), wsd_fp);
	    if (num_write) {
	        //std::cout << "Wrote " << num_write << " bytes to .cpct file" << std::endl;
	    } else {
	        std::cerr << "Error writing to .cpct file" << std::endl;
	        assert (0);
	    }
	    wsd_tap[i].caft_com_ds->totalBytes += num_write;
	    std::vector<std::uint8_t>().swap(wsd_tap[i].caft_fwd_diff_posns); //Free memory
    }
    for (std::uint32_t i = 0; i < wsd_tap[0].caft_in_args->threadCount; i++) {
	    my_size = wsd_tap[i].caft_bwd_diff_posns.size();
	    num_write = fwrite(&(my_size), 1, sizeof(std::size_t), wsd_fp);
	    if (num_write) {
	        //std::cout << "Bytes written : " << num_write << std::endl;
	    } else {
	        std::cerr << "Error writing to .cpct file" << std::endl;
	        assert (0);
	    }
	    wsd_tap[i].caft_com_ds->totalBytes += num_write;
	    num_write = fwrite(wsd_tap[i].caft_bwd_diff_posns.data(), 1, wsd_tap[i].caft_bwd_diff_posns.size(), wsd_fp);
	    if (num_write) {
	        //std::cout << "Wrote " << num_write << " bytes to .cpct file" << std::endl;
	    } else {
	        std::cerr << "Error writing to .cpct file" << std::endl;
	        assert (0);
	    }
	    wsd_tap[i].caft_com_ds->totalBytes += num_write;
	    std::vector<std::uint8_t>().swap(wsd_tap[i].caft_bwd_diff_posns); //Free memory
    }
    for (std::uint32_t i = 0; i < wsd_tap[0].caft_in_args->threadCount; i++) {
	    my_size = wsd_tap[i].caft_fwd_diff_values.size();
	    num_write = fwrite(&(my_size), 1, sizeof(std::size_t), wsd_fp);
	    if (num_write) {
	        //std::cout << "Bytes written : " << num_write << std::endl;
	    } else {
	        std::cerr << "Error writing to .cpct file" << std::endl;
	        assert (0);
	    }
	    wsd_tap[i].caft_com_ds->totalBytes += num_write;
	    num_write = fwrite(wsd_tap[i].caft_fwd_diff_values.data(), 1, wsd_tap[i].caft_fwd_diff_values.size(), wsd_fp);
	    if (num_write) {
	        //std::cout << "Wrote " << num_write << " bytes to .cpct file" << std::endl;
	    } else {
	        std::cerr << "Error writing to .cpct file" << std::endl;
	        assert (0);
	    }
	    wsd_tap[i].caft_com_ds->totalBytes += num_write;
	    std::vector<std::uint8_t>().swap(wsd_tap[i].caft_fwd_diff_values); //Free memory
    }
    for (std::uint32_t i = 0; i < wsd_tap[0].caft_in_args->threadCount; i++) {
	    my_size = wsd_tap[i].caft_bwd_diff_values.size();
	    num_write = fwrite(&(my_size), 1, sizeof(std::size_t), wsd_fp);
	    if (num_write) {
	        //std::cout << "Bytes written : " << num_write << std::endl;
	    } else {
	        std::cerr << "Error writing to .cpct file" << std::endl;
	        assert (0);
	    }
	    wsd_tap[i].caft_com_ds->totalBytes += num_write;
	    num_write = fwrite(wsd_tap[i].caft_bwd_diff_values.data(), 1, wsd_tap[i].caft_bwd_diff_values.size(), wsd_fp);
	    if (num_write) {
	        //std::cout << "Wrote " << num_write << " bytes to .cpct file" << std::endl;
	    } else {
	        std::cerr << "Error writing to .cpct file" << std::endl;
	        assert (0);
	    }
	    wsd_tap[i].caft_com_ds->totalBytes += num_write;
	    std::vector<std::uint8_t>().swap(wsd_tap[i].caft_bwd_diff_values); //Free memory
    }
}

void write_pe_data(FILE * wpd_fp, CompressArgsForThread * wpd_tap) {
    std::size_t num_write;
    std::size_t my_size;
    for (std::uint32_t i = 0; i < wpd_tap[0].caft_in_args->threadCount; i++) {
	    my_size = wpd_tap[i].caft_pe_rel_locns.size();
	    num_write = fwrite(&(my_size), 1, sizeof(std::size_t), wpd_fp);
	    if (num_write) {
	        //std::cout << "Bytes written : " << num_write << std::endl;
	    } else {
	        std::cerr << "Error writing to .cpct file" << std::endl;
	        assert (0);
	    }
        if (wpd_tap[0].caft_in_args->writeSepFiles) {
            wpd_tap[i].caft_com_ds->peBytes += num_write;
        } else {
	        wpd_tap[i].caft_com_ds->totalBytes += num_write;
        }
	    num_write = fwrite(wpd_tap[i].caft_pe_rel_locns.data(), 1, wpd_tap[i].caft_pe_rel_locns.size(), wpd_fp);
	    if (num_write) {
	        //std::cout << "Wrote " << num_write << " bytes to .cpct file" << std::endl;
	    } else {
	        std::cerr << "Error writing to .cpct file" << std::endl;
	        assert (0);
	    }
        if (wpd_tap[0].caft_in_args->writeSepFiles) {
            wpd_tap[i].caft_com_ds->peBytes += num_write;
        } else {
	        wpd_tap[i].caft_com_ds->totalBytes += num_write;
        }
	    std::vector<std::uint8_t>().swap(wpd_tap[i].caft_pe_rel_locns); //Free memory
    }
    for (std::uint32_t i = 0; i < wpd_tap[0].caft_in_args->threadCount; i++) {
	    my_size = wpd_tap[i].caft_pe_rel_posns.size();
	    num_write = fwrite(&(my_size), 1, sizeof(std::size_t), wpd_fp);
	    if (num_write) {
	        //std::cout << "Bytes written : " << num_write << std::endl;
	    } else {
	        std::cerr << "Error writing to .cpct file" << std::endl;
	        assert (0);
	    }
        if (wpd_tap[0].caft_in_args->writeSepFiles) {
            wpd_tap[i].caft_com_ds->peBytes += num_write;
        } else {
	        wpd_tap[i].caft_com_ds->totalBytes += num_write;
        }
	    num_write = fwrite(wpd_tap[i].caft_pe_rel_posns.data(), 1, wpd_tap[i].caft_pe_rel_posns.size(), wpd_fp);
	    if (num_write) {
	        //std::cout << "Wrote " << num_write << " bytes to .cpct file" << std::endl;
	    } else {
	        std::cerr << "Error writing to .cpct file" << std::endl;
	        assert (0);
	    }
        if (wpd_tap[0].caft_in_args->writeSepFiles) {
            wpd_tap[i].caft_com_ds->peBytes += num_write;
        } else {
	        wpd_tap[i].caft_com_ds->totalBytes += num_write;
        }
	    std::vector<std::uint8_t>().swap(wpd_tap[i].caft_pe_rel_posns); //Free memory
    }
}

__UTILS_UINT8_CHAR
void write_compact_unm_op(InputArgs& in_args, CompressionDataStructures& comDS) {
	FILE * jnk_fp_1 = fopen("./rd_1.jnk.compact", "w");
	if (jnk_fp_1 == NULL) {
	    std::cerr << "Cannnot open file for writing : " << "./rd_1.jnk.align" << std::endl;
	    assert (0);
	}
	FILE * jnk_fp_2 = fopen("./rd_2.jnk.compact", "w");
	if (jnk_fp_2 == NULL) {
	    std::cerr << "Cannnot open file for writing : " << "./rd_2.jnk.align" << std::endl;
	    assert (0);
	}
    
    int z, one_base, j;
    std::uint32_t ri, k;
    char read_1[READ_LEN_MAX], read_2[READ_LEN_MAX];
    
    for (std::size_t i = 0; i < comDS.unmapped_reads.size(); i++) {
        UnmappedRead & u1 = comDS.unmapped_reads[i];
        
        j = 0;
        for (ri = 0, k = 0; ri < in_args.rdLength; ri++) {
            one_base = 0;
            z = 1;
            if (TestBit(u1.fwd_read, k)) {
                one_base |= z;
            }
            k++;
            z = (z << 1);
            if (TestBit(u1.fwd_read, k)) {
                one_base |= z;
            }
            k++;
            z = (z << 1);
            if (TestBit(u1.fwd_read, k)) {
                one_base |= z;
            }
            k++;
            //z = (z << 1); //Not necessary
            read_1[j] = Uint8Tochar((uint8_t) one_base);
            j++;
        }
        read_1[j] = '\0';
        
        j = (int) in_args.rdLength;
        read_2[j] = '\0';
        j--;
        for (ri = 0, k = 0; ri < in_args.rdLength; ri++) {
            one_base = 0;
            z = 1;
            if (TestBit(u1.bwd_read, k)) {
                one_base |= z;
            }
            k++;
            z = (z << 1);
            if (TestBit(u1.bwd_read, k)) {
                one_base |= z;
            }
            k++;
            z = (z << 1);
            if (TestBit(u1.bwd_read, k)) {
                one_base |= z;
            }
            k++;
            //z = (z << 1); //Not necessary
            if (one_base != 4) {
                one_base = 3 - one_base;
            }
            assert (one_base < 5);
            read_2[j] = Uint8Tochar((uint8_t) one_base);
            j--;
        }
        assert (j == -1);
        
        fprintf(jnk_fp_1, ">\n");
        fprintf(jnk_fp_1, "%s\n", read_1);
        fprintf(jnk_fp_2, ">\n");
        fprintf(jnk_fp_2, "%s\n", read_2);
    }
    
    fclose(jnk_fp_1);
    fclose(jnk_fp_2);
    return;
}
void write_compact_map_op(InputArgs& in_args, CompressionDataStructures& comDS) {
	FILE * jnk_fp_1 = fopen("./rd_1.jnk.compact", "a");
	if (jnk_fp_1 == NULL) {
	    std::cerr << "Cannnot open file for writing : " << "./rd_1.jnk.align" << std::endl;
	    assert (0);
	}
	FILE * jnk_fp_2 = fopen("./rd_2.jnk.compact", "a");
	if (jnk_fp_2 == NULL) {
	    std::cerr << "Cannnot open file for writing : " << "./rd_2.jnk.align" << std::endl;
	    assert (0);
	}
    
    int z, one_base, j;
    std::uint32_t ri, k;
    char read_1[READ_LEN_MAX], read_2[READ_LEN_MAX];
    
    int my_diffs_1, my_diffs_2;
    for (std::size_t i = 0; i < comDS.fwd_mapped_reads.size(); i++) {
        MappedRead & mf = comDS.fwd_mapped_reads[i];
        MappedRead & mb = comDS.bwd_mapped_reads[i];
        
        my_diffs_1 = 0;
        k = 32*READ_LEN_U32 - 4;
        z = 1;
        if (TestBit(mf.read, k)) {
            my_diffs_1 |= z;
        }
        k++;
        z = (z << 1);
        if (TestBit(mf.read, k)) {
            my_diffs_1 |= z;
        }
        k++;
        z = (z << 1);
        if (TestBit(mf.read, k)) {
            my_diffs_1 |= z;
        }
        k++;
        z = (z << 1);
        if (TestBit(mf.read, k)) {
            my_diffs_1 |= z;
        }
        assert (my_diffs_1 < 16);
        
        j = 0;
        for (ri = 0, k = 0; ri < in_args.rdLength; ri++) {
            one_base = 0;
            z = 1;
            if (TestBit(mf.read, k)) {
                one_base |= z;
            }
            k++;
            z = (z << 1);
            if (TestBit(mf.read, k)) {
                one_base |= z;
            }
            k++;
            z = (z << 1);
            if (TestBit(mf.read, k)) {
                one_base |= z;
            }
            k++;
            //z = (z << 1); //Not necessary
            read_1[j] = Uint8Tochar((uint8_t) one_base);
            j++;
        }
        read_1[j] = '\0';
        
        my_diffs_2 = 0;
        k = 32*READ_LEN_U32 - 4;
        z = 1;
        if (TestBit(mf.read, k)) {
            my_diffs_2 |= z;
        }
        k++;
        z = (z << 1);
        if (TestBit(mf.read, k)) {
            my_diffs_2 |= z;
        }
        k++;
        z = (z << 1);
        if (TestBit(mf.read, k)) {
            my_diffs_2 |= z;
        }
        k++;
        z = (z << 1);
        if (TestBit(mf.read, k)) {
            my_diffs_2 |= z;
        }
        assert (my_diffs_2 < 16);
        
        j = (int) in_args.rdLength;
        read_2[j] = '\0';
        j--;
        for (ri = 0, k = 0; ri < in_args.rdLength; ri++) {
            one_base = 0;
            z = 1;
            if (TestBit(mb.read, k)) {
                one_base |= z;
            }
            k++;
            z = (z << 1);
            if (TestBit(mb.read, k)) {
                one_base |= z;
            }
            k++;
            z = (z << 1);
            if (TestBit(mb.read, k)) {
                one_base |= z;
            }
            k++;
            //z = (z << 1); //Not necessary
            if (one_base != 4) {
                one_base = 3 - one_base;
            }
            assert (one_base < 5);
            read_2[j] = Uint8Tochar((uint8_t) one_base);
            j--;
        }
        assert (j == -1);
        
        fprintf(jnk_fp_1, "> %u, %u, %d, %u, %s, %d\n", mf.location, mf.base_location, mf.diff_location, mf.id, mf.tmp_cigar, my_diffs_1);
        fprintf(jnk_fp_1, "%s\n", read_1);
        fprintf(jnk_fp_2, "> %u, %u, %d, %u, %s, %d\n", mb.location, mb.base_location, mb.diff_location, mb.id, mb.tmp_cigar, my_diffs_2);
        fprintf(jnk_fp_2, "%s\n", read_2);
    }
    
    fclose(jnk_fp_1);
    fclose(jnk_fp_2);
    return;
}

int perform_compaction(InputArgs& in_args, CompressionDataStructures& comDS) {
	std::size_t num_write;
    std::string unm_file_name(in_args.comFileName);
    unm_file_name.append(".unm");
	std::string cpct_file_name(in_args.comFileName);
	cpct_file_name.append(".cpct");
    std::string pe_file_name(in_args.comFileName);
    pe_file_name.append(".pe");
	FILE * unm_file_fp = fopen(unm_file_name.c_str(), "w");
	if (unm_file_fp == NULL) {
	    std::cerr << "Cannnot open file for writing : " << unm_file_name << std::endl;
	    assert (0);
	}
	FILE * cpct_file_fp = fopen(cpct_file_name.c_str(), "w");
	if (cpct_file_fp == NULL) {
	    std::cerr << "Cannnot open file for writing : " << cpct_file_name << std::endl;
	    assert (0);
	}
	FILE * pe_file_fp = fopen(pe_file_name.c_str(), "w");
	if (pe_file_fp == NULL) {
	    std::cerr << "Cannnot open file for writing : " << pe_file_name << std::endl;
	    assert (0);
	}
	
// #if !NDEBUG
//     in_args.threadCount = 2;
// #endif
// 	
    CompressionStatistics pc_cs;
    pc_cs.cs_pe_mapped_count = comDS.fwd_mapped_reads.size() * 2;
    CompressArgsForThread thread_args_array[in_args.threadCount];
    pthread_t CPUTaskHandle[in_args.threadCount];
    pthread_barrier_t my_barrier;
    std::uint32_t finished_thread_num; int err;
    
    pthread_barrier_init(&my_barrier, NULL, in_args.threadCount);
	std::uint8_t * save_unmapped_reads = (std::uint8_t *) 
	malloc (comDS.unmapped_reads.size() * (in_args.rdLength * 2) * sizeof(std::uint8_t));
	assert (save_unmapped_reads);
	comDS.fwd_pairing_info.resize(comDS.fwd_mapped_reads.size());
	comDS.bwd_pairing_info.resize(comDS.bwd_mapped_reads.size());
    for (std::uint32_t i = 0; i < in_args.threadCount; i++) {
        thread_args_array[i].caft_thread_id = i;
        thread_args_array[i].caft_in_args = &in_args;
        thread_args_array[i].caft_com_ds = &comDS;
        thread_args_array[i].caft_barrier = &my_barrier;
        thread_args_array[i].csp = &pc_cs;
        //thread_args_array[i].caft_fp = &cpct_file_fp;
        thread_args_array[i].caft_sur = save_unmapped_reads;
    }
    initialize_compression_stats(in_args, &pc_cs); //Also, allocates capacity
    
    std::cout << "Before unmapped sort" << std::endl;
    pss::parallel_stable_sort(comDS.unmapped_reads.begin(), comDS.unmapped_reads.end(), 
    [](const UnmappedRead& x, const UnmappedRead& y) {
        return (x.location < y.location);
    });
    std::cout << "After unmapped sort" << std::endl;
        
    finished_thread_num = 0; err = 0;
    for (std::uint32_t i = 0; i < in_args.threadCount; i++) {
        err += pthread_create(CPUTaskHandle + i, NULL, compactReads1Thread, thread_args_array + i);
    }
    if (err == 0) {
        std::cout << "Created threads for compact-reads-1 successfully" << std::endl;
    } else {
        assert (0);
    }
    
    for (std::uint32_t i = 0; i < in_args.threadCount; i++) {
        err += pthread_join(CPUTaskHandle[i], NULL);
        finished_thread_num++;
    }
    if (err == 0) {
        std::cout << "Threads for compact-reads-1 completed successfully" << std::endl;
    } else {
        assert (0);
    }
    assert (finished_thread_num == in_args.threadCount);
	
	num_write = fwrite(&(in_args.updateReference), 1, sizeof(std::uint32_t), cpct_file_fp);
	if (num_write) {
	    //std::cout << "Wrote update reference to " << cpct_file_name << std::endl;
	    //std::cout << "Bytes written : " << num_write << std::endl;
	} else {
	    std::cerr << "Error writing to " << cpct_file_name << std::endl;
	    assert (0);
	}
	comDS.totalBytes += num_write;
	num_write = fwrite(&(in_args.threadCount), 1, sizeof(std::uint32_t), cpct_file_fp);
	if (num_write) {
	    //std::cout << "Wrote thread count to " << cpct_file_name << std::endl;
	    //std::cout << "Bytes written : " << num_write << std::endl;
	} else {
	    std::cerr << "Error writing to " << cpct_file_name << std::endl;
	    assert (0);
	}
	comDS.totalBytes += num_write;
	num_write = fwrite(&(in_args.rdLength), 1, sizeof(std::uint32_t), cpct_file_fp);
	if (num_write) {
	    //std::cout << "Wrote read length to " << cpct_file_name << std::endl;
	    //std::cout << "Bytes written : " << num_write << std::endl;
	} else {
	    std::cerr << "Error writing to " << cpct_file_name << std::endl;
	    assert (0);
	}
	comDS.totalBytes += num_write;
	std::size_t my_size = comDS.fwd_mapped_reads.size();
	num_write = fwrite(&(my_size), 1, sizeof(std::size_t), cpct_file_fp);
	if (num_write) {
	    //std::cout << "Wrote mapped read count to " << cpct_file_name << std::endl;
	    //std::cout << "Bytes written : " << num_write << std::endl;
	} else {
	    std::cerr << "Error writing to " << cpct_file_name << std::endl;
	    assert (0);
	}
	comDS.totalBytes += num_write;
	my_size = comDS.unmapped_reads.size();
	num_write = fwrite(&(my_size), 1, sizeof(std::size_t), in_args.writeSepFiles ? unm_file_fp : cpct_file_fp);
	if (num_write) {
	    //std::cout << "Wrote unmapped read count to " << cpct_file_name << std::endl;
	    //std::cout << "Bytes written : " << num_write << std::endl;
	} else {
	    std::cerr << "Error writing to " << (in_args.writeSepFiles ? unm_file_name : cpct_file_name) << std::endl;
	    assert (0);
	}
    if (in_args.writeSepFiles) {
        comDS.unmBytes += num_write;
    } else {
	    comDS.totalBytes += num_write;
    }
	num_write = fwrite(save_unmapped_reads, 1, comDS.unmapped_reads.size() * (in_args.rdLength * 2), in_args.writeSepFiles ? unm_file_fp : cpct_file_fp);
	if (num_write) {
	    std::cout << "Wrote " << num_write << " bytes to " << (in_args.writeSepFiles ? unm_file_name : cpct_file_name) << std::endl;
	} else {
	    std::cerr << "Error writing to " << (in_args.writeSepFiles ? unm_file_name : cpct_file_name) << std::endl;
	    assert (0);
	}
//     
// #if !NDEBUG
//     write_compact_unm_op(in_args, comDS);
// #endif
// 
    if (in_args.writeSepFiles) {
        comDS.unmBytes += num_write;
    } else {
	    comDS.totalBytes += num_write;
    }
	free (save_unmapped_reads);
	std::vector<UnmappedRead>().swap(comDS.unmapped_reads); //Free memory
    
    std::cout << "Before mapped fwd sort" << std::endl;
    //TODO[LATER] : Use differences count during sorting
    pss::parallel_stable_sort(comDS.fwd_mapped_reads.begin(), comDS.fwd_mapped_reads.end(), 
    [](const MappedRead& x, const MappedRead& y) {
        return (x.location < y.location) ||
        ((x.location == y.location) && (x.id < y.id));
    });
    std::cout << "After mapped fwd sort" << std::endl;
    
    //double startTime = realtime();
    std::cout << "Before mapped bwd sort" << std::endl;
    //TODO[LATER] : Use differences count during sorting
    pss::parallel_stable_sort(comDS.bwd_mapped_reads.begin(), comDS.bwd_mapped_reads.end(), 
    [](const MappedRead& x, const MappedRead& y) {
        return (x.location < y.location) ||
        ((x.location == y.location) && (x.id < y.id));
    });
    std::cout << "After mapped bwd sort" << std::endl;
    //std::cout << "Sorted bwd reads in " << realtime() - startTime << " s." << std::endl;
    
    finished_thread_num = 0; err = 0;
    for (std::uint32_t i = 0; i < in_args.threadCount; i++) {
        err += pthread_create(CPUTaskHandle + i, NULL, compactReads2Thread, thread_args_array + i);
    }
    if (err == 0) {
        std::cout << "Created threads for compact-reads-2 successfully" << std::endl;
    } else {
        assert (0);
    }
    
    for (std::uint32_t i = 0; i < in_args.threadCount; i++) {
        err += pthread_join(CPUTaskHandle[i], NULL);
        finished_thread_num++;
    }
    if (err == 0) {
        std::cout << "Threads for compact-reads-2 completed successfully" << std::endl;
    } else {
        assert (0);
    }
    assert (finished_thread_num == in_args.threadCount);
    
    free(comDS.ref_bases); comDS.ref_length = 0; //Free memory
//     
// #if !NDEBUG
//     write_compact_map_op(in_args, comDS);
// #endif
// 
    std::vector<MappedRead>().swap(comDS.fwd_mapped_reads); //Free memory
    std::vector<MappedRead>().swap(comDS.bwd_mapped_reads); //Free memory
    write_se_data(cpct_file_fp, thread_args_array); //Also, free memory
    
    std::cout << "Before pairing fwd sort" << std::endl;
    //TODO[LATER] : Use differences count during sorting
    pss::parallel_stable_sort(comDS.fwd_pairing_info.begin(), comDS.fwd_pairing_info.end(), 
    [](const PairingInfoFwd& x, const PairingInfoFwd& y) {
        return (x.id < y.id);
    });
    std::cout << "After pairing fwd sort" << std::endl;
    finished_thread_num = 0; err = 0;
    for (std::uint32_t i = 0; i < in_args.threadCount; i++) {
        err += pthread_create(CPUTaskHandle + i, NULL, compactReads3Thread, thread_args_array + i);
    }
    if (err == 0) {
        std::cout << "Created threads for compact-reads-3 successfully" << std::endl;
    } else {
        assert (0);
    }
    
    for (std::uint32_t i = 0; i < in_args.threadCount; i++) {
        err += pthread_join(CPUTaskHandle[i], NULL);
        finished_thread_num++;
    }
    if (err == 0) {
        std::cout << "Threads for compact-reads-3 completed successfully" << std::endl;
    } else {
        assert (0);
    }
    assert (finished_thread_num == in_args.threadCount);
    std::cout << "Before pairing bwd sort" << std::endl;
    //TODO[LATER] : Use differences count during sorting
    pss::parallel_stable_sort(comDS.bwd_pairing_info.begin(), comDS.bwd_pairing_info.end(), 
    [](const PairingInfoBwd& x, const PairingInfoBwd& y) {
        return (x.my_location < y.my_location) ||
        ((x.my_location == y.my_location) && (x.id < y.id));
    });
    std::cout << "After pairing bwd sort" << std::endl;
    
    int64_t pe_rel_locn_sum = 0, pe_rel_locn_cnt = 0;
    for (std::uint32_t i = 0; i < in_args.threadCount; i++) {
        pe_rel_locn_sum += thread_args_array[i].caft_pe_rel_locn_sum;
        pe_rel_locn_cnt += thread_args_array[i].caft_pe_rel_locn_cnt;
    }
    int64_t pe_rel_locn_mean = pe_rel_locn_sum/pe_rel_locn_cnt;
    for (std::uint32_t i = 0; i < in_args.threadCount; i++) {
        thread_args_array[i].caft_pe_rel_locn_mean = pe_rel_locn_mean;
    }
#if !NDEBUG
    std::cout << "pe_rel_locn_mean : " << thread_args_array[0].caft_pe_rel_locn_mean << std::endl;
#endif
	num_write = fwrite(&(pe_rel_locn_mean), 1, sizeof(int64_t), in_args.writeSepFiles ? pe_file_fp : cpct_file_fp);
	if (num_write) {
	    //std::cout << "Wrote pe_rel_locn_mean to " << cpct_file_name << std::endl;
	    //std::cout << "Bytes written : " << num_write << std::endl;
	} else {
	    std::cerr << "Error writing to " << (in_args.writeSepFiles ? pe_file_name : cpct_file_name) << std::endl;
	    assert (0);
	}
    if (in_args.writeSepFiles) {
        comDS.peBytes += num_write;
    } else {
	    comDS.totalBytes += num_write;
    }
    
    finished_thread_num = 0; err = 0;
    for (std::uint32_t i = 0; i < in_args.threadCount; i++) {
        err += pthread_create(CPUTaskHandle + i, NULL, compactReads4Thread, thread_args_array + i);
    }
    if (err == 0) {
        std::cout << "Created threads for compact-reads-4 successfully" << std::endl;
    } else {
        assert (0);
    }
    
    for (std::uint32_t i = 0; i < in_args.threadCount; i++) {
        err += pthread_join(CPUTaskHandle[i], NULL);
        finished_thread_num++;
    }
    if (err == 0) {
        std::cout << "Threads for compact-reads-4 completed successfully" << std::endl;
    } else {
        assert (0);
    }
    assert (finished_thread_num == in_args.threadCount);
    
    std::vector<PairingInfoFwd>().swap(comDS.fwd_pairing_info); //Free memory
    std::vector<PairingInfoBwd>().swap(comDS.bwd_pairing_info); //Free memory
    write_pe_data(in_args.writeSepFiles ? pe_file_fp : cpct_file_fp, thread_args_array); //Also, free memory
    
    fclose(unm_file_fp);
    fclose(cpct_file_fp);
    fclose(pe_file_fp);
#if !NDEBUG
    display_compression_stats(in_args, &pc_cs);
#endif
	
	return 0;
}

int compress_reads(InputArgs& in_args, CompressionDataStructures& comDS) {
	double startTime;
    
    startTime = realtime();
    perform_compaction(in_args, comDS);
    comDS.totalTime += (realtime() - startTime);
    std::cout << "Compacted reads in " << realtime() - startTime << " s." << std::endl;
    std::cout << "CPU time : " << cputime() << " s." << std::endl;
    std::string my_run_cmd;
    if (in_args.writeSepFiles) {
        std::cout << "unmBytes after compaction:   " << comDS.unmBytes << std::endl;
        std::cout << "peBytes after compaction:    " << comDS.peBytes << std::endl;
    }
    std::cout << "totalBytes after compaction: " << comDS.totalBytes << std::endl;
    
	std::string unm_file_name(in_args.comFileName);
	unm_file_name.append(".unm");
    if (in_args.writeSepFiles) {
	    startTime = realtime();
        my_run_cmd = in_args.bscExecutable;
        my_run_cmd.append(" e ");
        my_run_cmd.append(unm_file_name);
        my_run_cmd.append(" ");
        my_run_cmd.append(in_args.comFileName);
        my_run_cmd.append(" -p -b");
        std::uint32_t myTC = in_args.threadCount;
        double tmpblockSize = ((double) comDS.unmBytes)/(myTC * 1000000);
        std::cout << "Initial block size for bsc: " << (std::uint32_t) tmpblockSize << std::endl;
        while (tmpblockSize > 1024) {
            myTC += in_args.threadCount;
            tmpblockSize = ((double) comDS.unmBytes)/(myTC * 1000000);
        }
        std::uint32_t blockSize = (std::uint32_t) tmpblockSize;
        if (blockSize == 0) blockSize = 1;
        std::cout << "Final block size for bsc  : " << blockSize << std::endl;
        my_run_cmd.append(std::to_string(blockSize));
        std::cout << my_run_cmd << std::endl;
        system((const char *) my_run_cmd.c_str());
        comDS.totalTime += (realtime() - startTime);
        std::cout << "Compressed unm using bsc in " << realtime() - startTime << " s." << std::endl;
        std::cout << "CPU time : " << cputime() << " s." << std::endl;
    }
    my_run_cmd = "ls -ltr ";
    my_run_cmd.append(unm_file_name);
#if !NDEBUG
    std::cout << my_run_cmd << std::endl;
    system((const char *) my_run_cmd.c_str());
#endif
    my_run_cmd = "rm -rf ";
    my_run_cmd.append(unm_file_name);
#if !NDEBUG
    std::cout << my_run_cmd << std::endl;
#endif
    system((const char *) my_run_cmd.c_str());
    
	std::string pe_file_name(in_args.comFileName);
	pe_file_name.append(".pe");
    if (in_args.writeSepFiles) {
	    startTime = realtime();
        my_run_cmd = in_args.bscExecutable;
        my_run_cmd.append(" e ");
        my_run_cmd.append(pe_file_name);
        my_run_cmd.append(" ");
        my_run_cmd.append(in_args.comFileName);
        my_run_cmd.append(" -p -b");
        std::uint32_t myTC = in_args.threadCount;
        double tmpblockSize = ((double) comDS.peBytes)/(myTC * 1000000);
        std::cout << "Initial block size for bsc: " << (std::uint32_t) tmpblockSize << std::endl;
        while (tmpblockSize > 1024) {
            myTC += in_args.threadCount;
            tmpblockSize = ((double) comDS.peBytes)/(myTC * 1000000);
        }
        std::uint32_t blockSize = (std::uint32_t) tmpblockSize;
        if (blockSize == 0) blockSize = 1;
        std::cout << "Final block size for bsc  : " << blockSize << std::endl;
        my_run_cmd.append(std::to_string(blockSize));
        std::cout << my_run_cmd << std::endl;
        system((const char *) my_run_cmd.c_str());
        comDS.totalTime += (realtime() - startTime);
        std::cout << "Compressed pe using bsc in " << realtime() - startTime << " s." << std::endl;
        std::cout << "CPU time : " << cputime() << " s." << std::endl;
    }
    my_run_cmd = "ls -ltr ";
    my_run_cmd.append(pe_file_name);
#if !NDEBUG
    std::cout << my_run_cmd << std::endl;
    system((const char *) my_run_cmd.c_str());
#endif
    my_run_cmd = "rm -rf ";
    my_run_cmd.append(pe_file_name);
#if !NDEBUG
    std::cout << my_run_cmd << std::endl;
#endif
    system((const char *) my_run_cmd.c_str());
    
	startTime = realtime();
	std::string cpct_file_name(in_args.comFileName);
	cpct_file_name.append(".cpct");
    my_run_cmd = in_args.bscExecutable;
    my_run_cmd.append(" e ");
    my_run_cmd.append(cpct_file_name);
    my_run_cmd.append(" ");
    my_run_cmd.append(in_args.comFileName);
    my_run_cmd.append(" -p -b");
    std::uint32_t myTC = in_args.threadCount;
    double tmpblockSize = ((double) comDS.totalBytes)/(myTC * 1000000);
    std::cout << "Initial block size for bsc: " << (std::uint32_t) tmpblockSize << std::endl;
    while (tmpblockSize > 1024) {
        myTC += in_args.threadCount;
        tmpblockSize = ((double) comDS.totalBytes)/(myTC * 1000000);
    }
    std::uint32_t blockSize = (std::uint32_t) tmpblockSize;
    if (blockSize == 0) blockSize = 1;
    std::cout << "Final block size for bsc  : " << blockSize << std::endl;
    my_run_cmd.append(std::to_string(blockSize));
    std::cout << my_run_cmd << std::endl;
    system((const char *) my_run_cmd.c_str());
    comDS.totalTime += (realtime() - startTime);
    std::cout << "Compressed using bsc in " << realtime() - startTime << " s." << std::endl;
    std::cout << "CPU time : " << cputime() << " s." << std::endl;
    my_run_cmd = "ls -ltr ";
    my_run_cmd.append(cpct_file_name);
#if !NDEBUG
    std::cout << my_run_cmd << std::endl;
    system((const char *) my_run_cmd.c_str());
#endif
    my_run_cmd = "rm -rf ";
    my_run_cmd.append(cpct_file_name);
#if !NDEBUG
    std::cout << my_run_cmd << std::endl;
#endif
    system((const char *) my_run_cmd.c_str());
    
	return 0;
}






















































