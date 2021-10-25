#include <omp.h>
#include "decompress.hpp"
//software.intel.com/en-us/articles/a-parallel-stable-sort-using-c11-for-tbb-cilk-plus-and-openmp
#include "parallel_stable_sort.h"

#define SetBit(A,k)   (A[(k/32)] |=  (1 << (k%32)))
#define TestBit(A,k)  (A[(k/32)] &   (1 << (k%32)))
#define ClearBit(A,k) (A[(k/32)] &= ~(1 << (k%32)))
KSEQ_INIT(int, read)
__UTILS_CHAR_UINT8
__UTILS_UINT8_CHAR

inline void compact_bwd_read (char * my_bases, std::uint32_t my_rl, std::uint32_t * my_compact_read) {
            std::uint32_t ri = 0, k = 0;
            for (; ri < my_rl; ri++) {
                if ((my_bases[ri] == 'A') || (my_bases[ri] == 'a')) {
                    k++;
                    k++;
                    k++;
                } else if ((my_bases[ri] == 'C') || (my_bases[ri] == 'c')) {
                    SetBit(my_compact_read, k);
                    k++;
                    k++;
                    k++;
                } else if ((my_bases[ri] == 'G') || (my_bases[ri] == 'g')) {
                    k++;
                    SetBit(my_compact_read, k);
                    k++;
                    k++;
                } else if ((my_bases[ri] == 'T') || (my_bases[ri] == 't')) {
                    SetBit(my_compact_read, k);
                    k++;
                    SetBit(my_compact_read, k);
                    k++;
                    k++;
                } else {
                    // N base
                    k++;
                    k++;
                    SetBit(my_compact_read, k);
                    k++;
                }
            }
}

void initialize_decompression_stats(InputArgs& in_args, DecompressionStatistics * ids_dsp) {
    (ids_dsp->ds_edit_count_befr).resize(in_args.threadCount, 0);
    (ids_dsp->ds_subs_count_befr).resize(in_args.threadCount, 0);
    (ids_dsp->ds_inss_count_befr).resize(in_args.threadCount, 0);
    (ids_dsp->ds_dels_count_befr).resize(in_args.threadCount, 0);
    (ids_dsp->ds_fwd_locns_count).resize(in_args.threadCount, 0);
    (ids_dsp->ds_bwd_locns_count).resize(in_args.threadCount, 0);
    (ids_dsp->ds_fwd_diff_counts_count).resize(in_args.threadCount, 0);
    (ids_dsp->ds_bwd_diff_counts_count).resize(in_args.threadCount, 0);
    (ids_dsp->ds_fwd_diff_posns_count).resize(in_args.threadCount, 0);
    (ids_dsp->ds_bwd_diff_posns_count).resize(in_args.threadCount, 0);
    (ids_dsp->ds_fwd_diff_values_count).resize(in_args.threadCount, 0);
    (ids_dsp->ds_bwd_diff_values_count).resize(in_args.threadCount, 0);
    (ids_dsp->ds_pe_rel_locns_count).resize(in_args.threadCount, 0);
    (ids_dsp->ds_pe_rel_posns_count).resize(in_args.threadCount, 0);
    
    (ids_dsp->ds_differences_count_befr).resize((in_args.threadCount * 16), 0);
}

void display_decompression_stats(InputArgs& in_args, DecompressionStatistics * dds_dsp) {
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
    for (std::uint32_t i = 0; i < in_args.threadCount; i++) {
        total_edit_count_befr += (dds_dsp->ds_edit_count_befr)[i];
        total_subs_count_befr += (dds_dsp->ds_subs_count_befr)[i];
        total_inss_count_befr += (dds_dsp->ds_inss_count_befr)[i];
        total_dels_count_befr += (dds_dsp->ds_dels_count_befr)[i];
        total_fwd_locns_count += (dds_dsp->ds_fwd_locns_count)[i];
        total_bwd_locns_count += (dds_dsp->ds_bwd_locns_count)[i];
        total_fwd_diff_counts_count += (dds_dsp->ds_fwd_diff_counts_count)[i];
        total_bwd_diff_counts_count += (dds_dsp->ds_bwd_diff_counts_count)[i];
        total_fwd_diff_posns_count += (dds_dsp->ds_fwd_diff_posns_count)[i];
        total_bwd_diff_posns_count += (dds_dsp->ds_bwd_diff_posns_count)[i];
        total_fwd_diff_values_count += (dds_dsp->ds_fwd_diff_values_count)[i];
        total_bwd_diff_values_count += (dds_dsp->ds_bwd_diff_values_count)[i];
        total_pe_rel_locns_count += (dds_dsp->ds_pe_rel_locns_count)[i];
        total_pe_rel_posns_count += (dds_dsp->ds_pe_rel_posns_count)[i];
    }
    std::size_t total_differences_count_befr[16] = {0};
    for (std::uint32_t i = 0; i < in_args.threadCount; i++) {
        for (std::uint32_t j = 0; j < 16; j++) {
            total_differences_count_befr[j] += (dds_dsp->ds_differences_count_befr)[i*16 + j];
        }
    }
    
/*
    fprintf(stderr, "total_edit_count_befr        : %lu, %.2lf.\n", total_edit_count_befr,
    (((double) total_edit_count_befr)/((double) dds_dsp->ds_pe_mapped_count))*100);
    fprintf(stderr, "total_subs_count_befr        : %lu, %.2lf.\n", total_subs_count_befr,
    (((double) total_subs_count_befr)/((double) total_edit_count_befr))*100);
    fprintf(stderr, "total_inss_count_befr        : %lu, %.2lf.\n", total_inss_count_befr,
    (((double) total_inss_count_befr)/((double) total_edit_count_befr))*100);
    fprintf(stderr, "total_dels_count_befr        : %lu, %.2lf.\n", total_dels_count_befr,
    (((double) total_dels_count_befr)/((double) total_edit_count_befr))*100);
*/
    fprintf(stderr, "total_fwd_locns_count        : %lu, %.2lf.\n", total_fwd_locns_count,
    (((double) total_fwd_locns_count)/((double) dds_dsp->ds_pe_mapped_count))*200);
    fprintf(stderr, "total_bwd_locns_count        : %lu, %.2lf.\n", total_bwd_locns_count,
    (((double) total_bwd_locns_count)/((double) dds_dsp->ds_pe_mapped_count))*200);
    fprintf(stderr, "total_fwd_diff_counts_count  : %lu, %.2lf.\n", total_fwd_diff_counts_count,
    (((double) total_fwd_diff_counts_count)/((double) dds_dsp->ds_pe_mapped_count))*200);
    fprintf(stderr, "total_bwd_diff_counts_count  : %lu, %.2lf.\n", total_bwd_diff_counts_count,
    (((double) total_bwd_diff_counts_count)/((double) dds_dsp->ds_pe_mapped_count))*200);
    fprintf(stderr, "total_fwd_diff_posns_count   : %lu, %.2lf.\n", total_fwd_diff_posns_count,
    (((double) total_fwd_diff_posns_count)/((double) dds_dsp->ds_pe_mapped_count))*200);
    fprintf(stderr, "total_bwd_diff_posns_count   : %lu, %.2lf.\n", total_bwd_diff_posns_count,
    (((double) total_bwd_diff_posns_count)/((double) dds_dsp->ds_pe_mapped_count))*200);
    fprintf(stderr, "total_fwd_diff_values_count  : %lu, %.2lf.\n", total_fwd_diff_values_count,
    (((double) total_fwd_diff_values_count)/((double) dds_dsp->ds_pe_mapped_count))*200);
    fprintf(stderr, "total_bwd_diff_values_count  : %lu, %.2lf.\n", total_bwd_diff_values_count,
    (((double) total_bwd_diff_values_count)/((double) dds_dsp->ds_pe_mapped_count))*200);
    fprintf(stderr, "total_pe_rel_locns_count     : %lu, %.2lf.\n", total_pe_rel_locns_count,
    (((double) total_pe_rel_locns_count)/((double) dds_dsp->ds_pe_mapped_count))*200);
    fprintf(stderr, "total_pe_rel_posns_count     : %lu, %.2lf.\n", total_pe_rel_posns_count,
    (((double) total_pe_rel_posns_count)/((double) dds_dsp->ds_pe_mapped_count))*200);

/*
    for (std::uint32_t j = 0; j < 16; j++) {
        fprintf(stdout, "total_differences_count_bef[%2u] : %lu, %.2lf.\n", j, total_differences_count_befr[j],
        (((double) total_differences_count_befr[j])/((double) dds_dsp->ds_pe_mapped_count))*100);
    }
*/
}

int parse_reference_decom(InputArgs& in_args, DecompressionDataStructures& decomDS) {
	assert (decomDS.ref_length == 0);
	assert (decomDS.ref_bases == NULL);
	//decomDS.ref_bases = (char *) aligned_alloc(64, sizeof(char) * REF_LEN_MAX);
	decomDS.ref_bases = (char *) malloc(sizeof(char) * REF_LEN_MAX);
	assert (decomDS.ref_bases);
	
	FILE * ref_file_fp = fopen(in_args.refFileName.c_str(), "r");
	if (ref_file_fp == NULL) {
	    std::cerr << "Cannnot open file for reading : " << in_args.refFileName << std::endl;
	    assert (0);
	}
	kseq_t * ref_seq;
	ref_seq = kseq_init(fileno(ref_file_fp));
	
	std::uint32_t scaffold_count  = 0;
	std::uint32_t all_bases_count = 0;
	std::uint32_t val_bases_count = 0;
	int l = kseq_read(ref_seq);
	while (l == 0) {
		l = kseq_read(ref_seq);
	}
	while (l > 0) {
        scaffold_count++;
		int i;
		for (i = 0; i < l; i++) {
			char c = ref_seq->seq.s[i];
			if (charToUint8(c) < 4) {
			    decomDS.ref_bases[val_bases_count] = c;
			    val_bases_count++;
			    if (val_bases_count == REF_LEN_MAX) {
			        std::cerr << "Reference sequence longer than currently supported value" << std::endl;
			        assert (0);
			    }
			}
			all_bases_count++;
		}
		l = kseq_read(ref_seq);
	}
    if (l == -1) {
		//End of file
	} else {
		std::cerr << "Truncated quality string" << std::endl;
        assert (0);
	}
	
	decomDS.ref_length = val_bases_count;
	std::cout << "Reference scaffold_count  : " << scaffold_count  << std::endl;
	std::cout << "Reference all_bases_count : " << all_bases_count << std::endl;
	std::cout << "Reference val_bases_count : " << val_bases_count << std::endl;
    
    kseq_destroy(ref_seq);
    fclose(ref_file_fp);
    
    return 0;
}
//process unmapped reads (0-align and 1-align)
int process_unmapped_reads(InputArgs& in_args, DecompressionDataStructures& decomDS) {
    std::size_t num_read, my_size;
	num_read = fread(&(my_size), 1, sizeof(std::size_t), *(decomDS.ip_fp));
	if (num_read) {
	    std::cout << "Unmapped read count : " << my_size << std::endl;
	    //std::cout << "Bytes read          : " << num_read << std::endl;
	} else {
	    std::cerr << "Error reading from " << *(decomDS.ip_fp) << std::endl;
	    assert (0);
	}
	decomDS.totalBytes += num_read;
	//std::uint8_t * uread_bytes = (std::uint8_t *) 
	//malloc (my_size * (in_args.rdLength * 2) * sizeof(std::uint8_t));
	//assert (uread_bytes);
	//std::uint8_t * save_uread_bytes = uread_bytes;
	decomDS.unmapped_count = my_size;
	decomDS.unm_bytes.resize(my_size * (in_args.rdLength * 2));
	num_read = fread(decomDS.unm_bytes.data(), 1, my_size * (in_args.rdLength * 2), *(decomDS.ip_fp));
	if (num_read) {
	    std::cout << "Read " << num_read << " bytes from " << *(decomDS.ip_fp) << std::endl;
	} else {
	    std::cerr << "Error reading from " << *(decomDS.ip_fp) << std::endl;
	    assert (0);
	}
	decomDS.totalBytes += num_read;
	decomDS.unm_fwd_reads.resize(decomDS.unmapped_count * (in_args.rdLength + 3));
	decomDS.unm_bwd_reads.resize(decomDS.unmapped_count * (in_args.rdLength + 3));
	
/*
	std::uint8_t my_rc_bases[READ_LEN_MAX];
	for (std::size_t i = 0; i < my_size; i++) {
	    fprintf(*(decomDS.o1_fp), ">%lu\n", (decomDS.r1_count + 1));
	    fprintf(*(decomDS.o2_fp), ">%lu\n", (decomDS.r2_count + 1));
	    for (std::uint32_t j = 0; j < in_args.rdLength; j++) {
	        fputc(Uint8Tochar(*uread_bytes), *(decomDS.o1_fp));
	        uread_bytes++;
	    }
	    for (std::uint32_t j = 0; j < in_args.rdLength; j++) {
			if ((*uread_bytes) != 4) {
				my_rc_bases[in_args.rdLength - 1 - j] = 3 - (*uread_bytes);
			} else {
				my_rc_bases[in_args.rdLength - 1 - j] = 4;
			}
	        uread_bytes++;
	    }
	    for (std::uint32_t j = 0; j < in_args.rdLength; j++) {
	        fputc(Uint8Tochar(my_rc_bases[j]), *(decomDS.o2_fp));
	    }
	    fputc('\n', *(decomDS.o1_fp));
	    fputc('\n', *(decomDS.o2_fp));
	    decomDS.r1_count += 1;
	    decomDS.r2_count += 1;
	}
	
	free (save_uread_bytes);
*/
    return 0;
}
//write unmapped data
int write_unm_data(InputArgs& in_args, DecompressionDataStructures& decomDS) {
    std::vector<std::uint8_t>().swap(decomDS.unm_bytes); //Free memory
    std::size_t num_write;
    
    num_write = fwrite(decomDS.unm_fwd_reads.data(), sizeof(char), decomDS.unm_fwd_reads.size(), *(decomDS.o1_fp));
    if (num_write) {
	    //std::cout << "Wrote " << num_write << " bytes to rd1 file" << std::endl;
#if !NDEBUG
	    std::cout << "To write " << decomDS.unm_fwd_reads.size() << " bytes to rd1 file" << std::endl;
        std::cout << "Wrote    " << num_write << " bytes to rd1 file" << std::endl;
#endif
    } else {
	    std::cerr << "Error writing to rd1 file" << std::endl;
	    assert (0);
    }
    std::vector<char>().swap(decomDS.unm_fwd_reads); //Free memory
    
    num_write = fwrite(decomDS.unm_bwd_reads.data(), sizeof(char), decomDS.unm_bwd_reads.size(), *(decomDS.o2_fp));
    if (num_write) {
	    //std::cout << "Wrote " << num_write << " bytes to rd2 file" << std::endl;
#if !NDEBUG
	    std::cout << "To write " << decomDS.unm_bwd_reads.size() << " bytes to rd2 file" << std::endl;
        std::cout << "Wrote    " << num_write << " bytes to rd2 file" << std::endl;
#endif
    } else {
	    std::cerr << "Error writing to rd2 file" << std::endl;
	    assert (0);
    }
    std::vector<char>().swap(decomDS.unm_bwd_reads); //Free memory
    
    return 0;
}

void *decompactReads0Thread(void *arg) {
    struct DecompressArgsForThread * tap;
    tap = (struct DecompressArgsForThread *) arg;
    std::uint32_t indices_per_thread = ceil(((double) tap->daft_decom_ds->unmapped_count)/((double) tap->daft_in_args->threadCount));
    std::uint32_t index_start = (tap->daft_thread_id) * indices_per_thread;
    std::uint32_t index_end = ((tap->daft_thread_id+1) * indices_per_thread) - 1;
    if ((tap->daft_decom_ds->unmapped_count - 1) < index_end) {
        index_end = tap->daft_decom_ds->unmapped_count - 1;
    }
    //std::cout << tap->daft_thread_id << " : " << index_start << " : " << index_end << std::endl;
    
    std::size_t ip_idx = ((std::size_t) index_start) * (tap->daft_in_args->rdLength * 2);
    std::size_t o1_idx = ((std::size_t) index_start) * (tap->daft_in_args->rdLength + 3);
    std::size_t o2_idx = ((std::size_t) index_start) * (tap->daft_in_args->rdLength + 3);
    std::uint8_t my_rc_bases[READ_LEN_MAX];
    for (std::uint32_t i = index_start; i <= index_end; i++) {
        tap->daft_decom_ds->unm_fwd_reads.data()[o1_idx] = '>'; o1_idx += 1;
        tap->daft_decom_ds->unm_fwd_reads.data()[o1_idx] = '\n'; o1_idx += 1;
        tap->daft_decom_ds->unm_bwd_reads.data()[o2_idx] = '>'; o2_idx += 1;
        tap->daft_decom_ds->unm_bwd_reads.data()[o2_idx] = '\n'; o2_idx += 1;
	    for (std::uint32_t j = 0; j < tap->daft_in_args->rdLength; j++) {
	        tap->daft_decom_ds->unm_fwd_reads.data()[o1_idx] = Uint8Tochar(tap->daft_decom_ds->unm_bytes.data()[ip_idx]);
	        ip_idx += 1; o1_idx += 1;
	    }
	    for (std::uint32_t j = 0; j < tap->daft_in_args->rdLength; j++) {
			if (tap->daft_decom_ds->unm_bytes.data()[ip_idx] != 4) {
				my_rc_bases[tap->daft_in_args->rdLength - 1 - j] = 3 - tap->daft_decom_ds->unm_bytes.data()[ip_idx];
			} else {
				my_rc_bases[tap->daft_in_args->rdLength - 1 - j] = 4;
			}
	        ip_idx += 1;
	    }
	    for (std::uint32_t j = 0; j < tap->daft_in_args->rdLength; j++) {
	        tap->daft_decom_ds->unm_bwd_reads.data()[o2_idx] = Uint8Tochar(my_rc_bases[j]);
	        o2_idx += 1;
	    }
        tap->daft_decom_ds->unm_fwd_reads.data()[o1_idx] = '\n'; o1_idx += 1;
        tap->daft_decom_ds->unm_bwd_reads.data()[o2_idx] = '\n'; o2_idx += 1;
    }
    
    return NULL;
}
//write to rd2 file (fastq 2 file)
int write_bwd_data(InputArgs& in_args, DecompressionDataStructures& decomDS) {
/*
	std::uint8_t my_bases[READ_LEN_MAX], my_rc_bases[READ_LEN_MAX];
	std::uint32_t ri, k;
	std::uint8_t z, one_base;
	
	for (std::size_t i = 0; i < decomDS.bwd_reads.size(); i++) {
        DecomBwdRead & u1 = decomDS.bwd_reads[i];
        for (ri = 0, k = 0; ri < in_args.rdLength; ri++) {
            one_base = 0;
            z = 1;
            if (TestBit(u1.read, k)) {
                one_base |= z;
            }
            k++;
            z = (z << 1);
            if (TestBit(u1.read, k)) {
                one_base |= z;
            }
            k++;
            z = (z << 1);
            if (TestBit(u1.read, k)) {
                one_base |= z;
            }
            k++;
            //z = (z << 1); //Not necessary
            my_bases[ri] = one_base;
        }
	    fprintf(*(decomDS.o2_fp), ">%lu\n", (decomDS.r2_count + 1));
	    for (std::uint32_t j = 0; j < in_args.rdLength; j++) {
			if (my_bases[j] != 4) {
				my_rc_bases[in_args.rdLength - 1 - j] = 3 - my_bases[j];
			} else {
				my_rc_bases[in_args.rdLength - 1 - j] = 4;
			}
	    }
	    for (std::uint32_t j = 0; j < in_args.rdLength; j++) {
	        fputc(Uint8Tochar(my_rc_bases[j]), *(decomDS.o2_fp));
	    }
	    fputc('\n', *(decomDS.o2_fp));
	    decomDS.r2_count += 1;
	}
*/
	
	std::vector<DecomBwdRead>().swap(decomDS.bwd_reads); //Free memory
	
    std::size_t num_write;
    num_write = fwrite(decomDS.bwd_print.data(), sizeof(char), decomDS.bwd_print.size(), *(decomDS.o2_fp));
    if (num_write) {
	    //std::cout << "Wrote " << num_write << " bytes to rd2 file" << std::endl;
#if !NDEBUG
	    std::cout << "To write " << decomDS.bwd_print.size() << " bytes to rd2 file" << std::endl;
        std::cout << "Wrote    " << num_write << " bytes to rd2 file" << std::endl;
#endif
    } else {
	    std::cerr << "Error writing to rd2 file" << std::endl;
	    assert (0);
    }
	std::vector<char>().swap(decomDS.bwd_print); //Free memory
    return 0;
}

void *decompactReads3Thread(void *arg) {
    struct DecompressArgsForThread * tap;
    tap = (struct DecompressArgsForThread *) arg;
    std::uint32_t indices_per_thread = ceil(((double) tap->daft_decom_ds->mapped_count)/((double) tap->daft_in_args->threadCount));
    std::uint32_t index_start = (tap->daft_thread_id) * indices_per_thread;
    std::uint32_t index_end = ((tap->daft_thread_id+1) * indices_per_thread) - 1;
    if ((tap->daft_decom_ds->mapped_count - 1) < index_end) {
        index_end = tap->daft_decom_ds->mapped_count - 1;
    }
    //std::cout << tap->daft_thread_id << " : " << index_start << " : " << index_end << std::endl;
    
	std::uint8_t my_bases[READ_LEN_MAX], my_rc_bases[READ_LEN_MAX];
	std::uint32_t ri, k;
	std::uint8_t z, one_base;
    std::size_t op_idx = ((std::size_t) index_start) * (tap->daft_in_args->rdLength + 3);
    for (std::uint32_t i = index_start; i <= index_end; i++) {
        DecomBwdRead & u1 = tap->daft_decom_ds->bwd_reads[i];
        for (ri = 0, k = 0; ri < tap->daft_in_args->rdLength; ri++) {
            one_base = 0;
            z = 1;
            if (TestBit(u1.read, k)) {
                one_base |= z;
            }
            k++;
            z = (z << 1);
            if (TestBit(u1.read, k)) {
                one_base |= z;
            }
            k++;
            z = (z << 1);
            if (TestBit(u1.read, k)) {
                one_base |= z;
            }
            k++;
            //z = (z << 1); //Not necessary
            my_bases[ri] = one_base;
        }
	    for (std::uint32_t j = 0; j < tap->daft_in_args->rdLength; j++) {
			if (my_bases[j] != 4) {
				my_rc_bases[tap->daft_in_args->rdLength - 1 - j] = 3 - my_bases[j];
			} else {
				my_rc_bases[tap->daft_in_args->rdLength - 1 - j] = 4;
			}
	    }
        tap->daft_decom_ds->bwd_print.data()[op_idx] = '>'; op_idx += 1;
        tap->daft_decom_ds->bwd_print.data()[op_idx] = '\n'; op_idx += 1;
	    for (std::uint32_t j = 0; j < tap->daft_in_args->rdLength; j++) {
	        tap->daft_decom_ds->bwd_print.data()[op_idx] = Uint8Tochar(my_rc_bases[j]);
	        op_idx += 1;
	    }
        tap->daft_decom_ds->bwd_print.data()[op_idx] = '\n'; op_idx += 1;
    }
    
    return NULL;
}
//write to rd1 file (fastq 1 file)
void write_fwd_data(DecompressArgsForThread * wfd_tap) {
	std::size_t num_write;
	FILE * wfd_fp = *(wfd_tap[0].daft_decom_ds->o1_fp);
    for (std::uint32_t i = 0; i < wfd_tap[0].daft_in_args->threadCount; i++) {
        /*
        std::size_t j;
        for (j = 0; j < wfd_tap[i].daft_fwd_reads.size(); j++) {
            if (j % wfd_tap[i].daft_in_args->rdLength == 0) {
                fprintf(wfd_fp, ">%lu\n", (wfd_tap[i].daft_decom_ds->r1_count + 1));
            }
            fputc(wfd_tap[i].daft_fwd_reads[j], wfd_fp);
            if (j % wfd_tap[i].daft_in_args->rdLength == (wfd_tap[i].daft_in_args->rdLength - 1)) {
                fputc('\n', wfd_fp);
                wfd_tap[i].daft_decom_ds->r1_count += 1;
            }
        }
        assert (j == wfd_tap[i].daft_fwd_reads.size());
        */
	    num_write = fwrite(wfd_tap[i].daft_fwd_reads.data(), sizeof(char), wfd_tap[i].daft_fwd_reads.size(), wfd_fp);
	    if (num_write) {
	        //std::cout << "Wrote " << num_write << " bytes to rd1 file" << std::endl;
#if !NDEBUG
	        std::cout << "To write " << wfd_tap[i].daft_fwd_reads.size() << " bytes to rd1 file" << std::endl;
            std::cout << "Wrote    " << num_write << " bytes to rd1 file" << std::endl;
#endif
	    } else {
	        std::cerr << "Error writing to rd1 file" << std::endl;
	        assert (0);
	    }
        std::vector<char>().swap(wfd_tap[i].daft_fwd_reads); //Free memory
    }
}
//read single end reads
void read_se_data(DecompressArgsForThread * rsd_tap) {
    std::size_t num_read, my_size;
	FILE * rsd_fp = *(rsd_tap[0].daft_decom_ds->ip_fp);
    for (std::uint32_t i = 0; i < rsd_tap[0].daft_in_args->threadCount; i++) {
	    num_read = fread(&(my_size), 1, sizeof(std::size_t), rsd_fp);
	    if (num_read) {
	        //std::cout << "Bytes read : " << num_read << std::endl;
	    } else {
	        std::cerr << "Error reading from .cpct file" << std::endl;
	        assert (0);
	    }
	    rsd_tap[i].daft_decom_ds->totalBytes += num_read;
	    rsd_tap[i].daft_fwd_locns.resize(my_size);
	    num_read = fread(rsd_tap[i].daft_fwd_locns.data(), 1, rsd_tap[i].daft_fwd_locns.size(), rsd_fp);
	    if (num_read) {
	        //std::cout << "Read " << num_read << " bytes from .cpct file" << std::endl;
	    } else {
	        std::cerr << "Error reading from .cpct file" << std::endl;
	        assert (0);
	    }
	    rsd_tap[i].daft_decom_ds->totalBytes += num_read;
	    (rsd_tap->dsp->ds_fwd_locns_count)[i] += rsd_tap[i].daft_fwd_locns.size();
    }
    for (std::uint32_t i = 0; i < rsd_tap[0].daft_in_args->threadCount; i++) {
	    num_read = fread(&(my_size), 1, sizeof(std::size_t), rsd_fp);
	    if (num_read) {
	        //std::cout << "Bytes read : " << num_read << std::endl;
	    } else {
	        std::cerr << "Error reading from .cpct file" << std::endl;
	        assert (0);
	    }
	    rsd_tap[i].daft_decom_ds->totalBytes += num_read;
	    rsd_tap[i].daft_bwd_locns.resize(my_size);
	    num_read = fread(rsd_tap[i].daft_bwd_locns.data(), 1, rsd_tap[i].daft_bwd_locns.size(), rsd_fp);
	    if (num_read) {
	        //std::cout << "Read " << num_read << " bytes from .cpct file" << std::endl;
	    } else {
	        std::cerr << "Error reading from .cpct file" << std::endl;
	        assert (0);
	    }
	    rsd_tap[i].daft_decom_ds->totalBytes += num_read;
	    (rsd_tap->dsp->ds_bwd_locns_count)[i] += rsd_tap[i].daft_bwd_locns.size();
    }
    for (std::uint32_t i = 0; i < rsd_tap[0].daft_in_args->threadCount; i++) {
	    num_read = fread(&(my_size), 1, sizeof(std::size_t), rsd_fp);
	    if (num_read) {
	        //std::cout << "Bytes read : " << num_read << std::endl;
	    } else {
	        std::cerr << "Error reading from .cpct file" << std::endl;
	        assert (0);
	    }
	    rsd_tap[i].daft_decom_ds->totalBytes += num_read;
	    rsd_tap[i].daft_fwd_diff_counts.resize(my_size);
	    num_read = fread(rsd_tap[i].daft_fwd_diff_counts.data(), 1, rsd_tap[i].daft_fwd_diff_counts.size(), rsd_fp);
	    if (num_read) {
	        //std::cout << "Read " << num_read << " bytes from .cpct file" << std::endl;
	    } else {
	        std::cerr << "Error reading from .cpct file" << std::endl;
	        assert (0);
	    }
	    rsd_tap[i].daft_decom_ds->totalBytes += num_read;
	    (rsd_tap->dsp->ds_fwd_diff_counts_count)[i] += rsd_tap[i].daft_fwd_diff_counts.size();
    }
    for (std::uint32_t i = 0; i < rsd_tap[0].daft_in_args->threadCount; i++) {
	    num_read = fread(&(my_size), 1, sizeof(std::size_t), rsd_fp);
	    if (num_read) {
	        //std::cout << "Bytes read : " << num_read << std::endl;
	    } else {
	        std::cerr << "Error reading from .cpct file" << std::endl;
	        assert (0);
	    }
	    rsd_tap[i].daft_decom_ds->totalBytes += num_read;
	    rsd_tap[i].daft_bwd_diff_counts.resize(my_size);
	    num_read = fread(rsd_tap[i].daft_bwd_diff_counts.data(), 1, rsd_tap[i].daft_bwd_diff_counts.size(), rsd_fp);
	    if (num_read) {
	        //std::cout << "Read " << num_read << " bytes from .cpct file" << std::endl;
	    } else {
	        std::cerr << "Error reading from .cpct file" << std::endl;
	        assert (0);
	    }
	    rsd_tap[i].daft_decom_ds->totalBytes += num_read;
	    (rsd_tap->dsp->ds_bwd_diff_counts_count)[i] += rsd_tap[i].daft_bwd_diff_counts.size();
    }
    for (std::uint32_t i = 0; i < rsd_tap[0].daft_in_args->threadCount; i++) {
	    num_read = fread(&(my_size), 1, sizeof(std::size_t), rsd_fp);
	    if (num_read) {
	        //std::cout << "Bytes read : " << num_read << std::endl;
	    } else {
	        std::cerr << "Error reading from .cpct file" << std::endl;
	        assert (0);
	    }
	    rsd_tap[i].daft_decom_ds->totalBytes += num_read;
	    rsd_tap[i].daft_fwd_diff_posns.resize(my_size);
	    num_read = fread(rsd_tap[i].daft_fwd_diff_posns.data(), 1, rsd_tap[i].daft_fwd_diff_posns.size(), rsd_fp);
	    if (num_read) {
	        //std::cout << "Read " << num_read << " bytes from .cpct file" << std::endl;
	    } else {
	        std::cerr << "Error reading from .cpct file" << std::endl;
	        assert (0);
	    }
	    rsd_tap[i].daft_decom_ds->totalBytes += num_read;
	    (rsd_tap->dsp->ds_fwd_diff_posns_count)[i] += rsd_tap[i].daft_fwd_diff_posns.size();
    }
    for (std::uint32_t i = 0; i < rsd_tap[0].daft_in_args->threadCount; i++) {
	    num_read = fread(&(my_size), 1, sizeof(std::size_t), rsd_fp);
	    if (num_read) {
	        //std::cout << "Bytes read : " << num_read << std::endl;
	    } else {
	        std::cerr << "Error reading from .cpct file" << std::endl;
	        assert (0);
	    }
	    rsd_tap[i].daft_decom_ds->totalBytes += num_read;
	    rsd_tap[i].daft_bwd_diff_posns.resize(my_size);
	    num_read = fread(rsd_tap[i].daft_bwd_diff_posns.data(), 1, rsd_tap[i].daft_bwd_diff_posns.size(), rsd_fp);
	    if (num_read) {
	        //std::cout << "Read " << num_read << " bytes from .cpct file" << std::endl;
	    } else {
	        std::cerr << "Error reading from .cpct file" << std::endl;
	        assert (0);
	    }
	    rsd_tap[i].daft_decom_ds->totalBytes += num_read;
	    (rsd_tap->dsp->ds_bwd_diff_posns_count)[i] += rsd_tap[i].daft_bwd_diff_posns.size();
    }
    for (std::uint32_t i = 0; i < rsd_tap[0].daft_in_args->threadCount; i++) {
	    num_read = fread(&(my_size), 1, sizeof(std::size_t), rsd_fp);
	    if (num_read) {
	        //std::cout << "Bytes read : " << num_read << std::endl;
	    } else {
	        std::cerr << "Error reading from .cpct file" << std::endl;
	        assert (0);
	    }
	    rsd_tap[i].daft_decom_ds->totalBytes += num_read;
	    rsd_tap[i].daft_fwd_diff_values.resize(my_size);
	    num_read = fread(rsd_tap[i].daft_fwd_diff_values.data(), 1, rsd_tap[i].daft_fwd_diff_values.size(), rsd_fp);
	    if (num_read) {
	        //std::cout << "Read " << num_read << " bytes from .cpct file" << std::endl;
	    } else {
	        std::cerr << "Error reading from .cpct file" << std::endl;
	        assert (0);
	    }
	    rsd_tap[i].daft_decom_ds->totalBytes += num_read;
	    (rsd_tap->dsp->ds_fwd_diff_values_count)[i] += rsd_tap[i].daft_fwd_diff_values.size();
    }
    for (std::uint32_t i = 0; i < rsd_tap[0].daft_in_args->threadCount; i++) {
	    num_read = fread(&(my_size), 1, sizeof(std::size_t), rsd_fp);
	    if (num_read) {
	        //std::cout << "Bytes read : " << num_read << std::endl;
	    } else {
	        std::cerr << "Error reading from .cpct file" << std::endl;
	        assert (0);
	    }
	    rsd_tap[i].daft_decom_ds->totalBytes += num_read;
	    rsd_tap[i].daft_bwd_diff_values.resize(my_size);
	    num_read = fread(rsd_tap[i].daft_bwd_diff_values.data(), 1, rsd_tap[i].daft_bwd_diff_values.size(), rsd_fp);
	    if (num_read) {
	        //std::cout << "Read " << num_read << " bytes from .cpct file" << std::endl;
	    } else {
	        std::cerr << "Error reading from .cpct file" << std::endl;
	        assert (0);
	    }
	    rsd_tap[i].daft_decom_ds->totalBytes += num_read;
	    (rsd_tap->dsp->ds_bwd_diff_values_count)[i] += rsd_tap[i].daft_bwd_diff_values.size();
    }
}
//read paired end data (using relative locations of other ends)
void read_pe_data(DecompressArgsForThread * rpd_tap) {
    std::size_t num_read, my_size;
	FILE * rpd_fp = *(rpd_tap[0].daft_decom_ds->ip_fp);
    for (std::uint32_t i = 0; i < rpd_tap[0].daft_in_args->threadCount; i++) {
	    num_read = fread(&(my_size), 1, sizeof(std::size_t), rpd_fp);
	    if (num_read) {
	        //std::cout << "Bytes read : " << num_read << std::endl;
	    } else {
	        std::cerr << "Error reading from .cpct file" << std::endl;
	        assert (0);
	    }
	    rpd_tap[i].daft_decom_ds->totalBytes += num_read;
	    rpd_tap[i].daft_pe_rel_locns.resize(my_size);
	    num_read = fread(rpd_tap[i].daft_pe_rel_locns.data(), 1, rpd_tap[i].daft_pe_rel_locns.size(), rpd_fp);
	    if (num_read) {
	        //std::cout << "Read " << num_read << " bytes from .cpct file" << std::endl;
	    } else {
	        std::cerr << "Error reading from .cpct file" << std::endl;
	        assert (0);
	    }
	    rpd_tap[i].daft_decom_ds->totalBytes += num_read;
	    (rpd_tap->dsp->ds_pe_rel_locns_count)[i] += rpd_tap[i].daft_pe_rel_locns.size();
    }
    for (std::uint32_t i = 0; i < rpd_tap[0].daft_in_args->threadCount; i++) {
	    num_read = fread(&(my_size), 1, sizeof(std::size_t), rpd_fp);
	    if (num_read) {
	        //std::cout << "Bytes read : " << num_read << std::endl;
	    } else {
	        std::cerr << "Error reading from .cpct file" << std::endl;
	        assert (0);
	    }
	    rpd_tap[i].daft_decom_ds->totalBytes += num_read;
	    rpd_tap[i].daft_pe_rel_posns.resize(my_size);
	    num_read = fread(rpd_tap[i].daft_pe_rel_posns.data(), 1, rpd_tap[i].daft_pe_rel_posns.size(), rpd_fp);
	    if (num_read) {
	        //std::cout << "Read " << num_read << " bytes from .cpct file" << std::endl;
	    } else {
	        std::cerr << "Error reading from .cpct file" << std::endl;
	        assert (0);
	    }
	    rpd_tap[i].daft_decom_ds->totalBytes += num_read;
	    (rpd_tap->dsp->ds_pe_rel_posns_count)[i] += rpd_tap[i].daft_pe_rel_posns.size();
    }
}
//decode list of starting locations
void *decompactReads1Thread(void *arg) {
    struct DecompressArgsForThread* tap;
    tap = (struct DecompressArgsForThread*)arg;
    std::uint32_t indices_per_thread = ceil(((double)tap->daft_decom_ds->mapped_count) / ((double)tap->daft_in_args->threadCount));
    std::uint32_t index_start = (tap->daft_thread_id) * indices_per_thread;
    std::uint32_t index_end = ((tap->daft_thread_id + 1) * indices_per_thread) - 1;
    if ((tap->daft_decom_ds->mapped_count - 1) < index_end) {
        index_end = tap->daft_decom_ds->mapped_count - 1;
    }
    //std::cout << tap->daft_thread_id << " : " << index_start << " : " << index_end << std::endl;

    // tap->daft_fwd_reads.resize(((std::size_t)(index_end - index_start + 1)) * (tap->daft_in_args->rdLength + 3));
    

    tap->daft_fwd_reads.resize(((std::size_t)(index_end - index_start + 1)) * (tap->daft_in_args->rdLength + 3));
    std::uint32_t loc_idx = 0, prev_locn = 0, curr_locn = 0;
    std::uint32_t dct_idx = 0, dpv_idx = 0 , mod_idx = 0;
    std::size_t fwd_idx = 0;
    std::uint32_t local_diffs_fwd[READ_LEN_MAX * 3][EDIT_DISTANCE][10];
    memset(local_diffs_fwd, 0, READ_LEN_MAX * 3 * EDIT_DISTANCE * 10 * sizeof(std::uint32_t));
    std::uint8_t  local_bases_fwd[READ_LEN_MAX * 3][EDIT_DISTANCE];
    memset(local_bases_fwd, 254, READ_LEN_MAX * 3 * EDIT_DISTANCE * sizeof(std::uint8_t));
    std::uint32_t local_count_fwd[READ_LEN_MAX * 3][EDIT_DISTANCE];
    memset(local_count_fwd, 0, READ_LEN_MAX * 3 * EDIT_DISTANCE * sizeof(std::uint32_t));
    std::uint32_t local_usage_fwd[READ_LEN_MAX * 3];
    memset(local_usage_fwd, 0, READ_LEN_MAX * 3 * sizeof(std::uint32_t));
    int temp_diff_cnt_fwd = 0;
    int prev_diff_loc_fwd = -1;
    std::uint32_t temp_ref_buf_pos_fwd = 0;

    std::uint32_t aftr_diff_locs[32];
    std::uint32_t aftr_diff_vals[32];
    std::uint32_t diff_locs[32];
    std::uint32_t diff_vals[32];

    //char my_bases_1[READ_LEN_MAX + EDIT_DISTANCE + 1]; //To account for indels
    //char my_bases_2[READ_LEN_MAX];
    bool is_mod_unused;


    if (tap->daft_in_args->updateReference) {
        for (std::uint32_t i = index_start; i <= index_end; i++) {
            if (tap->daft_fwd_locns[loc_idx] < 253) {
                curr_locn = (std::uint32_t)tap->daft_fwd_locns[loc_idx];
                loc_idx += 1;
            }
            else if (tap->daft_fwd_locns[loc_idx] == 253) {
                loc_idx += 1;
                curr_locn = (std::uint32_t)tap->daft_fwd_locns[loc_idx];
                loc_idx += 1;
                curr_locn |= (((std::uint32_t)tap->daft_fwd_locns[loc_idx]) << 8);
                loc_idx += 1;
            }
            else if (tap->daft_fwd_locns[loc_idx] == 254) {
                loc_idx += 1;
                curr_locn = (std::uint32_t)tap->daft_fwd_locns[loc_idx];
                loc_idx += 1;
                curr_locn |= (((std::uint32_t)tap->daft_fwd_locns[loc_idx]) << 8);
                loc_idx += 1;
                curr_locn |= (((std::uint32_t)tap->daft_fwd_locns[loc_idx]) << 16);
                loc_idx += 1;
            }
            else {
                loc_idx += 1;
                curr_locn = (std::uint32_t)tap->daft_fwd_locns[loc_idx];
                loc_idx += 1;
                curr_locn |= (((std::uint32_t)tap->daft_fwd_locns[loc_idx]) << 8);
                loc_idx += 1;
                curr_locn |= (((std::uint32_t)tap->daft_fwd_locns[loc_idx]) << 16);
                loc_idx += 1;
                curr_locn |= (((std::uint32_t)tap->daft_fwd_locns[loc_idx]) << 24);
                loc_idx += 1;
            }
            curr_locn += prev_locn;
#if !NDEBUG
            //if (curr_locn > (tap->daft_decom_ds->ref_length - tap->daft_in_args->rdLength)) {
            //    curr_locn = tap->daft_decom_ds->ref_length - tap->daft_in_args->rdLength;
            //}
            //if (curr_locn > (tap->daft_decom_ds->ref_length - tap->daft_in_args->rdLength + tap->daft_fwd_diff_counts[dct_idx])) {
            //    std::cerr << "tap->daft_thread_id : " << tap->daft_thread_id << std::endl;
            //    std::cerr << "prev_locn : " << prev_locn << std::endl;
            //    std::cerr << "curr_locn : " << curr_locn << std::endl;
            //    std::cerr << "i : " << i << std::endl;
            //    std::cerr << "diff_count : " << ((std::uint32_t) tap->daft_fwd_diff_counts[dct_idx]) << std::endl;
            //}
            assert(curr_locn <= (tap->daft_decom_ds->ref_length - tap->daft_in_args->rdLength));
#endif

            tap->daft_fwd_reads.data()[fwd_idx] = '>';
            fwd_idx += 1;
            tap->daft_fwd_reads.data()[fwd_idx] = '\n';
            fwd_idx += 1;
            if (tap->daft_fwd_diff_counts[dct_idx]) {
                //std::memcpy(my_bases_1, tap->daft_decom_ds->ref_bases + curr_locn, tap->daft_in_args->rdLength + EDIT_DISTANCE);
                //std::uint32_t mb1_idx = 0, mb2_idx = 0;
                //for (std::uint32_t j = 0; j < tap->daft_in_args->rdLength; j++) {
                //   my_bases_2[mb2_idx] = my_bases_1[mb1_idx];
                //   mb1_idx += 1; mb2_idx += 1;
                //}
                //dpv_idx += tap->daft_fwd_diff_counts[dct_idx];
                //std::memcpy(tap->daft_fwd_reads.data() + fwd_idx, my_bases_2, tap->daft_in_args->rdLength);
                //fwd_idx += tap->daft_in_args->rdLength;
                std::uint32_t ref_idx = 0, temp_mod_idx = 0;
                //mod_idx += tap->daft_fwd_diff_posns[dpv_idx];
                is_mod_unused = true;

                for (std::uint32_t k = 0; k < tap->daft_fwd_diff_counts[dct_idx]; k++) { //to set initial values in aftr_diff_locs and aftr_diff_vals
                    aftr_diff_vals[k] = tap->daft_bwd_diff_values[dpv_idx];
                    diff_vals[k] = tap->daft_bwd_diff_values[dpv_idx];
                    temp_mod_idx += tap->daft_fwd_diff_posns[dpv_idx];
                    aftr_diff_locs[k + 1] = temp_mod_idx;
                    diff_locs[k + 1] = temp_mod_idx;
                    dpv_idx += 1;
                }

                //tap->daft_fwd_diff_counts[dct_idx] -= 1;
                aftr_diff_locs[0] = 0;

                for (std::uint32_t i = prev_locn; i < curr_locn; i++) {
                    temp_diff_cnt_fwd = 0;
                    prev_diff_loc_fwd = -1;
                    if (prev_locn) {
                        for (std::uint32_t i = prev_locn; i < curr_locn; i++) {
                            temp_ref_buf_pos_fwd = i % (READ_LEN_MAX * 3);
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
                    for (int i = 1; i < tap->daft_fwd_diff_counts[dct_idx]; i++)
                    {
                        temp_ref_buf_pos_fwd = curr_locn + ref_idx;
                        if (diff_locs[i] == prev_diff_loc_fwd) {
                            temp_diff_cnt_fwd++;
                        }
                        else {
                            prev_diff_loc_fwd = diff_locs[i];
                            temp_diff_cnt_fwd = 1;
                        }
                        local_usage_fwd[temp_ref_buf_pos_fwd] = (std::uint32_t)temp_diff_cnt_fwd;

                        if (local_bases_fwd[temp_ref_buf_pos_fwd][temp_diff_cnt_fwd - 1] == ((std::uint8_t)diff_vals[i])) {
                            local_diffs_fwd[temp_ref_buf_pos_fwd][temp_diff_cnt_fwd - 1][diff_vals[i]] += 1;
                            local_count_fwd[temp_ref_buf_pos_fwd][temp_diff_cnt_fwd - 1] += 1;
                        }
                        else {
                            aftr_diff_locs[0] += 1;
                            aftr_diff_locs[aftr_diff_locs[0]] = diff_locs[i];
                            aftr_diff_vals[aftr_diff_locs[0]] = diff_vals[i];
                            local_diffs_fwd[temp_ref_buf_pos_fwd][temp_diff_cnt_fwd - 1][diff_vals[i]] += 1;
                            if (local_diffs_fwd[temp_ref_buf_pos_fwd][temp_diff_cnt_fwd - 1][diff_vals[i]] > (local_count_fwd[temp_ref_buf_pos_fwd][temp_diff_cnt_fwd - 1])) {
                                local_count_fwd[temp_ref_buf_pos_fwd][temp_diff_cnt_fwd - 1] = local_diffs_fwd[temp_ref_buf_pos_fwd][temp_diff_cnt_fwd - 1][diff_vals[i]];
                                local_bases_fwd[temp_ref_buf_pos_fwd][temp_diff_cnt_fwd - 1] = ((std::uint8_t)diff_vals[i]);
                            }
                        }
                    }



                }
                aftr_diff_locs[0] -= 1;
                
                for (std::uint32_t j = 0; j < tap->daft_in_args->rdLength; ) {
                    if (is_mod_unused && (ref_idx == aftr_diff_locs[mod_idx + 1])) {
                        if (aftr_diff_vals[mod_idx] < 4) {
                            //Substitution
                            if (aftr_diff_vals[mod_idx] < charToUint8(tap->daft_decom_ds->ref_bases[curr_locn + ref_idx])) {
                                tap->daft_fwd_reads.data()[fwd_idx] = Uint8Tochar(aftr_diff_vals[aftr_diff_locs[0]]);
                            }
                            else {
                                tap->daft_fwd_reads.data()[fwd_idx] = Uint8Tochar(aftr_diff_vals[aftr_diff_locs[0]] + 1);
                            }
                            fwd_idx += 1; j += 1;
                            ref_idx += 1;
                        }
                        else if (aftr_diff_vals[mod_idx] < 9) {
                            //Insertion
                            tap->daft_fwd_reads.data()[fwd_idx] = Uint8Tochar(aftr_diff_vals[aftr_diff_locs[0]] - 4);
                            fwd_idx += 1; j += 1;
                        }
                        else {
                            //Deletion
                            ref_idx += 1;
                        }
                        mod_idx += 1;
                        is_mod_unused = false;
                        if (aftr_diff_locs[0]) {
                            //mod_idx += tap->daft_fwd_diff_posns[dpv_idx];
                            is_mod_unused = true;
                            aftr_diff_locs[0] -= 1;
                        }
                    }
                    else {
                        tap->daft_fwd_reads.data()[fwd_idx] = tap->daft_decom_ds->ref_bases[curr_locn + ref_idx];
                        fwd_idx += 1; j += 1;
                        ref_idx += 1;
                    }
                }
            }
            else {
                //No differences
                std::memcpy(tap->daft_fwd_reads.data() + fwd_idx, tap->daft_decom_ds->ref_bases + curr_locn, tap->daft_in_args->rdLength);
                fwd_idx += tap->daft_in_args->rdLength;
            }
            tap->daft_fwd_reads.data()[fwd_idx] = '\n';
            fwd_idx += 1;


            dct_idx += 1;

            prev_locn = curr_locn;
        }
        assert(loc_idx == tap->daft_fwd_locns.size());
        std::vector<std::uint8_t>().swap(tap->daft_fwd_locns); //Free memory
        assert(dct_idx == tap->daft_fwd_diff_counts.size());
        std::vector<std::uint8_t>().swap(tap->daft_fwd_diff_counts); //Free memory
        assert(dpv_idx == tap->daft_fwd_diff_posns.size());
        std::vector<std::uint8_t>().swap(tap->daft_fwd_diff_posns); //Free memory
        std::vector<std::uint8_t>().swap(tap->daft_fwd_diff_values); //Free memory
        assert(fwd_idx == tap->daft_fwd_reads.size());

        return NULL;

        }
    else {
        for (std::uint32_t i = index_start; i <= index_end; i++) {
            if (tap->daft_fwd_locns[loc_idx] < 253) {
                curr_locn = (std::uint32_t)tap->daft_fwd_locns[loc_idx];
                loc_idx += 1;
            }
            else if (tap->daft_fwd_locns[loc_idx] == 253) {
                loc_idx += 1;
                curr_locn = (std::uint32_t)tap->daft_fwd_locns[loc_idx];
                loc_idx += 1;
                curr_locn |= (((std::uint32_t)tap->daft_fwd_locns[loc_idx]) << 8);
                loc_idx += 1;
            }
            else if (tap->daft_fwd_locns[loc_idx] == 254) {
                loc_idx += 1;
                curr_locn = (std::uint32_t)tap->daft_fwd_locns[loc_idx];
                loc_idx += 1;
                curr_locn |= (((std::uint32_t)tap->daft_fwd_locns[loc_idx]) << 8);
                loc_idx += 1;
                curr_locn |= (((std::uint32_t)tap->daft_fwd_locns[loc_idx]) << 16);
                loc_idx += 1;
            }
            else {
                loc_idx += 1;
                curr_locn = (std::uint32_t)tap->daft_fwd_locns[loc_idx];
                loc_idx += 1;
                curr_locn |= (((std::uint32_t)tap->daft_fwd_locns[loc_idx]) << 8);
                loc_idx += 1;
                curr_locn |= (((std::uint32_t)tap->daft_fwd_locns[loc_idx]) << 16);
                loc_idx += 1;
                curr_locn |= (((std::uint32_t)tap->daft_fwd_locns[loc_idx]) << 24);
                loc_idx += 1;
            }
            curr_locn += prev_locn;
#if !NDEBUG
            //if (curr_locn > (tap->daft_decom_ds->ref_length - tap->daft_in_args->rdLength)) {
            //    curr_locn = tap->daft_decom_ds->ref_length - tap->daft_in_args->rdLength;
            //}
            //if (curr_locn > (tap->daft_decom_ds->ref_length - tap->daft_in_args->rdLength + tap->daft_fwd_diff_counts[dct_idx])) {
            //    std::cerr << "tap->daft_thread_id : " << tap->daft_thread_id << std::endl;
            //    std::cerr << "prev_locn : " << prev_locn << std::endl;
            //    std::cerr << "curr_locn : " << curr_locn << std::endl;
            //    std::cerr << "i : " << i << std::endl;
            //    std::cerr << "diff_count : " << ((std::uint32_t) tap->daft_fwd_diff_counts[dct_idx]) << std::endl;
            //}
            assert(curr_locn <= (tap->daft_decom_ds->ref_length - tap->daft_in_args->rdLength));
#endif

            tap->daft_fwd_reads.data()[fwd_idx] = '>';
            fwd_idx += 1;
            tap->daft_fwd_reads.data()[fwd_idx] = '\n';
            fwd_idx += 1;
            if (tap->daft_fwd_diff_counts[dct_idx]) {
                //std::memcpy(my_bases_1, tap->daft_decom_ds->ref_bases + curr_locn, tap->daft_in_args->rdLength + EDIT_DISTANCE);
                //std::uint32_t mb1_idx = 0, mb2_idx = 0;
                //for (std::uint32_t j = 0; j < tap->daft_in_args->rdLength; j++) {
                //   my_bases_2[mb2_idx] = my_bases_1[mb1_idx];
                //   mb1_idx += 1; mb2_idx += 1;
                //}
                //dpv_idx += tap->daft_fwd_diff_counts[dct_idx];
                //std::memcpy(tap->daft_fwd_reads.data() + fwd_idx, my_bases_2, tap->daft_in_args->rdLength);
                //fwd_idx += tap->daft_in_args->rdLength;
                std::uint32_t ref_idx = 0, mod_idx = 0;
                mod_idx += tap->daft_fwd_diff_posns[dpv_idx];
                is_mod_unused = true;
                tap->daft_fwd_diff_counts[dct_idx] -= 1;
                for (std::uint32_t j = 0; j < tap->daft_in_args->rdLength; ) {
                    if (is_mod_unused && (ref_idx == mod_idx)) {
                        if (tap->daft_fwd_diff_values[dpv_idx] < 4) {
                            //Substitution
                            if (tap->daft_fwd_diff_values[dpv_idx] < charToUint8(tap->daft_decom_ds->ref_bases[curr_locn + ref_idx])) {
                                tap->daft_fwd_reads.data()[fwd_idx] = Uint8Tochar(tap->daft_fwd_diff_values[dpv_idx]);
                            }
                            else {
                                tap->daft_fwd_reads.data()[fwd_idx] = Uint8Tochar(tap->daft_fwd_diff_values[dpv_idx] + 1);
                            }
                            fwd_idx += 1; j += 1;
                            ref_idx += 1;
                        }
                        else if (tap->daft_fwd_diff_values[dpv_idx] < 9) {
                            //Insertion
                            tap->daft_fwd_reads.data()[fwd_idx] = Uint8Tochar(tap->daft_fwd_diff_values[dpv_idx] - 4);
                            fwd_idx += 1; j += 1;
                        }
                        else {
                            //Deletion
                            ref_idx += 1;
                        }
                        dpv_idx += 1;
                        is_mod_unused = false;
                        if (tap->daft_fwd_diff_counts[dct_idx]) {
                            mod_idx += tap->daft_fwd_diff_posns[dpv_idx];
                            is_mod_unused = true;
                            tap->daft_fwd_diff_counts[dct_idx] -= 1;
                        }
                    }
                    else {
                        tap->daft_fwd_reads.data()[fwd_idx] = tap->daft_decom_ds->ref_bases[curr_locn + ref_idx];
                        fwd_idx += 1; j += 1;
                        ref_idx += 1;
                    }
                }
            }
            else {
                //No differences
                std::memcpy(tap->daft_fwd_reads.data() + fwd_idx, tap->daft_decom_ds->ref_bases + curr_locn, tap->daft_in_args->rdLength);
                fwd_idx += tap->daft_in_args->rdLength;
            }
            tap->daft_fwd_reads.data()[fwd_idx] = '\n';
            fwd_idx += 1;
#if !NDEBUG
            assert(tap->daft_fwd_diff_counts[dct_idx] == 0);
#endif
            dct_idx += 1;

            prev_locn = curr_locn;
        }
        assert(loc_idx == tap->daft_fwd_locns.size());
        std::vector<std::uint8_t>().swap(tap->daft_fwd_locns); //Free memory
        assert(dct_idx == tap->daft_fwd_diff_counts.size());
        std::vector<std::uint8_t>().swap(tap->daft_fwd_diff_counts); //Free memory
        assert(dpv_idx == tap->daft_fwd_diff_posns.size());
        std::vector<std::uint8_t>().swap(tap->daft_fwd_diff_posns); //Free memory
        std::vector<std::uint8_t>().swap(tap->daft_fwd_diff_values); //Free memory
        assert(fwd_idx == tap->daft_fwd_reads.size());

        return NULL;
       }
               
    }
 

/* Populate list of forward aligned reads and reverse 
aligned reads
Identify and update differing bases
Compute reverse complement */
void *decompactReads2Thread(void *arg) {
    struct DecompressArgsForThread* tap;
    tap = (struct DecompressArgsForThread*)arg;
    std::uint32_t indices_per_thread = ceil(((double)tap->daft_decom_ds->mapped_count) / ((double)tap->daft_in_args->threadCount));
    std::uint32_t index_start = (tap->daft_thread_id) * indices_per_thread;
    std::uint32_t index_end = ((tap->daft_thread_id + 1) * indices_per_thread) - 1;
    if ((tap->daft_decom_ds->mapped_count - 1) < index_end) {
        index_end = tap->daft_decom_ds->mapped_count - 1;
    }
    //std::cout << tap->daft_thread_id << " : " << index_start << " : " << index_end << std::endl;

    assert(tap->daft_decom_ds->bwd_reads.size() == tap->daft_decom_ds->mapped_count);
    std::uint32_t loc_idx = 0, prev_locn = 0, curr_locn = 0;
    std::uint32_t dct_idx = 0, dpv_idx = 0, mod_idx = 0;
    std::uint32_t pel_idx = 0, pep_idx = 0;
    std::uint32_t pe_posn; int64_t pe_locn_1; std::uint64_t pe_locn_2;
    //char my_bases_1[READ_LEN_MAX + EDIT_DISTANCE + 1]; //To account for indels

    tap->daft_fwd_reads.resize(((std::size_t)(index_end - index_start + 1)) * (tap->daft_in_args->rdLength + 3));
    
    
    std::uint32_t local_diffs_bwd[READ_LEN_MAX * 3][EDIT_DISTANCE][10];
    memset(local_diffs_bwd, 0, READ_LEN_MAX * 3 * EDIT_DISTANCE * 10 * sizeof(std::uint32_t));
    std::uint8_t  local_bases_bwd[READ_LEN_MAX * 3][EDIT_DISTANCE];
    memset(local_bases_bwd, 254, READ_LEN_MAX * 3 * EDIT_DISTANCE * sizeof(std::uint8_t));
    std::uint32_t local_count_bwd[READ_LEN_MAX * 3][EDIT_DISTANCE];
    memset(local_count_bwd, 0, READ_LEN_MAX * 3 * EDIT_DISTANCE * sizeof(std::uint32_t));
    std::uint32_t local_usage_bwd[READ_LEN_MAX * 3];
    memset(local_usage_bwd, 0, READ_LEN_MAX * 3 * sizeof(std::uint32_t));
    int temp_diff_cnt_bwd = 0;
    int prev_diff_loc_bwd = -1;
    std::uint32_t temp_ref_buf_pos_bwd = 0;

    std::uint32_t aftr_diff_locs[32];
    std::uint32_t aftr_diff_vals[32];
    std::uint32_t diff_locs[32];
    std::uint32_t diff_vals[32];

    char my_bases_2[READ_LEN_MAX];
    bool is_mod_unused;
    
    if (tap->daft_in_args->updateReference)
    {
        for (std::uint32_t i = index_start; i <= index_end; i++)
        {
            if (tap->daft_bwd_locns[loc_idx] < 253) {
                curr_locn = (std::uint32_t)tap->daft_bwd_locns[loc_idx];
                loc_idx += 1;
            }
            else if (tap->daft_bwd_locns[loc_idx] == 253) {
                loc_idx += 1;
                curr_locn = (std::uint32_t)tap->daft_bwd_locns[loc_idx];
                loc_idx += 1;
                curr_locn |= (((std::uint32_t)tap->daft_bwd_locns[loc_idx]) << 8);
                loc_idx += 1;
            }
            else if (tap->daft_bwd_locns[loc_idx] == 254) {
                loc_idx += 1;
                curr_locn = (std::uint32_t)tap->daft_bwd_locns[loc_idx];
                loc_idx += 1;
                curr_locn |= (((std::uint32_t)tap->daft_bwd_locns[loc_idx]) << 8);
                loc_idx += 1;
                curr_locn |= (((std::uint32_t)tap->daft_bwd_locns[loc_idx]) << 16);
                loc_idx += 1;
            }
            else {
                loc_idx += 1;
                curr_locn = (std::uint32_t)tap->daft_bwd_locns[loc_idx];
                loc_idx += 1;
                curr_locn |= (((std::uint32_t)tap->daft_bwd_locns[loc_idx]) << 8);
                loc_idx += 1;
                curr_locn |= (((std::uint32_t)tap->daft_bwd_locns[loc_idx]) << 16);
                loc_idx += 1;
                curr_locn |= (((std::uint32_t)tap->daft_bwd_locns[loc_idx]) << 24);
                loc_idx += 1;
            }
            curr_locn += prev_locn;
#if !NDEBUG
            //if (curr_locn > (tap->daft_decom_ds->ref_length - tap->daft_in_args->rdLength)) {
            //    curr_locn = tap->daft_decom_ds->ref_length - tap->daft_in_args->rdLength;
            //}
            //if (curr_locn > (tap->daft_decom_ds->ref_length - tap->daft_in_args->rdLength + tap->daft_bwd_diff_counts[dct_idx])) {
            //    std::cerr << "tap->daft_thread_id : " << tap->daft_thread_id << std::endl;
            //    std::cerr << "prev_locn : " << prev_locn << std::endl;
            //    std::cerr << "curr_locn : " << curr_locn << std::endl;
            //    std::cerr << "i : " << i << std::endl;
            //    std::cerr << "diff_count : " << ((std::uint32_t) tap->daft_bwd_diff_counts[dct_idx]) << std::endl;
            //}
            assert(curr_locn <= (tap->daft_decom_ds->ref_length - tap->daft_in_args->rdLength));
#endif

            if (tap->daft_bwd_diff_counts[dct_idx]) {
                //std::memcpy(my_bases_1, tap->daft_decom_ds->ref_bases + curr_locn, tap->daft_in_args->rdLength + EDIT_DISTANCE);
                //std::uint32_t mb1_idx = 0, mb2_idx = 0;
                //for (std::uint32_t j = 0; j < tap->daft_in_args->rdLength; j++) {
                //   my_bases_2[mb2_idx] = my_bases_1[mb1_idx];
                //   mb1_idx += 1; mb2_idx += 1;
                //}
                //dpv_idx += tap->daft_bwd_diff_counts[dct_idx];
                std::uint32_t ref_idx = 0, temp_mod_idx = 0;
                // mod_idx += tap->daft_bwd_diff_posns[dpv_idx];
                is_mod_unused = true;

                for (int k = 0; k < tap->daft_bwd_diff_counts[dct_idx]; k++) {
                    aftr_diff_vals[k] = tap->daft_bwd_diff_values[dpv_idx];
                    diff_vals[k] = tap->daft_bwd_diff_values[dpv_idx];
                    temp_mod_idx += tap->daft_bwd_diff_posns[dpv_idx];
                    aftr_diff_locs[k + 1] = temp_mod_idx;
                    diff_locs[k + 1] = temp_mod_idx;
                    dpv_idx += 1;
                }
                aftr_diff_locs[0] = 0;

                //tap->daft_fwd_diff_counts[dct_idx] -= 1;


                for (std::uint32_t i = prev_locn; i < curr_locn; i++) {
                    temp_diff_cnt_bwd = 0;
                    prev_diff_loc_bwd = -1;
                    if (prev_locn) {
                        for (std::uint32_t i = prev_locn; i < curr_locn; i++) {
                            temp_ref_buf_pos_bwd = i % (READ_LEN_MAX * 3);
                            //                     assert (temp_ref_buf_pos_fwd < (READ_LEN_MAX*3));
                            //                     assert (local_usage_fwd[temp_ref_buf_pos_fwd] < EDIT_DISTANCE);
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
                    for (int i = 1; i < tap->daft_bwd_diff_counts[dct_idx]; i++)
                    {
                        temp_ref_buf_pos_bwd = curr_locn + ref_idx;
                        if (diff_locs[i] == prev_diff_loc_bwd) {
                            temp_diff_cnt_bwd++;
                        }
                        else {
                            prev_diff_loc_bwd = diff_locs[i];
                            temp_diff_cnt_bwd = 1;
                        }
                        local_usage_bwd[temp_ref_buf_pos_bwd] = (std::uint32_t)temp_diff_cnt_bwd;

                        if (local_bases_bwd[temp_ref_buf_pos_bwd][temp_diff_cnt_bwd - 1] == ((std::uint8_t)diff_vals[i])) {
                            local_diffs_bwd[temp_ref_buf_pos_bwd][temp_diff_cnt_bwd - 1][diff_vals[i]] += 1;
                            local_count_bwd[temp_ref_buf_pos_bwd][temp_diff_cnt_bwd - 1] += 1;
                        }
                        else {
                            aftr_diff_locs[0] += 1;
                            aftr_diff_locs[aftr_diff_locs[0]] = diff_locs[i];
                            aftr_diff_vals[aftr_diff_locs[0]] = diff_vals[i];
                            local_diffs_bwd[temp_ref_buf_pos_bwd][temp_diff_cnt_bwd - 1][diff_vals[i]] += 1;
                            if (local_diffs_bwd[temp_ref_buf_pos_bwd][temp_diff_cnt_bwd - 1][diff_vals[i]] > (local_count_bwd[temp_ref_buf_pos_bwd][temp_diff_cnt_bwd - 1])) {
                                local_count_bwd[temp_ref_buf_pos_bwd][temp_diff_cnt_bwd - 1] = local_diffs_bwd[temp_ref_buf_pos_bwd][temp_diff_cnt_bwd - 1][diff_vals[i]];
                                local_bases_bwd[temp_ref_buf_pos_bwd][temp_diff_cnt_bwd - 1] = ((std::uint8_t)diff_vals[i]);
                            }
                        }
                    }
                                        
                    }
                aftr_diff_locs[0] -= 1;

                for (std::uint32_t j = 0; j < tap->daft_in_args->rdLength; ) {
                    if (is_mod_unused && (ref_idx == aftr_diff_locs[mod_idx + 1])) {
                        if (aftr_diff_vals[mod_idx] < 4) {
                            //Substitution
                            if (aftr_diff_vals[mod_idx] < charToUint8(tap->daft_decom_ds->ref_bases[curr_locn + ref_idx])) {
                                my_bases_2[j] = Uint8Tochar(aftr_diff_vals[mod_idx]);
                            }
                            else {
                                my_bases_2[j] = Uint8Tochar(aftr_diff_vals[mod_idx]);
                            }
                            j += 1;
                            ref_idx += 1;
                        }
                        else if (aftr_diff_vals[mod_idx] < 9) {
                            //Insertion
                            my_bases_2[j] = Uint8Tochar(aftr_diff_vals[mod_idx] - 4);
                            j += 1;
                        }
                        else {
                            //Deletion
                            ref_idx += 1;
                        }
                        mod_idx += 1;
                        is_mod_unused = false;
                        if (aftr_diff_locs[0]) {
                            //mod_idx += tap->daft_bwd_diff_posns[dpv_idx];
                            is_mod_unused = true;
                            aftr_diff_locs[0] -= 1;
                        }
                    }
                    else {
                        my_bases_2[j] = tap->daft_decom_ds->ref_bases[curr_locn + ref_idx];
                        j += 1;
                        ref_idx += 1;
                    }
                }
                compact_bwd_read(my_bases_2, tap->daft_in_args->rdLength, tap->daft_decom_ds->bwd_reads[i].read);
            }
            else {
                //No differences
                compact_bwd_read(tap->daft_decom_ds->ref_bases + curr_locn, tap->daft_in_args->rdLength, tap->daft_decom_ds->bwd_reads[i].read);
            }


            dct_idx += 1;

            prev_locn = curr_locn;

            if (tap->daft_pe_rel_posns[pep_idx] < 253) {
                pe_posn = (std::uint32_t)tap->daft_pe_rel_posns[pep_idx];
                pep_idx += 1;
            }
            else if (tap->daft_pe_rel_posns[pep_idx] == 253) {
                pep_idx += 1;
                pe_posn = (std::uint32_t)tap->daft_pe_rel_posns[pep_idx];
                pep_idx += 1;
                pe_posn |= (((std::uint32_t)tap->daft_pe_rel_posns[pep_idx]) << 8);
                pep_idx += 1;
            }
            else if (tap->daft_pe_rel_posns[pep_idx] == 254) {
                pep_idx += 1;
                pe_posn = (std::uint32_t)tap->daft_pe_rel_posns[pep_idx];
                pep_idx += 1;
                pe_posn |= (((std::uint32_t)tap->daft_pe_rel_posns[pep_idx]) << 8);
                pep_idx += 1;
                pe_posn |= (((std::uint32_t)tap->daft_pe_rel_posns[pep_idx]) << 16);
                pep_idx += 1;
            }
            else {
                pep_idx += 1;
                pe_posn = (std::uint32_t)tap->daft_pe_rel_posns[pep_idx];
                pep_idx += 1;
                pe_posn |= (((std::uint32_t)tap->daft_pe_rel_posns[pep_idx]) << 8);
                pep_idx += 1;
                pe_posn |= (((std::uint32_t)tap->daft_pe_rel_posns[pep_idx]) << 16);
                pep_idx += 1;
                pe_posn |= (((std::uint32_t)tap->daft_pe_rel_posns[pep_idx]) << 24);
                pep_idx += 1;
            }
            tap->daft_decom_ds->bwd_reads[i].pe_rel_posn = pe_posn;

            if (tap->daft_pe_rel_locns[pel_idx] < 252) {
                pe_locn_2 = (std::uint64_t)tap->daft_pe_rel_locns[pel_idx];
                pel_idx += 1;
            }
            else if (tap->daft_pe_rel_locns[pel_idx] == 252) {
                pel_idx += 1;
                pe_locn_2 = (std::uint64_t)tap->daft_pe_rel_locns[pel_idx];
                pel_idx += 1;
                pe_locn_2 |= (((std::uint64_t)tap->daft_pe_rel_locns[pel_idx]) << 8);
                pel_idx += 1;
            }
            else if (tap->daft_pe_rel_locns[pel_idx] == 253) {
                pel_idx += 1;
                pe_locn_2 = (std::uint64_t)tap->daft_pe_rel_locns[pel_idx];
                pel_idx += 1;
                pe_locn_2 |= (((std::uint64_t)tap->daft_pe_rel_locns[pel_idx]) << 8);
                pel_idx += 1;
                pe_locn_2 |= (((std::uint64_t)tap->daft_pe_rel_locns[pel_idx]) << 16);
                pel_idx += 1;
            }
            else if (tap->daft_pe_rel_locns[pel_idx] == 254) {
                pel_idx += 1;
                pe_locn_2 = (std::uint64_t)tap->daft_pe_rel_locns[pel_idx];
                pel_idx += 1;
                pe_locn_2 |= (((std::uint64_t)tap->daft_pe_rel_locns[pel_idx]) << 8);
                pel_idx += 1;
                pe_locn_2 |= (((std::uint64_t)tap->daft_pe_rel_locns[pel_idx]) << 16);
                pel_idx += 1;
                pe_locn_2 |= (((std::uint64_t)tap->daft_pe_rel_locns[pel_idx]) << 24);
                pel_idx += 1;
            }
            else {
                pel_idx += 1;
                pe_locn_2 = (std::uint64_t)tap->daft_pe_rel_locns[pel_idx];
                pel_idx += 1;
                pe_locn_2 |= (((std::uint64_t)tap->daft_pe_rel_locns[pel_idx]) << 8);
                pel_idx += 1;
                pe_locn_2 |= (((std::uint64_t)tap->daft_pe_rel_locns[pel_idx]) << 16);
                pel_idx += 1;
                pe_locn_2 |= (((std::uint64_t)tap->daft_pe_rel_locns[pel_idx]) << 24);
                pel_idx += 1;
                pe_locn_2 |= (((std::uint64_t)tap->daft_pe_rel_locns[pel_idx]) << 32);
                pel_idx += 1;
            }
            pe_locn_1 = (int64_t)pe_locn_2;
            if (pe_locn_1 == 0) {
                //Do nothing
            }
            else if (pe_locn_1 % 2) {
                pe_locn_1 = (pe_locn_1 + 1) / 2;
            }
            else {
                pe_locn_1 = pe_locn_1 / (-2);
            }
            pe_locn_1 += tap->daft_decom_ds->pe_rel_locn_mean;
            pe_locn_1 = ((int64_t)curr_locn) - pe_locn_1;
#if !NDEBUG
            assert(pe_locn_1 >= 0);
#endif
            tap->daft_decom_ds->bwd_reads[i].pe_location = (std::uint32_t)pe_locn_1;
        }
     assert(loc_idx == tap->daft_bwd_locns.size());
    std::vector<std::uint8_t>().swap(tap->daft_bwd_locns); //Free memory
    assert(dct_idx == tap->daft_bwd_diff_counts.size());
    std::vector<std::uint8_t>().swap(tap->daft_bwd_diff_counts); //Free memory
    assert(dpv_idx == tap->daft_bwd_diff_posns.size());
    std::vector<std::uint8_t>().swap(tap->daft_bwd_diff_posns); //Free memory
    std::vector<std::uint8_t>().swap(tap->daft_bwd_diff_values); //Free memory
    assert(pel_idx == tap->daft_pe_rel_locns.size());
    std::vector<std::uint8_t>().swap(tap->daft_pe_rel_locns); //Free memory
    assert(pep_idx == tap->daft_pe_rel_posns.size());
    std::vector<std::uint8_t>().swap(tap->daft_pe_rel_posns); //Free memory

        return NULL;


        }
    

    else {
    for (std::uint32_t i = index_start; i <= index_end; i++) {
        if (tap->daft_bwd_locns[loc_idx] < 253) {
            curr_locn = (std::uint32_t)tap->daft_bwd_locns[loc_idx];
            loc_idx += 1;
        }
        else if (tap->daft_bwd_locns[loc_idx] == 253) {
            loc_idx += 1;
            curr_locn = (std::uint32_t)tap->daft_bwd_locns[loc_idx];
            loc_idx += 1;
            curr_locn |= (((std::uint32_t)tap->daft_bwd_locns[loc_idx]) << 8);
            loc_idx += 1;
        }
        else if (tap->daft_bwd_locns[loc_idx] == 254) {
            loc_idx += 1;
            curr_locn = (std::uint32_t)tap->daft_bwd_locns[loc_idx];
            loc_idx += 1;
            curr_locn |= (((std::uint32_t)tap->daft_bwd_locns[loc_idx]) << 8);
            loc_idx += 1;
            curr_locn |= (((std::uint32_t)tap->daft_bwd_locns[loc_idx]) << 16);
            loc_idx += 1;
        }
        else {
            loc_idx += 1;
            curr_locn = (std::uint32_t)tap->daft_bwd_locns[loc_idx];
            loc_idx += 1;
            curr_locn |= (((std::uint32_t)tap->daft_bwd_locns[loc_idx]) << 8);
            loc_idx += 1;
            curr_locn |= (((std::uint32_t)tap->daft_bwd_locns[loc_idx]) << 16);
            loc_idx += 1;
            curr_locn |= (((std::uint32_t)tap->daft_bwd_locns[loc_idx]) << 24);
            loc_idx += 1;
        }
        curr_locn += prev_locn;
#if !NDEBUG
        //if (curr_locn > (tap->daft_decom_ds->ref_length - tap->daft_in_args->rdLength)) {
        //    curr_locn = tap->daft_decom_ds->ref_length - tap->daft_in_args->rdLength;
        //}
        //if (curr_locn > (tap->daft_decom_ds->ref_length - tap->daft_in_args->rdLength + tap->daft_bwd_diff_counts[dct_idx])) {
        //    std::cerr << "tap->daft_thread_id : " << tap->daft_thread_id << std::endl;
        //    std::cerr << "prev_locn : " << prev_locn << std::endl;
        //    std::cerr << "curr_locn : " << curr_locn << std::endl;
        //    std::cerr << "i : " << i << std::endl;
        //    std::cerr << "diff_count : " << ((std::uint32_t) tap->daft_bwd_diff_counts[dct_idx]) << std::endl;
        //}
        assert(curr_locn <= (tap->daft_decom_ds->ref_length - tap->daft_in_args->rdLength));
#endif

        if (tap->daft_bwd_diff_counts[dct_idx]) {
            //std::memcpy(my_bases_1, tap->daft_decom_ds->ref_bases + curr_locn, tap->daft_in_args->rdLength + EDIT_DISTANCE);
            //std::uint32_t mb1_idx = 0, mb2_idx = 0;
            //for (std::uint32_t j = 0; j < tap->daft_in_args->rdLength; j++) {
            //   my_bases_2[mb2_idx] = my_bases_1[mb1_idx];
            //   mb1_idx += 1; mb2_idx += 1;
            //}
            //dpv_idx += tap->daft_bwd_diff_counts[dct_idx];
            std::uint32_t ref_idx = 0, mod_idx = 0;
            mod_idx += tap->daft_bwd_diff_posns[dpv_idx];
            is_mod_unused = true;
            tap->daft_bwd_diff_counts[dct_idx] -= 1;
            for (std::uint32_t j = 0; j < tap->daft_in_args->rdLength; ) {
                if (is_mod_unused && (ref_idx == mod_idx)) {
                    if (tap->daft_bwd_diff_values[dpv_idx] < 4) {
                        //Substitution
                        if (tap->daft_bwd_diff_values[dpv_idx] < charToUint8(tap->daft_decom_ds->ref_bases[curr_locn + ref_idx])) {
                            my_bases_2[j] = Uint8Tochar(tap->daft_bwd_diff_values[dpv_idx]);
                        }
                        else {
                            my_bases_2[j] = Uint8Tochar(tap->daft_bwd_diff_values[dpv_idx] + 1);
                        }
                        j += 1;
                        ref_idx += 1;
                    }
                    else if (tap->daft_bwd_diff_values[dpv_idx] < 9) {
                        //Insertion
                        my_bases_2[j] = Uint8Tochar(tap->daft_bwd_diff_values[dpv_idx] - 4);
                        j += 1;
                    }
                    else {
                        //Deletion
                        ref_idx += 1;
                    }
                    dpv_idx += 1;
                    is_mod_unused = false;
                    if (tap->daft_bwd_diff_counts[dct_idx]) {
                        mod_idx += tap->daft_bwd_diff_posns[dpv_idx];
                        is_mod_unused = true;
                        tap->daft_bwd_diff_counts[dct_idx] -= 1;
                    }
                }
                else {
                    my_bases_2[j] = tap->daft_decom_ds->ref_bases[curr_locn + ref_idx];
                    j += 1;
                    ref_idx += 1;
                }
            }
            compact_bwd_read(my_bases_2, tap->daft_in_args->rdLength, tap->daft_decom_ds->bwd_reads[i].read);
        }
        else {
            //No differences
            compact_bwd_read(tap->daft_decom_ds->ref_bases + curr_locn, tap->daft_in_args->rdLength, tap->daft_decom_ds->bwd_reads[i].read);
        }
#if !NDEBUG
        assert(tap->daft_bwd_diff_counts[dct_idx] == 0);
#endif
        dct_idx += 1;

        prev_locn = curr_locn;

        if (tap->daft_pe_rel_posns[pep_idx] < 253) {
            pe_posn = (std::uint32_t)tap->daft_pe_rel_posns[pep_idx];
            pep_idx += 1;
        }
        else if (tap->daft_pe_rel_posns[pep_idx] == 253) {
            pep_idx += 1;
            pe_posn = (std::uint32_t)tap->daft_pe_rel_posns[pep_idx];
            pep_idx += 1;
            pe_posn |= (((std::uint32_t)tap->daft_pe_rel_posns[pep_idx]) << 8);
            pep_idx += 1;
        }
        else if (tap->daft_pe_rel_posns[pep_idx] == 254) {
            pep_idx += 1;
            pe_posn = (std::uint32_t)tap->daft_pe_rel_posns[pep_idx];
            pep_idx += 1;
            pe_posn |= (((std::uint32_t)tap->daft_pe_rel_posns[pep_idx]) << 8);
            pep_idx += 1;
            pe_posn |= (((std::uint32_t)tap->daft_pe_rel_posns[pep_idx]) << 16);
            pep_idx += 1;
        }
        else {
            pep_idx += 1;
            pe_posn = (std::uint32_t)tap->daft_pe_rel_posns[pep_idx];
            pep_idx += 1;
            pe_posn |= (((std::uint32_t)tap->daft_pe_rel_posns[pep_idx]) << 8);
            pep_idx += 1;
            pe_posn |= (((std::uint32_t)tap->daft_pe_rel_posns[pep_idx]) << 16);
            pep_idx += 1;
            pe_posn |= (((std::uint32_t)tap->daft_pe_rel_posns[pep_idx]) << 24);
            pep_idx += 1;
        }
        tap->daft_decom_ds->bwd_reads[i].pe_rel_posn = pe_posn;

        if (tap->daft_pe_rel_locns[pel_idx] < 252) {
            pe_locn_2 = (std::uint64_t)tap->daft_pe_rel_locns[pel_idx];
            pel_idx += 1;
        }
        else if (tap->daft_pe_rel_locns[pel_idx] == 252) {
            pel_idx += 1;
            pe_locn_2 = (std::uint64_t)tap->daft_pe_rel_locns[pel_idx];
            pel_idx += 1;
            pe_locn_2 |= (((std::uint64_t)tap->daft_pe_rel_locns[pel_idx]) << 8);
            pel_idx += 1;
        }
        else if (tap->daft_pe_rel_locns[pel_idx] == 253) {
            pel_idx += 1;
            pe_locn_2 = (std::uint64_t)tap->daft_pe_rel_locns[pel_idx];
            pel_idx += 1;
            pe_locn_2 |= (((std::uint64_t)tap->daft_pe_rel_locns[pel_idx]) << 8);
            pel_idx += 1;
            pe_locn_2 |= (((std::uint64_t)tap->daft_pe_rel_locns[pel_idx]) << 16);
            pel_idx += 1;
        }
        else if (tap->daft_pe_rel_locns[pel_idx] == 254) {
            pel_idx += 1;
            pe_locn_2 = (std::uint64_t)tap->daft_pe_rel_locns[pel_idx];
            pel_idx += 1;
            pe_locn_2 |= (((std::uint64_t)tap->daft_pe_rel_locns[pel_idx]) << 8);
            pel_idx += 1;
            pe_locn_2 |= (((std::uint64_t)tap->daft_pe_rel_locns[pel_idx]) << 16);
            pel_idx += 1;
            pe_locn_2 |= (((std::uint64_t)tap->daft_pe_rel_locns[pel_idx]) << 24);
            pel_idx += 1;
        }
        else {
            pel_idx += 1;
            pe_locn_2 = (std::uint64_t)tap->daft_pe_rel_locns[pel_idx];
            pel_idx += 1;
            pe_locn_2 |= (((std::uint64_t)tap->daft_pe_rel_locns[pel_idx]) << 8);
            pel_idx += 1;
            pe_locn_2 |= (((std::uint64_t)tap->daft_pe_rel_locns[pel_idx]) << 16);
            pel_idx += 1;
            pe_locn_2 |= (((std::uint64_t)tap->daft_pe_rel_locns[pel_idx]) << 24);
            pel_idx += 1;
            pe_locn_2 |= (((std::uint64_t)tap->daft_pe_rel_locns[pel_idx]) << 32);
            pel_idx += 1;
        }
        pe_locn_1 = (int64_t)pe_locn_2;
        if (pe_locn_1 == 0) {
            //Do nothing
        }
        else if (pe_locn_1 % 2) {
            pe_locn_1 = (pe_locn_1 + 1) / 2;
        }
        else {
            pe_locn_1 = pe_locn_1 / (-2);
        }
        pe_locn_1 += tap->daft_decom_ds->pe_rel_locn_mean;
        pe_locn_1 = ((int64_t)curr_locn) - pe_locn_1;
#if !NDEBUG
        assert(pe_locn_1 >= 0);
#endif
        tap->daft_decom_ds->bwd_reads[i].pe_location = (std::uint32_t)pe_locn_1;
    }

    }
    assert(loc_idx == tap->daft_bwd_locns.size());
    std::vector<std::uint8_t>().swap(tap->daft_bwd_locns); //Free memory
    assert(dct_idx == tap->daft_bwd_diff_counts.size());
    std::vector<std::uint8_t>().swap(tap->daft_bwd_diff_counts); //Free memory
    assert(dpv_idx == tap->daft_bwd_diff_posns.size());
    std::vector<std::uint8_t>().swap(tap->daft_bwd_diff_posns); //Free memory
    std::vector<std::uint8_t>().swap(tap->daft_bwd_diff_values); //Free memory
    assert(pel_idx == tap->daft_pe_rel_locns.size());
    std::vector<std::uint8_t>().swap(tap->daft_pe_rel_locns); //Free memory
    assert(pep_idx == tap->daft_pe_rel_posns.size());
    std::vector<std::uint8_t>().swap(tap->daft_pe_rel_posns); //Free memory

    return NULL;


    }

    


int perform_decompaction(InputArgs& in_args, DecompressionDataStructures& decomDS) {
    //File format : Thread count, Read length, Mapped PE read count, Unmapped PE read count
    //Unmapped PE reads, write_se_data, pe_rel_locn_mean, write_pe_data
	double startTime;
    std::size_t num_read;
	std::string cpct_file_name(in_args.comFileName);
	cpct_file_name.append(".cpct");
	FILE * cpct_file_fp = fopen(cpct_file_name.c_str(), "r");
	if (cpct_file_fp == NULL) {
	    std::cerr << "Cannnot open file for reading : " << cpct_file_name << std::endl;
	    assert (0);
	}
	decomDS.ip_fp = &cpct_file_fp;
	FILE * read_seq_fp1 = fopen(in_args.rd1FileName.c_str(), "w");
	if (read_seq_fp1 == NULL) {
	    std::cerr << "Cannnot open file for writing : " << in_args.rd1FileName << std::endl;
	    assert (0);
	}
	decomDS.o1_fp = &read_seq_fp1;
	FILE * read_seq_fp2 = fopen(in_args.rd2FileName.c_str(), "w");
	if (read_seq_fp2 == NULL) {
	    std::cerr << "Cannnot open file for writing : " << in_args.rd2FileName << std::endl;
	    assert (0);
	}
	decomDS.o2_fp = &read_seq_fp2;
    
	num_read = fread(&(in_args.updateReference), 1, sizeof(std::uint32_t), cpct_file_fp);
	if (num_read) {
	    std::cout << "Update reference : " << in_args.updateReference << std::endl;
	    //std::cout << "Bytes read   : " << num_read << std::endl;
	} else {
	    std::cerr << "Error reading from " << cpct_file_name << std::endl;
	    assert (0);
	}
	decomDS.totalBytes += num_read;
	num_read = fread(&(in_args.threadCount), 1, sizeof(std::uint32_t), cpct_file_fp);
	if (num_read) {
	    std::cout << "Thread count : " << in_args.threadCount << std::endl;
	    //std::cout << "Bytes read   : " << num_read << std::endl;
	} else {
	    std::cerr << "Error reading from " << cpct_file_name << std::endl;
	    assert (0);
	}
	decomDS.totalBytes += num_read;
    omp_set_num_threads(in_args.threadCount);
	num_read = fread(&(in_args.rdLength), 1, sizeof(std::uint32_t), cpct_file_fp);
	in_args.rd1Length = in_args.rdLength; in_args.rd2Length = in_args.rdLength;
	if (num_read) {
	    std::cout << "Read length : " << in_args.rdLength << std::endl;
	    //std::cout << "Bytes read  : " << num_read << std::endl;
	} else {
	    std::cerr << "Error reading from " << cpct_file_name << std::endl;
	    assert (0);
	}
	decomDS.totalBytes += num_read;
	num_read = fread(&(decomDS.mapped_count), 1, sizeof(std::size_t), cpct_file_fp);
	if (num_read) {
	    std::cout << "Mapped read count : " << decomDS.mapped_count << std::endl;
	    //std::cout << "Bytes read        : " << num_read << std::endl;
	} else {
	    std::cerr << "Error reading from " << cpct_file_name << std::endl;
	    assert (0);
	}
	decomDS.totalBytes += num_read;
    
    DecompressionStatistics pd_ds;
    pd_ds.ds_pe_mapped_count = decomDS.mapped_count * 2;
    DecompressArgsForThread thread_args_array[in_args.threadCount];
    pthread_t CPUTaskHandle[in_args.threadCount];
    pthread_barrier_t my_barrier;
    std::uint32_t finished_thread_num; int err;
    
    pthread_barrier_init(&my_barrier, NULL, in_args.threadCount);
    for (std::uint32_t i = 0; i < in_args.threadCount; i++) {
        thread_args_array[i].daft_thread_id = i;
        thread_args_array[i].daft_in_args = &in_args;
        thread_args_array[i].daft_decom_ds = &decomDS;
        thread_args_array[i].daft_barrier = &my_barrier;
        thread_args_array[i].dsp = &pd_ds;
    }
    initialize_decompression_stats(in_args, &pd_ds); //Also, allocates capacity
    
    startTime = realtime();
    process_unmapped_reads(in_args, decomDS);
    finished_thread_num = 0; err = 0;
    for (std::uint32_t i = 0; i < in_args.threadCount; i++) {
        err += pthread_create(CPUTaskHandle + i, NULL, decompactReads0Thread, thread_args_array + i);
    }
    if (err == 0) {
        std::cout << "Created threads for decompact-reads-0 successfully" << std::endl;
    } else {
        assert (0);
    }
    
    for (std::uint32_t i = 0; i < in_args.threadCount; i++) {
        err += pthread_join(CPUTaskHandle[i], NULL);
        finished_thread_num++;
    }
    if (err == 0) {
        std::cout << "Threads for decompact-reads-0 completed successfully" << std::endl;
    } else {
        assert (0);
    }
    assert (finished_thread_num == in_args.threadCount);
    write_unm_data(in_args, decomDS);
    read_se_data(thread_args_array);
    decomDS.fileIOTime += (realtime() - startTime);
	
    finished_thread_num = 0; err = 0;
    for (std::uint32_t i = 0; i < in_args.threadCount; i++) {
        err += pthread_create(CPUTaskHandle + i, NULL, decompactReads1Thread, thread_args_array + i);
    }
    if (err == 0) {
        std::cout << "Created threads for decompact-reads-1 successfully" << std::endl;
    } else {
        assert (0);
    }
    
    for (std::uint32_t i = 0; i < in_args.threadCount; i++) {
        err += pthread_join(CPUTaskHandle[i], NULL);
        finished_thread_num++;
    }
    if (err == 0) {
        std::cout << "Threads for decompact-reads-1 completed successfully" << std::endl;
    } else {
        assert (0);
    }
    assert (finished_thread_num == in_args.threadCount);
	num_read = fread(&(decomDS.pe_rel_locn_mean), 1, sizeof(int64_t), cpct_file_fp);
	if (num_read) {
	    std::cout << "pe_rel_locn_mean : " << decomDS.pe_rel_locn_mean << std::endl;
	    //std::cout << "Bytes read        : " << num_read << std::endl;
	} else {
	    std::cerr << "Error reading from " << cpct_file_name << std::endl;
	    assert (0);
	}
	decomDS.totalBytes += num_read;
    startTime = realtime();
    write_fwd_data(thread_args_array);
    decomDS.bwd_reads.resize(decomDS.mapped_count);
    read_pe_data(thread_args_array);
    decomDS.fileIOTime += (realtime() - startTime);
	
    finished_thread_num = 0; err = 0;
    for (std::uint32_t i = 0; i < in_args.threadCount; i++) {
        err += pthread_create(CPUTaskHandle + i, NULL, decompactReads2Thread, thread_args_array + i);
    }
    if (err == 0) {
        std::cout << "Created threads for decompact-reads-2 successfully" << std::endl;
    } else {
        assert (0);
    }
    
    for (std::uint32_t i = 0; i < in_args.threadCount; i++) {
        err += pthread_join(CPUTaskHandle[i], NULL);
        finished_thread_num++;
    }
    if (err == 0) {
        std::cout << "Threads for decompact-reads-2 completed successfully" << std::endl;
    } else {
        assert (0);
    }
    assert (finished_thread_num == in_args.threadCount);
    
    std::cout << "Before bwd sort" << std::endl;
    //TODO[LATER] : Use differences count during sorting
    pss::parallel_stable_sort(decomDS.bwd_reads.begin(), decomDS.bwd_reads.end(), 
    [](const DecomBwdRead& x, const DecomBwdRead& y) {
        return (x.pe_location < y.pe_location) ||
        ((x.pe_location == y.pe_location) && (x.pe_rel_posn < y.pe_rel_posn));
    });
    std::cout << "After bwd sort" << std::endl;
    startTime = realtime();
    decomDS.bwd_print.resize(decomDS.mapped_count * (in_args.rdLength + 3));
    finished_thread_num = 0; err = 0;
    for (std::uint32_t i = 0; i < in_args.threadCount; i++) {
        err += pthread_create(CPUTaskHandle + i, NULL, decompactReads3Thread, thread_args_array + i);
    }
    if (err == 0) {
        std::cout << "Created threads for decompact-reads-3 successfully" << std::endl;
    } else {
        assert (0);
    }
    
    for (std::uint32_t i = 0; i < in_args.threadCount; i++) {
        err += pthread_join(CPUTaskHandle[i], NULL);
        finished_thread_num++;
    }
    if (err == 0) {
        std::cout << "Threads for decompact-reads-3 completed successfully" << std::endl;
    } else {
        assert (0);
    }
    assert (finished_thread_num == in_args.threadCount);
    write_bwd_data(in_args, decomDS);
    decomDS.fileIOTime += (realtime() - startTime);
    
	fclose(cpct_file_fp);
	fclose(read_seq_fp1);
	fclose(read_seq_fp2);
#if !NDEBUG
    display_decompression_stats(in_args, &pd_ds);
#endif
	//std::cout << "r1_count : " << decomDS.r1_count << std::endl;
	//std::cout << "r2_count : " << decomDS.r2_count << std::endl;
	
	return 0;
}

int decompress_reads(InputArgs& in_args, DecompressionDataStructures& decomDS) {
	double startTime;
    
	startTime = realtime();
	std::string cpct_file_name(in_args.comFileName);
	cpct_file_name.append(".cpct");
    std::string my_run_cmd(in_args.bscExecutable);
    my_run_cmd.append(" d ");
    my_run_cmd.append(in_args.comFileName);
    my_run_cmd.append(" ");
    my_run_cmd.append(cpct_file_name);
    //std::cout << my_run_cmd << std::endl;
    system((const char *) my_run_cmd.c_str());
    decomDS.totalTime += (realtime() - startTime);
    std::cout << "Decompressed using bsc in " << realtime() - startTime << " s." << std::endl;
    std::cout << "CPU time : " << cputime() << " s." << std::endl;
    
    startTime = realtime();
    parse_reference_decom(in_args, decomDS);
    decomDS.totalTime += (realtime() - startTime);
    decomDS.fileIOTime += (realtime() - startTime);
    std::cout << "Loaded reference into memory in " << realtime() - startTime << " s." << std::endl;
    std::cout << "CPU time : " << cputime() << " s." << std::endl;
    
    startTime = realtime();
    perform_decompaction(in_args, decomDS);
    decomDS.totalTime += (realtime() - startTime);
    std::cout << "Uncompacted reads in " << realtime() - startTime << " s." << std::endl;
    std::cout << "CPU time : " << cputime() << " s." << std::endl;
    std::cout << "totalBytes for uncompaction: " << decomDS.totalBytes << std::endl;
    
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
    
    my_run_cmd = "ls -ltr ";
    my_run_cmd.append(in_args.rd1FileName.c_str());
#if !NDEBUG
    std::cout << my_run_cmd << std::endl;
    system((const char *) my_run_cmd.c_str());
#endif
//     my_run_cmd = "wc ";
//     my_run_cmd.append(in_args.rd1FileName.c_str());
// #if !NDEBUG
//     std::cout << my_run_cmd << std::endl;
//     system((const char *) my_run_cmd.c_str());
// #endif
    
    my_run_cmd = "ls -ltr ";
    my_run_cmd.append(in_args.rd2FileName.c_str());
#if !NDEBUG
    std::cout << my_run_cmd << std::endl;
    system((const char *) my_run_cmd.c_str());
#endif
//     my_run_cmd = "wc ";
//     my_run_cmd.append(in_args.rd2FileName.c_str());
// #if !NDEBUG
//     std::cout << my_run_cmd << std::endl;
//     system((const char *) my_run_cmd.c_str());
// #endif

#if !NDEBUG
    printf("size of char: %lu\n", sizeof(char));
    printf("size of unsigned int: %lu\n", sizeof(unsigned int));
    printf("size of size_t: %lu\n", sizeof(std::size_t));
#endif
    
	return 0;
}

