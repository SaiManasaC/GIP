#include "aligner.hpp"

//KSEQ_INIT(int, read)
__UTILS_CHAR_UINT8
//__UTILS_HASH_VALUE

int getSingleRead(Read * read, kseq_t * read_seq, std::uint32_t& rdLength) {
	int l = kseq_read(read_seq);
	while (l == 0) {
		l = kseq_read(read_seq);
		//printf("!");
	}
	if (l > 0) {
		if (rdLength == 0) {
		    rdLength = (std::uint32_t) l;
		}
		if (l != ((int) rdLength)) {assert (0);}
		read->length = l;
		//memset(read->name, '\0', sizeof(char) * READ_NAME_LEN_MAX);
		//Redundant memset(read->charBases, '\0', sizeof(char) * READ_LEN_MAX);
		//strcpy(read->name, read_seq->name.s);
		strcpy(read->charBases, read_seq->seq.s);
		int i;
		for (i = 0; i < l; ++i) {
			char c = read_seq->seq.s[i];
			uint8_t base = charToUint8(c);
			read->bases[i] = base;
			if (base != 4) {
				read->rc_bases[l - 1 - i] = 3 - base;
			} else {
				read->rc_bases[l - 1 - i] = 4;
			}
		}
		return 1;
	} else if (l == -1) {
		/*end of file*/
		return -1;
	} else {
        std::cerr << "Truncated quality string" << std::endl;
        assert (0);
	}
	/*never come to here*/
	assert (0);
	return 0;
}

bool fetchReads(std::queue<Read>& fr_queue_1, std::queue<Read>& fr_queue_2, AlignArgsForThread * fr_tap) {
    kseq_t * fr_seq_1 = *(fr_tap->aaft_read_seq1);
    kseq_t * fr_seq_2 = *(fr_tap->aaft_read_seq2);
    Read fr_read_1, fr_read_2;
    int fr_flag_1, fr_flag_2;
    
    for (std::uint32_t i = 0; i < READ_BLOCK_SIZE; i++) {
        fr_flag_1 = getSingleRead(&fr_read_1, fr_seq_1, fr_tap->aaft_in_args->rd1Length);
        fr_flag_2 = getSingleRead(&fr_read_2, fr_seq_2, fr_tap->aaft_in_args->rd2Length);
        if (fr_tap->aaft_in_args->rdLength == 0) {
            assert (fr_tap->aaft_in_args->rd1Length == fr_tap->aaft_in_args->rd2Length);
            fr_tap->aaft_in_args->rdLength = fr_tap->aaft_in_args->rd1Length;
        }
        if ((fr_flag_1 == -1) && (fr_flag_2 == -1)) {
            break;
        } else if ((fr_flag_1 == -1) || (fr_flag_2 == -1)) {
            assert (0);
        }
        fr_queue_1.push(fr_read_1);
        fr_queue_2.push(fr_read_2);
    }
    
    if ((fr_flag_1 == -1) && (fr_flag_2 == -1)) {
        return true;
    } else {
        return false;
    }
}

void populateResults(std::queue<MappedRead>& pr_f_read_queue, std::queue<MappedRead>& pr_b_read_queue, 
std::queue<UnmappedRead>& pr_u_read_queue, AlignArgsForThread * pr_tap) {
    UnmappedRead pr_unm_read;
    MappedRead pr_fwd_read, pr_bwd_read;

    while (!pr_f_read_queue.empty()) {
        pr_fwd_read = pr_f_read_queue.front();
        pr_f_read_queue.pop();
        pr_fwd_read.id = (std::uint32_t) (pr_tap->aaft_com_ds->fwd_mapped_reads).size();
        (pr_tap->aaft_com_ds->fwd_mapped_reads).push_back(pr_fwd_read);
        pr_bwd_read = pr_b_read_queue.front();
        pr_b_read_queue.pop();
        pr_bwd_read.id = (std::uint32_t) (pr_tap->aaft_com_ds->bwd_mapped_reads).size();
        (pr_tap->aaft_com_ds->bwd_mapped_reads).push_back(pr_bwd_read);
    }
#if !NDEBUG
    assert (pr_b_read_queue.empty());
#endif

    while (!pr_u_read_queue.empty()) {
        pr_unm_read = pr_u_read_queue.front();
        pr_u_read_queue.pop();
        (pr_tap->aaft_com_ds->unmapped_reads).push_back(pr_unm_read);
    }
}

void *alignReads1Thread(void *arg) {
    struct AlignArgsForThread * tap;
    tap = (struct AlignArgsForThread *) arg;
    //std::uint32_t my_id = tap->aaft_thread_id;
    //InputArgs * my_in_args = tap->aaft_in_args;
    //CompressionDataStructures * my_com_ds = tap->aaft_com_ds;
    pthread_mutex_t * my_mutex = tap->aaft_mutex;
    bool is_reading_complete = false;
    bool are_both_reads_mapped = false;
    
    std::queue<Read> one_read_queue;
    std::queue<Read> two_read_queue;
    std::queue<MappedRead> fwd_read_queue;
    std::queue<MappedRead> bwd_read_queue;
    std::queue<UnmappedRead> unm_read_queue;
    VerifyArgsForThread vaft;
    vaft.vaft_thread_id = tap->aaft_thread_id;
    vaft.vaft_in_args = tap->aaft_in_args;
    vaft.vaft_com_ds = tap->aaft_com_ds;
    vaft.asp = tap->asp;
    vaft.read_1 = &(tap->read_1);
    vaft.read_2 = &(tap->read_2);
    vaft.fwd_read = &(tap->fwd_read);
    vaft.bwd_read = &(tap->bwd_read);
    vaft.unm_read = &(tap->unm_read);
    
    while (1) {
        while (!one_read_queue.empty()) {
            tap->read_1 = one_read_queue.front();
            one_read_queue.pop();
            tap->read_2 = two_read_queue.front();
            two_read_queue.pop();
            are_both_reads_mapped = NJ_CPUVBM(&vaft); //Align reads
            if (are_both_reads_mapped) {
                fwd_read_queue.push(tap->fwd_read);
                bwd_read_queue.push(tap->bwd_read);
            } else {
                unm_read_queue.push(tap->unm_read);
            }
            if (one_read_queue.size() <= 1024) {
                if (one_read_queue.size() % 128 == 0) {
                    //Try lock
                    if (!pthread_mutex_trylock(my_mutex)) {
                        //Fetch reads and populate results
                        is_reading_complete = fetchReads(one_read_queue, two_read_queue, tap);
                        populateResults(fwd_read_queue, bwd_read_queue, unm_read_queue, tap);
                        pthread_mutex_unlock(my_mutex);
                    }
                }
            }
            if (is_reading_complete) {break;}
        }
        if (is_reading_complete) {break;}
        pthread_mutex_lock(my_mutex);
        //Fetch reads and populate results
        is_reading_complete = fetchReads(one_read_queue, two_read_queue, tap);
        populateResults(fwd_read_queue, bwd_read_queue, unm_read_queue, tap);
        pthread_mutex_unlock(my_mutex);
        if (is_reading_complete) {break;}
    }
    
    while (!one_read_queue.empty()) {
            tap->read_1 = one_read_queue.front();
            one_read_queue.pop();
            tap->read_2 = two_read_queue.front();
            two_read_queue.pop();
            are_both_reads_mapped = NJ_CPUVBM(&vaft); //Align reads
            if (are_both_reads_mapped) {
                fwd_read_queue.push(tap->fwd_read);
                bwd_read_queue.push(tap->bwd_read);
            } else {
                unm_read_queue.push(tap->unm_read);
            }
    }
    pthread_mutex_lock(my_mutex);
    //Populate results
    assert (is_reading_complete); //Done fetching reads
    populateResults(fwd_read_queue, bwd_read_queue, unm_read_queue, tap);
    pthread_mutex_unlock(my_mutex);
    
    return NULL;
}

void initialize_alignment_stats(InputArgs& in_args, AlignmentStatistics * ias_asp) {
    (ias_asp->as_unmapped_reads_n_count).resize(in_args.threadCount, 0);
    (ias_asp->as_unmapped_diffs_n_count).resize(in_args.threadCount, 0);
    (ias_asp->as_paired_read_count).resize(in_args.threadCount, 0);
    (ias_asp->as_paired_mapped_count).resize(in_args.threadCount, 0);
    (ias_asp->as_single_mapped_count).resize(in_args.threadCount, 0);
    (ias_asp->as_both_unmapped_count).resize(in_args.threadCount, 0);
    (ias_asp->as_can_unmapped_count).resize(in_args.threadCount, 0);
    (ias_asp->as_pla_unmapped_count).resize(in_args.threadCount, 0);
    (ias_asp->as_nil_unmapped_count).resize(in_args.threadCount, 0);
    
    (ias_asp->as_differences_count).resize((in_args.threadCount * 16), 0);
    
    (ias_asp->as_less_than_length_count).resize(in_args.threadCount, 0);
    (ias_asp->as_equal_to_length_count).resize(in_args.threadCount, 0);
    (ias_asp->as_more_than_length_count).resize(in_args.threadCount, 0);
    (ias_asp->as_neg_diff_count_bef).resize(in_args.threadCount, 0);
    (ias_asp->as_neg_diff_count_aft).resize(in_args.threadCount, 0);
    (ias_asp->as_cigar_unequal_count).resize(in_args.threadCount, 0);
    (ias_asp->as_location_unequal_count).resize(in_args.threadCount, 0);
    (ias_asp->as_soft_clipped_count).resize(in_args.threadCount, 0);
    
    (ias_asp->as_paired_read_order).resize((in_args.threadCount * 8), 0);
}

void display_alignment_stats(InputArgs& in_args, AlignmentStatistics * das_asp) {
    std::size_t total_paired_read_count = 0;
    std::size_t total_paired_mapped_count = 0;
    std::size_t total_single_mapped_count = 0;
    std::size_t total_both_unmapped_count = 0;
    std::size_t total_can_unmapped_count = 0;
    std::size_t total_pla_unmapped_count = 0;
    std::size_t total_nil_unmapped_count = 0;
    std::size_t total_less_than_length_count = 0;
    std::size_t total_equal_to_length_count = 0;
    std::size_t total_more_than_length_count = 0;
    std::size_t total_location_unequal_count = 0;
    std::size_t total_cigar_unequal_count = 0;
    std::size_t total_neg_diff_count_bef = 0;
    std::size_t total_neg_diff_count_aft = 0;
    std::size_t total_soft_clipped_count = 0;
    std::size_t total_unmapped_reads_n_count = 0;
    std::size_t total_unmapped_diffs_n_count = 0;
    for (std::uint32_t i = 0; i < in_args.threadCount; i++) {
        total_paired_read_count += (das_asp->as_paired_read_count)[i];
        total_paired_mapped_count += (das_asp->as_paired_mapped_count)[i];
        total_single_mapped_count += (das_asp->as_single_mapped_count)[i];
        total_both_unmapped_count += (das_asp->as_both_unmapped_count)[i];
        total_can_unmapped_count += (das_asp->as_can_unmapped_count)[i];
        total_pla_unmapped_count += (das_asp->as_pla_unmapped_count)[i];
        total_nil_unmapped_count += (das_asp->as_nil_unmapped_count)[i];
        total_less_than_length_count += (das_asp->as_less_than_length_count)[i];
        total_equal_to_length_count += (das_asp->as_equal_to_length_count)[i];
        total_more_than_length_count += (das_asp->as_more_than_length_count)[i];
        total_location_unequal_count += (das_asp->as_location_unequal_count)[i];
        total_cigar_unequal_count += (das_asp->as_cigar_unequal_count)[i];
        total_neg_diff_count_bef += (das_asp->as_neg_diff_count_bef)[i];
        total_neg_diff_count_aft += (das_asp->as_neg_diff_count_aft)[i];
        total_soft_clipped_count += (das_asp->as_soft_clipped_count)[i];
        total_unmapped_reads_n_count += (das_asp->as_unmapped_reads_n_count)[i];
        total_unmapped_diffs_n_count += (das_asp->as_unmapped_diffs_n_count)[i];
    }
    std::size_t total_differences_count[16] = {0};
    std::size_t all_differences_count = 0;
    for (std::uint32_t i = 0; i < in_args.threadCount; i++) {
        for (std::uint32_t j = 0; j < 16; j++) {
            total_differences_count[j] += (das_asp->as_differences_count)[i*16 + j];
        }
    }
    std::size_t total_paired_read_order[8] = {0};
    for (std::uint32_t i = 0; i < in_args.threadCount; i++) {
        for (std::uint32_t j = 0; j < 8; j++) {
            total_paired_read_order[j] += (das_asp->as_paired_read_order)[i*8 + j];
        }
    }
    fprintf(stdout, "total_paired_read_count : %lu.\n", total_paired_read_count);
    fprintf(stdout, "total_paired_mapped_count : %lu, %.2lf.\n", total_paired_mapped_count,
    (((double) total_paired_mapped_count)/((double) total_paired_read_count))*100);
    fprintf(stdout, "total_single_mapped_count : %lu, %.2lf.\n", total_single_mapped_count,
    (((double) total_single_mapped_count)/((double) total_paired_read_count))*100);
    fprintf(stdout, "total_both_unmapped_count : %lu, %.2lf.\n", total_both_unmapped_count,
    (((double) total_both_unmapped_count)/((double) total_paired_read_count))*100);
    fprintf(stdout, "total_can_unmapped_count  :   %lu, %.2lf.\n", total_can_unmapped_count,
    (((double) total_can_unmapped_count)/((double) total_paired_read_count))*100);
    fprintf(stdout, "total_pla_unmapped_count  :   %lu, %.2lf.\n", total_pla_unmapped_count,
    (((double) total_pla_unmapped_count)/((double) total_paired_read_count))*100);
    fprintf(stdout, "total_nil_unmapped_count  :   %lu, %.2lf.\n", total_nil_unmapped_count,
    (((double) total_nil_unmapped_count)/((double) total_paired_read_count))*100);
    fprintf(stdout, "total_less_than_length_count :  %lu.\n", total_less_than_length_count);
    fprintf(stdout, "total_equal_to_length_count  :  %lu.\n", total_equal_to_length_count);
    fprintf(stdout, "total_more_than_length_count :  %lu.\n", total_more_than_length_count);
    fprintf(stdout, "total_location_unequal_count :  %lu.\n", total_location_unequal_count);
    fprintf(stdout, "total_cigar_unequal_count    :  %lu.\n", total_cigar_unequal_count);
    fprintf(stdout, "total_neg_diff_count_bef     :  %lu.\n", total_neg_diff_count_bef);
    fprintf(stdout, "total_neg_diff_count_aft     :  %lu.\n", total_neg_diff_count_aft);
    fprintf(stdout, "total_soft_clipped_count     :  %lu.\n", total_soft_clipped_count);
    std::size_t total_single_read_count = total_paired_read_count * 2;
    fprintf(stdout, "total_unmapped_reads_n_count : %lu, %.2lf.\n", total_unmapped_reads_n_count,
    (((double) total_unmapped_reads_n_count)/((double) total_single_read_count))*100);
    fprintf(stdout, "total_unmapped_diffs_n_count : %lu.\n", total_unmapped_diffs_n_count);
    for (std::uint32_t j = 0; j < 16; j++) {
        fprintf(stdout, "total_differences_count[%2u] : %lu, %.2lf.\n", j, total_differences_count[j],
        (((double) total_differences_count[j])/((double) total_single_read_count))*100);
        all_differences_count += total_differences_count[j];
    }
    fprintf(stdout, "all_differences_count : %lu, %.2lf.\n", all_differences_count,
    (((double) all_differences_count)/((double) total_single_read_count))*100);
    for (std::uint32_t j = 0; j < 8; j++) {
        fprintf(stdout, "total_paired_read_order[%u] : %lu, %.2lf.\n", j, total_paired_read_order[j],
        (((double) total_paired_read_order[j])/((double) total_paired_mapped_count))*100);
    }
}

int perform_alignment(InputArgs& in_args, CompressionDataStructures& comDS) {
    kseq_t * read_seq1;
    kseq_t * read_seq2;
	FILE * read_seq_fp1 = fopen(in_args.rd1FileName.c_str(), "r");
	if (read_seq_fp1 == NULL) {
	    std::cerr << "Cannnot open file for reading : " << in_args.rd1FileName << std::endl;
	    assert (0);
	} else {
		read_seq1 = kseq_init(fileno(read_seq_fp1));
	}
	FILE * read_seq_fp2 = fopen(in_args.rd2FileName.c_str(), "r");
	if (read_seq_fp2 == NULL) {
	    std::cerr << "Cannnot open file for reading : " << in_args.rd2FileName << std::endl;
	    assert (0);
	} else {
		read_seq2 = kseq_init(fileno(read_seq_fp2));
	}
    
    AlignmentStatistics pa_as;
    AlignArgsForThread thread_args_array[in_args.threadCount];
    pthread_t CPUTaskHandle[in_args.threadCount];
    pthread_mutex_t alignment_mutex;
    std::uint32_t finished_thread_num; int err;
    
    if (pthread_mutex_init(&alignment_mutex, NULL)) {assert (0);}
    for (std::uint32_t i = 0; i < in_args.threadCount; i++) {
        thread_args_array[i].aaft_thread_id = i;
        thread_args_array[i].aaft_in_args = &in_args;
        thread_args_array[i].aaft_com_ds = &comDS;
        thread_args_array[i].aaft_mutex = &alignment_mutex;
        thread_args_array[i].aaft_read_seq1 = &read_seq1;
        thread_args_array[i].aaft_read_seq2 = &read_seq2;
        thread_args_array[i].asp = &pa_as;
    }
    initialize_alignment_stats(in_args, &pa_as); //Also, allocates capacity
    
    finished_thread_num = 0; err = 0;
    for (std::uint32_t i = 0; i < in_args.threadCount; i++) {
        err += pthread_create(CPUTaskHandle + i, NULL, alignReads1Thread, thread_args_array + i);
    }
    if (err == 0) {
        std::cout << "Created threads for align-reads-1 successfully" << std::endl;
    } else {
        assert (0);
    }
    
    for (std::uint32_t i = 0; i < in_args.threadCount; i++) {
        err += pthread_join(CPUTaskHandle[i], NULL);
        finished_thread_num++;
    }
    if (err == 0) {
        std::cout << "Threads for align-reads-1 completed successfully" << std::endl;
    } else {
        assert (0);
    }
    assert (finished_thread_num == in_args.threadCount);
    
#if !NDEBUG
    display_alignment_stats(in_args, &pa_as);
#endif
    
    if (pthread_mutex_destroy(&alignment_mutex)) {assert (0);}
	kseq_destroy(read_seq1);
	fclose(read_seq_fp1);
	kseq_destroy(read_seq2);
	fclose(read_seq_fp2);
    
    std::vector<std::uint32_t>().swap(comDS.lookup_table); //Free memory
    std::vector<std::uint32_t>().swap(comDS.occurrence_table); //Free memory
    
    std::vector<MappedRead>(comDS.fwd_mapped_reads).swap(comDS.fwd_mapped_reads);
    std::vector<MappedRead>(comDS.bwd_mapped_reads).swap(comDS.bwd_mapped_reads);
    std::vector<UnmappedRead>(comDS.unmapped_reads).swap(comDS.unmapped_reads);
	
	return 0;
}

int align_reads(InputArgs& in_args, CompressionDataStructures& comDS) {
	double startTime;
    
    startTime = realtime();
    perform_alignment(in_args, comDS);
    comDS.totalTime += (realtime() - startTime);
    std::cout << "Aligned reads in " << realtime() - startTime << " s." << std::endl;
    std::cout << "CPU time : " << cputime() << " s." << std::endl;
    
	return 0;
}

















































