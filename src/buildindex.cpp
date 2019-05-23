#include "buildindex.hpp"
//software.intel.com/en-us/articles/a-parallel-stable-sort-using-c11-for-tbb-cilk-plus-and-openmp
#include "parallel_stable_sort.h"

//using namespace std;

std::vector<HashLocationPair> tmpOccurrenceTable;
std::vector<int> last_hash_values;
std::vector<std::uint32_t> last_hash_counts;

KSEQ_INIT(int, read)
__UTILS_CHAR_UINT8
__UTILS_HASH_VALUE

int parse_reference_com(InputArgs& in_args, CompressionDataStructures& comDS) {
	assert (comDS.ref_length == 0);
	assert (comDS.ref_bases == NULL);
	//comDS.ref_bases = (std::uint8_t *) aligned_alloc(64, sizeof(std::uint8_t) * REF_LEN_MAX);
	comDS.ref_bases = (std::uint8_t *) malloc(sizeof(std::uint8_t) * REF_LEN_MAX);
	assert (comDS.ref_bases);
	
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
			    comDS.ref_bases[val_bases_count] = charToUint8(c);
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
	
	comDS.ref_length = val_bases_count;
	std::cout << "Reference scaffold_count  : " << scaffold_count  << std::endl;
	std::cout << "Reference all_bases_count : " << all_bases_count << std::endl;
	std::cout << "Reference val_bases_count : " << val_bases_count << std::endl;
    
    kseq_destroy(ref_seq);
    fclose(ref_file_fp);
    
    return 0;
}

void *constructIndex1Thread(void *arg) {
    struct IndexArgsForThread * tap;
    tap = (struct IndexArgsForThread *) arg;
    std::uint32_t indices_per_thread = ceil(((double) tmpOccurrenceTable.size())/((double) tap->iaft_in_args->threadCount));
    std::uint32_t index_start = (tap->iaft_thread_id) * indices_per_thread;
    std::uint32_t index_end = ((tap->iaft_thread_id+1) * indices_per_thread) - 1;
    if ((tmpOccurrenceTable.size() - 1) < index_end) {
        index_end = tmpOccurrenceTable.size() - 1;
    }
    //std::cout << tap->iaft_thread_id << " : " << index_start << " : " << index_end << std::endl;
    
    //Compute hash value and location
    for (std::uint32_t i = index_start; i <= index_end; i++) {
        tmpOccurrenceTable[i].location = i;
        int hashVal = hashValue(tap->iaft_com_ds->ref_bases + i, KMER_LENGTH);
#if !NDEBUG
        assert (hashVal >= 0);
#endif
        tmpOccurrenceTable[i].hashValue = hashVal;
    }
    
    return NULL;
}

void *constructIndex2Thread(void *arg) {
    struct IndexArgsForThread * tap;
    tap = (struct IndexArgsForThread *) arg;
    std::uint32_t indices_per_thread = ceil(((double) tap->iaft_com_ds->lookup_table.size())/((double) tap->iaft_in_args->threadCount));
    std::uint32_t index_start = (tap->iaft_thread_id) * indices_per_thread;
    std::uint32_t index_end = ((tap->iaft_thread_id+1) * indices_per_thread) - 1;
    if ((tap->iaft_com_ds->lookup_table.size() - 1) < index_end) {
        index_end = tap->iaft_com_ds->lookup_table.size() - 1;
    }
    //std::cout << tap->iaft_thread_id << " : " << index_start << " : " << index_end << std::endl;
    for (std::uint32_t i = index_start; i <= index_end; i++) {
        tap->iaft_com_ds->lookup_table[i] = 0;
    }
    
    pthread_barrier_wait(tap->iaft_barrier);
    
    indices_per_thread = ceil(((double) tmpOccurrenceTable.size())/((double) tap->iaft_in_args->threadCount));
    index_start = (tap->iaft_thread_id) * indices_per_thread;
    index_end = ((tap->iaft_thread_id+1) * indices_per_thread) - 1;
    if ((tmpOccurrenceTable.size() - 1) < index_end) {
        index_end = tmpOccurrenceTable.size() - 1;
    }
    //std::cout << tap->iaft_thread_id << " : " << index_start << " : " << index_end << std::endl;
    
    std::uint32_t i_prev = index_start;
    std::uint32_t i_curr = index_start;
    tap->iaft_com_ds->occurrence_table[i_curr] = tmpOccurrenceTable[i_curr].location;
    std::uint32_t first_lt_count, curr_lt_count = 1;
    int first_hash_val = tmpOccurrenceTable[i_curr].hashValue;
    i_curr += 1;
    for (; i_curr <= index_end; i_curr++) {
        if (tmpOccurrenceTable[i_curr].hashValue == tmpOccurrenceTable[i_prev].hashValue) {
            tap->iaft_com_ds->occurrence_table[i_curr] = tmpOccurrenceTable[i_curr].location;
            curr_lt_count++;
        } else {
            break;
        }
        i_prev = i_curr;
    }
    first_lt_count = curr_lt_count;
    tap->iaft_com_ds->occurrence_table[i_curr] = tmpOccurrenceTable[i_curr].location;
    curr_lt_count = 1;
    i_prev = i_curr;
    i_curr += 1;
    for (; i_curr <= index_end; i_curr++) {
#if !NDEBUG
        assert (tmpOccurrenceTable[i_prev].hashValue <= tmpOccurrenceTable[i_curr].hashValue);
        if (tmpOccurrenceTable[i_prev].hashValue == tmpOccurrenceTable[i_curr].hashValue) {
            assert (tmpOccurrenceTable[i_prev].location < tmpOccurrenceTable[i_curr].location);
        }
#endif
        
        tap->iaft_com_ds->occurrence_table[i_curr] = tmpOccurrenceTable[i_curr].location;
        if (tmpOccurrenceTable[i_curr].hashValue == tmpOccurrenceTable[i_prev].hashValue) {
            curr_lt_count++;
        } else {
            tap->iaft_com_ds->lookup_table[tmpOccurrenceTable[i_prev].hashValue + 1] = curr_lt_count;
            curr_lt_count = 1;
        }
        i_prev = i_curr;
    }
    
    i_curr -= 1;
    last_hash_values[tap->iaft_thread_id] = tmpOccurrenceTable[i_curr].hashValue;
    last_hash_counts[tap->iaft_thread_id] = curr_lt_count;
    
    pthread_barrier_wait(tap->iaft_barrier);
    
    if (tap->iaft_thread_id == 0) {
        //First thread. Process first value
        tap->iaft_com_ds->lookup_table[first_hash_val + 1] = first_lt_count;
    } else if (tap->iaft_thread_id == (tap->iaft_in_args->threadCount - 1)) {
        //Last thread. Process previous last, first, and last values
        if (last_hash_values[tap->iaft_thread_id - 1] == first_hash_val) {
            tap->iaft_com_ds->lookup_table[first_hash_val + 1] = last_hash_counts[tap->iaft_thread_id - 1] + first_lt_count;
        } else {
            tap->iaft_com_ds->lookup_table[first_hash_val + 1] = first_lt_count;
            tap->iaft_com_ds->lookup_table[last_hash_values[tap->iaft_thread_id - 1] + 1] = last_hash_counts[tap->iaft_thread_id - 1];
        }
        tap->iaft_com_ds->lookup_table[last_hash_values[tap->iaft_thread_id] + 1] = last_hash_counts[tap->iaft_thread_id];
    } else {
        //In-between thread. Process previous last and first values
        if (last_hash_values[tap->iaft_thread_id - 1] == first_hash_val) {
            tap->iaft_com_ds->lookup_table[first_hash_val + 1] = last_hash_counts[tap->iaft_thread_id - 1] + first_lt_count;
        } else {
            tap->iaft_com_ds->lookup_table[first_hash_val + 1] = first_lt_count;
            tap->iaft_com_ds->lookup_table[last_hash_values[tap->iaft_thread_id - 1] + 1] = last_hash_counts[tap->iaft_thread_id - 1];
        }
    }
    
    pthread_barrier_wait(tap->iaft_barrier);
    
    indices_per_thread = ceil(((double) tap->iaft_com_ds->lookup_table.size())/((double) tap->iaft_in_args->threadCount));
    index_start = (tap->iaft_thread_id) * indices_per_thread;
    index_end = ((tap->iaft_thread_id+1) * indices_per_thread) - 1;
    if ((tap->iaft_com_ds->lookup_table.size() - 1) < index_end) {
        index_end = tap->iaft_com_ds->lookup_table.size() - 1;
    }
    //tap->iaft_lt_sum = 0;
    std::uint32_t sumVar = 0;
    for (std::uint32_t i = index_start; i <= index_end; i++) {
        //tap->iaft_lt_sum += tap->iaft_com_ds->lookup_table[i];
        sumVar += tap->iaft_com_ds->lookup_table[i];
        tap->iaft_com_ds->lookup_table[i] = sumVar;
    }
    last_hash_counts[tap->iaft_thread_id] = sumVar;
    
    //std::cout << tap->iaft_thread_id << " : " << tap->iaft_com_ds->lookup_table[index_start] << 
    //" : " << tap->iaft_com_ds->lookup_table[index_end] << " : " << sumVar << std::endl;
    pthread_barrier_wait(tap->iaft_barrier);
    
    sumVar = 0;
    for (std::uint32_t i = 0; i < tap->iaft_thread_id; i++) {
        sumVar += last_hash_counts[i];
    }
    for (std::uint32_t i = index_start; i <= index_end; i++) {
        tap->iaft_com_ds->lookup_table[i] += sumVar;
    }
    //std::cout << tap->iaft_thread_id << " : " << tap->iaft_com_ds->lookup_table[index_start] << 
    //" : " << tap->iaft_com_ds->lookup_table[index_end] << " : " << sumVar << std::endl;
    
    return NULL;
}

/*
void *constructIndex3Thread(void *arg) {
    struct IndexArgsForThread * tap;
    tap = (struct IndexArgsForThread *) arg;
    std::uint32_t indices_per_thread = ceil(((double) tap->iaft_com_ds->lookup_table.size())/((double) tap->iaft_in_args->threadCount));
    std::uint32_t index_start = (tap->iaft_thread_id) * indices_per_thread;
    std::uint32_t index_end = ((tap->iaft_thread_id+1) * indices_per_thread) - 1;
    if ((tap->iaft_com_ds->lookup_table.size() - 1) < index_end) {
        index_end = tap->iaft_com_ds->lookup_table.size() - 1;
    }
    //std::cout << tap->iaft_thread_id << " : " << index_start << " : " << index_end << std::endl;
    for (std::uint32_t i = index_start; i <= index_end; i++) {
        tap->iaft_com_ds->tmp_lookup_table[i] = 0;
    }
    
    pthread_barrier_wait(tap->iaft_barrier);
    
    indices_per_thread = ceil(((double) tmpOccurrenceTable.size())/((double) tap->iaft_in_args->threadCount));
    index_start = (tap->iaft_thread_id) * indices_per_thread;
    index_end = ((tap->iaft_thread_id+1) * indices_per_thread) - 1;
    if ((tmpOccurrenceTable.size() - 1) < index_end) {
        index_end = tmpOccurrenceTable.size() - 1;
    }
    last_hash_counts[tap->iaft_thread_id] = 0;
    for (std::uint32_t i = index_start; i <= index_end; i++) {
        tap->iaft_com_ds->tmp_occurrence_table[i] = tmpOccurrenceTable[i].location;
        if (tap->iaft_com_ds->tmp_occurrence_table[i] != tap->iaft_com_ds->occurrence_table[i]) {
            last_hash_counts[tap->iaft_thread_id] += 1;
        }
    }
    
    return NULL;
}
*/

int construct_index(InputArgs& in_args, CompressionDataStructures& comDS) {
    std::uint32_t hash_table_size = comDS.ref_length - KMER_LENGTH + 1;
    int max_hash_value = ((int) std::pow(4, KMER_LENGTH)) - 1;
    
    comDS.lookup_table.resize(max_hash_value+2);
    comDS.occurrence_table.resize(hash_table_size);
    assert (comDS.lookup_table.size() == ((std::size_t) (max_hash_value+2)));
    assert (comDS.occurrence_table.size() == hash_table_size);
    
    //comDS.tmp_lookup_table.resize(max_hash_value+2);
    //comDS.tmp_occurrence_table.resize(hash_table_size);
    //assert (comDS.tmp_lookup_table.size() == ((std::size_t) (max_hash_value+2)));
    //assert (comDS.tmp_occurrence_table.size() == hash_table_size);
    
#if !NDEBUG
	std::cout << "lookup_table size is     : " << comDS.lookup_table.size() << std::endl;
	std::cout << "occurrence_table size is : " << comDS.occurrence_table.size() << std::endl;
#endif
    
    tmpOccurrenceTable.resize(hash_table_size);
    assert (tmpOccurrenceTable.size() == hash_table_size);
    
    IndexArgsForThread thread_args_array[in_args.threadCount];
    pthread_t CPUTaskHandle[in_args.threadCount];
    pthread_barrier_t my_barrier;
    std::uint32_t finished_thread_num;
    int err;
    
    pthread_barrier_init(&my_barrier, NULL, in_args.threadCount);
    for (std::uint32_t i = 0; i < in_args.threadCount; i++) {
        thread_args_array[i].iaft_barrier = &my_barrier;
        thread_args_array[i].iaft_thread_id = i;
        thread_args_array[i].iaft_in_args = &in_args;
        thread_args_array[i].iaft_com_ds = &comDS;
    }
    last_hash_values.resize(in_args.threadCount);
    last_hash_counts.resize(in_args.threadCount);
    
    finished_thread_num = 0;
    err = 0;
    for (std::uint32_t i = 0; i < in_args.threadCount; i++) {
        err += pthread_create(CPUTaskHandle + i, NULL, constructIndex1Thread, thread_args_array + i);
    }
    if (err == 0) {
        std::cout << "Created threads for construct-index-1 successfully" << std::endl;
    } else {
        assert (0);
    }
    
    for (std::uint32_t i = 0; i < in_args.threadCount; i++) {
        err += pthread_join(CPUTaskHandle[i], NULL);
        finished_thread_num++;
    }
    if (err == 0) {
        std::cout << "Threads for construct-index-1 completed successfully" << std::endl;
    } else {
        assert (0);
    }
    assert (finished_thread_num == in_args.threadCount);
    
    //std::cout << "tmpOccurrenceTable.front().hashValue : " << tmpOccurrenceTable.front().hashValue << std::endl;
    //std::cout << "tmpOccurrenceTable.front().location  : " << tmpOccurrenceTable.front().location << std::endl;
    //std::cout << "tmpOccurrenceTable.back().hashValue  : " << tmpOccurrenceTable.back().hashValue << std::endl;
    //std::cout << "tmpOccurrenceTable.back().location   : " << tmpOccurrenceTable.back().location << std::endl;
    
    std::cout << "Before parallel sort" << std::endl;
    //std::sort(std::execution::par, tmpOccurrenceTable.begin(), tmpOccurrenceTable.end(), 
    pss::parallel_stable_sort(tmpOccurrenceTable.begin(), tmpOccurrenceTable.end(), 
    [](const HashLocationPair& x, const HashLocationPair& y) {
        return (x.hashValue < y.hashValue) ||
        ((x.hashValue == y.hashValue) && (x.location < y.location));
    });
    std::cout << "After parallel sort" << std::endl;
    
    //std::cout << "tmpOccurrenceTable.front().hashValue : " << tmpOccurrenceTable.front().hashValue << std::endl;
    //std::cout << "tmpOccurrenceTable.front().location  : " << tmpOccurrenceTable.front().location << std::endl;
    //std::cout << "tmpOccurrenceTable.back().hashValue  : " << tmpOccurrenceTable.back().hashValue << std::endl;
    //std::cout << "tmpOccurrenceTable.back().location   : " << tmpOccurrenceTable.back().location << std::endl;
    
    finished_thread_num = 0;
    err = 0;
    for (std::uint32_t i = 0; i < in_args.threadCount; i++) {
        err += pthread_create(CPUTaskHandle + i, NULL, constructIndex2Thread, thread_args_array + i);
    }
    if (err == 0) {
        std::cout << "Created threads for construct-index-2 successfully" << std::endl;
    } else {
        assert (0);
    }
    
    for (std::uint32_t i = 0; i < in_args.threadCount; i++) {
        err += pthread_join(CPUTaskHandle[i], NULL);
        finished_thread_num++;
    }
    if (err == 0) {
        std::cout << "Threads for construct-index-2 completed successfully" << std::endl;
    } else {
        assert (0);
    }
    assert (finished_thread_num == in_args.threadCount);
    
/*
    finished_thread_num = 0;
    err = 0;
    for (std::uint32_t i = 0; i < in_args.threadCount; i++) {
        err += pthread_create(CPUTaskHandle + i, NULL, constructIndex3Thread, thread_args_array + i);
    }
    if (err == 0) {
        std::cout << "Created threads for construct-index-3 successfully" << std::endl;
    } else {
        assert (0);
    }
    
    for (std::uint32_t i = 0; i < in_args.threadCount; i++) {
        err += pthread_join(CPUTaskHandle[i], NULL);
        finished_thread_num++;
    }
    if (err == 0) {
        std::cout << "Threads for construct-index-3 completed successfully" << std::endl;
    } else {
        assert (0);
    }
    assert (finished_thread_num == in_args.threadCount);
    
    int currentHashVal;
    for (std::uint32_t i = 0; i < tmpOccurrenceTable.size(); i++) {
        currentHashVal = tmpOccurrenceTable[i].hashValue + 1;
        assert (currentHashVal <= (max_hash_value + 1));
        (comDS.tmp_lookup_table[currentHashVal])++;
    }
    
    std::uint32_t my_sum = 0;
    for (int i = 0; i <= (max_hash_value + 1); i++) {
        my_sum += comDS.tmp_lookup_table[i];
        comDS.tmp_lookup_table[i] = my_sum;
    }
    
    std::uint32_t my_mismatches = 0;
    for (int i = 0; i <= (max_hash_value + 1); i++) {
        if (comDS.tmp_lookup_table[i] != comDS.lookup_table[i]) {
            my_mismatches++;
            std::cout << "mismatch at               : " << i << std::endl;
            std::cout << "comDS.tmp_lookup_table[i] : " << comDS.tmp_lookup_table[i] << std::endl;
            std::cout << "comDS.lookup_table[i]     : " << comDS.lookup_table[i] << std::endl;
        }
    }
    std::cout << "my_mismatches : " << my_mismatches << std::endl;
*/
    
    std::vector<HashLocationPair>().swap(tmpOccurrenceTable); //Free memory
    std::vector<int>().swap(last_hash_values); //Free memory
    std::vector<std::uint32_t>().swap(last_hash_counts); //Free memory
    
    //std::uint32_t total_lt_sum = 0;
    //for (std::uint32_t i = 0; i < in_args.threadCount; i++) {
    //    total_lt_sum += thread_args_array[i].iaft_lt_sum;
    //    //std::cout << "i : " << i << ", " << thread_args_array[i].iaft_lt_sum << std::endl;
    //}
    //std::cout << "total_lt_sum : " << total_lt_sum << std::endl;
    
    //std::cout << "Before inclusive_scan" << std::endl;
    //std::inclusive_scan(std::execution::par, comDS.lookup_table.begin(), comDS.lookup_table.end(), comDS.lookup_table.begin());
    //std::cout << "After inclusive_scan" << std::endl;
    
    //std::cout << "lookup_table.front() : " << comDS.lookup_table.front() << std::endl;
    //std::cout << "lookup_table.back()  : " << comDS.lookup_table.back() << std::endl;
    //std::cout << "lookup_table[ 0] : " << comDS.tmp_lookup_table[0] << std::endl;
    //std::cout << "lookup_table[ 1] : " << comDS.tmp_lookup_table[1] << std::endl;
    //std::cout << "lookup_table[-2] : " << comDS.tmp_lookup_table[max_hash_value] << std::endl;
    //std::cout << "lookup_table[-1] : " << comDS.tmp_lookup_table[max_hash_value+1] << std::endl;
#if !NDEBUG
    std::cout << "lookup_table[ 0] : " << comDS.lookup_table[0] << std::endl;
    std::cout << "lookup_table[ 1] : " << comDS.lookup_table[1] << std::endl;
    std::cout << "lookup_table[-2] : " << comDS.lookup_table[max_hash_value] << std::endl;
    std::cout << "lookup_table[-1] : " << comDS.lookup_table[max_hash_value+1] << std::endl;
#endif
    
    return 0;
}

int build_index(InputArgs& in_args, CompressionDataStructures& comDS) {
	double startTime;
    
    startTime = realtime();
    parse_reference_com(in_args, comDS);
    comDS.totalTime += (realtime() - startTime);
    std::cout << "Loaded reference into memory in " << realtime() - startTime << " s." << std::endl;
    std::cout << "CPU time : " << cputime() << " s." << std::endl;
    
    startTime = realtime();
    construct_index(in_args, comDS);
    comDS.totalTime += (realtime() - startTime);
    std::cout << "Built index in " << realtime() - startTime << " s." << std::endl;
    std::cout << "CPU time : " << cputime() << " s." << std::endl;
    
	return 0;
}























































