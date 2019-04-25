#ifndef ALIGNER_HPP_
#define ALIGNER_HPP_

#include "utils.hpp"
#include "verifier.hpp"
#include "kseq.h"

KSEQ_INIT(int, read)

struct AlignArgsForThread{
    std::uint32_t aaft_thread_id;
    InputArgs * aaft_in_args;
    CompressionDataStructures * aaft_com_ds;
    pthread_mutex_t * aaft_mutex;
    kseq_t ** aaft_read_seq1;
    kseq_t ** aaft_read_seq2;
    AlignmentStatistics * asp;
    Read read_1;
    Read read_2;
    MappedRead fwd_read;
    MappedRead bwd_read;
    UnmappedRead unm_read;
};

int perform_alignment(InputArgs&, CompressionDataStructures&);
int align_reads(InputArgs&, CompressionDataStructures&);

#endif /* ALIGNER_HPP_ */

