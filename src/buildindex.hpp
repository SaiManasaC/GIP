#ifndef BUILDINDEX_HPP_
#define BUILDINDEX_HPP_

#include "utils.hpp"
#include "kseq.h"

struct HashLocationPair {
    int hashValue;
    std::uint32_t location;
};

struct IndexArgsForThread{
    pthread_barrier_t * iaft_barrier;
    std::uint32_t iaft_thread_id;
    InputArgs * iaft_in_args;
    CompressionDataStructures * iaft_com_ds;
    std::uint32_t iaft_lt_sum;
};

int parse_reference_com(InputArgs&, CompressionDataStructures&);
int construct_index(InputArgs&, CompressionDataStructures&);
int build_index(InputArgs&, CompressionDataStructures&);

#endif /* BUILDINDEX_HPP_ */

