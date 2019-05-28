#ifndef REFCOM_HPP_
#define REFCOM_HPP_

#include "utils.hpp"
#include "argvparser.hpp"
#include "buildindex.hpp"
#include "aligner.hpp"
#include "compress.hpp"
#include "decompress.hpp"

int parse_args(int, char**, InputArgs&);
int refcom_compression(InputArgs&, CompressionDataStructures&);
int refcom_decompression(InputArgs&, DecompressionDataStructures&);

#endif // REFCOM_HPP_
























































