Parallel Reference-based Compression of Paired-end Genomics Read Datasets
=========================================================================

## Dependencies

- `g++` with support for c++11
- `libbsc` software (included in this project)

## Compilation

Inside refcom directory,

```sh
cd libbsc
make
cd ..
make
```

## Execution

### Compression 

Inside refcom directory,

```sh
./REFCOM --operation compression     \
--ref_file_name <REFERENCE_FILE>     \
--rd1_file_name <INPUT_FASTQ_FILE_1> \
--rd2_file_name <INPUT_FASTQ_FILE_2> \
--com_file_name <COMPRESSION_OUTPUT> \
--thread_count <THREAD_COUNT>
```

### Decompression 

Inside refcom directory,

```sh
./REFCOM --operation decompression    \
--ref_file_name <REFERENCE_FILE>      \
--com_file_name <COMPRESSION_OUTPUT>  \
--rd1_file_name <OUTPUT_FASTQ_FILE_1> \
--rd2_file_name <OUTPUT_FASTQ_FILE_2>
```

> <REFERENCE_FILE\> and <COMPRESSION_OUTPUT\> are the same for both compression and decompression

> Decompression uses the same number of threads as compression
