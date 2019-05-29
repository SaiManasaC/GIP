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
./REFCOM --operation compression          \
--ref\_file\_name <REFERENCE\_FILE>       \
--rd1\_file\_name <INPUT\_FASTQ\_FILE\_1> \
--rd2\_file\_name <INPUT\_FASTQ\_FILE\_2> \
--com\_file\_name <COMPRESSION\_OUTPUT>   \
--thread\_count <THREAD\_COUNT>
```

### Decompression 

Inside refcom directory,

```sh
./REFCOM --operation decompression         \
--ref\_file\_name <REFERENCE\_FILE>        \
--com\_file\_name <COMPRESSION\_OUTPUT>    \
--rd1\_file\_name <OUTPUT\_FASTQ\_FILE\_1> \
--rd2\_file\_name <OUTPUT\_FASTQ\_FILE\_2>
```

> <REFERENCE\_FILE\> and <COMPRESSION\_OUTPUT\> are the same for compression and decompression

> Decompression uses the same number of threads as compression
