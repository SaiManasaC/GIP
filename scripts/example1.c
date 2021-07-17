#include <stdio.h>
#include <zlib.h>
#include <string.h>
#include <errno.h>
#include <assert.h>
#include "bwamem.h"
#include "kseq.h" // for the FASTA/Q parser
#include "utils.h"

KSEQ_DECLARE(gzFile)

typedef struct {
	kseq_t *ks, *ks2;
	mem_opt_t *opt; //match score and other option for given subsequece/bases
	mem_pestat_t *pes0; //data regarding read pair
	int64_t n_processed;
	int copy_comment, actual_chunk_size;
	bwaidx_t *idx; //FM index, reference sequence, 
} ktp_aux_t;

typedef struct {
	ktp_aux_t *aux;
	int n_seqs;
	bseq1_t *seqs;
} ktp_data_t;

static void *process(void *shared, int step, void *_data)
{
	ktp_aux_t *aux = (ktp_aux_t*)shared;
	ktp_data_t *data = (ktp_data_t*)_data;
	int i;
	if (step == 0) {
		ktp_data_t *ret;
		int64_t size = 0;
		ret = calloc(1, sizeof(ktp_data_t));
		ret->seqs = bseq_read(aux->actual_chunk_size, &ret->n_seqs, aux->ks, aux->ks2); // from kseq to bseq 
		if (ret->seqs == 0) {
			free(ret);
			return 0;
		}
		if (!aux->copy_comment)
			for (i = 0; i < ret->n_seqs; ++i) {
				free(ret->seqs[i].comment);
				ret->seqs[i].comment = 0;
			}
		for (i = 0; i < ret->n_seqs; ++i) size += ret->seqs[i].l_seq;
		if (bwa_verbose >= 3)
			fprintf(stderr, "[M::%s] read %d sequences (%ld bp)...\n", __func__, ret->n_seqs, (long)size);
		return ret;
	} else if (step == 1) {
		const mem_opt_t *opt = aux->opt;
		const bwaidx_t *idx = aux->idx;
		//check flag output
		if (opt->flag & MEM_F_SMARTPE) {
			bseq1_t *sep[2];
			int n_sep[2];
			mem_opt_t tmp_opt = *opt;
			bseq_classify(data->n_seqs, data->seqs, n_sep, sep);
			if (bwa_verbose >= 3)
				fprintf(stderr, "[M::%s] %d single-end sequences; %d paired-end sequences\n", __func__, n_sep[0], n_sep[1]);
			if (n_sep[0]) {
				tmp_opt.flag &= ~MEM_F_PE;
				mem_process_seqs(&tmp_opt, idx->bwt, idx->bns, idx->pac, aux->n_processed, n_sep[0], sep[0], 0);
				for (i = 0; i < n_sep[0]; ++i)
					data->seqs[sep[0][i].id].sam = sep[0][i].sam;
			}
			if (n_sep[1]) {
				tmp_opt.flag |= MEM_F_PE;
				mem_process_seqs(&tmp_opt, idx->bwt, idx->bns, idx->pac, aux->n_processed + n_sep[0], n_sep[1], sep[1], aux->pes0);
				for (i = 0; i < n_sep[1]; ++i)
					data->seqs[sep[1][i].id].sam = sep[1][i].sam;
			}
			free(sep[0]); free(sep[1]);
		} else mem_process_seqs(opt, idx->bwt, idx->bns, idx->pac, aux->n_processed, data->n_seqs, data->seqs, aux->pes0);
		aux->n_processed += data->n_seqs;
		return data;
	} else if (step == 2) {
		for (i = 0; i < data->n_seqs; ++i) {
			if (data->seqs[i].sam) err_fputs(data->seqs[i].sam, stdout);
			free(data->seqs[i].name); free(data->seqs[i].comment);
			free(data->seqs[i].seq); free(data->seqs[i].qual); free(data->seqs[i].sam);
		}
		free(data->seqs); free(data);
		return 0;
	}
	return 0;
}

static void update_a(mem_opt_t *opt, const mem_opt_t *opt0)
{
	if (opt0->a) { // matching score is changed
		if (!opt0->b) opt->b *= opt->a;
		if (!opt0->T) opt->T *= opt->a;
		if (!opt0->o_del) opt->o_del *= opt->a;
		if (!opt0->e_del) opt->e_del *= opt->a;
		if (!opt0->o_ins) opt->o_ins *= opt->a;
		if (!opt0->e_ins) opt->e_ins *= opt->a;
		if (!opt0->zdrop) opt->zdrop *= opt->a;
		if (!opt0->pen_clip5) opt->pen_clip5 *= opt->a;
		if (!opt0->pen_clip3) opt->pen_clip3 *= opt->a;
		if (!opt0->pen_unpaired) opt->pen_unpaired *= opt->a;
	}
}

int main(int argc, char *argv[])
{
    bwaidx_t *idx;
	gzFile fp, fp2;
	kseq_t *ks;
	mem_opt_t *opt, opt0;
	int fd, fd2;
	void *ko = 0, *ko2 = 0;
	ktp_aux_t aux;
	int no_mt_io = 0;
	char *p, *rg_line = 0, *hdr_line = 0;

	int fixed_chunk_size = -1;

	memset(&aux, 0, sizeof(ktp_aux_t));

	if (argc < 4) {
		fprintf(stderr, "Usage: bwamem  <idx.base> <read1s.fq> <read2s.fq> \n");
		return 1;
	}

    ks = kseq_init(fp); // initialize the FASTA/Q parser
	aux.opt = opt = mem_opt_init(); // initialize the BWA-MEM parameters to the default values
	memset(&opt0, 0, sizeof(mem_opt_t));
	// if we want to set, a, b and n_threads value
	/*
	opt->a = 3; opt0.a = 1;
	opt->b = 5; opt0.b = 1;
	update_a(opt, &opt0);
	opt->n_threads = 4; // to set number of threads
	*/
	bwa_fill_scmat(opt->a, opt->b, opt->mat);
	fprintf(stderr, "ref file %s", argv[1]);

	aux.idx = bwa_idx_load_from_shm(argv[1]);
	if (aux.idx == 0) {
		if ((aux.idx = bwa_idx_load(argv[1], 0x7)) == 0) return 1;
	} else if (bwa_verbose >= 3)
		fprintf(stderr, "[M::%s] load the bwa index from shared memory\n", __func__);

	fprintf(stderr, "ref loaded");

	fprintf(stderr, "main file1 %s", argv[2]);
	ko = kopen(argv[2], &fd);
	if (ko == 0) {
		if (bwa_verbose >= 1) fprintf(stderr, "[E::%s] fail to open file `%s'.\n", __func__, argv[2]);
		return 1;
	}
	fp = gzdopen(fd, "r");
	aux.ks = kseq_init(fp);

	fprintf(stderr, "main file1 loaded");
	fprintf(stderr, "main file2 %s", argv[3]);
	ko2 = kopen(argv[3], &fd2);
	if (ko2 == 0) {
		if (bwa_verbose >= 1) fprintf(stderr, "[E::%s] fail to open file `%s'.\n", __func__, argv[3]);
		return 1;
	}
	fp2 = gzdopen(fd2, "r");
	aux.ks2 = kseq_init(fp2);

	fprintf(stderr, "main file2 loaded");
	opt->flag |= MEM_F_PE;

	bwa_print_sam_hdr(aux.idx->bns, hdr_line);
	aux.actual_chunk_size = fixed_chunk_size > 0? fixed_chunk_size : opt->chunk_size * opt->n_threads;
	fprintf(stderr, "process called");
	kt_pipeline(no_mt_io? 1 : 2, process, &aux, 3);
	fprintf(stderr, "process finished");
	free(hdr_line);
	free(opt);
	bwa_idx_destroy(aux.idx);
	kseq_destroy(aux.ks);
	err_gzclose(fp); kclose(ko);
	if (aux.ks2) {
		kseq_destroy(aux.ks2);
		err_gzclose(fp2); kclose(ko2);
	}
	return 0;
}