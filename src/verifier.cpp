#include "verifier.hpp"

__UTILS_HASH_VALUE

inline void permute(__m128i *patterns) {
	__m128i tp[V_CPU];
	tp[0] = _mm_unpacklo_epi8(patterns[0], patterns[1]);
	tp[1] = _mm_unpackhi_epi8(patterns[0], patterns[1]);
	tp[2] = _mm_unpacklo_epi8(patterns[2], patterns[3]);
	tp[3] = _mm_unpackhi_epi8(patterns[2], patterns[3]);
	tp[4] = _mm_unpacklo_epi8(patterns[4], patterns[5]);
	tp[5] = _mm_unpackhi_epi8(patterns[4], patterns[5]);
	tp[6] = _mm_unpacklo_epi8(patterns[6], patterns[7]);
	tp[7] = _mm_unpackhi_epi8(patterns[6], patterns[7]);
	patterns[0] = _mm_unpacklo_epi16(tp[0], tp[2]);
	patterns[1] = _mm_unpacklo_epi16(tp[1], tp[3]);
	patterns[2] = _mm_unpackhi_epi16(tp[0], tp[2]);
	patterns[3] = _mm_unpackhi_epi16(tp[1], tp[3]);
	patterns[4] = _mm_unpacklo_epi16(tp[4], tp[6]);
	patterns[5] = _mm_unpacklo_epi16(tp[5], tp[7]);
	patterns[6] = _mm_unpackhi_epi16(tp[4], tp[6]);
	patterns[7] = _mm_unpackhi_epi16(tp[5], tp[7]);
	tp[0] = _mm_unpacklo_epi32(patterns[0], patterns[4]);
	tp[1] = _mm_unpackhi_epi32(patterns[0], patterns[4]);
	tp[2] = _mm_unpacklo_epi32(patterns[2], patterns[6]);
	tp[3] = _mm_unpackhi_epi32(patterns[2], patterns[6]);
	tp[4] = _mm_unpacklo_epi32(patterns[1], patterns[5]);
	tp[5] = _mm_unpackhi_epi32(patterns[1], patterns[5]);
	tp[6] = _mm_unpacklo_epi32(patterns[3], patterns[7]);
	tp[7] = _mm_unpackhi_epi32(patterns[3], patterns[7]);
	for (int pi = 0; pi < V_CPU; ++pi) {
		_mm_storeu_si128(patterns + pi, tp[pi]);
	}
}

void check_compact_mapped_read (uint8_t * my_bases, int my_rl, char * my_char_bases, uint32_t * my_compact_read) {
    int ri = 0, k = 0;//, l = 0, m = 0, n = 0;
    for (; ri < my_rl; ri++) {
        //l = k; m = k+1; n = k+2;
        //fprintf (stderr, "k3: %u, k2: %u, k1: %u.\n",
        //TestBit(my_compact_read, n), TestBit(my_compact_read, m), TestBit(my_compact_read, l));
        //fprintf (stderr, "ri : %d, k: %d, base: %c-%u.\n", ri , k, my_char_bases[ri], (uint32_t) my_bases[ri]);
        if (TestBit(my_compact_read,k)) {
            assert ((my_char_bases[ri] == 'C') || (my_char_bases[ri] == 'T'));
        } else {
            assert ((my_char_bases[ri] == 'A') || (my_char_bases[ri] == 'G') || (my_char_bases[ri] == 'N'));
        }
        k++;
        //fprintf (stderr, "ri : %d, k: %d, base: %c-%u.\n", ri , k, my_char_bases[ri], (uint32_t) my_bases[ri]);
        if (TestBit(my_compact_read,k)) {
            assert ((my_char_bases[ri] == 'G') || (my_char_bases[ri] == 'T'));
        } else {
            assert ((my_char_bases[ri] == 'A') || (my_char_bases[ri] == 'C') || (my_char_bases[ri] == 'N'));
        }
        k++;
        //fprintf (stderr, "ri : %d, k: %d, base: %c-%u.\n", ri , k, my_char_bases[ri], (uint32_t) my_bases[ri]);
        if (TestBit(my_compact_read,k)) {
            assert ((my_char_bases[ri] == 'N'));
        } else {
            assert ((my_char_bases[ri] == 'A') || (my_char_bases[ri] == 'C') || (my_char_bases[ri] == 'G') || (my_char_bases[ri] == 'T'));
        }
        k++;
    }
}

void check_compact_unmapped_read (uint8_t * my_bases, int my_rl, char * my_char_bases, uint32_t * my_compact_read) {
    int ri = 0, k = 0;//, l = 0, m = 0, n = 0;
    for (; ri < my_rl; ri++) {
        if (TestBit(my_compact_read,k)) {
            assert ((my_char_bases[ri] == 'C') || (my_char_bases[ri] == 'T'));
        } else {
            assert ((my_char_bases[ri] == 'A') || (my_char_bases[ri] == 'G') || (my_char_bases[ri] == 'N'));
        }
        k++;
        if (TestBit(my_compact_read,k)) {
            assert ((my_char_bases[ri] == 'G') || (my_char_bases[ri] == 'T'));
        } else {
            assert ((my_char_bases[ri] == 'A') || (my_char_bases[ri] == 'C') || (my_char_bases[ri] == 'N'));
        }
        k++;
        if (TestBit(my_compact_read,k)) {
            assert ((my_char_bases[ri] == 'N'));
        } else {
            assert ((my_char_bases[ri] == 'A') || (my_char_bases[ri] == 'C') || (my_char_bases[ri] == 'G') || (my_char_bases[ri] == 'T'));
        }
        k++;
    }
}

void compact_mapped_read (uint8_t * my_bases, int my_rl, uint32_t * my_compact_read, int my_best_error) {
            int ri = 0, k = 0;//, l = 0, m = 0, n = 0;
            for (; ri < my_rl; ri++) {
                //l = k; m = k+1; n = k+2;
                //if (ri < 4) {
                //    fprintf (stderr, "ri : %d, base: %u, k: %d.\n", ri, (uint32_t) my_bases[ri], k);
                //}
                if (my_bases[ri] == 0) {
                    k++;
                    k++;
                    k++;
                } else if (my_bases[ri] == 1) {
                    SetBit(my_compact_read, k);
                    k++;
                    k++;
                    k++;
                } else if (my_bases[ri] == 2) {
                    k++;
                    SetBit(my_compact_read, k);
                    k++;
                    k++;
                } else if (my_bases[ri] == 3) {
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
                //if (ri < 4) {
                //    fprintf (stderr, "k3: %u, k2: %u, k1: %u.\n",
                //    TestBit(my_compact_read, n), TestBit(my_compact_read, m), TestBit(my_compact_read, l));
                //    fprintf (stderr, "5: %u, 4: %u, 3: %u.\n",
                //    TestBit(my_compact_read, 5), TestBit(my_compact_read, 4), TestBit(my_compact_read, 3));
                //}
                //if (TestBit(my_compact_read, 5)) {
                //    fprintf (stderr, "[5] is being set. ri : %d, base: %u, k: %d.\n", ri, (uint32_t) my_bases[ri], k);
                //    ClearBit(my_compact_read, 5);
                //}
                //if (TestBit(my_compact_read, 11)) {
                //    fprintf (stderr, "[11] is being set. ri : %d, base: %u, k: %d.\n", ri, (uint32_t) my_bases[ri], k);
                //    ClearBit(my_compact_read, 11);
                //}
            }
            //l = 0; m = 1; n = 2;
            //fprintf (stderr, "k3: %u, k2: %u, k1: %u.\n",
            //TestBit(my_compact_read, n), TestBit(my_compact_read, m), TestBit(my_compact_read, l));
            //l += 3; m += 3; n += 3;
            //fprintf (stderr, "k3: %u, k2: %u, k1: %u.\n",
            //TestBit(my_compact_read, n), TestBit(my_compact_read, m), TestBit(my_compact_read, l));
            //l += 3; m += 3; n += 3;
            //fprintf (stderr, "k3: %u, k2: %u, k1: %u.\n",
            //TestBit(my_compact_read, n), TestBit(my_compact_read, m), TestBit(my_compact_read, l));
            //l += 3; m += 3; n += 3;
            //fprintf (stderr, "k3: %u, k2: %u, k1: %u.\n",
            //TestBit(my_compact_read, n), TestBit(my_compact_read, m), TestBit(my_compact_read, l));
            
#if !NDEBUG
            assert (my_best_error < 16);
            assert ((my_rl*3 + 4) <= (32 * READ_LEN_U32));
#endif
            k = 32 * READ_LEN_U32 - 4;
            int z = 1;
            
            if (my_best_error & z)
                SetBit(my_compact_read, k);
            z = (z << 1);
            k++; //-3
            if (my_best_error & z)
                SetBit(my_compact_read, k);
            z = (z << 1);
            k++; //-2
            if (my_best_error & z)
                SetBit(my_compact_read, k);
            z = (z << 1);
            k++; //-1
            if (my_best_error & z)
                SetBit(my_compact_read, k);
            z = (z << 1); //Not necessary
            k++; //Not necessary
}

std::size_t compact_unmapped_read (uint8_t * my_bases, int my_rl, uint32_t * my_compact_read) {//, int threadID) {
            int ri = 0, k = 0;
            //bool read_has_n = false;
            std::size_t my_n_count = 0;
            for (; ri < my_rl; ri++) {
/*
                if (my_bases[ri] == 0) {
                    k++;
                    k++;
                } else if (my_bases[ri] == 1) {
                    SetBit(my_compact_read, k);
                    k++;
                    k++;
                } else if (my_bases[ri] == 2) {
                    k++;
                    SetBit(my_compact_read, k);
                    k++;
                } else if (my_bases[ri] == 3) {
                    SetBit(my_compact_read, k);
                    k++;
                    SetBit(my_compact_read, k);
                    k++;
                } else {
                    // N base. Interpret as A
                    k++;
                    k++;
                    read_has_n = true;
                    unmapped_diffs_n_count[threadID] += 1;
                }
*/
                if (my_bases[ri] == 0) {
                    k++;
                    k++;
                    k++;
                } else if (my_bases[ri] == 1) {
                    SetBit(my_compact_read, k);
                    k++;
                    k++;
                    k++;
                } else if (my_bases[ri] == 2) {
                    k++;
                    SetBit(my_compact_read, k);
                    k++;
                    k++;
                } else if (my_bases[ri] == 3) {
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
                    //read_has_n = true;
                    //unmapped_diffs_n_count[threadID] += 1;
                    my_n_count++;
                }

            }
            //if (read_has_n) {
            //    unmapped_reads_n_count[threadID] += 1;
            //}
            return my_n_count;
}

int compare_uints(const void *a, const void *b) {
    uint32_t ta = *((uint32_t *) a);
    uint32_t tb = *((uint32_t *) b);
    if (ta < tb) {
        return -1;
    } else {
        return 1;
    }
}

int compare_init_candis(const void *a, const void *b) {
    TwoTuple * t1 = (TwoTuple *) a;
    TwoTuple * t2 = (TwoTuple *) b;
	//Sort in descending order of counts
	if (t1->b > t2->b) {
		return -1;
	} else if (t1->b == t2->b) {
		if (t1->a < t2->a)
			return -1;
		else
			return 1;
	} else
		return 1;
}

inline int NJ_generateCandidates (VerifyArgsForThread * vaftp, uint8_t * my_bases, int my_rl, uint32_t * my_places, 
uint32_t * my_candis, TwoTuple * my_init_candis, uint32_t * plaCnt, int my_edit_distance) {
    std::vector<std::uint32_t>& my_lookup_table = vaftp->vaft_com_ds->lookup_table;
    std::vector<std::uint32_t>& my_occurrence_table = vaftp->vaft_com_ds->occurrence_table;
    int hashVal;
    uint32_t occIdx;
    uint32_t occCnt;
    uint32_t my_plaCnt;
    uint32_t my_canCnt;
    
    my_plaCnt = 0;
    //int nCount1 = 0; //Not used
    //int qCount1 = 0; //Not used
    //TODO[LATER]: What if my_rl > 100? See assert below
    for (int ri = 5; ri < (my_rl - KMER_LENGTH); ri += KMER_LENGTH) {
        hashVal = hashValue(my_bases + ri, KMER_LENGTH);
        if (hashVal == -1) {
            //nCount1++;
        } else {
            occCnt = my_lookup_table[hashVal+1] - my_lookup_table[hashVal];
            for (uint32_t j=0; j<(CANDIDATES_COUNT/2); j++) {
                if (j == occCnt) {
                    break;
                }
                if (my_plaCnt == PLACES_COUNT) {
                    break;
                }
                occIdx = my_lookup_table[hashVal] + j;
                if (my_occurrence_table[occIdx] > ((uint32_t) ri)) {
                    my_places[my_plaCnt] = my_occurrence_table[occIdx] - ((uint32_t) ri);
                    my_plaCnt++;
                }
            }
            //qCount1++;
        }
        if (my_plaCnt == PLACES_COUNT) {
            break;
        }
    }
    
    //assert ((qCount1 + nCount1) == 6);
#if !NDEBUG
    assert (my_plaCnt <= PLACES_COUNT);
#endif
    *plaCnt = my_plaCnt;
    my_canCnt = 0;
    if (my_plaCnt > 1) {
        qsort(my_places, my_plaCnt, sizeof(uint32_t), compare_uints);
        uint32_t ta, tb;
        uint32_t currCan = my_places[0];
        int runCnt = 1;
        for (uint32_t j=1; j<my_plaCnt; j++) {
            //TODO[LATER]: Is below conditional block ncessary?
            if (currCan <= my_places[j]) {
                ta = currCan;
                tb = my_places[j];
            } else {
                ta = my_places[j];
                tb = currCan;
            }
            if ((tb - ta) <= ((uint32_t) my_edit_distance)) {//Tolerate up to my_edit_distance indels
                runCnt++;
            } else {
                if (runCnt > 1) {
                    my_init_candis[my_canCnt].a = currCan;
                    my_init_candis[my_canCnt].b = runCnt;
                    my_canCnt++;
                }
                currCan = my_places[j];
                runCnt = 1;
            }
            if (my_canCnt == (CANDIDATES_COUNT)) {//Maximum candidates found
                break;
            }
        }
        if (my_canCnt < (CANDIDATES_COUNT)) {
            if (runCnt > 1) {//Process last run
                my_init_candis[my_canCnt].a = currCan;
                my_init_candis[my_canCnt].b = runCnt;
                my_canCnt++;
            }
        }
    }
    
    qsort(my_init_candis, my_canCnt, sizeof(TwoTuple), compare_init_candis);
    for (uint32_t j=0; j<my_canCnt; j++) {
        my_candis[j] = my_init_candis[j].a;
    }
    
    return my_canCnt;
}

bool NJ_CPUVBM (VerifyArgsForThread * vaftp) {
    std::uint32_t threadID = vaftp->vaft_thread_id;
    std::uint8_t * my_ref_bases = vaftp->vaft_com_ds->ref_bases;
    
    Read * read1 = vaftp->read_1;
    Read * read2 = vaftp->read_2;
    MappedRead * m_read_fwd = vaftp->fwd_read;
    MappedRead * m_read_bwd = vaftp->bwd_read;
    UnmappedRead * u_read = vaftp->unm_read;
    std::uint32_t my_places[PLACES_COUNT];
    TwoTuple my_init_candis[CANDIDATES_COUNT];
    
    int differences = EDIT_DISTANCE;
    int strict_differences = EDIT_DISTANCE-6;
    int candis_distance = EDIT_DISTANCE-4;
    int r1_fwd_diffs = 1024, r1_bwd_diffs = 1024; //big number
    int r2_fwd_diffs = 1024, r2_bwd_diffs = 1024; //big number
    uint32_t canCnt, r1_canFwdCnt = 0, r1_canBwdCnt = 0, r2_canFwdCnt = 0, r2_canBwdCnt = 0;
    uint32_t r1_plaFwdCnt = 0, r1_plaBwdCnt = 0, r2_plaFwdCnt = 0, r2_plaBwdCnt = 0;
    bool r1_fwd = false, r1_map = false, r2_fwd = false, r2_map = false;
    uint32_t r1_fwd_candis[CANDIDATES_COUNT];
    uint32_t r1_bwd_candis[CANDIDATES_COUNT];
    int16_t r1_locs[V_CPU], r1_errs[V_CPU];
    uint32_t r2_fwd_candis[CANDIDATES_COUNT];
    uint32_t r2_bwd_candis[CANDIDATES_COUNT];
    int16_t r2_locs[V_CPU], r2_errs[V_CPU];
    uint32_t regNum;
    
    int cnt_r1_fle = 0, cnt_r1_ble = 0, cnt_r2_fle = 0, cnt_r2_ble = 0;
    int best_r1_fwd_err = 2048, best_r1_bwd_err = 2048, best_r2_fwd_err = 2048, best_r2_bwd_err = 2048;
    uint32_t r1_fwd_loc, r1_bwd_loc, r2_fwd_loc, r2_bwd_loc;
    uint32_t best_r1_fwd_loc = (2^32) - 2, best_r1_bwd_loc = (2^32) - 2;
    uint32_t best_r2_fwd_loc = (2^32) - 2, best_r2_bwd_loc = (2^32) - 2;
    uint32_t reg_r1_fle = 8096, reg_r1_ble = 8096, reg_r2_fle = 8096, reg_r2_ble = 8096;
    int r1_canFwdMax = 0, r1_canBwdMax = 0, r2_canFwdMax = 0, r2_canBwdMax = 0;
    r1_fwd_candis[0] = (2^32) - 2; r1_bwd_candis[0] = (2^32) - 2;
    r2_fwd_candis[0] = (2^32) - 2; r2_bwd_candis[0] = (2^32) - 2;
    uint32_t best_r1_fwd_pla = (2^32) - 2, best_r1_bwd_pla = (2^32) - 2;
    uint32_t best_r2_fwd_pla = (2^32) - 2, best_r2_bwd_pla = (2^32) - 2;
    int diff_r1_fwd_loc = 0, diff_r1_bwd_loc = 0, diff_r2_fwd_loc = 0, diff_r2_bwd_loc = 0;
    
    
    r1_canFwdCnt = NJ_generateCandidates(vaftp, read1->bases, read1->length, my_places, r1_fwd_candis,
    my_init_candis, &r1_plaFwdCnt, candis_distance);
//    fprintf (stderr, "canCnt read1->bases : %u.\n", canCnt);
    if (r1_plaFwdCnt) {best_r1_fwd_pla = my_places[0];}
    if (r1_canFwdCnt) {
        r1_canFwdMax = my_init_candis[0].b;
        regNum = r1_canFwdCnt / V_CPU;
        if (r1_canFwdCnt % V_CPU) { regNum += 1; }
        for (uint32_t j = 0; j < regNum; j++) {
            CPUVMBMKernel(vaftp, j, read1->bases, read1->length, r1_fwd_candis, r1_canFwdCnt, r1_errs, r1_locs);
            if ((j == (regNum - 1)) && (r1_canFwdCnt % V_CPU)) {
                canCnt = r1_canFwdCnt % V_CPU;
            } else {
                canCnt = V_CPU;
            }
            for (uint32_t mi = 0; mi < canCnt; mi++) {
                r1_fwd_diffs = (int) r1_errs[mi];
                if (r1_fwd_diffs <= differences) {
                    cnt_r1_fle++;
                    r1_map = true;
                    if (r1_fwd_diffs <= strict_differences) {
                        r1_fwd = true;
                    }
                    
                    r1_fwd_loc = r1_fwd_candis[j * V_CPU + mi];
                    if (r1_fwd_diffs < best_r1_fwd_err) {
                        best_r1_fwd_err = r1_fwd_diffs;
                        best_r1_fwd_loc = r1_fwd_loc;
                        diff_r1_fwd_loc = (int) r1_locs[mi];
                    }
                    if (r1_fwd_diffs == best_r1_fwd_err) {
                        if (r1_fwd_loc < best_r1_fwd_loc) {
                            best_r1_fwd_loc = r1_fwd_loc;
                            diff_r1_fwd_loc = (int) r1_locs[mi];
                        }
                    }
                    if (reg_r1_fle == 8096) {
                        reg_r1_fle = j;
                    }
                    if (cnt_r1_fle == SAVE_LOCS_COUNT) { break; }
                }
            }
            if (cnt_r1_fle == SAVE_LOCS_COUNT) { break; }
            if (j == (reg_r1_fle+3)) { break; }
        }
    }
    
    if (!r1_fwd) {//r1_fwd mapping failed
        r1_canBwdCnt = NJ_generateCandidates(vaftp, read1->rc_bases, read1->length, my_places, r1_bwd_candis,
        my_init_candis, &r1_plaBwdCnt, candis_distance);
//        fprintf (stderr, "canCnt read1->rc_bases : %u.\n", canCnt);
        if (r1_plaBwdCnt) {best_r1_bwd_pla = my_places[0];}
        if (r1_canBwdCnt) {
            r1_canBwdMax = my_init_candis[0].b;
            regNum = r1_canBwdCnt / V_CPU;
            if (r1_canBwdCnt % V_CPU) { regNum += 1; }
            for (uint32_t j = 0; j < regNum; j++) {
                CPUVMBMKernel(vaftp, j, read1->rc_bases, read1->length, r1_bwd_candis, r1_canBwdCnt, r1_errs, r1_locs);
                if ((j == (regNum - 1)) && (r1_canBwdCnt % V_CPU)) {
                    canCnt = r1_canBwdCnt % V_CPU;
                } else {
                    canCnt = V_CPU;
                }
                for (uint32_t mi = 0; mi < canCnt; mi++) {
                    r1_bwd_diffs = (int) r1_errs[mi];
                    if (r1_bwd_diffs <= differences) {
                        cnt_r1_ble++;
                        r1_map = true;
                        
                        r1_bwd_loc = r1_bwd_candis[j * V_CPU + mi];
                        if (r1_bwd_diffs < best_r1_bwd_err) {
                            best_r1_bwd_err = r1_bwd_diffs;
                            best_r1_bwd_loc = r1_bwd_loc;
                            diff_r1_bwd_loc = (int) r1_locs[mi];
                        }
                        if (r1_bwd_diffs == best_r1_bwd_err) {
                            if (r1_bwd_loc < best_r1_bwd_loc) {
                                best_r1_bwd_loc = r1_bwd_loc;
                                diff_r1_bwd_loc = (int) r1_locs[mi];
                            }
                        }
                        if (reg_r1_ble == 8096) {
                            reg_r1_ble = j;
                        }
                        if (cnt_r1_ble == SAVE_LOCS_COUNT) { break; }
                    }
                }
                if (cnt_r1_ble == SAVE_LOCS_COUNT) { break; }
                if (j == (reg_r1_ble+3)) { break; }
            }
        }
    }
    
    if (r1_map) {
        if (best_r1_fwd_err <= best_r1_bwd_err) {
            r1_fwd = true;
            (vaftp->asp->as_differences_count[threadID*16 + best_r1_fwd_err])++;
        } else {
            r1_fwd = false;
            (vaftp->asp->as_differences_count[threadID*16 + best_r1_bwd_err])++;
        }
        
        if (r1_fwd) {
            r2_canBwdCnt = NJ_generateCandidates(vaftp, read2->rc_bases, read2->length, my_places, r2_bwd_candis,
            my_init_candis, &r2_plaBwdCnt, candis_distance);
            if (r2_plaBwdCnt) {best_r2_bwd_pla = my_places[0];}
            if (r2_canBwdCnt) {
                r2_canBwdMax = my_init_candis[0].b;
                regNum = r2_canBwdCnt / V_CPU;
                if (r2_canBwdCnt % V_CPU) { regNum += 1; }
                for (uint32_t j = 0; j < regNum; j++) {
                    CPUVMBMKernel(vaftp, j, read2->rc_bases, read2->length, r2_bwd_candis, r2_canBwdCnt, r2_errs, r2_locs);
                    if ((j == (regNum - 1)) && (r2_canBwdCnt % V_CPU)) {
                        canCnt = r2_canBwdCnt % V_CPU;
                    } else {
                        canCnt = V_CPU;
                    }
                    for (uint32_t mi = 0; mi < canCnt; mi++) {
                        r2_bwd_diffs = (int) r2_errs[mi];
                        if (r2_bwd_diffs <= differences) {
                            cnt_r2_ble++;
                            r2_map = true;
                            r2_bwd_loc = r2_bwd_candis[j * V_CPU + mi];
                            if (r2_bwd_diffs < best_r2_bwd_err) {
                                best_r2_bwd_err = r2_bwd_diffs;
                                best_r2_bwd_loc = r2_bwd_loc;
                                diff_r2_bwd_loc = (int) r2_locs[mi];
                            }
                            if (r2_bwd_diffs == best_r2_bwd_err) {
                                if (r2_bwd_loc <= best_r1_fwd_loc) {
                                    if (best_r2_bwd_loc <= best_r1_fwd_loc) {
                                        if ((best_r1_fwd_loc - r2_bwd_loc) < (best_r1_fwd_loc - best_r2_bwd_loc)) {
                                            best_r2_bwd_loc = r2_bwd_loc;
                                            diff_r2_bwd_loc = (int) r2_locs[mi];
                                        }
                                    } else {
                                        if ((best_r1_fwd_loc - r2_bwd_loc) < (best_r2_bwd_loc - best_r1_fwd_loc)) {
                                            best_r2_bwd_loc = r2_bwd_loc;
                                            diff_r2_bwd_loc = (int) r2_locs[mi];
                                        }
                                    }
                                } else {
                                    if (best_r2_bwd_loc <= best_r1_fwd_loc) {
                                        if ((r2_bwd_loc - best_r1_fwd_loc) < (best_r1_fwd_loc - best_r2_bwd_loc)) {
                                            best_r2_bwd_loc = r2_bwd_loc;
                                            diff_r2_bwd_loc = (int) r2_locs[mi];
                                        }
                                    } else {
                                        if ((r2_bwd_loc - best_r1_fwd_loc) < (best_r2_bwd_loc - best_r1_fwd_loc)) {
                                            best_r2_bwd_loc = r2_bwd_loc;
                                            diff_r2_bwd_loc = (int) r2_locs[mi];
                                        }
                                    }
                                }
                            }
                            if (reg_r2_ble == 8096) {
                                reg_r2_ble = j;
                            }
                            if (cnt_r2_ble == SAVE_LOCS_COUNT) { break; }
                        }
                    }
                    if (cnt_r2_ble == SAVE_LOCS_COUNT) { break; }
                    if (j == (reg_r2_ble+3)) { break; }
                }
            }
        } else {//r2_bwd mapping failed
            r2_canFwdCnt = NJ_generateCandidates(vaftp, read2->bases, read2->length, my_places, r2_fwd_candis,
            my_init_candis, &r2_plaFwdCnt, candis_distance);
            if (r2_plaFwdCnt) {best_r2_fwd_pla = my_places[0];}
            if (r2_canFwdCnt) {
                r2_canFwdMax = my_init_candis[0].b;
                regNum = r2_canFwdCnt / V_CPU;
                if (r2_canFwdCnt % V_CPU) { regNum += 1; }
                for (uint32_t j = 0; j < regNum; j++) {
                    CPUVMBMKernel(vaftp, j, read2->bases, read2->length, r2_fwd_candis, r2_canFwdCnt, r2_errs, r2_locs);
                    if ((j == (regNum - 1)) && (r2_canFwdCnt % V_CPU)) {
                        canCnt = r2_canFwdCnt % V_CPU;
                    } else {
                        canCnt = V_CPU;
                    }
                    for (uint32_t mi = 0; mi < canCnt; mi++) {
                        r2_fwd_diffs = (int) r2_errs[mi];
                        if (r2_fwd_diffs <= differences) {
                            cnt_r2_fle++;
                            r2_map = true;
                            r2_fwd_loc = r2_fwd_candis[j * V_CPU + mi];
                            if (r2_fwd_diffs < best_r2_fwd_err) {
                                best_r2_fwd_err = r2_fwd_diffs;
                                best_r2_fwd_loc = r2_fwd_loc;
                                diff_r2_fwd_loc = (int) r2_locs[mi];
                            }
                            if (r2_fwd_diffs == best_r2_fwd_err) {
                                if (r2_fwd_loc <= best_r1_bwd_loc) {
                                    if (best_r2_fwd_loc <= best_r1_bwd_loc) {
                                        if ((best_r1_bwd_loc - r2_fwd_loc) < (best_r1_bwd_loc - best_r2_fwd_loc)) {
                                            best_r2_fwd_loc = r2_fwd_loc;
                                            diff_r2_fwd_loc = (int) r2_locs[mi];
                                        }
                                    } else {
                                        if ((best_r1_bwd_loc - r2_fwd_loc) < (best_r2_fwd_loc - best_r1_bwd_loc)) {
                                            best_r2_fwd_loc = r2_fwd_loc;
                                            diff_r2_fwd_loc = (int) r2_locs[mi];
                                        }
                                    }
                                } else {
                                    if (best_r2_fwd_loc <= best_r1_bwd_loc) {
                                        if ((r2_fwd_loc - best_r1_bwd_loc) < (best_r1_bwd_loc - best_r2_fwd_loc)) {
                                            best_r2_fwd_loc = r2_fwd_loc;
                                            diff_r2_fwd_loc = (int) r2_locs[mi];
                                        }
                                    } else {
                                        if ((r2_fwd_loc - best_r1_bwd_loc) < (best_r2_fwd_loc - best_r1_bwd_loc)) {
                                            best_r2_fwd_loc = r2_fwd_loc;
                                            diff_r2_fwd_loc = (int) r2_locs[mi];
                                        }
                                    }
                                }
                            }
                            if (reg_r2_fle == 8096) {
                                reg_r2_fle = j;
                            }
                            if (cnt_r2_fle == SAVE_LOCS_COUNT) { break; }
                        }
                    }
                    if (cnt_r2_fle == SAVE_LOCS_COUNT) { break; }
                    if (j == (reg_r2_fle+3)) { break; }
                }
            }
        }
        
        if (r2_map) {
            if (r1_fwd) {
                r2_fwd = false;
                (vaftp->asp->as_differences_count[threadID*16 + best_r2_bwd_err])++;
            } else {
                r2_fwd = true;
                (vaftp->asp->as_differences_count[threadID*16 + best_r2_fwd_err])++;
            }
        }
    } else {
        //r1 (both fwd and bwd) failed to map
        //Generating r2 mapping is not useful for reference-based compression
        //but is helpful for positioning/ordering and orienting the PE read
        r2_canFwdCnt = NJ_generateCandidates(vaftp, read2->bases, read2->length, my_places, r2_fwd_candis,
        my_init_candis, &r2_plaFwdCnt, candis_distance);
        if (r2_plaFwdCnt) {best_r2_fwd_pla = my_places[0];}
        if (r2_canFwdCnt) {
            r2_canFwdMax = my_init_candis[0].b;
            regNum = r2_canFwdCnt / V_CPU;
            if (r2_canFwdCnt % V_CPU) { regNum += 1; }
            for (uint32_t j = 0; j < regNum; j++) {
                CPUVMBMKernel(vaftp, j, read2->bases, read2->length, r2_fwd_candis, r2_canFwdCnt, r2_errs, r2_locs);
                if ((j == (regNum - 1)) && (r2_canFwdCnt % V_CPU)) {
                    canCnt = r2_canFwdCnt % V_CPU;
                } else {
                    canCnt = V_CPU;
                }
                for (uint32_t mi = 0; mi < canCnt; mi++) {
                    r2_fwd_diffs = (int) r2_errs[mi];
                    if (r2_fwd_diffs <= differences) {
                        cnt_r2_fle++;
                        r2_map = true;
                        if (r2_fwd_diffs <= strict_differences) {
                            r2_fwd = true;
                        }
                        r2_fwd_loc = r2_fwd_candis[j * V_CPU + mi];
                        if (r2_fwd_diffs < best_r2_fwd_err) {
                            best_r2_fwd_err = r2_fwd_diffs;
                            best_r2_fwd_loc = r2_fwd_loc;
                            diff_r2_fwd_loc = (int) r2_locs[mi];
                        }
                        if (r2_fwd_diffs == best_r2_fwd_err) {
                            if (r2_fwd_loc < best_r2_fwd_loc) {
                                best_r2_fwd_loc = r2_fwd_loc;
                                diff_r2_fwd_loc = (int) r2_locs[mi];
                            }
                        }
                        if (reg_r2_fle == 8096) {
                            reg_r2_fle = j;
                        }
                        if (cnt_r2_fle == SAVE_LOCS_COUNT) { break; }
                    }
                }
                if (cnt_r2_fle == SAVE_LOCS_COUNT) { break; }
                if (j == (reg_r2_fle+3)) { break; }
            }
        }
        
        if (!r2_fwd) {//r2_fwd mapping failed
            r2_canBwdCnt = NJ_generateCandidates(vaftp, read2->rc_bases, read2->length, my_places, r2_bwd_candis,
            my_init_candis, &r2_plaBwdCnt, candis_distance);
            if (r2_plaBwdCnt) {best_r2_bwd_pla = my_places[0];}
            if (r2_canBwdCnt) {
                r2_canBwdMax = my_init_candis[0].b;
                regNum = r2_canBwdCnt / V_CPU;
                if (r2_canBwdCnt % V_CPU) { regNum += 1; }
                for (uint32_t j = 0; j < regNum; j++) {
                    CPUVMBMKernel(vaftp, j, read2->rc_bases, read2->length, r2_bwd_candis, r2_canBwdCnt, r2_errs, r2_locs);
                    if ((j == (regNum - 1)) && (r2_canBwdCnt % V_CPU)) {
                        canCnt = r2_canBwdCnt % V_CPU;
                    } else {
                        canCnt = V_CPU;
                    }
                    for (uint32_t mi = 0; mi < canCnt; mi++) {
                        r2_bwd_diffs = (int) r2_errs[mi];
                        if (r2_bwd_diffs <= differences) {
                            cnt_r2_ble++;
                            r2_map = true;
                            r2_bwd_loc = r2_bwd_candis[j * V_CPU + mi];
                            if (r2_bwd_diffs < best_r2_bwd_err) {
                                best_r2_bwd_err = r2_bwd_diffs;
                                best_r2_bwd_loc = r2_bwd_loc;
                                diff_r2_bwd_loc = (int) r2_locs[mi];
                            }
                            if (r2_bwd_diffs == best_r2_bwd_err) {
                                if (r2_bwd_loc < best_r2_bwd_loc) {
                                    best_r2_bwd_loc = r2_bwd_loc;
                                    diff_r2_bwd_loc = (int) r2_locs[mi];
                                }
                            }
                            if (reg_r2_ble == 8096) {
                                reg_r2_ble = j;
                            }
                            if (cnt_r2_ble == SAVE_LOCS_COUNT) { break; }
                        }
                    }
                    if (cnt_r2_ble == SAVE_LOCS_COUNT) { break; }
                    if (j == (reg_r2_ble+3)) { break; }
                }
            }
        }
        
        if (r2_map) {
            if (best_r2_fwd_err <= best_r2_bwd_err) {
                r2_fwd = true;
                (vaftp->asp->as_differences_count[threadID*16 + best_r2_fwd_err])++;
            } else {
                r2_fwd = false;
                (vaftp->asp->as_differences_count[threadID*16 + best_r2_bwd_err])++;
            }
        }
    }
    
    vaftp->asp->as_paired_read_count[threadID] += 1;
    if (r1_map && r2_map) {
        //r1 and r2 both map
        vaftp->asp->as_paired_mapped_count[threadID] += 1;
        
        for (uint32_t j = 0; j < READ_LEN_U32; j++) {
            m_read_fwd->read[j] = 0;
            m_read_bwd->read[j] = 0;
        }
        char cigar[READ_LEN_MAX*2];
        char copy_cigar[READ_LEN_MAX*2];
#if !NDEBUG
        int copy_diff_location;
        assert (r1_fwd != r2_fwd);
#endif
        //read length is not accounted for below
        if (r1_fwd && (!r2_fwd)) {
            compact_mapped_read(read1->bases, read1->length, m_read_fwd->read, best_r1_fwd_err);
            if (best_r1_fwd_err) {
#if !NDEBUG
                if (diff_r1_fwd_loc < (read1->length - 1)) {
                    vaftp->asp->as_less_than_length_count[threadID] += 1;
                } else if (diff_r1_fwd_loc == (read1->length - 1)) {
                    vaftp->asp->as_equal_to_length_count[threadID] += 1;
                } else {
                    vaftp->asp->as_more_than_length_count[threadID] += 1;
                }
                if ((diff_r1_fwd_loc - read1->length + 1) < 0) {
                    vaftp->asp->as_neg_diff_count_bef[threadID] += 1;
                }
#endif
                m_read_fwd->base_location = best_r1_fwd_loc;
                m_read_fwd->diff_location = diff_r1_fwd_loc;
                diff_r1_fwd_loc = 
                CPUAlignKernel((my_ref_bases + best_r1_fwd_loc), read1->bases, diff_r1_fwd_loc, cigar, read1->length, best_r1_fwd_err);
#if !NDEBUG
                if (diff_r1_fwd_loc < 0) {
                    vaftp->asp->as_neg_diff_count_aft[threadID] += 1;
                    int z1 = -diff_r1_fwd_loc;
                    assert (best_r1_fwd_loc > ((uint32_t) z1));
                }
#endif
                m_read_fwd->location = best_r1_fwd_loc + (uint32_t) diff_r1_fwd_loc;
#if !NDEBUG
                copy_diff_location = m_read_fwd->diff_location;
                copy_diff_location = 
                CPUAlignKernel((my_ref_bases + best_r1_fwd_loc), read1->bases, copy_diff_location, copy_cigar, read1->length, best_r1_fwd_err);
                if (strcmp(copy_cigar, cigar)) {
                    vaftp->asp->as_cigar_unequal_count[threadID] += 1;
                    fprintf (stderr, "base : %u, length : %d, copy_cigar : %s, cigar : %s.\n",
                    m_read_fwd->base_location, read1->length, copy_cigar, cigar);
                    fprintf (stderr, "Read string : %s.\n", read1->charBases);
                }
                if (copy_diff_location != diff_r1_fwd_loc) {
                    vaftp->asp->as_location_unequal_count[threadID] += 1;
                    fprintf (stderr, "base : %u, length : %d, diff : %d, diff_r1_fwd_loc : %d, copy_diff : %d.\n",
                    m_read_fwd->base_location, read1->length, m_read_fwd->diff_location, diff_r1_fwd_loc, copy_diff_location);
                    fprintf (stderr, "Read string : %s.\n", read1->charBases);
                }
#endif
                if (strchr(cigar, 'S')) {
                    vaftp->asp->as_soft_clipped_count[threadID] += 1;
                }
                if (strlen(cigar) < 16) {
                    strcpy(m_read_fwd->tmp_cigar, cigar);
                } else {
                    strcpy(m_read_fwd->tmp_cigar, "Y");
                }
            } else {
                //zero errors
                m_read_fwd->location = best_r1_fwd_loc;
                strcpy(m_read_fwd->tmp_cigar, "Z");
            }
            compact_mapped_read(read2->rc_bases, read2->length, m_read_bwd->read, best_r2_bwd_err);
            if (best_r2_bwd_err) {
#if !NDEBUG
                if (diff_r2_bwd_loc < (read2->length - 1)) {
                    vaftp->asp->as_less_than_length_count[threadID] += 1;
                } else if (diff_r2_bwd_loc == (read2->length - 1)) {
                    vaftp->asp->as_equal_to_length_count[threadID] += 1;
                } else {
                    vaftp->asp->as_more_than_length_count[threadID] += 1;
                }
                if ((diff_r2_bwd_loc - read2->length + 1) < 0) {
                    vaftp->asp->as_neg_diff_count_bef[threadID] += 1;
                }
#endif
                m_read_bwd->base_location = best_r2_bwd_loc;
                m_read_bwd->diff_location = diff_r2_bwd_loc;
                diff_r2_bwd_loc = 
                CPUAlignKernel((my_ref_bases + best_r2_bwd_loc), read2->rc_bases, diff_r2_bwd_loc, cigar, read2->length, best_r2_bwd_err);
#if !NDEBUG
                if (diff_r2_bwd_loc < 0) {
                    vaftp->asp->as_neg_diff_count_aft[threadID] += 1;
                    int z1 = -diff_r2_bwd_loc;
                    assert (best_r2_bwd_loc > ((uint32_t) z1));
                }
#endif
                m_read_bwd->location = best_r2_bwd_loc + (uint32_t) diff_r2_bwd_loc;
#if !NDEBUG
                copy_diff_location = m_read_bwd->diff_location;
                copy_diff_location = 
                CPUAlignKernel((my_ref_bases + best_r2_bwd_loc), read2->rc_bases, copy_diff_location, copy_cigar, read2->length, best_r2_bwd_err);
                if (strcmp(copy_cigar, cigar)) {
                    vaftp->asp->as_cigar_unequal_count[threadID] += 1;
                    fprintf (stderr, "base : %u, length : %d, copy_cigar : %s, cigar : %s.\n",
                    m_read_bwd->base_location, read2->length, copy_cigar, cigar);
                }
                if (copy_diff_location != diff_r2_bwd_loc) {
                    vaftp->asp->as_location_unequal_count[threadID] += 1;
                    fprintf (stderr, "base : %u, length : %d, diff : %d, diff_r2_bwd_loc : %d, copy_diff : %d.\n",
                    m_read_bwd->base_location, read2->length, m_read_bwd->diff_location, diff_r2_bwd_loc, copy_diff_location);
                }
#endif
                if (strchr(cigar, 'S')) {
                    vaftp->asp->as_soft_clipped_count[threadID] += 1;
                }
                if (strlen(cigar) < 16) {
                    strcpy(m_read_bwd->tmp_cigar, cigar);
                } else {
                    strcpy(m_read_bwd->tmp_cigar, "Y");
                }
            } else {
                //zero errors
                m_read_bwd->location = best_r2_bwd_loc;
                strcpy(m_read_bwd->tmp_cigar, "Z");
            }
            
#if !NDEBUG
		    check_compact_mapped_read(read1->bases, read1->length, read1->charBases, m_read_fwd->read);
            
            if (best_r1_fwd_loc < best_r2_bwd_loc) {
                vaftp->asp->as_paired_read_order[threadID*8 + 0] += 1;
            } else if (best_r1_fwd_loc == best_r2_bwd_loc) {
                vaftp->asp->as_paired_read_order[threadID*8 + 1] += 1;
            } else {
                vaftp->asp->as_paired_read_order[threadID*8 + 2] += 1;
            }
#endif
        } else {
            compact_mapped_read(read2->bases, read2->length, m_read_fwd->read, best_r2_fwd_err);
            if (best_r2_fwd_err) {
#if !NDEBUG
                if (diff_r2_fwd_loc < (read2->length - 1)) {
                    vaftp->asp->as_less_than_length_count[threadID] += 1;
                } else if (diff_r2_fwd_loc == (read2->length - 1)) {
                    vaftp->asp->as_equal_to_length_count[threadID] += 1;
                } else {
                    vaftp->asp->as_more_than_length_count[threadID] += 1;
                }
                if ((diff_r2_fwd_loc - read2->length + 1) < 0) {
                    vaftp->asp->as_neg_diff_count_bef[threadID] += 1;
                }
#endif
                m_read_fwd->base_location = best_r2_fwd_loc;
                m_read_fwd->diff_location = diff_r2_fwd_loc;
                diff_r2_fwd_loc = 
                CPUAlignKernel((my_ref_bases + best_r2_fwd_loc), read2->bases, diff_r2_fwd_loc, cigar, read2->length, best_r2_fwd_err);
#if !NDEBUG
                if (diff_r2_fwd_loc < 0) {
                    vaftp->asp->as_neg_diff_count_aft[threadID] += 1;
                    int z1 = -diff_r2_fwd_loc;
                    assert (best_r2_fwd_loc > ((uint32_t) z1));
                }
#endif
                m_read_fwd->location = best_r2_fwd_loc + (uint32_t) diff_r2_fwd_loc;
#if !NDEBUG
                copy_diff_location = m_read_fwd->diff_location;
                copy_diff_location = 
                CPUAlignKernel((my_ref_bases + best_r2_fwd_loc), read2->bases, copy_diff_location, copy_cigar, read2->length, best_r2_fwd_err);
                if (strcmp(copy_cigar, cigar)) {
                    vaftp->asp->as_cigar_unequal_count[threadID] += 1;
                    fprintf (stderr, "base : %u, length : %d, copy_cigar : %s, cigar : %s.\n",
                    m_read_fwd->base_location, read2->length, copy_cigar, cigar);
                    fprintf (stderr, "Read string : %s.\n", read2->charBases);
                }
                if (copy_diff_location != diff_r2_fwd_loc) {
                    vaftp->asp->as_location_unequal_count[threadID] += 1;
                    fprintf (stderr, "base : %u, length : %d, diff : %d, diff_r2_fwd_loc : %d, copy_diff : %d.\n",
                    m_read_fwd->base_location, read2->length, m_read_fwd->diff_location, diff_r2_fwd_loc, copy_diff_location);
                    fprintf (stderr, "Read string : %s.\n", read2->charBases);
                }
#endif
                if (strchr(cigar, 'S')) {
                    vaftp->asp->as_soft_clipped_count[threadID] += 1;
                }
                if (strlen(cigar) < 16) {
                    strcpy(m_read_fwd->tmp_cigar, cigar);
                } else {
                    strcpy(m_read_fwd->tmp_cigar, "Y");
                }
            } else {
                //zero errors
                m_read_fwd->location = best_r2_fwd_loc;
                strcpy(m_read_fwd->tmp_cigar, "Z");
            }
            compact_mapped_read(read1->rc_bases, read1->length, m_read_bwd->read, best_r1_bwd_err);
            if (best_r1_bwd_err) {
#if !NDEBUG
                if (diff_r1_bwd_loc < (read1->length - 1)) {
                    vaftp->asp->as_less_than_length_count[threadID] += 1;
                } else if (diff_r1_bwd_loc == (read1->length - 1)) {
                    vaftp->asp->as_equal_to_length_count[threadID] += 1;
                } else {
                    vaftp->asp->as_more_than_length_count[threadID] += 1;
                }
                if ((diff_r1_bwd_loc - read1->length + 1) < 0) {
                    vaftp->asp->as_neg_diff_count_bef[threadID] += 1;
                }
#endif
                m_read_bwd->base_location = best_r1_bwd_loc;
                m_read_bwd->diff_location = diff_r1_bwd_loc;
                diff_r1_bwd_loc = 
                CPUAlignKernel((my_ref_bases + best_r1_bwd_loc), read1->rc_bases, diff_r1_bwd_loc, cigar, read1->length, best_r1_bwd_err);
#if !NDEBUG
                if (diff_r1_bwd_loc < 0) {
                    vaftp->asp->as_neg_diff_count_aft[threadID] += 1;
                    int z1 = -diff_r1_bwd_loc;
                    assert (best_r1_bwd_loc > ((uint32_t) z1));
                }
#endif
                m_read_bwd->location = best_r1_bwd_loc + (uint32_t) diff_r1_bwd_loc;
#if !NDEBUG
                copy_diff_location = m_read_bwd->diff_location;
                copy_diff_location = 
                CPUAlignKernel((my_ref_bases + best_r1_bwd_loc), read1->rc_bases, copy_diff_location, copy_cigar, read1->length, best_r1_bwd_err);
                if (strcmp(copy_cigar, cigar)) {
                    vaftp->asp->as_cigar_unequal_count[threadID] += 1;
                    fprintf (stderr, "base : %u, length : %d, copy_cigar : %s, cigar : %s.\n",
                    m_read_bwd->base_location, read1->length, copy_cigar, cigar);
                }
                if (copy_diff_location != diff_r1_bwd_loc) {
                    vaftp->asp->as_location_unequal_count[threadID] += 1;
                    fprintf (stderr, "base : %u, length : %d, diff : %d, diff_r1_bwd_loc : %d, copy_diff : %d.\n",
                    m_read_bwd->base_location, read1->length, m_read_bwd->diff_location, diff_r1_bwd_loc, copy_diff_location);
                    
                }
#endif
                if (strchr(cigar, 'S')) {
                    vaftp->asp->as_soft_clipped_count[threadID] += 1;
                }
                if (strlen(cigar) < 16) {
                    strcpy(m_read_bwd->tmp_cigar, cigar);
                } else {
                    strcpy(m_read_bwd->tmp_cigar, "Y");
                }
            } else {
                //zero errors
                m_read_bwd->location = best_r1_bwd_loc;
                strcpy(m_read_bwd->tmp_cigar, "Z");
            }
            
#if !NDEBUG
		    check_compact_mapped_read(read2->bases, read2->length, read2->charBases, m_read_fwd->read);
            
            if (best_r2_fwd_loc < best_r1_bwd_loc) {
                vaftp->asp->as_paired_read_order[threadID*8 + 4] += 1;
            } else if (best_r2_fwd_loc == best_r1_bwd_loc) {
                vaftp->asp->as_paired_read_order[threadID*8 + 5] += 1;
            } else {
                vaftp->asp->as_paired_read_order[threadID*8 + 6] += 1;
            }
#endif
        }
    } else if (r1_map && !r2_map) {
        //r1 maps but r2 does not
        vaftp->asp->as_single_mapped_count[threadID] += 1;
        std::size_t my_unm_n_count;
        
        for (uint32_t j = 0; j < READ_LEN_U32; j++) {
            u_read->fwd_read[j] = 0;
            u_read->bwd_read[j] = 0;
        }
        if (r1_fwd) {
            my_unm_n_count = compact_unmapped_read (read1->bases, read1->length, u_read->fwd_read);//, threadID);
            if (my_unm_n_count) {
                vaftp->asp->as_unmapped_reads_n_count[threadID] += 1;
                vaftp->asp->as_unmapped_diffs_n_count[threadID] += my_unm_n_count;
            }
            my_unm_n_count = compact_unmapped_read (read2->rc_bases, read2->length, u_read->bwd_read);//, threadID);
            if (my_unm_n_count) {
                vaftp->asp->as_unmapped_reads_n_count[threadID] += 1;
                vaftp->asp->as_unmapped_diffs_n_count[threadID] += my_unm_n_count;
            }
            u_read->location = best_r1_fwd_loc;
#if !NDEBUG
            check_compact_unmapped_read(read1->bases, read1->length, read1->charBases, u_read->fwd_read);
#endif
        } else {
            my_unm_n_count = compact_unmapped_read (read1->rc_bases, read1->length, u_read->bwd_read);//, threadID);
            if (my_unm_n_count) {
                vaftp->asp->as_unmapped_reads_n_count[threadID] += 1;
                vaftp->asp->as_unmapped_diffs_n_count[threadID] += my_unm_n_count;
            }
            my_unm_n_count = compact_unmapped_read (read2->bases, read2->length, u_read->fwd_read);//, threadID);
            if (my_unm_n_count) {
                vaftp->asp->as_unmapped_reads_n_count[threadID] += 1;
                vaftp->asp->as_unmapped_diffs_n_count[threadID] += my_unm_n_count;
            }
            u_read->location = best_r1_bwd_loc;
#if !NDEBUG
            check_compact_unmapped_read(read2->bases, read2->length, read2->charBases, u_read->fwd_read);
#endif
        }
    } else if (!r1_map && r2_map) {
        //r1 does not map but r2 does
        vaftp->asp->as_single_mapped_count[threadID] += 1;
        std::size_t my_unm_n_count;
        
        for (uint32_t j = 0; j < READ_LEN_U32; j++) {
            u_read->fwd_read[j] = 0;
            u_read->bwd_read[j] = 0;
        }
        if (r2_fwd) {
            my_unm_n_count = compact_unmapped_read (read2->bases, read2->length, u_read->fwd_read);//, threadID);
            if (my_unm_n_count) {
                vaftp->asp->as_unmapped_reads_n_count[threadID] += 1;
                vaftp->asp->as_unmapped_diffs_n_count[threadID] += my_unm_n_count;
            }
            my_unm_n_count = compact_unmapped_read (read1->rc_bases, read1->length, u_read->bwd_read);//, threadID);
            if (my_unm_n_count) {
                vaftp->asp->as_unmapped_reads_n_count[threadID] += 1;
                vaftp->asp->as_unmapped_diffs_n_count[threadID] += my_unm_n_count;
            }
            u_read->location = best_r2_fwd_loc;
#if !NDEBUG
            check_compact_unmapped_read(read2->bases, read2->length, read2->charBases, u_read->fwd_read);
#endif
        } else {
            my_unm_n_count = compact_unmapped_read (read2->rc_bases, read2->length, u_read->bwd_read);//, threadID);
            if (my_unm_n_count) {
                vaftp->asp->as_unmapped_reads_n_count[threadID] += 1;
                vaftp->asp->as_unmapped_diffs_n_count[threadID] += my_unm_n_count;
            }
            my_unm_n_count = compact_unmapped_read (read1->bases, read1->length, u_read->fwd_read);//, threadID);
            if (my_unm_n_count) {
                vaftp->asp->as_unmapped_reads_n_count[threadID] += 1;
                vaftp->asp->as_unmapped_diffs_n_count[threadID] += my_unm_n_count;
            }
            u_read->location = best_r2_bwd_loc;
#if !NDEBUG
            check_compact_unmapped_read(read1->bases, read1->length, read1->charBases, u_read->fwd_read);
#endif
        }
    } else {
        //r1 and r2 both do not map
        vaftp->asp->as_both_unmapped_count[threadID] += 1;
        std::size_t my_unm_n_count;
        
        for (uint32_t j = 0; j < READ_LEN_U32; j++) {
            u_read->fwd_read[j] = 0;
            u_read->bwd_read[j] = 0;
        }
        if (r1_canFwdCnt || r1_canBwdCnt || r2_canFwdCnt || r2_canBwdCnt) {
            vaftp->asp->as_can_unmapped_count[threadID] += 1;
            
            r1_fwd = true;
            r2_fwd = false;
            if ((r1_canFwdMax + r2_canBwdMax) < (r1_canBwdMax + r2_canFwdMax)) {
                r1_fwd = false;
                r2_fwd = true;
            }
#if !NDEBUG
            assert (r1_fwd != r2_fwd);
#endif
            if (r1_fwd && (!r2_fwd)) {
                my_unm_n_count = compact_unmapped_read (read1->bases, read1->length, u_read->fwd_read);//, threadID);
            if (my_unm_n_count) {
                vaftp->asp->as_unmapped_reads_n_count[threadID] += 1;
                vaftp->asp->as_unmapped_diffs_n_count[threadID] += my_unm_n_count;
            }
                my_unm_n_count = compact_unmapped_read (read2->rc_bases, read2->length, u_read->bwd_read);//, threadID);
            if (my_unm_n_count) {
                vaftp->asp->as_unmapped_reads_n_count[threadID] += 1;
                vaftp->asp->as_unmapped_diffs_n_count[threadID] += my_unm_n_count;
            }
                if (r2_bwd_candis[0] < r1_fwd_candis[0]) {
                    u_read->location = r2_bwd_candis[0];
                } else {
                    u_read->location = r1_fwd_candis[0];
                }
                //check_compact_unmapped_read(read1->bases, read1->length, read1->charBases, u_read->fwd_read);
            } else {
                my_unm_n_count = compact_unmapped_read (read1->rc_bases, read1->length, u_read->bwd_read);//, threadID);
            if (my_unm_n_count) {
                vaftp->asp->as_unmapped_reads_n_count[threadID] += 1;
                vaftp->asp->as_unmapped_diffs_n_count[threadID] += my_unm_n_count;
            }
                my_unm_n_count = compact_unmapped_read (read2->bases, read2->length, u_read->fwd_read);//, threadID);
            if (my_unm_n_count) {
                vaftp->asp->as_unmapped_reads_n_count[threadID] += 1;
                vaftp->asp->as_unmapped_diffs_n_count[threadID] += my_unm_n_count;
            }
                if (r1_bwd_candis[0] < r2_fwd_candis[0]) {
                    u_read->location = r1_bwd_candis[0];
                } else {
                    u_read->location = r2_fwd_candis[0];
                }
                //check_compact_unmapped_read(read2->bases, read2->length, read2->charBases, u_read->fwd_read);
            }
        } else if (r1_plaFwdCnt || r1_plaBwdCnt || r2_plaFwdCnt || r2_plaBwdCnt) {
            vaftp->asp->as_pla_unmapped_count[threadID] += 1;
            
            r1_fwd = true;
            r2_fwd = false;
            if ((r1_plaFwdCnt + r2_plaBwdCnt) < (r1_plaBwdCnt + r2_plaFwdCnt)) {
                r1_fwd = false;
                r2_fwd = true;
            }
#if !NDEBUG
            assert (r1_fwd != r2_fwd);
#endif
            if (r1_fwd && (!r2_fwd)) {
                my_unm_n_count = compact_unmapped_read (read1->bases, read1->length, u_read->fwd_read);//, threadID);
            if (my_unm_n_count) {
                vaftp->asp->as_unmapped_reads_n_count[threadID] += 1;
                vaftp->asp->as_unmapped_diffs_n_count[threadID] += my_unm_n_count;
            }
                my_unm_n_count = compact_unmapped_read (read2->rc_bases, read2->length, u_read->bwd_read);//, threadID);
            if (my_unm_n_count) {
                vaftp->asp->as_unmapped_reads_n_count[threadID] += 1;
                vaftp->asp->as_unmapped_diffs_n_count[threadID] += my_unm_n_count;
            }
                if (best_r2_bwd_pla < best_r1_fwd_pla) {
                    u_read->location = best_r2_bwd_pla;
                } else {
                    u_read->location = best_r1_fwd_pla;
                }
                //check_compact_unmapped_read(read1->bases, read1->length, read1->charBases, u_read->fwd_read);
            } else {
                my_unm_n_count = compact_unmapped_read (read1->rc_bases, read1->length, u_read->bwd_read);//, threadID);
            if (my_unm_n_count) {
                vaftp->asp->as_unmapped_reads_n_count[threadID] += 1;
                vaftp->asp->as_unmapped_diffs_n_count[threadID] += my_unm_n_count;
            }
                my_unm_n_count = compact_unmapped_read (read2->bases, read2->length, u_read->fwd_read);//, threadID);
            if (my_unm_n_count) {
                vaftp->asp->as_unmapped_reads_n_count[threadID] += 1;
                vaftp->asp->as_unmapped_diffs_n_count[threadID] += my_unm_n_count;
            }
                if (best_r1_bwd_pla < best_r2_fwd_pla) {
                    u_read->location = best_r1_bwd_pla;
                } else {
                    u_read->location = best_r2_fwd_pla;
                }
                //check_compact_unmapped_read(read2->bases, read2->length, read2->charBases, u_read->fwd_read);
            }
        } else {
            vaftp->asp->as_nil_unmapped_count[threadID] += 1;

            my_unm_n_count = compact_unmapped_read (read1->bases, read1->length, u_read->fwd_read);//, threadID);
            if (my_unm_n_count) {
                vaftp->asp->as_unmapped_reads_n_count[threadID] += 1;
                vaftp->asp->as_unmapped_diffs_n_count[threadID] += my_unm_n_count;
            }
            my_unm_n_count = compact_unmapped_read (read2->rc_bases, read2->length, u_read->bwd_read);//, threadID);
            if (my_unm_n_count) {
                vaftp->asp->as_unmapped_reads_n_count[threadID] += 1;
                vaftp->asp->as_unmapped_diffs_n_count[threadID] += my_unm_n_count;
            }
            //check_compact_unmapped_read(read1->bases, read1->length, read1->charBases, u_read->fwd_read);
#if !NDEBUG
            assert (best_r1_fwd_loc == ((2^32) - 2));
#endif
            u_read->location = best_r1_fwd_loc;
        }
    }
    
    return (r1_map && r2_map);
}

#pragma GCC push_options
#pragma GCC optimize ("unroll-loops")
inline void CPUVMBMKernel(VerifyArgsForThread * vaftp, const uint32_t regIndex, const uint8_t* text, const int readLen, 
const uint32_t *candidates, const uint32_t candiNum, int16_t *errors, int16_t *locations) {
	std::uint8_t * my_ref_bases = vaftp->vaft_com_ds->ref_bases;
	int refLen = readLen + 2 * EDIT_DISTANCE;
	int band_length = 2 * EDIT_DISTANCE + 1;
	int bandLenMask = 0;
	//1 on least significant bits
	for (int i = 0; i < band_length; ++i) {
		bandLenMask <<= 1;
		++bandLenMask;
	}
	__m128i Peq[ALPHABET_SIZE];
	for (int ai = 0; ai < ALPHABET_SIZE; ai++) {
		Peq[ai] = _mm_set1_epi16(0);
	}
	__m128i patterns[V_CPU];
//#pragma unroll (V_CPU)
	//for (int pi = 0; pi < V_CPU; ++pi) {
	//	if (regIndex * V_CPU + pi >= candiNum) {
//			patterns[pi] = _mm_set1_epi8(4);
//		} else {
//			patterns[pi] =
//					_mm_loadu_si128(
//							(__m128i *) (refs[refIndex].bases
//									+ candidates[regIndex * V_CPU + pi]));
//		}
//		__m128i result;
//		int tmp = 0;
//#pragma unroll (ALPHABET_SIZE-1)
//		for (int ai = 0; ai < ALPHABET_SIZE - 1; ++ai) {
//			__m128i letter = _mm_set1_epi8((uint8_t) ai);
//			result = _mm_cmpeq_epi8(letter, patterns[pi]);
//			tmp = _mm_movemask_epi8(result);
//			tmp &= bandLenMask;
//			Peq[ai] = _mm_insert_epi16(Peq[ai], tmp, pi);
//		}
//	}

    if (regIndex * V_CPU + 0 >= candiNum) {
        patterns[0] = _mm_set1_epi8(4);
    } else {    
        patterns[0] =  _mm_loadu_si128( (__m128i *) (my_ref_bases + candidates[regIndex * V_CPU + 0]));
    }
    __m128i result;
    int tmp = 0;
    for (int ai = 0; ai < ALPHABET_SIZE - 1; ++ai) {
        __m128i letter = _mm_set1_epi8((uint8_t) ai);
        result = _mm_cmpeq_epi8(letter, patterns[0]);
        tmp = _mm_movemask_epi8(result);
        tmp &= bandLenMask;
        Peq[ai] = _mm_insert_epi16(Peq[ai], tmp, 0);
    }

    if (regIndex * V_CPU + 1 >= candiNum) {
        patterns[1] = _mm_set1_epi8(4);
    } else {                                   
        patterns[1] =  _mm_loadu_si128( (__m128i *) (my_ref_bases + candidates[regIndex * V_CPU + 1]));
    }                                                    
    for (int ai = 0; ai < ALPHABET_SIZE - 1; ++ai) {
        __m128i letter = _mm_set1_epi8((uint8_t) ai);
        result = _mm_cmpeq_epi8(letter, patterns[1]);
        tmp = _mm_movemask_epi8(result);            
        tmp &= bandLenMask;                                      
        Peq[ai] = _mm_insert_epi16(Peq[ai], tmp, 1);
    }


    if (regIndex * V_CPU + 2 >= candiNum) {
        patterns[2] = _mm_set1_epi8(4);
    } else {                                   
        patterns[2] =  _mm_loadu_si128( (__m128i *) (my_ref_bases + candidates[regIndex * V_CPU + 2]));
    }                                                    
    for (int ai = 0; ai < ALPHABET_SIZE - 1; ++ai) {
        __m128i letter = _mm_set1_epi8((uint8_t) ai);
        result = _mm_cmpeq_epi8(letter, patterns[2]);
        tmp = _mm_movemask_epi8(result);            
        tmp &= bandLenMask;                                      
        Peq[ai] = _mm_insert_epi16(Peq[ai], tmp, 2);
    }


    if (regIndex * V_CPU + 3 >= candiNum) {
        patterns[3] = _mm_set1_epi8(4);
    } else {                                   
        patterns[3] =  _mm_loadu_si128( (__m128i *) (my_ref_bases + candidates[regIndex * V_CPU + 3]));
    }                                                    
    for (int ai = 0; ai < ALPHABET_SIZE - 1; ++ai) {
        __m128i letter = _mm_set1_epi8((uint8_t) ai);
        result = _mm_cmpeq_epi8(letter, patterns[3]);
        tmp = _mm_movemask_epi8(result);            
        tmp &= bandLenMask;                                      
        Peq[ai] = _mm_insert_epi16(Peq[ai], tmp, 3);
    }


    if (regIndex * V_CPU + 4 >= candiNum) {
        patterns[4] = _mm_set1_epi8(4);
    } else {                                   
        patterns[4] =  _mm_loadu_si128( (__m128i *) (my_ref_bases + candidates[regIndex * V_CPU + 4]));
    }                                                    
    for (int ai = 0; ai < ALPHABET_SIZE - 1; ++ai) {
        __m128i letter = _mm_set1_epi8((uint8_t) ai);
        result = _mm_cmpeq_epi8(letter, patterns[4]);
        tmp = _mm_movemask_epi8(result);            
        tmp &= bandLenMask;                                      
        Peq[ai] = _mm_insert_epi16(Peq[ai], tmp, 4);
    }


    if (regIndex * V_CPU + 5 >= candiNum) {
        patterns[5] = _mm_set1_epi8(4);
    } else {                                   
        patterns[5] =  _mm_loadu_si128( (__m128i *) (my_ref_bases + candidates[regIndex * V_CPU + 5]));
    }                                                    
    for (int ai = 0; ai < ALPHABET_SIZE - 1; ++ai) {
        __m128i letter = _mm_set1_epi8((uint8_t) ai);
        result = _mm_cmpeq_epi8(letter, patterns[5]);
        tmp = _mm_movemask_epi8(result);            
        tmp &= bandLenMask;                                      
        Peq[ai] = _mm_insert_epi16(Peq[ai], tmp, 5);
    }


    if (regIndex * V_CPU + 6 >= candiNum) {
        patterns[6] = _mm_set1_epi8(4);
    } else {                                   
        patterns[6] =  _mm_loadu_si128( (__m128i *) (my_ref_bases + candidates[regIndex * V_CPU + 6]));
    }                                                    
    for (int ai = 0; ai < ALPHABET_SIZE - 1; ++ai) {
        __m128i letter = _mm_set1_epi8((uint8_t) ai);
        result = _mm_cmpeq_epi8(letter, patterns[6]);
        tmp = _mm_movemask_epi8(result);            
        tmp &= bandLenMask;                                      
        Peq[ai] = _mm_insert_epi16(Peq[ai], tmp, 6);
    }


    if (regIndex * V_CPU + 7 >= candiNum) {
        patterns[7] = _mm_set1_epi8(4);
    } else {                                   
        patterns[7] =  _mm_loadu_si128( (__m128i *) (my_ref_bases + candidates[regIndex * V_CPU + 7]));
    }                                                    
    for (int ai = 0; ai < ALPHABET_SIZE - 1; ++ai) {
        __m128i letter = _mm_set1_epi8((uint8_t) ai);
        result = _mm_cmpeq_epi8(letter, patterns[7]);
        tmp = _mm_movemask_epi8(result);            
        tmp &= bandLenMask;                                      
        Peq[ai] = _mm_insert_epi16(Peq[ai], tmp, 7);
    }



//a0,a1,a2...a7,a8,a9...->a0,b0,c0...a1,b1,c1...
	permute(patterns);
	uint16_t Mask = (uint16_t) 1 << (band_length - 1);
	__m128i PeqMask = _mm_set1_epi16(Mask);
	__m128i VP = _mm_set1_epi16(0);
	__m128i VN = _mm_set1_epi16(0);
	__m128i X = _mm_set1_epi16(0);
	__m128i D0 = _mm_set1_epi16(0);
	__m128i HN = _mm_set1_epi16(0);
	__m128i HP = _mm_set1_epi16(0);
	__m128i maxMask = _mm_set1_epi16(0xffff);
	__m128i errMask = _mm_set1_epi16(1);
	__m128i errs = _mm_set1_epi16(0);
	int i_bd = 2 * EDIT_DISTANCE;
	int last_high = 2 * EDIT_DISTANCE;
	__m128i threshold = _mm_set1_epi16(EDIT_DISTANCE * 3);//_mm_set1_epi16(EDIT_DISTANCE + 13);
	for (int i = 0; i < readLen; i++) {
		X = _mm_or_si128(Peq[text[i]], VN);
		D0 = _mm_and_si128(X, VP);
		D0 = _mm_add_epi16(D0, VP);
		D0 = _mm_xor_si128(D0, VP);
		D0 = _mm_or_si128(D0, X);
		HN = _mm_and_si128(VP, D0);
		HP = _mm_or_si128(VP, D0);
		HP = _mm_xor_si128(HP, maxMask);
		HP = _mm_or_si128(HP, VN);
		X = _mm_srli_epi16(D0, 1);
		VN = _mm_and_si128(X, HP);
		VP = _mm_or_si128(X, HP);
		VP = _mm_xor_si128(VP, maxMask);
		VP = _mm_or_si128(VP, HN);
		__m128i E = _mm_and_si128(D0, errMask);
		E = _mm_xor_si128(E, errMask);
		errs = _mm_add_epi16(errs, E);
		__m128i earlyEnd = _mm_cmpgt_epi16(errs, threshold);
		int tmp = _mm_movemask_epi8(earlyEnd);
		if (tmp == 0xffff) {
			_mm_store_si128((__m128i *) errors, errs);
			return;
		}
		++i_bd;
		if (i_bd % 16 == 0) {
//#pragma unroll (V_CPU)
			for (int pi = 0; pi < V_CPU; ++pi) {
				if (regIndex * V_CPU + pi >= candiNum) {
					patterns[pi] = _mm_set1_epi8(4);
				} else {
					patterns[pi] = _mm_loadu_si128( (__m128i *) (my_ref_bases + candidates[regIndex * V_CPU + pi] + i_bd));
				}
			}
			permute(patterns);
		}
		int pi = (i_bd % 16) / 2;
		if ((i_bd % 16) % 2 != 0) {
			patterns[pi] = _mm_unpackhi_epi64(patterns[pi], patterns[pi]);
		}
//#pragma unroll (ALPHABET_SIZE-1)
		for (int ai = 0; ai < ALPHABET_SIZE - 1; ++ai) {
			Peq[ai] = _mm_srli_epi16(Peq[ai], 1);
			__m128i letter = _mm_set1_epi16((int16_t) ai);
			__m128i tmpPat = _mm_cvtepu8_epi16(patterns[pi]);
			__m128i mask = _mm_cmpeq_epi16(letter, tmpPat);
			mask = _mm_and_si128(mask, PeqMask);
			Peq[ai] = _mm_or_si128(Peq[ai], mask);
		}
	}
	int site = refLen - last_high - 1;
	__m128i tmpLocs = _mm_set1_epi16(site);
	__m128i tmpErrs = _mm_set1_epi16(EDIT_DISTANCE + 1);//_mm_set1_epi16(EDIT_DISTANCE + 13);
	tmpErrs = _mm_min_epu16(tmpErrs, errs);
	threshold = _mm_set1_epi16(EDIT_DISTANCE + 1);//_mm_set1_epi16(EDIT_DISTANCE + 13);
	for (int i = 0; i < last_high; i++) {
		__m128i tmpVP = _mm_and_si128(VP, errMask);
		__m128i tmpVN = _mm_and_si128(VN, errMask);
		errs = _mm_add_epi16(errs, tmpVP);
		errs = _mm_sub_epi16(errs, tmpVN);
		__m128i tmpMask1 = _mm_cmplt_epi16(errs, threshold);
		__m128i tmpMask2 = _mm_cmplt_epi16(errs, tmpErrs);
        //__m128i tmpMask3 = _mm_cmpeq_epi16(errs,tmpErrs);
        //tmpMask2 =  _mm_or_si128(tmpMask2, tmpMask3);
		tmpErrs = _mm_min_epu16(tmpErrs, errs);
		tmpMask1 = _mm_and_si128(tmpMask1, tmpMask2);
		int tmp = _mm_movemask_epi8(tmpMask1);
//#pragma unroll (V_CPU)
		for (int pi = 0; pi < V_CPU; ++pi) {
			if (tmp & (1 << (pi * 2))) {
			   tmpLocs = _mm_insert_epi16(tmpLocs, site + i + 1, pi);
			}
		}
		VP = _mm_srli_epi16(VP, 1);
		VN = _mm_srli_epi16(VN, 1);
	}
	_mm_store_si128((__m128i *) errors, tmpErrs);
	_mm_store_si128((__m128i *) locations, tmpLocs);
}
#pragma GCC pop_options

int CPUAlignKernel(const uint8_t* pattern, const uint8_t* text, int match_site,
    char* cigar, const int readLen, int best_error) {
	int refLen = readLen + 2 * EDIT_DISTANCE;
	int band_down = 2 * EDIT_DISTANCE;
	int band_length = 2 * EDIT_DISTANCE + 1;

	cigar[0] = 0;

	int start_location = match_site - readLen + 1;
	int tmp_err = 0;

	for (int i = 0; i < readLen; i++)
		if (text[i] != pattern[i + start_location])
			tmp_err++;

	if (tmp_err <= best_error) {
		sprintf(cigar + strlen(cigar), "%d%c", readLen, 'M');
		return start_location;
	}
	//assert (tmp_err > best_error);
	//Encountered atleast one read for which the above assertion failed

	uint32_t D0_arry_64[READ_LEN_MAX];
	uint32_t HP_arry_64[READ_LEN_MAX];
	int Route_Size_Whole[READ_LEN_MAX];
	char Route_Char_Whole[READ_LEN_MAX];

	uint32_t Peq[ALPHABET_SIZE];

	for (int ai = 0; ai < ALPHABET_SIZE; ai++) {
		Peq[ai] = (uint32_t) 0;
	}

	uint32_t tmp = (uint32_t) 1;
	for (int i = 0; i < band_length; i++) {
		Peq[pattern[i]] = Peq[pattern[i]] | tmp;
		tmp = tmp << 1;
	}

	uint32_t Mask = (uint32_t) 1 << (band_length - 1);
	uint32_t VP = 0;
	uint32_t VN = 0;
	uint32_t X = 0;
	uint32_t HN = 0;

	int j = 0;
	int i_bd = band_down;
	int last_high = band_length - readLen + refLen - band_down - 1;

	uint32_t tmp_D0 = 0;
	uint32_t tmp_HP = 0;
	int i = 0;
	for (i = 0; i < readLen; i++) {
		X = Peq[text[i]] | VN;
		tmp_D0 = ((VP + (X & VP)) ^ VP) | X;
		HN = VP & tmp_D0;
		tmp_HP = VN | ~(VP | tmp_D0);
		X = tmp_D0 >> 1;
		VN = X & tmp_HP;
		VP = HN | ~(X | tmp_HP);
		D0_arry_64[i] = tmp_D0;
		HP_arry_64[i] = tmp_HP;

		for (int ai = 0; ai < ALPHABET_SIZE - 1; ai++) {
			Peq[ai] >>= 1;
		}

		++i_bd;
		if ((i_bd) < refLen)
			Peq[pattern[i_bd]] = Peq[pattern[i_bd]] | Mask;
	}

	int site = refLen - last_high - 1;
	int search_site = match_site - site;
	int pre_size = 1;
	char pre_char = 'N';
	uint32_t Mask_1 = (uint32_t) 1;
	i = readLen - 1;
	int sum_err = 0;

	j = 1;
	if (((D0_arry_64[i] >> search_site) & Mask_1)
			&& (pattern[match_site] == text[i])) {
		i--;
		pre_size = 1;
		pre_char = 'M';
		match_site--;
	} else if (!((D0_arry_64[i] >> search_site) & Mask_1)
			&& (pattern[match_site] != text[i])) {
		i--;
		pre_size = 1;
		pre_char = 'S';
		match_site--;
		sum_err++;
	} else if ((HP_arry_64[i] >> search_site) & Mask_1) {
		i--;
		search_site++;
		pre_size = 1;

		pre_char = 'I';
		sum_err++;
		start_location++;
	} else {
		search_site--;
		pre_size = 1;
		pre_char = 'D';

		match_site--;
		sum_err++;
		start_location--;
	}

	while (i >= 0) {
		if (sum_err == EDIT_DISTANCE)
			break;

		if (((D0_arry_64[i] >> search_site) & Mask_1)
				&& (pattern[match_site] == text[i])) {
			i--;
			match_site--;

			if (pre_char != 'M') {
				Route_Size_Whole[j] = pre_size;
				Route_Char_Whole[j++] = pre_char;
				pre_size = 1;
				pre_char = 'M';
			} else
				pre_size++;
		} else if (!((D0_arry_64[i] >> search_site) & Mask_1)
				&& (pattern[match_site] != text[i])) {
			i--;
			match_site--;
			sum_err++;

			if (pre_char != 'S') {
				Route_Size_Whole[j] = pre_size;
				Route_Char_Whole[j++] = pre_char;
				pre_size = 1;
				pre_char = 'S';
			} else
				pre_size++;
		} else if ((HP_arry_64[i] >> search_site) & Mask_1) {
			i--;
			search_site++;
			sum_err++;
			if (pre_char != 'I') {
				Route_Size_Whole[j] = pre_size;
				Route_Char_Whole[j++] = pre_char;
				pre_size = 1;
				pre_char = 'I';
			} else
				pre_size++;
			start_location++;
		} else {
			search_site--;
			match_site--;
			sum_err++;
			if (pre_char != 'D') {
				Route_Size_Whole[j] = pre_size;
				Route_Char_Whole[j++] = pre_char;
				pre_size = 1;
				pre_char = 'D';
			} else
				pre_size++;
			start_location--;
		}
	}

	Route_Size_Whole[j] = pre_size;
	Route_Char_Whole[j++] = pre_char;
	if (i >= 0) {
		Route_Size_Whole[j] = i + 1;
		Route_Char_Whole[j++] = 'M';
	}

	int size_SM = 0;
	for (j = j - 1; j > 0; j--) {
		if (Route_Char_Whole[j] == 'M' || Route_Char_Whole[j] == 'S') {
			size_SM = 0;
			while (j > 0
					&& (Route_Char_Whole[j] == 'M' || Route_Char_Whole[j] == 'S')) {
				size_SM = size_SM + Route_Size_Whole[j];
				j--;
			}
			j++;
			sprintf(cigar + strlen(cigar), "%d%c", size_SM, 'M');
		} else
			sprintf(cigar + strlen(cigar), "%d%c", Route_Size_Whole[j],
					Route_Char_Whole[j]);
	}
	
	//if (sum_err > best_error) {
	//    fprintf(stderr, "sum_err: %d, best_error: %d.\n", sum_err, best_error);
	//    fprintf(stderr, "cigar: %s.\n", cigar);
	//    fprintf(stderr, "mapped_read_queue_size: %d.\n", mapped_read_queue_size);
	//}
	//assert (sum_err <= best_error);
	assert (sum_err <= EDIT_DISTANCE);

	return start_location;

}

void CPUMDTag(const uint8_t* pattern, const uint8_t* text, uint32_t startLocation, char* cigar, int* my_diff_locs, int* my_diff_vals) {
	int cigarLen = strlen(cigar);
	char num[READ_LEN_MAX];
	int numIndex = 0;
	//int matchNum = 0;
	const uint8_t *read = text;
	const uint8_t *reference = pattern + startLocation;
	int readIndex = 0;
	int refIndex = 0;
	my_diff_locs[0] = 0;
	my_diff_vals[0] = 0;

	for (int ci = 0; ci < cigarLen; ci++) {
		char c = cigar[ci];
		if (c >= '0' && c <= '9') {
			num[numIndex++] = c;
		} else {
			num[numIndex] = '\0';
			int opNum = atoi(num);
			if (c == 'M') {
				for (int opi = 0; opi < opNum; opi++) {
					if (reference[refIndex] != read[readIndex]) {
						my_diff_locs[0] += 1;
						my_diff_locs[my_diff_locs[0]] = refIndex;
						if (read[readIndex] < reference[refIndex]) {
						    my_diff_vals[my_diff_locs[0]] = read[readIndex];
						} else {
						    my_diff_vals[my_diff_locs[0]] = read[readIndex] - 1;
						}
						//if (matchNum != 0) {
						//	matchNum = 0;
						//}
					//} else {
					//	++matchNum;
					}
					++refIndex;
					++readIndex;
				}
			} else if (c == 'I') {
				for (int opi = 0; opi < opNum; opi++) {
				    my_diff_locs[0] += 1;
				    my_diff_locs[my_diff_locs[0]] = refIndex;
				    my_diff_vals[my_diff_locs[0]] = 4 + read[readIndex];
				    readIndex++;
				}
			} else if (c == 'D') {
				//if (matchNum != 0) {
				//	matchNum = 0;
				//}
				for (int opi = 0; opi < opNum; opi++) {
				    my_diff_locs[0] += 1;
				    my_diff_locs[my_diff_locs[0]] = refIndex;
				    my_diff_vals[my_diff_locs[0]] = 9;
					refIndex++;
				}
			}
			numIndex = 0;
		}
	}
}

void uncompact_mapped_read (uint8_t * my_bases, uint32_t * my_compact_read, uint32_t my_read_length) {
    uint32_t ri = 0, k = 0;
    int z, one_base;
    for (; ri < my_read_length; ri++) {
        one_base = 0;
        z = 1;
        if (TestBit(my_compact_read, k)) {
            one_base |= z;
        }
        k++;
        z = (z << 1);
        if (TestBit(my_compact_read, k)) {
            one_base |= z;
        }
        k++;
        z = (z << 1);
        if (TestBit(my_compact_read, k)) {
            one_base |= z;
        }
        k++;
        //z = (z << 1); //Not necessary
        my_bases[ri] = (uint8_t) one_base;
    }
}


