/*
 *   quaqc: QUick Atac-seq Quality Control
 *   Copyright (C) 2024-2025  Benjamin Jean-Marie Tremblay
 *
 *   This program is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   This program is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with this program.  If not, see <https://www.gnu.org/licenses/>.
 *
 */

#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <math.h>  // sqrt(), lround()
#include <limits.h>
#include <assert.h>
#include "htslib/hts.h"
#include "htslib/sam.h"
#include "bedidx.h"
#include "quaqc.h"
#include "hash.h"
#include "depth.h"
#include "output.h"
#include "zlib.h"

// Chimeric reads generate multiple lines in a SAM. A random fragment of the
// read will be assigned as a primary alignment, whereas the remaining fragments
// will be assigned as supplementary alignments. All fragments will have
// bits 0x40 and 0x80 set.
//
// Secondary alignments are when a read aligns to multiple locations. The first/best
// alignment will be marked as primary, and additional ones as secondary.

#define is_1st(x) (((x)->core.flag & (BAM_FSECONDARY|BAM_FSUPPLEMENTARY)) == 0)
#define is_2nd(x) (((x)->core.flag & BAM_FSECONDARY) != 0)
#define is_sup(x) (((x)->core.flag & BAM_FSUPPLEMENTARY) != 0)
#define is_proper_mate(x) (((x)->core.flag & BAM_FPROPER_PAIR) != 0 && \
    ((x)->core.flag & BAM_FMUNMAP) == 0 && (x)->core.tid == (x)->core.mtid)
#define is_frag_r1(x) (((x)->core.flag & (BAM_FPROPER_PAIR|BAM_FREAD1)) == (BAM_FPROPER_PAIR|BAM_FREAD1))
#define is_pos_strand(x) (((x)->core.flag & BAM_FREVERSE) == 0)

// bed utilities ----------------------------------------------------------------------

static inline hts_pos_t actual_size(const char *seq, const hts_pos_t size, const params_t *params) {
  if (params->tlist != NULL) {
    return bed_total(params->tlist, seq, size);
  } else if (params->blist != NULL) {
    return size - bed_total(params->blist, seq, size);
  } else {
    return size;
  }
}

static inline bool whitelisted(const void *tlist, const char *seq, const bam1_t *aln, const hts_pos_t qend) {
  return tlist == NULL || bed_overlap_within(tlist, seq, aln->core.pos, qend);
}

static inline bool blacklisted(const void *blist, const char *seq, const bam1_t *aln, const hts_pos_t qend) {
  return blist != NULL && bed_overlap(blist, seq, aln->core.pos, qend);
}

// aln utilities ----------------------------------------------------------------------

static inline bool use_seq(const char *seq, void *tseqs) {
  return tseqs == NULL || str_hash_exists(tseqs, seq);
}

static inline bool frag_filter(const bam1_t *aln, const hts_pos_t qlen, const params_t *params) {
  const bool proper_mate = is_proper_mate(aln);
  const hts_pos_t flen = llabs(aln->core.isize);
  // Pass if: (1) is SE, or
  //          (2) is PE, has a passing frag size, no dovetailing, and a propely mapped mate (|-><-|), or
  //          (3) is PE and --use-nomate is set
  return
/*1*/ ((aln->core.flag & BAM_FPAIRED) == 0)       ||
/*2*/ (proper_mate
       && flen <= params->flen_max
       && flen >= params->flen_min
       && (params->use_dovetail || flen >= qlen)) ||
/*3*/ (!proper_mate && params->use_nomate)
  ;
}

static inline bool check_rg(const bam1_t *aln, const params_t *params) {
  uint8_t *rg_p = bam_aux_get(aln, params->rg_tag);
  if (rg_p == NULL) return false;
  char *rg = bam_aux2Z(rg_p);
  if (rg == NULL) return false;
  return str_hash_exists(params->trg, rg);
}

static inline bool use_rg(const bam1_t *aln, const params_t *params) {
  return params->trg == NULL || check_rg(aln, params);
}

static inline bool read_is_visible(const bam1_t *aln, const char *seq, const hts_pos_t qend, const params_t *params) {
  return whitelisted(params->tlist, seq, aln, qend) && !blacklisted(params->blist, seq, aln, qend) && use_rg(aln, params);
}

static inline bool passes_filters(const bam1_t *aln, const hts_pos_t qlen, const params_t *params) {
  return params->use_all ||
    // -t and -b are not considered filters, but rather "visibility" options.
    (
      !(params->no_se && (aln->core.flag & BAM_FPAIRED) == 0)                               &&
      aln->core.qual >= params->mapq                                                        &&
      qlen <= params->qlen_max                                                              &&
      qlen >= params->qlen_min                                                              &&
      (params->use_dups || ((aln->core.flag & BAM_FDUP) == 0))                              &&
      (is_1st(aln) || (params->use_2nd && is_2nd(aln)) || (params->use_chi && is_sup(aln))) &&
      frag_filter(aln, qlen, params)
    )
  ;
}

static hts_pos_t bases_covered(const bam1_t *b) {
  int n_cigar = b->core.n_cigar;
  const uint32_t *cigar = bam_get_cigar(b);
  hts_pos_t l;
  if (b->core.l_qseq) {
    // Does l_qseq exclude intra splices?
    l = b->core.l_qseq;
    int kl, kr;
    for (kl = 0; kl < n_cigar; kl++) {
      if (bam_cigar_op(cigar[kl]) == BAM_CSOFT_CLIP) {
        l -= bam_cigar_oplen(cigar[kl]);
      } else {
        break;
      }
    }
    for (kr = n_cigar-1; kr > kl; kr--) {
      if (bam_cigar_op(cigar[kr]) == BAM_CSOFT_CLIP) {
        l -= bam_cigar_oplen(cigar[kr]);
      } else {
        break;
      }
    }
  } else {
    static const int query[16] = {
    //M I D N  S H P =  X B ? ?  ? ? ? ?
      1,1,0,0, 0,0,0,1, 1,0,0,0, 0,0,0,0
    };
    for (int k = l = 0; k < n_cigar; k++) {
      if (query[bam_cigar_op(cigar[k])]) {
        l += bam_cigar_oplen(cigar[k]);
      }
    }
  }
  return l;
}

static void read_gc(const bam1_t *aln, int *at, int *gc) {
  const uint32_t *cigar = bam_get_cigar(aln);
  static const int query[16] = {
  //M I D N  S H P =  X B ? ?  ? ? ? ?
    1,1,0,0, 0,0,0,1, 1,0,0,0, 0,0,0,0
  };                         // A C G T N
  static const int is_at[5] = { 1,0,0,1,0 };
  static const int is_gc[5] = { 0,1,1,0,0 };
  for (int k = 0, r = 0, b; k < aln->core.n_cigar; k++) {
    b = bam_cigar_oplen(cigar[k]);
    if (query[bam_cigar_op(cigar[k])]) {
      for (int i = r; i < r + b; i++) {
        *at += is_at[seq_nt16_int[bam_seqi(bam_get_seq(aln), i)]];
        *gc += is_gc[seq_nt16_int[bam_seqi(bam_get_seq(aln), i)]];
      }
    }
    r += b;
  }
}

static void add_read_to_tss(int32_t *tss, const int tss_size, int tss_offset, int qlen) {
  if (tss_offset < 0) {
    qlen -= (-1 * tss_offset);
    tss_offset = 0;
  }
  if (qlen + tss_offset > tss_size) {
    qlen = tss_size - tss_offset;
  }
  for (int i = tss_offset; i < tss_offset + qlen; i++) {
    tss[i]++;
  }
  if (tss_offset > tss_size) abort();
}

static void calc_flag_stats(const bam1_t *aln, stats_t *stats) {
  stats->dups_n     += (aln->core.flag & BAM_FDUP) != 0;
  stats->dups_pri_n += ((aln->core.flag & (BAM_FSECONDARY|BAM_FSUPPLEMENTARY)) == 0) && ((aln->core.flag & BAM_FDUP) != 0);;
  stats->se_n       += (aln->core.flag & BAM_FPAIRED) == 0;
  stats->pe_n       += (aln->core.flag & BAM_FPAIRED) != 0;
  stats->pri_n      += (aln->core.flag & (BAM_FSECONDARY|BAM_FSUPPLEMENTARY)) == 0;
  stats->sec_n      += (aln->core.flag & BAM_FSECONDARY) != 0;
  stats->sup_n      += (aln->core.flag & BAM_FSUPPLEMENTARY) != 0;
  stats->mated_n    += (aln->core.flag & (BAM_FPAIRED|BAM_FPROPER_PAIR)) == (BAM_FPAIRED|BAM_FPROPER_PAIR);
}

// stats_t handling ----------------------------------------------------------------------

static stats_t *init_stats(const int i, const bam_hdr_t *hdr, const params_t *params) {
  stats_t *stats = alloc(sizeof(stats_t));
  stats->id = i;
  stats->type = SEQ_NUCL;
  stats->size = hdr->target_len[i];
  stats->peaks_n = bed_n_in_chr(params->peaks, hdr->target_name[i]);
  stats->tss_n = bed_n_in_chr(params->tss, hdr->target_name[i]);
  stats->blist_n = bed_n_in_chr(params->blist, hdr->target_name[i]);
  stats->tlist_n = bed_n_in_chr(params->tlist, hdr->target_name[i]);
  stats->peaks_total = bed_total(params->peaks, hdr->target_name[i], hdr->target_len[i]);
  stats->tss_total = bed_total(params->tss, hdr->target_name[i], hdr->target_len[i]);
  if (use_seq(hdr->target_name[i], params->tseqs)) {
    stats->actual = actual_size(hdr->target_name[i], hdr->target_len[i], params);
  }
  if (params->peaks != NULL) {
    // TODO: Restrict this to peaks overlapping target list/outside blacklist.
    stats->peak_total = bed_total(params->peaks, hdr->target_name[i], hdr->target_len[i]);
  }
  if (str_hash_exists(params->mito, hdr->target_name[i])) {
    stats->type = SEQ_MITO;
  }
  if (str_hash_exists(params->pltd, hdr->target_name[i])) {
    if (stats->type == SEQ_MITO) {
      quit("Error: Found matching names in --mitochondria and --plastids.");
    }
    stats->type = SEQ_PLTD;
  }
  return stats;
}

static void init_stats_list(results_t *results, const bam_hdr_t *hdr, const params_t *params) {
  results->seqs = init_stats(0, hdr, params);
  stats_t *prev = results->seqs;
  for (int i = 1; i < hdr->n_targets; i++) {
    prev->next = init_stats(i, hdr, params);
    prev = prev->next;
  }
}

static void calc_summary_stats(stats_t *stats, const params_t *params) {
  if (stats->actual == 0) return;
  if ((stats->at + stats->gc + stats->n) > 0) {
    stats->gc_pct = (double) stats->gc / (double) (stats->at + stats->gc + stats->n);
  }
  if (stats->filt_n > 0) {
    stats->avg_qlen = (double) stats->qlen_total / (double) stats->filt_n;
    stats->avg_mapq = (double) stats->mapq_total / (double) stats->filt_n;
  }
  if (stats->frags_n > 0) {
    stats->avg_flen = (double) stats->flen_total / (double) stats->frags_n;
  }
  if (params->peaks != NULL && stats->filt_n > 0) {
    stats->frip = (double) stats->rip_n / (double) stats->filt_n;
  }
  if (stats->pri_n > 0) {
    stats->nrf = (double) (stats->pri_n - stats->dups_pri_n) / (double) stats->pri_n;
  }
  stats->avg_cov = (double) stats->filt_n * stats->avg_qlen / (double) stats->actual;
  // TODO: Try to fix this, would be nice to do even when --omit-depth is used.
  // Not sure why it's not giving the right numbers sometimes (!).
  stats->naive_cov = (double) stats->cov / (double) stats->actual;
}

static void calc_nucl_shared_stats_core(const int32_t *arr, const int arr_max, const double arr_avg, int64_t *mem_pct1, int64_t *mem_pct99, int64_t *mem_max, int64_t *mem_min, double *mem_avg, double *mem_sd, const bool use_arr_avg) {
  int64_t tmp_total = 0, tmp_n = 0, tmp_max, tmp_min = 0;
  for (int i = 0; i < arr_max + 1; i++) {
    if (arr[i]) {
      tmp_min = i;
      break;
    }
  }
  // Be careful, very easy to run into integer overflow issues with regular ints here.
  for (int64_t i = 0; i < arr_max + 1; i++) {
    tmp_total += (int64_t) arr[i] * i;
    tmp_n += arr[i];
    if (arr[i]) tmp_max = i;
  }
  if (tmp_n == 0) return;
  int64_t tmp_i1 = tmp_total / 100, tmp_i1n = 0;
  for (int64_t i = 0; i < arr_max + 1; i++) {
    tmp_i1n += (int64_t) arr[i] * i;
    if (tmp_i1n > tmp_i1) {
      *mem_pct1 = i;
      break;
    }
  }
  int64_t tmp_i99 = (tmp_total * 99) / 100, tmp_i99n = 0;
  for (int64_t i = 0; i < arr_max + 1; i++) {
    tmp_i99n +=  (int64_t) arr[i] * i;
    if (tmp_i99n > tmp_i99) {
      *mem_pct99 = i;
      break;
    }
  }
  *mem_max = tmp_max;
  *mem_min = tmp_min;
  double tmp_avg = (double) tmp_total / (double) tmp_n;
  if (use_arr_avg) {
    tmp_avg = arr_avg;
  } else {
    if (mem_avg != NULL) {
      *mem_avg = tmp_avg;
    }
  }
  double var = 0, square_sum = 0;
  for (int i = 0; i < arr_max + 1; i++) {
    var = (double) i - tmp_avg;
    var = var * var;
    square_sum += var * (double) arr[i];
  }
  square_sum /= (double) tmp_n;
  *mem_sd = sqrt(square_sum);
}

static void calc_nucl_shared_stats(stats_t *stats, const globals_t *nucl_shared, const params_t *params) {
  if (params->tss != NULL) {
    int tes_bkg_n = params->tss_size / 4;
    int tes_tss_n = params->tss_size / 10;
    int64_t tes_bkg_sum = 0, tes_tss_max = 0;
    for (int i = 0; i < tes_bkg_n; i++) {
      tes_bkg_sum += nucl_shared->tss[i];
    }
    int tes_tss_n_i = params->tss_size / 2 - tes_tss_n / 2;
    for (int i = tes_tss_n_i; i < tes_tss_n_i + tes_tss_n; i++) {
      tes_tss_max = max(tes_tss_max, nucl_shared->tss[i]);
    }
    if (tes_tss_max == 0) {
      stats->tes = 0.0;
    } else if (tes_bkg_sum == 0) {
      stats->tes = (double) tes_tss_max;
    } else {
      double tes_bkg_avg = (double) tes_bkg_sum / (double) tes_bkg_n;
      stats->tes = (double)tes_tss_max / tes_bkg_avg;
    }
  }
  // This can be done because quaqc requires the depths struct to at minimum contain depths 0 and 1.
  // The current approach adds potentially quite bad accuracy when --target-list is used, since reads
  // that overlap outside the target regions will contribute extra depths.
  if (!params->omit_depth) {
    int64_t depths_total = 0;
    for (int i = 0; i < params->depth_max + 1; i++) {
      depths_total += nucl_shared->depths[i];
    }
    nucl_shared->depths[0] -= (depths_total - stats->actual);
    stats->genom_cov = (double) (stats->actual - nucl_shared->depths[0]) / (double) stats->actual;
    if (stats->genom_cov > 1.0) stats->genom_cov = 1.0;
    calc_nucl_shared_stats_core(nucl_shared->depths, params->depth_max, 0,
        &(stats->pct1_depth), &(stats->pct99_depth), &(stats->max_depth),
        &(stats->min_depth), &(stats->avg_depth), &(stats->sd_depth), false);
  }
  if (!params->omit_gc) {
    calc_nucl_shared_stats_core(nucl_shared->read_gc, 100, 100.0 * stats->gc_pct,
        &(stats->pct1_gc), &(stats->pct99_gc), &(stats->max_gc),
        &(stats->min_gc), NULL, &(stats->sd_gc), true);
  }
  calc_nucl_shared_stats_core(nucl_shared->mapqs, 254, 0,
      &(stats->pct1_mapq), &(stats->pct99_mapq), &(stats->max_mapq),
      &(stats->min_mapq), &(stats->avg_mapq2), &(stats->sd_mapq), false);
  calc_nucl_shared_stats_core(nucl_shared->read_sizes, params->qhist_max, 0,
      &(stats->pct1_qlen), &(stats->pct99_qlen), &(stats->max_qlen),
      &(stats->min_qlen), NULL, &(stats->sd_qlen), false);
  calc_nucl_shared_stats_core(nucl_shared->frag_sizes, params->fhist_max, 0,
      &(stats->pct1_flen), &(stats->pct99_flen), &(stats->max_flen),
      &(stats->min_flen), NULL, &(stats->sd_flen), false);
}

static void merge_two_stats(stats_t *s, const stats_t *s0) {
#ifndef NO_STAT_SIZE_CHECK  // Since the check is not portable...
#ifndef STAT_SIZE
#define STAT_SIZE 552
#endif
// If this assert triggers, make sure to update this function!!!
  assert(sizeof(stats_t) == STAT_SIZE);
#endif
  s->size        += s0->size;
  s->actual      += s0->actual;
  s->peak_total  += s0->peak_total;
  s->peaks_total += s0->peaks_total;
  s->peaks_n     += s0->peaks_n;
  s->tss_n       += s0->tss_n;
  s->blist_n     += s0->blist_n;
  s->tlist_n     += s0->tlist_n;
  s->reads_n     += s0->reads_n;
  s->dups_n      += s0->dups_n;
  s->dups_pri_n  += s0->dups_pri_n;
  s->se_n        += s0->se_n;
  s->pe_n        += s0->pe_n;
  s->un_n        += s0->un_n;
  s->pri_n       += s0->pri_n;
  s->sec_n       += s0->sec_n;
  s->sup_n       += s0->sup_n;
  s->mated_n     += s0->mated_n;
  s->bl_n        += s0->bl_n;
  s->filt_n      += s0->filt_n;
  s->rip_n       += s0->rip_n;
  s->frags_n     += s0->frags_n;
  s->at          += s0->at;
  s->gc          += s0->gc;
  s->n           += s0->n;
  s->cov         += s0->cov;
  s->qlen_total  += s0->qlen_total;
  s->flen_total  += s0->flen_total;
  s->mapq_total  += s0->mapq_total;
}

static stats_t *copy_stats(stats_t *stats) {
  stats_t *copy = alloc(sizeof(stats_t));
  *copy = *stats;
  copy->id = 0;
  copy->next = NULL;
  return copy;
}

static stats_t *merge_all_stats_by_type(stats_t *seqs, const int n, enum seq_type type, const params_t *params) {
  if (n == 0) {
    return NULL;
  } else if (n == 1) {
    stats_t *tmp;
    while (seqs != NULL) {
      tmp = seqs;
      seqs = seqs->next;
      if (tmp->actual > 0 && tmp->type == type) {
        return copy_stats(tmp);
      }
    }
    quit("Internal error: could not find any sequence matching requested type!");
  } else {
    stats_t *merged = alloc(sizeof(stats_t));
    merged->type = type;
    stats_t *tmp;
    while (seqs != NULL) {
      tmp = seqs;
      seqs = seqs->next;
      if (tmp->actual > 0 && tmp->type == type) {
        merge_two_stats(merged, tmp);
      }
    }
    calc_summary_stats(merged, params);
    return merged;
  }
}

// quaqc_run ----------------------------------------------------------------------

void quaqc_run(htsFile *bam, results_t *results, const params_t *params) {

  results->time_start = time(NULL);
  time_t time_start = time(NULL);
  bam1_t *aln = NULL;
  bam_hdr_t *hdr = NULL;
  hts_itr_t *itr = NULL;
  void *depths = NULL;
  void *bedGraph = NULL;
  samFile *filt_bam = NULL;
  depths = init_depths();
  gzFile bedGraph_f = NULL;

  if (params->bedGraph) {
    bedGraph = init_bedGraph(params);
    bedGraph_f = init_bedGraph_f(bam->fn, params);
  }

  hdr = sam_hdr_read(bam);
  if (hdr == NULL) {
    error(params->qerr, "Failed to open header: %s", bam->fn);
    goto run_quaqc_end;
  }
  // Check read groups if necessary.
  if (params->trg != NULL && strcmp("RG", params->rg_tag) == 0) {
    int rg_n = sam_hdr_count_lines(hdr, "RG");
    if (rg_n < 1) {
      error(params->qerr, "No read group info in header for file '%s'.", bam->fn);
      goto run_quaqc_end;
    }
    kstring_t hdr_line = { 0, 0, NULL };
    int rg_trg_n = 0;
    for (unsigned int i = 0; i != str_hash_end(params->trg); i++) {
      if (str_ind_exists(params->trg, i) &&
          sam_hdr_find_line_id(hdr, "RG", "ID", str_hash_key(params->trg, i), &hdr_line) == 0) {
        rg_trg_n++;
      }
    }
    ks_free(&hdr_line);
    if (rg_trg_n == 0) {
      error(params->qerr, "No target read groups in header for file '%s'.", bam->fn);
      goto run_quaqc_end;
    }
  }
  aln = bam_init1();

  if (!hdr->n_targets) {
    if (!bam->fn_aux) {
      error(params->qerr, "Missing header for file '%s'.", bam->fn);
    } else {
      error(params->qerr, "Missing header or bad chrom.sizes for file '%s'.", bam->fn);
    }
    goto run_quaqc_end;
  }

  if (params->save) {
    filt_bam = init_filtered_bam(hdr, bam->fn, params);
    if (filt_bam == NULL) {
      warn("Giving up on creating new BAM for this sample, continuing with QC.");
    }
  }

  results->seq_n = hdr->n_targets;
  results->r_unmapped += hts_idx_get_n_no_coor(bam->idx);
  for (uint64_t i = 0, m, u; i < hdr->n_targets; i++) {
    hts_idx_get_stat(bam->idx, i, &m, &u);
    results->r_mapped   += (int64_t) m;
    results->r_unmapped += (int64_t) u;
  }
  results->r_total = results->r_mapped + results->r_unmapped;

  init_stats_list(results, hdr, params);

  const hts_pos_t tn5_fwd = params->tn5_shift ? TN5_FOWARD_SHIFT : 0;
  const hts_pos_t tn5_rev = params->tn5_shift ? TN5_REVERSE_SHIFT : 0;
  const hts_pos_t bg_tn5_fwd = params->bg_tn5 ? TN5_FOWARD_SHIFT : 0;
  const hts_pos_t bg_tn5_rev = params->bg_tn5 ? TN5_REVERSE_SHIFT : 0;
  int itr_ret, at, gc, n, nucl_n = 0, pltd_n = 0, mito_n = 0, tss_offset;
  int64_t proc_n = 0;
  hts_pos_t qend, qlen, flen, last_start, last_end, tss_qbeg, tss_qend, bg_qbeg0, bg_qbeg, bg_qend;
  stats_t *stats = results->seqs;
  for (int64_t i = 0; i < hdr->n_targets; i++) {
    itr = sam_itr_querys(bam->idx, hdr, hdr->target_name[i]);
    if (itr == NULL) {
      error(params->qerr, "Failed to iterate over input: %s", bam->fn);
      goto run_quaqc_end;
    }
    last_start = last_end = 0;

    results->seq_sum += stats->size;
    switch (stats->type) {
      case SEQ_NUCL:
        results->nuc_sum += stats->size;
        results->nuc_n++;
        break;
      case SEQ_MITO:
        results->mito_sum += stats->size;
        results->mito_n++;
        break;
      case SEQ_PLTD:
        results->pltd_sum += stats->size;
        results->pltd_n++;
        break;
    }

    if (use_seq(hdr->target_name[i], params->tseqs)) {
      results->seq_act += stats->actual;
      if (stats->actual > 0) {
        switch (stats->type) {
          case SEQ_NUCL:
            results->nuc_act += stats->actual;
            nucl_n++;
            break;
          case SEQ_PLTD:
            results->pltd_act += stats->actual;
            pltd_n++;
            break;
          case SEQ_MITO:
            results->mito_act += stats->actual;
            mito_n++;
            break;
        }
      }
      while ((itr_ret = sam_itr_next(bam, itr, aln)) >= 0) {
        if (unlikely((aln->core.flag & BAM_FUNMAP) != 0)) continue;
        check_rg(aln, params);
        proc_n++; stats->reads_n++; results->r_seen++;

        if (proc_n % 10000000 == 0) {
          vvmsg("... %s: processed %'"PRId64" M reads, %.1f%% passing filters\n", bam->fn, proc_n / 1000000, 100.0 * (double)results->r_filt / (double)(proc_n - 1));
        }

        qend = bam_endpos(aln); // Half open!!!
        qlen = bases_covered(aln);

        calc_flag_stats(aln, stats);

        // Post-filter stats

        if (likely(read_is_visible(aln, hdr->target_name[i], qend, params))) {
          if (likely(passes_filters(aln, qlen, params))) {
            stats->filt_n++; results->r_filt++;

            // Stats for all sequence types

            if (aln->core.pos > last_end) {
              stats->cov += (int64_t) (last_end - last_start);
              last_start = aln->core.pos;
              last_end = qend;
            } else {
              last_end = qend;
            }
            flen = llabs(aln->core.isize);

            // TODO: this 255 score thing doesn't always apply to all mappers
            stats->mapq_total += aln->core.qual < 255 ? (int64_t) aln->core.qual : 0;
            stats->qlen_total += qlen;

            at = gc = n = 0;
            if (!params->omit_gc) {
              read_gc(aln, &at, &gc);
              n = qlen - at - gc;
              stats->at += at;
              stats->gc += gc;
              stats->n  += n;
            }

            if (params->peaks != NULL) {
              stats->rip_n += bed_overlap(params->peaks, hdr->target_name[i], aln->core.pos, qend);
            }

            // Read stats for nuclear sequences only (nucl_shared)

            if (stats->type == SEQ_NUCL) {
              if (params->tss != NULL) {
                if (params->tss_qlen == 0) {
                  tss_qbeg = aln->core.pos + tn5_fwd;
                  tss_qend = max(qend - tn5_rev, tss_qbeg + 1);
                } else {
                  if (is_pos_strand(aln)) {
                    tss_qbeg = aln->core.pos + tn5_fwd;
                    tss_qend = tss_qbeg + params->tss_qlen / 2 + 1;
                    tss_qbeg = tss_qbeg - params->tss_qlen / 2;
                  } else {
                    tss_qend = qend - tn5_rev;
                    tss_qbeg = tss_qend - (params->tss_qlen / 2 + 1);
                    tss_qend = tss_qend + params->tss_qlen / 2;
                  }
                  if (params->tss_qlen % 2 == 0) {
                    if (is_pos_strand(aln)) {
                      tss_qend--;
                    } else {
                      tss_qbeg++;
                    }
                  }
                }
                tss_qbeg = max(0, tss_qbeg);
                tss_qend = min(tss_qend, hdr->target_len[i]);
                tss_offset = bed_overlap_offset(params->tss, hdr->target_name[i], tss_qbeg, tss_qend);
                if (tss_offset != INT_MIN) {
                  add_read_to_tss(results->nucl_shared->tss, params->tss_size, tss_offset, tss_qend - tss_qbeg);
                }
              }
              results->nucl_shared->mapqs[aln->core.qual]++;
              if ((at + gc + n) > 0) {
                const double gc_frac = (double)(100 * gc) / (double)(at + gc + n); 
                results->nucl_shared->read_gc[lround(gc_frac)]++;
              }
              results->nucl_shared->read_sizes[min(qlen, params->qhist_max)]++;
              if (!params->omit_depth) {
                add_read_to_depths(aln, qend, depths, results->nucl_shared->depths, params->depth_max);
              }
              if (bedGraph_f != NULL) {
                if (params->bg_qlen == 0) {
                  bg_qbeg = aln->core.pos + bg_tn5_fwd;
                  bg_qbeg0 = bg_qbeg;
                  bg_qend = max(qend - bg_tn5_rev, bg_qbeg + 1);
                } else {
                  if (is_pos_strand(aln)) {
                    bg_qbeg = aln->core.pos + bg_tn5_fwd;
                    bg_qend = bg_qbeg + params->bg_qlen / 2 + 1;
                    bg_qbeg = bg_qbeg - params->bg_qlen / 2;
                    bg_qbeg0 = bg_qbeg;
                  } else {
                    bg_qend = qend - bg_tn5_rev;
                    bg_qbeg = bg_qend - (params->bg_qlen / 2 + 1);
                    bg_qend = bg_qend + params->bg_qlen / 2;
                    bg_qbeg0 = min(aln->core.pos, bg_qbeg);
                  }
                  if (params->bg_qlen % 2 == 0) {
                    if (is_pos_strand(aln)) {
                      bg_qend--;
                    } else {
                      bg_qbeg++;
                    }
                  }
                }
                bg_qbeg0 = max(0, bg_qbeg0);
                bg_qbeg = max(0, bg_qbeg);
                bg_qend = min(bg_qend, hdr->target_len[i]);
                add_read_to_bedGraph(bedGraph_f, bedGraph, bg_qbeg0, bg_qbeg, bg_qend, hdr->target_name[i]);
              }
              if (filt_bam != NULL && sam_write1(filt_bam, hdr, aln) == -1) {
                error(params->qerr, "Failed to write to new BAM '%s'.", filt_bam->fn);
                warn("Giving up on creating new BAM for this sample, continuing with QC.");
                hts_close(filt_bam);
                filt_bam = NULL;
              }
            } // END: if (stats->type == SEQ_NUCL)

            // Fragment stats

            if (is_frag_r1(aln)) {
              stats->frags_n++;
              stats->flen_total += flen;
              if (stats->type == SEQ_NUCL) {
                results->nucl_shared->frag_sizes[min(flen, params->fhist_max)]++;
              }
            }

          } // END: if (passes_filters(aln, hdr->target_name[i], qend, params))
        } else {
          stats->bl_n++; results->r_mapped_bl++;
        } // END: if (likely(read_is_visible(aln, hdr->target_name[i], qend, params))) {
      } // END: while ((itr_ret = sam_itr_next(bam, itr, aln)) >= 0)
      if (itr_ret < -1) {
        error(params->qerr, "Encountered error reading sample '%s'.", bam->fn);
        goto run_quaqc_end;
      }

      // Once-per-sequence stats

      stats->cov += (int64_t) (last_end - last_start);
      calc_summary_stats(stats, params);

      // Clean up depths struct and add remaining positions to results->nucl_shared

      purge_and_reset_depths(depths, results->nucl_shared->depths, params->depth_max);
      if (bedGraph_f != NULL) {
        purge_and_reset_bedGraph(bedGraph_f, bedGraph, hdr->target_name[i]);
      }

    } else {
      uint64_t m, u;
      hts_idx_get_stat(bam->idx, i, &m, &u);
      stats->bl_n += (int64_t) m;
      // TODO: The stats->bl_n value never gets merged with other sequences of the same type...
      results->r_mapped_bl += (int64_t) m;
    } // END: if (use_seq(hdr->target_name[i], params->tseqs))

#ifdef DEBUG0
    print_stats(hdr->target_name[i], stats);
#endif

    sam_itr_destroy(itr);
    stats = stats->next;

  } // END: for (int64_t i = 0; i < hdr->n_targets; i++)

  // Combine stats per sequence type

  if (params->tseqs != NULL && (nucl_n + pltd_n + mito_n) == 0) {
    error(params->qerr, "Could not find any of the specified targets in '%s'.", bam->fn);
    goto run_quaqc_end;
  }

  results->nucl = merge_all_stats_by_type(results->seqs, nucl_n, SEQ_NUCL, params);
  results->pltd = merge_all_stats_by_type(results->seqs, pltd_n, SEQ_PLTD, params);
  results->mito = merge_all_stats_by_type(results->seqs, mito_n, SEQ_MITO, params);

  // Calculate extra nuclear sequence-only stats

  if (results->nucl != NULL) {
    if (results->nucl_shared->depths[params->depth_max]) {
      msg("Warning: Some depths likely have been truncated for '%s' (--max-depth)\n", bam->fn);
    }
    if (results->nucl_shared->frag_sizes[params->fhist_max]) {
      msg("Warning: Some fragment sizes likely have been truncated for '%s' (--max-fhist)\n", bam->fn);
    }
    if (results->nucl_shared->read_sizes[params->qhist_max]) {
      msg("Warning: Some read sizes likely have been truncated for '%s' (--max-qhist)\n", bam->fn);
    }
    calc_nucl_shared_stats(results->nucl, results->nucl_shared, params);
  }

  if (params->v) {
    time_t time_end = time(NULL);
    double reads_p = proc_n == 0 ? 0.0 : 100.0 * (double) results->r_filt / (double) proc_n;
    msg("Finished file: %s\n", bam->fn);
    if (proc_n > (1 << 30)) {
      msg("|--> Processed %'.1f B reads, ", (((double) proc_n / 1000.0) / 1000.0) / 1000.0);
    } else if (proc_n > (1 << 20)) {
      msg("|--> Processed %'.1f M reads, ", ((double) proc_n / 1000.0) / 1000.0);
    } else if (proc_n > (1 << 10)) {
      msg("|--> Processed %'.1f K reads, ", (double) proc_n / 1000.0);
    } else {
      msg("|--> Processed %'"PRId64" read%s, ", proc_n, pluralize(proc_n));
    }
    msg("%.1f%% passing filters.\n|--> Time to process: ", reads_p);
    print_time(difftime(time_end, time_start));
    if (filt_bam != NULL) {
      msg("|--> Created filtered BAM: %s\n", filt_bam->fn);
    }
  }

  results->success = true;
  results->bam_out = filt_bam != NULL;

run_quaqc_end:

  if (depths != NULL) destroy_depths(depths);
  if (hdr != NULL) sam_hdr_destroy(hdr);
  if (aln != NULL) bam_destroy1(aln);
  if (filt_bam != NULL) hts_close(filt_bam);
  if (bedGraph != NULL) {
    destroy_depths(bedGraph);
    if (bedGraph_f != NULL && gzclose(bedGraph_f) != Z_OK) {
      int e;
      error(params->qerr, "Failed to close bedGraph file: %s", gzerror(bedGraph_f, &e));
    }
  }

  results->time_end = time(NULL);

}

