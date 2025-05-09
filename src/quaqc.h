#ifndef QUAQC_H
#define QUAQC_H

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

#include <stdint.h>
#include <stdbool.h>
#include <time.h>
#include "htslib/hts.h"

#define QUAQC_VERSION  "1.3c"
#define QUAQC_YEAR    "2025"

// Command line & other defaults -------------------------------------------------

#define DEFAULT_MAPQ                  30
#define DEFAULT_MIN_QLEN              15
#define DEFAULT_MIN_FLEN              15
#define DEFAULT_MAX_QLEN             250
#define DEFAULT_MAX_FLEN            2000
#define DEFAULT_MAX_DEPTH         100000
#define DEFAULT_TSS_SIZE            2001
#define DEFAULT_TSS_QLEN             100
#define DEFAULT_BG_QLEN              100
#define DEFAULT_THREADS                1

#define DEFAULT_OUT_EXT     ".quaqc.txt"
#define DEFAULT_BAM_EXT      ".filt.bam"
#define DEFAULT_BG_EXT    ".bedGraph.gz"

#define DEFAULT_MITO      "chrM,ChrM,Mt,MT,MtDNA,mit,Mito,mitochondria,mitochondrion"
#define DEFAULT_PLTD      "chrC,ChrC,Pt,PT,Pltd,Chloro,chloroplast"

#define DEFAULT_RG_TAG              "RG"

#define TN5_FOWARD_SHIFT               4
#define TN5_REVERSE_SHIFT              5

// Functions, macros ------------------------------------------------------------

#define I32_ALLOC_WARNING_TRIGGER 26214400  // x4 = 100 MB
#define check_for_large_i32(size, need_for) do { \
  if ((size) >= I32_ALLOC_WARNING_TRIGGER) { \
    fprintf(stderr, "Warning: Large memory request due to %s=%d (%.2f MB)\n", \
      need_for, (size), (double)((size) * 4) / 1024.0 / 1024.0); \
  } } while (0) 

// Maximum number of elements allowed for hash table construction when parsing
// -p, -m, -n, -r, and -R.
#define MAX_HASH_SIZE    1000000000

#define   likely(cond) __builtin_expect((cond), 1)
#define unlikely(cond) __builtin_expect((cond), 0)

#define min(x, y) ((x) < (y) ? (x) : (y))
#define max(x, y) ((x) > (y) ? (x) : (y))

#define error(do_exit, msg, ...) do { \
    fprintf(stderr, "[E::%s] " msg "\n", __func__, ##__VA_ARGS__); \
    if (do_exit) { \
      fputs("Encountered fatal error, exiting. Run quaqc -h for usage.\n", stderr); \
      exit(EXIT_FAILURE); \
    } } while (0) 
#define quit(msg, ...) error(true, msg, ##__VA_ARGS__)
#define warn(msg, ...) error(false, msg, ##__VA_ARGS__)
#define TODO(m) quit("Feature '%s' not yet implemented.", m)
#define msg(...) if (params->v) fprintf(stderr, __VA_ARGS__)
#define vvmsg(...) if (params->vv) fprintf(stderr, __VA_ARGS__)

#define pluralize(C) ((C) == 1 ? "" : "s")

void *calloc_or_die(size_t size, const char *func_name);
#define alloc(size) calloc_or_die((size), __func__)

void print_time(const time_t s);

long get_mem(void);

// Enums, structs ------------------------------------------------------------

enum seq_type {
  SEQ_NUCL,
  SEQ_PLTD,
  SEQ_MITO
};

typedef struct params_t {
  char **argv; int argc, flag_n;
  void *mito, *pltd, *tseqs;
  char *title, *out_dir, *out_ext, *keep_dir, *keep_ext, *json;
  char *bg_ext, *bg_dir, *rg_tag;
  void *peaks, *tss, *blist, *tlist, *trg;
  int peaks_n, tss_n, blist_n, tlist_n, trg_n, tn5_fwd, tn5_rev;
  uint8_t mapq;
  hts_pos_t qlen_min, flen_min, qlen_max, flen_max;
  int threads, depth_max, qhist_max, fhist_max, tss_size, tss_qlen, bg_qlen;
  bool use_2nd, use_chi, use_nomate, use_dups, use_all;
  bool lenient, strict, nfr, nbr, use_dovetail;
  bool no_se, no_out, save, bedGraph, bg_tn5, tn5_shift, omit_gc, omit_depth;
  bool fast, qerr, v, vv, footprint, low_mem, chip;
} params_t;

// Use one instance of this struct for all nuclear chromosomes.
typedef struct globals_t {
  // Keep adding to these and process at the very end of the BAM:
  // (Note: 255 == MAPQ is unavailable)
  int32_t mapqs[256], read_gc[101];
  int32_t *read_sizes, *frag_sizes, *tss;
  int32_t *depths;
} globals_t;

// One per sequence.
typedef struct stats_t {
  struct stats_t *next;
  int id;
  enum seq_type type;
  // --- Start of what matters for merge_two_stats() ---
  int64_t size, actual, peak_total;
  int peaks_n, tss_n, blist_n, tlist_n;
  int64_t peaks_total, tss_total;
  // All reads stats
  int64_t reads_n, dups_n, dups_pri_n;
  int64_t se_n, pe_n, un_n, pri_n, sec_n, sup_n, mated_n;
  int64_t bl_n;
  // Post-filter reads stats
  int64_t filt_n, rip_n, frags_n, at, gc, n, cov;
  int64_t qlen_total, flen_total, mapq_total;
  // --- End of what matters for merge_two_stats() ---
  // calc_summary_stats()
  double nrf, gc_pct, avg_flen, avg_qlen, avg_mapq;
  double avg_cov, genom_cov, naive_cov, frip;
  // calc_nucl_shared_stats()
  double tes, avg_depth, sd_depth, avg_mapq2, sd_mapq;
  double sd_qlen, sd_flen, sd_gc;
  int64_t min_gc, max_gc, pct1_gc, pct99_gc;
  int64_t max_qlen, max_flen, min_qlen, min_flen;
  int64_t pct1_depth, pct1_mapq, pct1_qlen, pct1_flen;
  int64_t pct99_depth, pct99_mapq, pct99_qlen, pct99_flen;
  int64_t min_depth, max_depth, med_depth, iqr_depth, max_mapq, min_mapq;
} stats_t;

typedef struct results_t {
  bool success, bam_out;
  int64_t seq_n, seq_sum, seq_act;
  int64_t nuc_n, nuc_sum, nuc_act;
  int64_t mito_act, mito_sum, mito_n;
  int64_t pltd_act, pltd_sum, pltd_n;
  int64_t r_mapped, r_mapped_bl, r_unmapped;
  int64_t r_total, r_seen, r_filt;
  globals_t *nucl_shared;
  stats_t *seqs;
  stats_t *nucl;
  stats_t *pltd;
  stats_t *mito;
  time_t time_start, time_end;
} results_t;

#endif
