/*
 *   quaqc: QUick Atac-seq Quality Control
 *   Copyright (C) 2024  Benjamin Jean-Marie Tremblay
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

// TODO: Update with new options (--use-all, --lenient, --nfr, --nbr)
// TODO: Change pileup title when --chip is used

#include <stdlib.h>
#include <string.h>
#include <libgen.h>
#include <time.h>
#include <ctype.h>
#include <errno.h>
#include <wchar.h>
#include "htslib/sam.h"
#include "quaqc.h"
#include "zlib.h"
#include "hash.h"

#define zero_to_one(N0) ((N0) == 0 ? 1 : (N0))
#define calc_pct(NUM, DEN) (100.0 * (double) (NUM) / (double) (zero_to_one(DEN)))
#define null_get(OBJ, MEM) (((OBJ) != NULL) ? ((OBJ)->MEM) : 0)
#define null_get_all(RES, MEM) (null_get((RES)->nucl, MEM) + null_get((RES)->mito, MEM) + null_get((RES)->pltd, MEM))

// TODO: Use static buffers instead of dynamic memory

// Utilities ----------------------------------------------------------------------

static char *paste_strings(char *strs[], int size) {
  int total = 0;
  for (int i = 0; i < size; i++) total += strlen(strs[i]);
  char *out = alloc(total + size + 1);
  for (int j = 0, i = 0; i < size; i++) {
    for (int k = 0; k < strlen(strs[i]); k++) {
      out[j++] = strs[i][k] == '\t' ? ' ' : strs[i][k];
    }
    out[j++] = ' ';
  }
  return out;
}

static char *lower(const char *str) {
  char *lwr = alloc(strlen(str) + 1);
  for (int i = 0; i < strlen(str); i++) {
    lwr[i] = tolower(str[i]);
  }
  return lwr;
}

static bool check_if_matching_fn(const char *old, const char *new) {
  bool answer = false;
  char *old_real = realpath(old, NULL);
  char *new_real = realpath(new, NULL);
  if (old_real != NULL && new_real != NULL && strcmp(old_real, new_real) == 0) {
    answer = true;
  }
  free(old_real);
  free(new_real);
  return answer;
}

static char *make_new_fn(const char *old_fn, char *new_dir, char *new_ext, const bool strip) {
  char *old_fn1 = strdup(old_fn), *old_fn2 = strdup(old_fn);
  char *base_fn = basename(old_fn1);
  if (strip) {
    int bnlen = strlen(base_fn);
    char *ext_lwr = NULL;
    if (bnlen > 3) {
      ext_lwr = lower(base_fn + bnlen - 4);
      if (strcmp(ext_lwr, ".bam") == 0) base_fn[bnlen - 4] = '\0';
    } else if (bnlen > 4) {
      ext_lwr = lower(base_fn + bnlen - 5);
      if (strcmp(ext_lwr, ".cram") == 0) base_fn[bnlen - 5] = '\0';
      base_fn[bnlen - 5] = '\0';
    }
    if (ext_lwr != NULL) free(ext_lwr);
  }
  char *dir_fn = new_dir != NULL ? new_dir : dirname(old_fn2);
  int new_fn_len = strlen(dir_fn) + 1 + strlen(base_fn) + strlen(new_ext) + 1;
  char *new_fn = alloc(new_fn_len + 1);
  strncpy(new_fn, dir_fn, new_fn_len - 1);
  new_fn[strlen(dir_fn)] = '/';
  new_fn_len -= strlen(dir_fn);
  strncpy(new_fn + strlen(dir_fn) + 1, base_fn, new_fn_len - 1);
  new_fn_len -= strlen(base_fn);
  strncpy(new_fn + strlen(dir_fn) + 1 + strlen(base_fn), new_ext, new_fn_len - 1);
  free(old_fn1);
  free(old_fn2);
  return new_fn;
}

static void get_time(time_t time, char *str) {
  struct tm *curr_time_struct = localtime(&time);
  strftime(str, 128, "%F %T %Z", curr_time_struct);
}

// --keep ----------------------------------------------------------------------

samFile *init_filtered_bam(bam_hdr_t *hdr, const char *fn, const params_t *params) {
  char *new_fn = make_new_fn(fn, params->keep_dir, params->keep_ext, true);
  if (check_if_matching_fn(fn, new_fn)) {
    error(params->qerr, "New BAM cannot have an identical path as input '%s'.", fn);
    free(new_fn);
    return NULL;
  }
  samFile *new_bam = sam_open(new_fn, "wb");
  char *CMD = paste_strings(params->argv, params->argc);
  if (sam_hdr_add_pg(hdr, "quaqc", "VN", QUAQC_VERSION, "CL", CMD, NULL) == -1) {
    goto init_filtered_bam_fail;
  }
  if (sam_hdr_write(new_bam, hdr) == -1) {
    goto init_filtered_bam_fail;
  }
  free(new_fn);
  free(CMD);
  return new_bam;
init_filtered_bam_fail:
  error(params->qerr, "Failed to write to new BAM '%s'.", new_fn);
  free(new_fn);
  free(CMD);
  hts_close(new_bam);
  return NULL;
}

// --json ----------------------------------------------------------------------

static bool use_gz = false;
static union fjson_t {
  FILE   *f;
  gzFile  gz;
} fjson;

#define fjson_write(CON, MSG, ...) do { \
  if (use_gz) { \
    int ret = gzprintf(CON.gz, MSG, ##__VA_ARGS__); \
    if (ret <= 0 && strlen(MSG) > 0) { \
      int e; \
      warn("Could not write JSON: %s", gzerror(CON.gz, &e)); \
      return 1; \
    } \
  } else { \
    int ret = fprintf(CON.f, MSG, ##__VA_ARGS__); \
    if (ret < 0) { \
      warn("Could not write JSON: %s", strerror(errno)); \
      return 1; \
    } \
  } } while (0)

#define strbool(COND) ((COND) ? "true" : "false")

int init_json(const params_t *params) {
  char time_start_str[128];
  get_time(time(NULL), time_start_str);
  char *CMD = paste_strings(params->argv + 1, params->flag_n);
  for (int i = 0; i < strlen(CMD); i++) {
    if (CMD[i] == '"') CMD[i] = '\'';
  }
  if (strlen(params->json) > 2 && strcmp(params->json + (strlen(params->json) - 3), ".gz") == 0) {
    use_gz = true;
  }
  if (use_gz) {
    fjson.gz = gzopen(params->json, "wb");
  } else if (strcmp(params->json, "-") == 0) {
    fjson.f = fdopen(1, "w");
  } else {
    fjson.f = fopen(params->json, "w");
  }
  if (use_gz && fjson.gz == NULL) {
    int e;
    warn("Cannot create file '%s': %s", params->json, gzerror(fjson.gz, &e));
    return 1;
  } else if (!use_gz && fjson.f == NULL) {
    warn("Cannot create file '%s': %s", params->json, strerror(errno));
    return 1;
  }

  fjson_write(fjson, "{\n");
  fjson_write(fjson, "  \"quaqc_version\": \"" QUAQC_VERSION "\",\n");
  fjson_write(fjson, "  \"quaqc_args\": \"%s\",\n", CMD);
  fjson_write(fjson, "  \"quaqc_time_start\": \"%s\",\n", time_start_str);

  fjson_write(fjson, "  \"quaqc_params\": {\n");
  fjson_write(fjson, "    \"sample_n\": %d,\n", params->argc - 1 - params->flag_n);
  fjson_write(fjson, "    \"target_seqs\": %s,\n", strbool(params->tseqs != NULL));
  fjson_write(fjson, "    \"target_seqs_n\": %d,\n", str_hash_size(params->tseqs));
  fjson_write(fjson, "    \"peak_bed\": %s,\n", strbool(params->peaks != NULL));
  fjson_write(fjson, "    \"peak_bed_n\": %d,\n", params->peaks_n);
  fjson_write(fjson, "    \"tss_bed\": %s,\n", strbool(params->tss != NULL));
  fjson_write(fjson, "    \"tss_bed_n\": %d,\n", params->tss_n);
  fjson_write(fjson, "    \"target_list_bed\": %s,\n", strbool(params->tlist != NULL));
  fjson_write(fjson, "    \"target_list_bed_n\": %d,\n", params->tlist_n);
  fjson_write(fjson, "    \"blacklist_bed\": %s,\n", strbool(params->blist != NULL));
  fjson_write(fjson, "    \"blacklist_bed_n\": %d,\n", params->blist_n);
  fjson_write(fjson, "    \"mapq_min\": %d,\n", (int) params->mapq);
  fjson_write(fjson, "    \"alignment_size_min\": %d,\n", (int) params->qlen_min);
  fjson_write(fjson, "    \"alignment_size_max\": %d,\n", (int) params->qlen_max);
  fjson_write(fjson, "    \"fragment_size_min\": %d,\n", (int) params->flen_min);
  fjson_write(fjson, "    \"fragment_size_max\": %d,\n", (int) params->flen_max);
  fjson_write(fjson, "    \"alignment_histogram_max\": %d,\n", params->qhist_max);
  fjson_write(fjson, "    \"fragment_histogram_max\": %d,\n", params->fhist_max);
  fjson_write(fjson, "    \"depth_histogram_max\": %d,\n", params->depth_max);
  fjson_write(fjson, "    \"tss_histogram_size\": %d,\n", params->tss_size);
  fjson_write(fjson, "    \"tss_read_size\": %d,\n", params->tss_qlen);
  fjson_write(fjson, "    \"tss_tn5_shift\": %s,\n", strbool(params->tn5_shift));
  fjson_write(fjson, "    \"use_secondary_alignments\": %s,\n", strbool(params->use_2nd));
  fjson_write(fjson, "    \"use_supplementary_alignments\": %s,\n", strbool(params->use_chi));
  fjson_write(fjson, "    \"use_improper_mates\": %s,\n", strbool(params->use_nomate));
  fjson_write(fjson, "    \"use_duplicates\": %s,\n", strbool(params->use_dups));
  fjson_write(fjson, "    \"no_se\": %s,\n", strbool(params->no_se));
  fjson_write(fjson, "    \"use_dovetails\": %s,\n", strbool(params->use_dovetail));
  fjson_write(fjson, "    \"use_all\": %s,\n", strbool(params->use_all));
  fjson_write(fjson, "    \"no_quaqc_reports\": %s,\n", strbool(params->no_out));
  fjson_write(fjson, "    \"save_as_json\": %s,\n", strbool(params->save));
  fjson_write(fjson, "    \"omit_gc_stats\": %s,\n", strbool(params->omit_gc));
  fjson_write(fjson, "    \"omit_depth_stats\": %s,\n", strbool(params->omit_depth));
  fjson_write(fjson, "    \"fast_mode\": %s,\n", strbool(params->fast));
  fjson_write(fjson, "    \"low_mem_mode\": %s,\n", strbool(params->low_mem));
  fjson_write(fjson, "    \"lenient_mode\": %s,\n", strbool(params->lenient));
  fjson_write(fjson, "    \"nfr_mode\": %s,\n", strbool(params->nfr));
  fjson_write(fjson, "    \"nbr_mode\": %s,\n", strbool(params->nbr));
  fjson_write(fjson, "    \"footprint_mode\": %s,\n", strbool(params->footprint));
  fjson_write(fjson, "    \"chip_mode\": %s,\n", strbool(params->chip));
  fjson_write(fjson, "    \"quit_on_sample_error\": %s,\n", strbool(params->qerr));
  fjson_write(fjson, "    \"verbose\": %s,\n", strbool(params->v));
  fjson_write(fjson, "    \"very_verbose\": %s,\n", strbool(params->vv));
  fjson_write(fjson, "    \"threads_n\": %d\n", params->threads);
  fjson_write(fjson, "  },\n");

  fjson_write(fjson, "  \"quaqc_reports\": [");

  free(CMD);
  return 0;
}

// Need to be able to keep track to know if a sample has to be printed first
// or not, which changes the formatting. The array_first is for the array
// printing macros. All the json functions are behind a mutex in quaqc.c.
static bool json_first = true, array_first = true;

int append_json_fail(char *fn) {
  if (json_first) {
    fjson_write(fjson, "\n    {\n");
    json_first = false;
  } else {
    fjson_write(fjson, ",\n    {\n");
  }
  fjson_write(fjson, "      \"sample\": \"%s\",\n", fn);
  fjson_write(fjson, "      \"status_success\": false,\n");
  fjson_write(fjson, "      \"report\": null\n");
  fjson_write(fjson, "    }");
  return 0;
}

#define fjson_fmt_int_if_first(conn, ii) \
  if (array_first) { \
      fjson_write(conn, " %d", ii); \
      array_first = false; \
    } else { \
      fjson_write(conn, ", %d", ii); \
    } \
  }

#define fjson_array_x(con, array, start, end) do { if (!(start == 0 && end == 0)) { \
  array_first = true; \
  for (int j = 0, i = start; i <= end; j++, i++) { \
    if (array[j] > 0) { \
      fjson_fmt_int_if_first(con, i); \
  } } } while (0)

#define fjson_array_y(con, array, start, end) do { if (!(start == 0 && end == 0)) { \
  array_first = true; \
  for (int j = 0, i = start; i <= end; j++, i++) { \
    if (array[j] > 0) { \
      fjson_fmt_int_if_first(con, array[j]); \
  } } } while (0)

int append_json_result(char *fn, results_t *results, const params_t *params) {
  char time_start_str[128], time_end_str[128];
  get_time(results->time_start, time_start_str);
  get_time(results->time_end, time_end_str);
  int eff_total = 0, eff_nucl = 0, eff_mito = 0, eff_pltd = 0;
  stats_t *tmp1, *tmp2 = results->seqs;
  while (tmp2 != NULL) {
    tmp1 = tmp2;
    tmp2 = tmp2->next;
    if (tmp1->actual > 0) {
      eff_total++;
      switch (tmp1->type) {
        case SEQ_NUCL: eff_nucl++; break;
        case SEQ_MITO: eff_mito++; break;
        case SEQ_PLTD: eff_pltd++; break;
      }
    }
  }

  if (json_first) {
    fjson_write(fjson, "\n    {\n");
    json_first = false;
  } else {
    fjson_write(fjson, ",\n    {\n");
  }

  fjson_write(fjson, "      \"sample\": \"%s\",\n", fn);
  fjson_write(fjson, "      \"status_success\": true,\n");
  fjson_write(fjson, "      \"report\": {\n");
  fjson_write(fjson, "        \"time_start\": \"%s\",\n", time_start_str);
  fjson_write(fjson, "        \"time_end\": \"%s\",\n", time_end_str);
  fjson_write(fjson, "        \"genome_stats\": {\n");
  fjson_write(fjson, "          \"total\": {\n");
  fjson_write(fjson, "            \"genome_wide\": {\n");
  fjson_write(fjson, "               \"seq_n\": %lld,\n", results->seq_n);
  fjson_write(fjson, "               \"seq_size\": %lld\n", results->seq_sum);
  fjson_write(fjson, "            },\n");
  fjson_write(fjson, "            \"nuclear\": {\n");
  fjson_write(fjson, "              \"seq_n\": %lld,\n", results->nuc_n);
  fjson_write(fjson, "              \"seq_size\": %lld\n", results->nuc_sum);
  fjson_write(fjson, "            },\n");
  fjson_write(fjson, "            \"mitochondrial\": {\n");
  fjson_write(fjson, "              \"seq_n\": %lld,\n", results->mito_n);
  fjson_write(fjson, "              \"seq_size\": %lld\n", results->mito_sum);
  fjson_write(fjson, "            },\n");
  fjson_write(fjson, "            \"plastidic\": {\n");
  fjson_write(fjson, "              \"seq_n\": %lld,\n", results->pltd_n);
  fjson_write(fjson, "              \"seq_size\": %lld\n", results->pltd_sum);
  fjson_write(fjson, "            }\n");
  fjson_write(fjson, "          },\n");
  fjson_write(fjson, "          \"effective\": {\n");
  fjson_write(fjson, "            \"genome_wide\": {\n");
  fjson_write(fjson, "              \"seq_n\": %d,\n", eff_total);
  fjson_write(fjson, "              \"seq_size\": %lld\n", results->seq_act);
  fjson_write(fjson, "            },\n");
  fjson_write(fjson, "            \"nuclear\": {\n");
  fjson_write(fjson, "              \"seq_n\": %d,\n", eff_nucl);
  fjson_write(fjson, "              \"seq_size\": %lld\n", results->nuc_act);
  fjson_write(fjson, "            },\n");
  fjson_write(fjson, "            \"mitochondrial\": {\n");
  fjson_write(fjson, "              \"seq_n\": %d,\n", eff_mito);
  fjson_write(fjson, "              \"seq_size\": %lld\n", results->mito_act);
  fjson_write(fjson, "            },\n");
  fjson_write(fjson, "            \"plastidic\": {\n");
  fjson_write(fjson, "              \"seq_n\": %d,\n", eff_pltd);
  fjson_write(fjson, "              \"seq_size\": %lld\n", results->pltd_act);
  fjson_write(fjson, "            }\n");
  fjson_write(fjson, "          }\n");
  fjson_write(fjson, "        },\n");
  fjson_write(fjson, "        \"unfiltered_read_stats\": {\n");
  fjson_write(fjson, "          \"total\": {\n");
  fjson_write(fjson, "            \"reads\": %lld,\n", results->r_total);
  fjson_write(fjson, "            \"unmapped\": %lld,\n", results->r_unmapped);
  fjson_write(fjson, "            \"mapped\": %lld,\n", results->r_mapped);
  fjson_write(fjson, "            \"blacklisted\": %lld\n", results->r_mapped_bl);
  fjson_write(fjson, "          },\n");
  fjson_write(fjson, "          \"effective\": {\n");
  fjson_write(fjson, "            \"mapped_n\": %lld,\n", results->r_seen);
  fjson_write(fjson, "            \"nuclear\": {\n");
  fjson_write(fjson, "              \"n\": %lld,\n", null_get(results->nucl, reads_n));
  fjson_write(fjson, "              \"duplicated\": %lld,\n", null_get(results->nucl, dups_n));
  fjson_write(fjson, "              \"single_end\": %lld,\n", null_get(results->nucl, se_n));
  fjson_write(fjson, "              \"paired_end\": %lld,\n", null_get(results->nucl, pe_n));
  fjson_write(fjson, "              \"properly_mated\": %lld,\n", null_get(results->nucl, mated_n));
  fjson_write(fjson, "              \"primary_alignments\": %lld,\n", null_get(results->nucl, pri_n));
  fjson_write(fjson, "              \"duplicated_primary_alignments\": %lld,\n", null_get(results->nucl, dups_pri_n));
  double nucl_nfr = (double) (null_get(results->nucl, pri_n) - null_get(results->nucl, dups_pri_n));
  if (null_get(results->nucl, pri_n) != 0) nucl_nfr /= (double) null_get(results->nucl, pri_n);
  fjson_write(fjson, "              \"non_redundant_fraction\": %f,\n", nucl_nfr);
  fjson_write(fjson, "              \"secondary_alignments\": %lld,\n", null_get(results->nucl, sec_n));
  fjson_write(fjson, "              \"supplementary_alignments\": %lld\n", null_get(results->nucl, sup_n));
  fjson_write(fjson, "            },\n");
  fjson_write(fjson, "            \"mitochondrial\": {\n");
  fjson_write(fjson, "              \"n\": %lld,\n", null_get(results->mito, reads_n));
  fjson_write(fjson, "              \"duplicated\": %lld,\n", null_get(results->mito, dups_n));
  fjson_write(fjson, "              \"single_end\": %lld,\n", null_get(results->mito, se_n));
  fjson_write(fjson, "              \"paired_end\": %lld,\n", null_get(results->mito, pe_n));
  fjson_write(fjson, "              \"properly_mated\": %lld,\n", null_get(results->mito, mated_n));
  fjson_write(fjson, "              \"primary_alignments\": %lld,\n", null_get(results->mito, pri_n));
  fjson_write(fjson, "              \"duplicated_primary_alignments\": %lld,\n", null_get(results->mito, dups_pri_n));
  double mito_nfr = (double) (null_get(results->mito, pri_n) - null_get(results->mito, dups_pri_n));
  if (null_get(results->mito, pri_n) != 0) mito_nfr /= (double) null_get(results->mito, pri_n);
  fjson_write(fjson, "              \"non_redundant_fraction\": %f,\n", mito_nfr);
  fjson_write(fjson, "              \"secondary_alignments\": %lld,\n", null_get(results->mito, sec_n));
  fjson_write(fjson, "              \"supplementary_alignments\": %lld\n", null_get(results->mito, sup_n));
  fjson_write(fjson, "            },\n");
  fjson_write(fjson, "            \"plastidic\": {\n");
  fjson_write(fjson, "              \"n\": %lld,\n", null_get(results->pltd, reads_n));
  fjson_write(fjson, "              \"duplicated\": %lld,\n", null_get(results->pltd, dups_n));
  fjson_write(fjson, "              \"single_end\": %lld,\n", null_get(results->pltd, se_n));
  fjson_write(fjson, "              \"paired_end\": %lld,\n", null_get(results->pltd, pe_n));
  fjson_write(fjson, "              \"properly_mated\": %lld,\n", null_get(results->pltd, mated_n));
  fjson_write(fjson, "              \"primary_alignments\": %lld,\n", null_get(results->pltd, pri_n));
  fjson_write(fjson, "              \"duplicated_primary_alignments\": %lld,\n", null_get(results->pltd, dups_pri_n));
  double pltd_nfr = (double) (null_get(results->pltd, pri_n) - null_get(results->pltd, dups_pri_n));
  if (null_get(results->pltd, pri_n) != 0) pltd_nfr /= (double) null_get(results->pltd, pri_n);
  fjson_write(fjson, "              \"non_redundant_fraction\": %f,\n", pltd_nfr);
  fjson_write(fjson, "              \"secondary_alignments\": %lld,\n", null_get(results->pltd, sec_n));
  fjson_write(fjson, "              \"supplementary_alignments\": %lld\n", null_get(results->pltd, sup_n));
  fjson_write(fjson, "            }\n");
  fjson_write(fjson, "          }\n");
  fjson_write(fjson, "        },\n");
  fjson_write(fjson, "        \"filtered_read_stats\": {\n");
  /* fjson_write(fjson, "          \"mapped_n\": %lld,\n", results->r_filt); */  // Same as filt_n?
  fjson_write(fjson, "          \"nuclear\": {\n");
  fjson_write(fjson, "            \"alignment\": {\n");
  bool warn_qlen = params->qhist_max < params->qlen_max && !!results->nucl_shared->read_sizes[params->qhist_max];
  fjson_write(fjson, "              \"passing_filters\": %lld,\n", null_get(results->nucl, filt_n));
  fjson_write(fjson, "              \"size_average\": %f,\n", null_get(results->nucl, avg_qlen));
  fjson_write(fjson, "              \"addn_stats_are_suspect\": %s,\n", strbool(warn_qlen));
  fjson_write(fjson, "              \"size_min\": %lld,\n", null_get(results->nucl, min_qlen));
  fjson_write(fjson, "              \"size_1st_pctile\": %lld,\n", null_get(results->nucl, pct1_qlen));
  fjson_write(fjson, "              \"size_sd\": %f,\n", null_get(results->nucl, sd_qlen));
  fjson_write(fjson, "              \"size_99th_pctile\": %lld,\n", null_get(results->nucl, pct99_qlen));
  fjson_write(fjson, "              \"size_max\": %lld,\n", null_get(results->nucl, max_qlen));
  fjson_write(fjson, "              \"size_histogram\": {\n");
  fjson_write(fjson, "                \"range\": [ %d, %d ],\n", 0, params->qhist_max);
  fjson_write(fjson, "                \"x\": [");
  fjson_array_x(fjson, results->nucl_shared->read_sizes, 0, params->qhist_max);
  fjson_write(fjson, " ],\n");
  fjson_write(fjson, "                \"y\": [");
  fjson_array_y(fjson, results->nucl_shared->read_sizes, 0, params->qhist_max);
  fjson_write(fjson, " ]\n");
  fjson_write(fjson, "              }\n");
  fjson_write(fjson, "            },\n");
  fjson_write(fjson, "            \"fragment\": {\n");
  bool warn_flen = params->fhist_max < params->flen_max && !!!results->nucl_shared->frag_sizes[params->fhist_max];
  fjson_write(fjson, "              \"passing_filters\": %lld,\n", null_get(results->nucl, frags_n));
  fjson_write(fjson, "              \"size_average\": %f,\n", null_get(results->nucl, avg_flen));
  fjson_write(fjson, "              \"addn_stats_are_suspect\": %s,\n", strbool(warn_flen));
  fjson_write(fjson, "              \"size_min\": %lld,\n", null_get(results->nucl, min_flen));
  fjson_write(fjson, "              \"size_1st_pctile\": %lld,\n", null_get(results->nucl, pct1_flen));
  fjson_write(fjson, "              \"size_sd\": %f,\n", null_get(results->nucl, sd_flen));
  fjson_write(fjson, "              \"size_99th_pctile\": %lld,\n", null_get(results->nucl, pct99_flen));
  fjson_write(fjson, "              \"size_max\": %lld,\n", null_get(results->nucl, max_flen));
  fjson_write(fjson, "              \"size_histogram\": {\n");
  fjson_write(fjson, "                \"range\": [ %d, %d ],\n", 0, params->fhist_max);
  fjson_write(fjson, "                \"x\": [");
  fjson_array_x(fjson, results->nucl_shared->frag_sizes, 0, params->fhist_max);
  fjson_write(fjson, " ],\n");
  fjson_write(fjson, "                \"y\": [");
  fjson_array_y(fjson, results->nucl_shared->frag_sizes, 0, params->fhist_max);
  fjson_write(fjson, " ]\n");
  fjson_write(fjson, "              }\n");
  fjson_write(fjson, "            },\n");
  fjson_write(fjson, "            \"mapq\": {\n");
  fjson_write(fjson, "              \"no_mapq_n\": %d,\n", results->nucl_shared->mapqs[255]);
  fjson_write(fjson, "              \"score_average\": %f,\n", null_get(results->nucl, avg_mapq2));
  fjson_write(fjson, "              \"score_min\": %lld,\n", null_get(results->nucl, min_mapq));
  fjson_write(fjson, "              \"score_1st_pctile\": %lld,\n", null_get(results->nucl, pct1_mapq));
  fjson_write(fjson, "              \"score_sd\": %f,\n", null_get(results->nucl, sd_mapq));
  fjson_write(fjson, "              \"score_99th_pctile\": %lld,\n", null_get(results->nucl, pct99_mapq));
  fjson_write(fjson, "              \"score_max\": %lld,\n", null_get(results->nucl, max_mapq));
  fjson_write(fjson, "              \"score_histogram\": {\n");
  fjson_write(fjson, "                \"range\": [ 0, 255 ],\n");
  fjson_write(fjson, "                \"x\": [");
  fjson_array_x(fjson, results->nucl_shared->mapqs, 0, 255);
  fjson_write(fjson, " ],\n");
  fjson_write(fjson, "                \"y\": [");
  fjson_array_y(fjson, results->nucl_shared->mapqs, 0, 255);
  fjson_write(fjson, " ]\n");
  fjson_write(fjson, "              }\n");
  fjson_write(fjson, "            },\n");
  double LN_G = (double) null_get(results->nucl, avg_qlen) * (double) null_get(results->nucl, filt_n);
  if (results->nuc_act) LN_G /= results->nuc_act;
  if (params->omit_depth) { // TODO: I could still print LN_G if I wanted to...
    fjson_write(fjson, "            \"depth\": null,\n");
  } else {
    fjson_write(fjson, "            \"depth\": {\n");
    bool warn_depth = !!results->nucl_shared->depths[params->depth_max];
    fjson_write(fjson, "              \"genome_cov\": %f,\n", null_get(results->nucl, genom_cov));
    if (warn_depth) {
      fjson_write(fjson, "              \"read_depth_average\": %f,\n", LN_G);
    } else {
      fjson_write(fjson, "              \"read_depth_average\": %f,\n", null_get(results->nucl, avg_depth));
    }
    fjson_write(fjson, "              \"addn_stats_are_suspect\": %s,\n", strbool(warn_depth));
    fjson_write(fjson, "              \"depths_min\": %lld,\n", null_get(results->nucl, min_depth));
    fjson_write(fjson, "              \"depths_1st_pctile\": %lld,\n", null_get(results->nucl, pct1_depth));
    fjson_write(fjson, "              \"depths_sd\": %f,\n", null_get(results->nucl, sd_depth));
    fjson_write(fjson, "              \"depths_99th_pctile\": %lld,\n", null_get(results->nucl, pct99_depth));
    fjson_write(fjson, "              \"depths_max\": %lld,\n", null_get(results->nucl, max_depth));
    fjson_write(fjson, "              \"depths_histogram\": {\n");
    fjson_write(fjson, "                \"range\": [ %d, %d ],\n", 0, params->fhist_max);
    fjson_write(fjson, "                \"x\": [");
    fjson_array_x(fjson, results->nucl_shared->depths, 0, params->depth_max);
    fjson_write(fjson, " ],\n");
    fjson_write(fjson, "                \"y\": [");
    fjson_array_y(fjson, results->nucl_shared->depths, 0, params->depth_max);
    fjson_write(fjson, " ]\n");
    fjson_write(fjson, "              }\n");
    fjson_write(fjson, "            },\n");
  }
  if (params->omit_gc) {
    fjson_write(fjson, "            \"gc\": null,\n");
  } else {
    fjson_write(fjson, "            \"gc\": {\n");
    fjson_write(fjson, "              \"pct_average\": %f,\n", 100.0 * null_get(results->nucl, gc_pct));
    fjson_write(fjson, "              \"pct_min\": %lld,\n", null_get(results->nucl, min_gc));
    fjson_write(fjson, "              \"pct_1st_pctile\": %lld,\n", null_get(results->nucl, pct1_gc));
    fjson_write(fjson, "              \"pct_sd\": %f,\n", null_get(results->nucl, sd_gc));
    fjson_write(fjson, "              \"pct_99th_pctile\": %lld,\n", null_get(results->nucl, pct99_gc));
    fjson_write(fjson, "              \"pct_max\": %lld,\n", null_get(results->nucl, max_gc));
    fjson_write(fjson, "              \"pct_histogram\": {\n");
    fjson_write(fjson, "                \"range\": [ 0, 100 ],\n");
    fjson_write(fjson, "                \"x\": [");
    fjson_array_x(fjson, results->nucl_shared->read_gc, 0, 100);
    fjson_write(fjson, " ],\n");
    fjson_write(fjson, "                \"y\": [");
    fjson_array_y(fjson, results->nucl_shared->read_gc, 0, 100);
    fjson_write(fjson, " ]\n");
    fjson_write(fjson, "              }\n");
    fjson_write(fjson, "            },\n");
  }
  if (params->peaks == NULL) {
    fjson_write(fjson, "            \"peaks\": null,\n");
  } else {
    double peaks_cov = calc_pct(null_get(results->nucl, peaks_total), null_get(results->nucl, actual));
    fjson_write(fjson, "            \"peaks\": {\n");
    fjson_write(fjson, "              \"n\": %d,\n", null_get(results->nucl, peaks_n));
    fjson_write(fjson, "              \"coverage\": %f,\n", peaks_cov / 100.0);
    fjson_write(fjson, "              \"fraction_of_reads_in_peaks\": %f\n", null_get(results->nucl, frip));
    fjson_write(fjson, "            },\n");
  }
  if (params->tss == NULL) {
    fjson_write(fjson, "            \"tss\": null,\n");
  } else {
    fjson_write(fjson, "            \"tss\": {\n");
    fjson_write(fjson, "              \"n\": %d,\n", null_get(results->nucl, tss_n));
    fjson_write(fjson, "              \"tss_enrichment_score\": %f,\n", null_get(results->nucl, tes));
    fjson_write(fjson, "              \"tss_pileup\": {\n");
    int tss_min = params->tss_size / 2;
    int tss_max = tss_min;
    if (params->tss_size % 2 == 0) tss_min--;
    fjson_write(fjson, "                \"range\": [ %d, %d ],\n", -(tss_min), tss_max);
    fjson_write(fjson, "                \"x\": [");
    fjson_array_x(fjson, results->nucl_shared->tss, -(tss_min), tss_max);
    fjson_write(fjson, " ],\n");
    fjson_write(fjson, "                \"y\": [");
    fjson_array_y(fjson, results->nucl_shared->tss, -(tss_min), tss_max);
    fjson_write(fjson, " ]\n");
    fjson_write(fjson, "              }\n");
    fjson_write(fjson, "            }\n");
    fjson_write(fjson, "          },\n");
  }
  // So, what can I print for these guys:
  // - Basics: reads_n, dups_n, dups_pri_n, se_n, pe_n, un_n, pri_n, sec_n, sup_n, mated_n, filt_n, frags_n
  // - Advanced: nrf, gc_pct, avg_flen, avg_qlen, avg_cov, naive_cov, avg_mapq
  // Just be aware that naive_cov may not be working correctly...
  fjson_write(fjson, "          \"mitochondrial\": {\n");
  fjson_write(fjson, "            \"alignments_passing_filters\": %lld,\n", null_get(results->mito, filt_n));
  fjson_write(fjson, "            \"alignment_size_average\": %f,\n", null_get(results->mito, avg_qlen));
  fjson_write(fjson, "            \"read_depth_average\": %f,\n", null_get(results->mito, avg_cov));
  fjson_write(fjson, "            \"fragments_passing_filters\": %lld,\n", null_get(results->mito, frags_n));
  fjson_write(fjson, "            \"fragment_size_average\": %f,\n", null_get(results->mito, avg_flen));
  fjson_write(fjson, "            \"mapq_score_average\": %f,\n", null_get(results->mito, avg_mapq));
  fjson_write(fjson, "            \"gc_pct_average\": %f\n", 100.0 * null_get(results->mito, gc_pct));
  fjson_write(fjson, "          },\n");
  fjson_write(fjson, "          \"plastidic\": {\n");
  fjson_write(fjson, "            \"alignments_passing_filters\": %lld,\n", null_get(results->pltd, filt_n));
  fjson_write(fjson, "            \"alignment_size_average\": %f,\n", null_get(results->pltd, avg_qlen));
  fjson_write(fjson, "            \"read_depth_average\": %f,\n", null_get(results->pltd, avg_cov));
  fjson_write(fjson, "            \"fragments_passing_filters\": %lld,\n", null_get(results->pltd, frags_n));
  fjson_write(fjson, "            \"fragment_size_average\": %f,\n", null_get(results->pltd, avg_flen));
  fjson_write(fjson, "            \"mapq_score_average\": %f,\n", null_get(results->pltd, avg_mapq));
  fjson_write(fjson, "            \"gc_pct_average\": %f\n", 100.0 * null_get(results->pltd, gc_pct));
  fjson_write(fjson, "          }\n");
  fjson_write(fjson, "        }\n");
  fjson_write(fjson, "      }\n");
  fjson_write(fjson, "    }");
  return 0;
}

int finish_json(void) {
  char time_end_str[128];
  get_time(time(NULL), time_end_str);
  fjson_write(fjson, "\n  ],\n  \"quaqc_time_end\": \"%s\"\n}\n", time_end_str);
  if (use_gz && gzclose(fjson.gz) != Z_OK) {
    int e;
    warn("Failed to close JSON file: %s", gzerror(fjson.gz, &e));
    return 1;
  } else if (!use_gz && fclose(fjson.f) != 0) {
    warn("Failed to close JSON file: %s", strerror(errno));
    return 1;
  }
  return 0;
}

// QC report ----------------------------------------------------------------------

static void print_centered(char *text, int width, FILE *con) {
  int tlen = strlen(text) - 1;
  int pad_left = (width - 1 - tlen) / 2;
  int pad_right = pad_left + (width - 1 - tlen) % 2;
  fprintf(con, "%*.*s%s%*.s", pad_left, pad_left, " ", text, pad_right, " ");
}

static void repeat_wchar(wchar_t wc, int width, FILE *con) {
  for (int i = 0; i < width; i++) {
    fputwc(wc, con);
  }
}

static void print_in_sbox(char *text, int width, FILE *con) {
  fwprintf(con, L"%lc", (wchar_t) 0x250C);
  repeat_wchar((wchar_t) 0x2500, width - 2, con);
  fwprintf(con, L"%lc", (wchar_t) 0x2510);
  fputc('\n', con);
  fwprintf(con, L"%lc", (wchar_t) 0x2502);
  print_centered(text, width - 2, con);
  fwprintf(con, L"%lc", (wchar_t) 0x2502);
  fputc('\n', con);
  fwprintf(con, L"%lc", (wchar_t) 0x2514);
  repeat_wchar((wchar_t) 0x2500, width - 2, con);
  fwprintf(con, L"%lc", (wchar_t) 0x2518);
  fputc('\n', con);
}

void fprint_time(const time_t s, FILE *con) {
  if (s > 7200) {
    fprintf(con, "%'.2f hours\n", ((double) s / 60.0) / 60.0);
  } else if (s > 120) {
    fprintf(con, "%'.2f minutes\n", (double) s / 60.0);
  } else {
    fprintf(con, "%'ld second%s\n", s, pluralize(s));
  }
}

#define bool2str(COND) (COND ? "yes" : "no")

int print_results(char *fn, results_t *results, const params_t *params) {
  char time_start_str[128], time_end_str[128];
  time_t time_diff = difftime(results->time_end, results->time_start);
  get_time(results->time_start, time_start_str);
  get_time(results->time_end, time_end_str);
  // Make the strip a user option?
  char *new_fn = make_new_fn(fn, params->out_dir, params->out_ext, true); 
  if (check_if_matching_fn(fn, new_fn)) {
    error(params->qerr, "QC output cannot have an identical path as input '%s'.", fn);
    free(new_fn);
    return 1;
  }
  char *CMD = paste_strings(params->argv + 1, params->flag_n);
  FILE *fout = fopen(new_fn, "w");
  if (fout == NULL) {
    char errbuf[1024];
    strerror_r(errno, errbuf, 1024);
    warn("Failed to create QC output '%s'.", new_fn);
    warn("%s.", errbuf);
    free(new_fn);
    free(CMD);
    return 1;
  }

  print_centered("quaqc v"QUAQC_VERSION, 80, fout);
  fputc('\n', fout);
  repeat_wchar((wchar_t) 0x2550, 80, fout);
  fputc('\n', fout);
  fputc('\n', fout);

  fprintf(fout, "Input file for this report: %s\n", fn);
  fputc('\n', fout);
  fprintf(fout, "Command arguments used: %s\n", CMD);
  fputc('\n', fout);
  fprintf(fout, "Time of run start: %s\n", time_start_str);
  fprintf(fout, "Time of run end:   %s\n", time_end_str);
  fprintf(fout, "Time ellapsed:     ");
  fprint_time(time_diff, fout);
  fputc('\n', fout);
  fprintf(fout, "Target names specified:                    %s\n", bool2str(params->tseqs != NULL));
  fprintf(fout, "Target BED provided:                       %s\n", bool2str(params->tlist != NULL));
  fprintf(fout, "Blacklist BED provided:                    %s\n", bool2str(params->blist != NULL));
  fprintf(fout, "Include secondary alignments:              %s\n", bool2str(params->use_2nd));
  fprintf(fout, "Include chimeric read alignments:          %s\n", bool2str(params->use_chi));
  fprintf(fout, "Include PE alignments with improper mate:  %s\n", bool2str(params->use_nomate));
  fprintf(fout, "Include duplicate alignments:              %s\n", bool2str(params->use_dups));
  fputc('\n', fout);
  fprintf(fout, "Minimum MAPQ:               %'10d\n", params->mapq);
  fprintf(fout, "Minimum alignment size:     %'10lld\n", params->qlen_min);
  fprintf(fout, "Maximum alignment size:     %'10lld\n", params->qlen_max);
  fprintf(fout, "Minimum fragment size:      %'10lld\n", params->flen_min);
  fprintf(fout, "Maximum fragment size:      %'10lld\n", params->flen_max);
  /* fprintf(fout, "Maximum read depth:         %'10d\n", params->depth_max); */
  fputc('\n', fout);
  fprintf(fout, "Peak count in BED:          %'10d\n", params->peaks_n);
  fprintf(fout, "TSS count in BED:           %'10d\n", params->tss_n);
  fprintf(fout, "Target list count in BED:   %'10d\n", params->tlist_n);
  fprintf(fout, "Blacklist count in BED:     %'10d\n", params->blist_n);
  fputc('\n', fout);
  if (params->qhist_max < params->qlen_max) {
    fprintf(fout, "Warning: Alignment size stats may be incorrect ([--max-qhist] < [--max-qlen])\n");
  }
  if (params->fhist_max < params->flen_max) {
    fprintf(fout, "Warning: Fragment size stats may be incorrect ([--max-fhist] < [--max-flen])\n");
  }
  if (params->qhist_max < params->qlen_max || params->fhist_max < params->flen_max) {
    fputc('\n', fout);
  }

  print_in_sbox("Target genome stats", 80, fout);
  fputc('\n', fout);

  fprintf(fout, "Total sequences:               %'14lld\n", results->seq_n);
  fprintf(fout, "- Nuclear sequences:           %'14lld\n", results->nuc_n);
  fprintf(fout, "- Mitochondrial sequences:     %'14lld\n", results->mito_n);
  fprintf(fout, "- Plastid sequences:           %'14lld\n", results->pltd_n);
  fputc('\n', fout);
  fprintf(fout, "Total genome size:             %'14lld\n", results->seq_sum);
  double nucl_pct = calc_pct(results->nuc_sum, results->seq_sum);
  fprintf(fout, "- Nuclear genome size:         %'14lld (%'.1f%%)\n", results->nuc_sum, nucl_pct);
  double mito_pct = calc_pct(results->mito_sum, results->seq_sum);
  fprintf(fout, "- Mitochondrial genome size:   %'14lld (%'.1f%%)\n", results->mito_sum, mito_pct);
  double pltd_pct = calc_pct(results->pltd_sum, results->seq_sum);
  fprintf(fout, "- Plastid genome size:         %'14lld (%'.1f%%)\n", results->pltd_sum, pltd_pct);
  fputc('\n', fout);
  double nucl_act_pct = calc_pct(results->nuc_act, results->seq_sum);
  fprintf(fout, "Effective nuclear genome size: %'14lld (%'.1f%%)\n", results->nuc_act, nucl_act_pct);
  fputc('\n', fout);

  print_in_sbox("Total read stats", 80, fout);
  fputc('\n', fout);

  fprintf(fout, "Total reads:                   %'14lld\n", results->r_total);
  double map_pct = calc_pct(results->r_mapped, results->r_total);
  fprintf(fout, "- Mapped:                      %'14lld (%'.1f%%)\n", results->r_mapped, map_pct);
  double unm_pct = calc_pct(results->r_unmapped, results->r_total);
  fprintf(fout, "- Unmapped:                    %'14lld (%'.1f%%)\n", results->r_unmapped, unm_pct);
  fputc('\n', fout);
  double eff_n_pct = calc_pct(results->r_seen, results->r_total);
  fprintf(fout, "Effective reads:               %'14lld (%'.1f%%)\n", results->r_seen, eff_n_pct);
  int64_t nucl_n = null_get(results->nucl, reads_n);
  double nucl_n_pct = calc_pct(nucl_n, results->r_seen);
  fprintf(fout, "- Nuclear:                     %'14lld - (%'.1f%%)\n", nucl_n, nucl_n_pct);
  int64_t nucl_dup = null_get(results->nucl, dups_n);
  double nucl_dup_pct = calc_pct(nucl_dup, null_get(results->nucl, reads_n));
  fprintf(fout, "--- Duplicated:                %'14lld --- (%'.1f%%)\n", nucl_dup, nucl_dup_pct);
  int64_t mito_n = null_get(results->mito, reads_n);
  double mito_n_pct = calc_pct(mito_n, results->r_seen);
  fprintf(fout, "- Mitochondrial:               %'14lld - (%'.1f%%)\n", mito_n, mito_n_pct);
  int64_t mito_dup = null_get(results->mito, dups_n);
  double mito_dup_pct = calc_pct(mito_dup, null_get(results->mito, reads_n));
  fprintf(fout, "--- Duplicated:                %'14lld --- (%'.1f%%)\n", mito_dup, mito_dup_pct);
  int64_t pltd_n = null_get(results->pltd, reads_n);
  double pltd_n_pct = calc_pct(pltd_n, results->r_seen);
  fprintf(fout, "- Plastid:                     %'14lld - (%'.1f%%)\n", pltd_n, pltd_n_pct);
  int64_t pltd_dup = null_get(results->pltd, dups_n);
  double pltd_dup_pct = calc_pct(pltd_dup, null_get(results->pltd, reads_n));
  fprintf(fout, "--- Duplicated:                %'14lld --- (%'.1f%%)\n", pltd_dup, pltd_dup_pct);
  int64_t se_total = null_get_all(results, se_n);
  double se_pct = calc_pct(se_total, results->r_seen);
  fprintf(fout, "- SE reads:                    %'14lld - (%'.1f%%)\n", se_total, se_pct);
  int64_t pe_total = null_get_all(results, pe_n);
  double pe_pct = calc_pct(pe_total, results->r_seen);
  fprintf(fout, "- PE reads:                    %'14lld - (%'.1f%%)\n", pe_total, pe_pct);
  int64_t pe_mated = null_get_all(results, mated_n);
  double pe_mated_pct = calc_pct(pe_mated, results->r_seen);
  fprintf(fout, "--- Properly mated:            %'14lld --- (%'.1f%%)\n", pe_mated, pe_mated_pct);
  int64_t pri_total = null_get_all(results, pri_n);
  double pri_pct = calc_pct(pri_total, results->r_seen);
  fprintf(fout, "- Primary alignments:          %'14lld - (%'.1f%%)\n", pri_total, pri_pct);
  int64_t pri_dedup_total = pri_total - null_get_all(results, dups_pri_n);
  double pri_dedup_pct = calc_pct(pri_dedup_total, pri_total);
  fprintf(fout, "--- Non redundant fraction:    %'14lld --- (%'.1f%%)\n", pri_dedup_total, pri_dedup_pct);
  int64_t sec_total = null_get_all(results, sec_n);
  double sec_pct = calc_pct(se_total, results->r_seen);
  fprintf(fout, "- Secondary alignments:        %'14lld - (%'.1f%%)\n", sec_total, sec_pct);
  int64_t sup_total = null_get_all(results, sup_n);
  double sup_pct = calc_pct(sup_total, results->r_seen);
  fprintf(fout, "- Supplementary alignments:    %'14lld - (%'.1f%%)\n", sup_total, sup_pct);
  fputc('\n', fout);

  print_in_sbox("Nuclear read stats after filtering", 80, fout);
  fputc('\n', fout);

  double filt_pct = calc_pct(null_get(results->nucl, filt_n), results->r_seen);
  fprintf(fout, "Reads passing filters:         %'14lld (%'.1f%%)\n", null_get(results->nucl, filt_n), filt_pct);
  fputc('\n', fout);
  bool warn_qlen = params->qhist_max < params->qlen_max && !!results->nucl_shared->read_sizes[params->qhist_max];
  fprintf(fout, "Alignment size min:            %'14lld%s\n", null_get(results->nucl, min_qlen), warn_qlen ? "     (!)" : "");
  fprintf(fout, "Alignment size 1st pctile:     %'14lld%s\n", null_get(results->nucl, pct1_qlen), warn_qlen ? "     (!)" : "");
  fprintf(fout, "Alignment size average:        %'18.3f\n", null_get(results->nucl, avg_qlen));
  fprintf(fout, "Alignment size SD:             %'18.3f%s\n", null_get(results->nucl, sd_qlen), warn_qlen ? " (!)" : "");
  fprintf(fout, "Alignment size 99th pctile:    %'14lld%s\n", null_get(results->nucl, pct99_qlen), warn_qlen ? "     (!)" : "");
  fprintf(fout, "Alignment size max:            %'14lld%s\n", null_get(results->nucl, max_qlen), warn_qlen ? "     (!)" : "");
  fputc('\n', fout);
  if (warn_qlen) {
    fprintf(fout, "(!) Note: Alignment size stats likely have been truncated.\n");
    fprintf(fout, "(!)       Try increasing --max-qhist.\n");
    fprintf(fout, "--> Max filter size:           %'14lld (--max-qlen)\n", params->qlen_max);
    fprintf(fout, "--> Max histogram size:        %'14d (--max-qhist)\n", params->qhist_max);
    fprintf(fout, "----> Max recorded size:       %'14lld\n", null_get(results->nucl, max_qlen));
    fputc('\n', fout);
  }
  fprintf(fout, "Fragments passing filters:     %'14lld\n", null_get(results->nucl, frags_n));
  fputc('\n', fout);
  bool warn_flen = params->fhist_max < params->flen_max && !!!results->nucl_shared->frag_sizes[params->fhist_max];
  fprintf(fout, "Fragment size min:             %'14lld%s\n", null_get(results->nucl, min_flen), warn_flen ? "     (!)" : "");
  fprintf(fout, "Fragment size 1st pctile:      %'14lld%s\n", null_get(results->nucl, pct1_flen), warn_flen ? "     (!)" : "");
  fprintf(fout, "Fragment size average:         %'18.3f\n", null_get(results->nucl, avg_flen));
  fprintf(fout, "Fragment size SD:              %'18.3f%s\n", null_get(results->nucl, sd_flen), warn_flen ? " (!)" : "");
  fprintf(fout, "Fragment size 99th pctile:     %'14lld%s\n", null_get(results->nucl, pct99_flen), warn_flen ? "     (!)" : "");
  fprintf(fout, "Fragment size max:             %'14lld%s\n", null_get(results->nucl, max_flen), warn_flen ? "     (!)" : "");
  fputc('\n', fout);
  if (warn_flen) {
    fprintf(fout, "(!) Note: Fragment size stats likely have been truncated.\n");
    fprintf(fout, "(!)       Try increasing --max-fhist.\n");
    fprintf(fout, "--> Max filter size:           %'14lld (--max-flen)\n", params->flen_max);
    fprintf(fout, "--> Max histogram size:        %'14d (--max-fhist)\n", params->fhist_max);
    fprintf(fout, "----> Max recorded size:       %'14lld\n", null_get(results->nucl, max_flen));
    fputc('\n', fout);
  }
  if (params->omit_depth) {
    // TODO: Fix the cov numbers and use that instead when --omit-depth is used!
    fprintf(fout, "Genome coverage:                          ---\n");
  } else {
    fprintf(fout, "Genome coverage:               %'18.3f%%\n", 100.0 * null_get(results->nucl, genom_cov));
  }
  /* fprintf(fout, "Naive coverage:                %'18.3f%%\n", 100.0 * null_get(results->nucl, naive_cov)); */
  fputc('\n', fout);
  double LN_G = (double) null_get(results->nucl, avg_qlen) * (double) null_get(results->nucl, filt_n);
  if (results->nuc_act) LN_G /= results->nuc_act;
  if (params->omit_depth) {
    fprintf(fout, "Read depth min:                           ---\n");
    fprintf(fout, "Read depth 1st pctile:                    ---\n");
    fprintf(fout, "Read depth average:                       ---\n");
    fprintf(fout, "Read depth average:            %'18.3f\n", LN_G);
    fprintf(fout, "Read depth SD:                            ---\n");
    fprintf(fout, "Read depth 99th pctile:                   ---\n");
    fprintf(fout, "Read depth max:                           ---\n");
    fputc('\n', fout);
    fprintf(fout, "(!) Note: Read depth stats omitted. (--omit-depth)\n");
  } else {
    bool warn_depth = !!results->nucl_shared->depths[params->depth_max];
    fprintf(fout, "Read depth min:                %'14lld%s\n", null_get(results->nucl, min_depth), warn_depth ? "     (!)" : "");
    fprintf(fout, "Read depth 1st pctile:         %'14lld%s\n", null_get(results->nucl, pct1_depth), warn_depth ? "     (!)" : "");
    // I think I should still prefer avg_depth over LN/G when I have it,
    // since LN/G is probably less robust to --target-list problems.
    if (warn_depth) {
      fprintf(fout, "Read depth average:            %'18.3f\n", LN_G);
    } else {
      fprintf(fout, "Read depth average:            %'18.3f\n", null_get(results->nucl, avg_depth));
    }
    fprintf(fout, "Read depth SD:                 %'18.3f%s\n", null_get(results->nucl, sd_depth), warn_depth ? " (!)" : "");
    fprintf(fout, "Read depth 99th pctile:        %'14lld%s\n", null_get(results->nucl, pct99_depth), warn_depth ? "     (!)" : "");
    fprintf(fout, "Read depth max:                %'14lld%s\n", null_get(results->nucl, max_depth), warn_depth ? "     (!)" : "");
    if (warn_depth) {
      fputc('\n', fout);
      fprintf(fout, "(!) Note: Read depth stats likely have been truncated.\n");
      fprintf(fout, "(!)       Try increasing --max-depth.\n");
      fprintf(fout, "--> Max allowed depth:         %'14d (--max-depth)\n", params->depth_max);
      fprintf(fout, "--> Max recorded depth:        %'14lld\n", null_get(results->nucl, max_depth));
    }
  }
  fputc('\n', fout);
  int no_mapq = results->nucl_shared->mapqs[255];
  double no_mapq_pct = calc_pct(no_mapq, results->r_seen);
  fprintf(fout, "Alignment MAPQ min:            %'14lld\n", null_get(results->nucl, min_mapq));
  fprintf(fout, "Alignment MAPQ 1st pctile:     %'14lld\n", null_get(results->nucl, pct1_mapq));
  fprintf(fout, "Alignment MAPQ average:        %'18.3f\n", null_get(results->nucl, avg_mapq2));
  fprintf(fout, "Alignment MAPQ SD:             %'18.3f\n", null_get(results->nucl, sd_mapq));
  fprintf(fout, "Alignment MAPQ 99th pctile:    %'14lld\n", null_get(results->nucl, pct99_mapq));
  fprintf(fout, "Alignment MAPQ max:            %'14lld\n", null_get(results->nucl, max_mapq));
  fprintf(fout, "Alignments without MAPQ:       %'14d (%'.1f%%)\n", no_mapq, no_mapq_pct);
  fputc('\n', fout);
  if (params->omit_gc) {
    fprintf(fout, "Alignment GC min:                         ---\n");
    fprintf(fout, "Alignment GC 1st pctile:                  ---\n");
    fprintf(fout, "Alignment GC average:                     ---\n");
    fprintf(fout, "Alignment GC SD:                          ---\n");
    fprintf(fout, "Alignment GC 99th pctile:                 ---\n");
    fprintf(fout, "Alignment GC max:                         ---\n");
    fputc('\n', fout);
    fprintf(fout, "(!) Note: GC stats omitted. (--omit-gc)\n");
  } else {
    fprintf(fout, "Alignment GC min:              %'14lld%%\n", null_get(results->nucl, min_gc));
    fprintf(fout, "Alignment GC 1st pctile:       %'14lld%%\n", null_get(results->nucl, pct1_gc));
    fprintf(fout, "Alignment GC average:          %'18.3f%%\n", 100.0 * null_get(results->nucl, gc_pct));
    fprintf(fout, "Alignment GC SD:               %'18.3f%%\n", null_get(results->nucl, sd_gc));
    fprintf(fout, "Alignment GC 99th pctile:      %'14lld%%\n", null_get(results->nucl, pct99_gc));
    fprintf(fout, "Alignment GC max:              %'14lld%%\n", null_get(results->nucl, max_gc));
  }
  fputc('\n', fout);
  int peaks_n = params->peaks_n;
  int peaks_nuc_n = null_get(results->nucl, peaks_n);
  double peaks_nuc_pct = calc_pct(peaks_nuc_n, peaks_n);
  fprintf(fout, "Number of peaks:               %'14d (%'.1f%%)\n", peaks_nuc_n, peaks_nuc_pct);
  double peaks_cov = calc_pct(null_get(results->nucl, peaks_total), null_get(results->nucl, actual));
  fprintf(fout, "Peak coverage:                 %18.3f%%\n", peaks_cov);
  if (params->peaks == NULL) {
    fprintf(fout, "Fraction of reads in peaks:               ---\n");
  } else {
    fprintf(fout, "Fraction of reads in peaks:    %'18.3f%%\n", 100.0 * null_get(results->nucl, frip));
  }
  fputc('\n', fout);
  int tss_n = params->tss_n;
  int tss_nuc_n = null_get(results->nucl, tss_n);
  double tss_nuc_pct = calc_pct(tss_nuc_n, tss_n);
  fprintf(fout, "Number of TSSs:                %'14d (%'.1f%%)\n", tss_nuc_n, tss_nuc_pct);
  if (params->tss == NULL) {
    fprintf(fout, "TSS enrichment score:                     ---\n");
  } else {
    fprintf(fout, "TSS enrichment score:          %'18.3f\n", null_get(results->nucl, tes));
  }
  fputc('\n', fout);
  fprintf(fout, "Note: Min/max and percentile stats are always printed as integers.\n");
  fputc('\n', fout);

  /*
  print_in_sbox("Alignment size histogram", 80, fout);
  fputc('\n', fout);

  print_in_sbox("Fragment size histogram", 80, fout);
  fputc('\n', fout);

  print_in_sbox("Alignment MAPQ histogram", 80, fout);
  fputc('\n', fout);

  print_in_sbox("GC content histogram", 80, fout);
  fputc('\n', fout);

  print_in_sbox("Read depth histogram", 80, fout);
  fputc('\n', fout);

  print_in_sbox("TSS read pileup", 80, fout);
  fputc('\n', fout);
  */

  repeat_wchar((wchar_t) 0x2550, 80, fout);
  fputc('\n', fout);
  print_centered("End of quaqc report", 80, fout);
  fputc('\n', fout);

  fclose(fout);
  if (params->v) fprintf(stderr, "|--> Created QC output: %s\n", new_fn);
  free(new_fn);
  free(CMD);
  return 0;
}

