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

// Consider getting rid of --use-chimeric?
//
// > One nucleosome = (145-)147 bp
// > Linker DNA between nucleosomes = 5-60 bp (avg 33 bp) (other sources say 20-90?)
// > Therefore, one nucleosome + linker = 180 bp (360, 540, 720, 900)
//
// A note on coordinates:
// - BED is 0-based, half open [0, 1)
// - SAM is 1-based, closed [1, 1]
// - BAM is 0-based, half open [0, 1) (at least that's what bam_endpos() returns)
//
// Additional metrics ideas.
// - Base composition along reads
// - Average base qualities along reads
// - C->T/G->A conversion rates along reads in WGBS, EMS experiments
//   + Or just any variant really

#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>
#include <stdint.h>
#include <getopt.h>
#include <locale.h>
#include <time.h>
#include <libgen.h>
#include <sys/resource.h>
#include <sys/stat.h>
#include <pthread.h>
#include "htslib/hts.h"
#include "htslib/sam.h"
#include "htslib/kstring.h"
#include "htslib/kseq.h"
#include "bedidx.h"
#include "quaqc.h"
#include "bam.h"
#include "output.h"
#include "hash.h"
#include "zlib.h"

KSTREAM_INIT(gzFile, gzread, 8192)

// Thread stuff ----------------------------------------------------------------------

static int                sample_count        = 0;
static int                sample_success      = 0;
static int               *threads_ind;
static int               *argv_ind;
static pthread_t         *threads;
static pthread_mutex_t    json_lock           = PTHREAD_MUTEX_INITIALIZER;
static pthread_mutex_t    sample_success_lock = PTHREAD_MUTEX_INITIALIZER;

// Command line arguments ------------------------------------------------------------

enum opts_enum {
  MITOCHONDRIA  = 'm',
  PLASTIDS      = 'p',
  TARGET_NAMES  = 'n',
  TARGET_LIST   = 't',
  BLACKLIST     = 'b',
  USE_SECONDARY = '2',
  MAPQ          = 'q',
  OUTPUT_DIR    = 'o',
  FAST          = 'f',
  HELP          = 'h',
  CONTINUE      = 'c',
  VERBOSE       = 'v',
  KEEP_FILTERED = 'S',
  KEEP_DIR      = 'k',
  KEEP_EXT      = 'K',
  THREADS       = 'j',
  JSON          = 'J',
  NO_OUTPUT     = '0',
  PEAKS         = 'P',
  TSS           = 'T',
  USE_NOMATE    = 'N',
  USE_DUPS      = 'd',
  OUTPUT_EXT    = 'O',
  USE_ALL       = 'A',
  USE_DOVETAILS = 'D',
  LOW_MEM       = 'L',
  TITLE         = 'i',
  RG_NAMES      = 'r',
  RG_LIST       = 'R',

  MIN_QLEN      = 1000,
  MIN_FLEN,
  MAX_QLEN,
  MAX_FLEN,
  MAX_DEPTH,
  MAX_QHIST,
  MAX_FHIST,
  USE_CHIMERIC,
  OMIT_GC,
  OMIT_DEPTH,
  TSS_SIZE,
  TSS_QLEN,
  TN5_SHIFT,
  LENIENT,
  NFR,
  NBR,
  NO_SE,
  FOOTPRINT,
  CHIP,
  VERSION,
};

static const char *opts_short = "m:p:n:t:b:2q:o:t:Sk:K:j:fhcvJ:0P:T:NDO:AdLi:r:R:";

static struct option opts_long[] = {
  { "mitochondria",  required_argument, 0, MITOCHONDRIA  },
  { "plastids",      required_argument, 0, PLASTIDS      },
  { "peaks",         required_argument, 0, PEAKS         },
  { "tss",           required_argument, 0, TSS           },
  { "target-names",  required_argument, 0, TARGET_NAMES  },
  { "target-list",   required_argument, 0, TARGET_LIST   },
  { "blacklist",     required_argument, 0, BLACKLIST     },
  { "use-secondary", no_argument,       0, USE_SECONDARY },
  { "use-chimeric",  no_argument,       0, USE_CHIMERIC  },
  { "use-nomate",    no_argument,       0, USE_NOMATE    },
  { "use-dups",      no_argument,       0, USE_DUPS      },
  { "mapq",          required_argument, 0, MAPQ          },
  { "min-qlen",      required_argument, 0, MIN_QLEN      },
  { "min-flen",      required_argument, 0, MIN_FLEN      },
  { "max-qlen",      required_argument, 0, MAX_QLEN      },
  { "max-flen",      required_argument, 0, MAX_FLEN      },
  { "max-depth",     required_argument, 0, MAX_DEPTH     },
  { "max-qhist",     required_argument, 0, MAX_QHIST     },
  { "max-fhist",     required_argument, 0, MAX_FHIST     },
  { "rg-names",      required_argument, 0, RG_NAMES      },
  { "rg-list",       required_argument, 0, RG_LIST       },
  { "tss-size",      required_argument, 0, TSS_SIZE      },
  { "tss-qlen",      required_argument, 0, TSS_QLEN      },
  { "tss-tn5",       no_argument,       0, TN5_SHIFT     },
  { "use-all",       no_argument,       0, USE_ALL       },
  { "lenient",       no_argument,       0, LENIENT       },
  { "nfr",           no_argument,       0, NFR           },
  { "nbr",           no_argument,       0, NBR           },
  { "footprint",     no_argument,       0, FOOTPRINT     },
  { "no-se",         no_argument,       0, NO_SE         },
  { "use-dovetails", no_argument,       0, USE_DOVETAILS },
  { "output-dir",    required_argument, 0, OUTPUT_DIR    },
  { "output-ext",    required_argument, 0, OUTPUT_EXT    },
  { "no-output",     no_argument,       0, NO_OUTPUT     },
  { "json",          required_argument, 0, JSON          },
  { "keep",          no_argument,       0, KEEP_FILTERED },
  { "keep-dir",      required_argument, 0, KEEP_DIR      },
  { "keep-ext",      required_argument, 0, KEEP_EXT      },
  { "omit-gc",       no_argument,       0, OMIT_GC       },
  { "omit-depth",    no_argument,       0, OMIT_DEPTH    },
  { "threads",       required_argument, 0, THREADS       },
  { "fast",          no_argument,       0, FAST          },
  { "low-mem",       no_argument,       0, LOW_MEM       },
  { "chip",          no_argument,       0, CHIP          },
  { "help",          no_argument,       0, HELP          },
  { "title",         required_argument, 0, TITLE         },
  { "continue",      no_argument,       0, CONTINUE      },
  { "verbose",       no_argument,       0, VERBOSE       },
  { "version",       no_argument,       0, VERSION       },
  { 0,               0,                 0, 0             }
};

static void help(void) {
  printf(
    "quaqc v%s  Copyright (C) %s  Benjamin Jean-Marie Tremblay\n"
    "\n"
    "Usage:  quaqc [options] file1.bam [file2.bam [...]]\n"
    "\n"
    "*  Recommended way to generate BAMs for quaqc during read alignment:   *\n"
    "*    [bowtie2/etc  ... ]  | samtools fixmate -m -u - - \\               *\n"
    "*    | samtools sort -u - | samtools markdup --write-index - file.bam  *\n"
    "\n"
    "Input options:\n"
    " -m, --mitochondria   STR   Mitochondria references names. [%s]\n"
    " -p, --plastids       STR   Plastid reference names. [%s]\n"
    " -P, --peaks          FILE  Peak coordinates in a BED file for FRIP calculation.\n"
    " -T, --tss            FILE  TSS coordinates in a BED file for TSS pileup.\n"
    "\n"
    "Visible read options:\n"
    " -n, --target-names   STR   Only consider reads on target sequences.\n"
    " -t, --target-list    FILE  Only consider reads overlapping ranges in a BED file.\n"
    " -b, --blacklist      FILE  Only consider reads outside ranges in a BED file.\n"
    " -r, --rg-names       STR   Only consider reads with these read groups (RG).\n"
    " -R, --rg-list        FILE  Only consider reads with read groups (RG) in a file.\n"
    "\n"
    "Filter options:\n"
    " -2, --use-secondary        Allow secondary alignments.\n"
    " -N, --use-nomate           Allow PE reads when the mate does not align properly.\n"
    " -d, --use-dups             Allow duplicate reads.\n"
    "     --use-chimeric         Allow supplemental/chimeric alignments.\n"
    " -D, --use-dovetails        Allow dovetailing PE reads.\n"
    "     --no-se                Discard SE reads.\n"
    " -q, --mapq           INT   Min read MAPQ score. [%d]\n"
    "     --min-qlen       INT   Min alignment length. [%d]\n"
    "     --min-flen       INT   Min fragment length. [%d]\n"
    "     --max-qlen       INT   Max alignment length. [%d]\n"
    "     --max-flen       INT   Max fragment length. [%d]\n"
    " -A, --use-all              Do not apply any filters to mapped reads.\n"
    "\n"
    "Stat options:\n"
    "     --max-depth      INT   Max base depth for read depth histogram. [%d]\n"
    "     --max-qhist      INT   Max alignment length for histogram. [--max-qlen]\n"
    "     --max-fhist      INT   Max fragment length for histogram. [--max-flen]\n"
    "     --tss-size       INT   Size of the TSS region for pileup. [%d]\n"
    "     --tss-qlen       INT   Resize reads (centered on the 5p end) for pileup. [%d]\n"
    "     --tss-tn5              Shift 5-prime end coordinates +4/-5 bases for pileup.\n"
    "     --omit-gc              Omit calculation of read GC content.\n"
    "     --omit-depth           Omit calculation of read depths.\n"
    "\n"
    "Presets:\n"
    " -f, --fast                 --omit-gc --omit-depth (~15%% shorter runtime)\n"
    /* " -L, --low-mem              --max-depth=1 --max-qhist=1 --max-fhist=1 --tss-size=1\n" */ // hide for now
    "     --lenient              --use-nomate --use-dups --use-dovetails --mapq=10\n"
    "     --nfr                  --no-se --max-flen=120 --tss-tn5\n"
    "     --nbr                  --no-se --min-flen=150 --max-flen=1000 --tss-qlen=0\n"
    "     --footprint            --tss-qlen=1 --tss-size=501 --tss-tn5\n"
    "     --chip                 --tss-qlen=0 --tss-size=5001 (uses --peaks for pileup)\n"
    "\n"
    "Output options:\n"
    " -o, --output-dir     DIR   Directory to save QC report if not that of input.\n"
    " -O, --output-ext     STR   Filename extension for output files. [%s]\n"
    " -0, --no-output            Suppress creation of output QC reports.\n"
    " -J, --json           FILE  Save combined QC as a JSON file. Use '-' for stdout.\n"
    " -S, --keep                 Save passing nuclear reads in a new BAM file.\n"
    " -k, --keep-dir       DIR   Directory to create new BAM in if not that of input.\n"
    " -K, --keep-ext       STR   Extension of new BAM. [%s]\n"
    "\n"
    "Program options:\n"
    " -j, --threads        INT   Number of worker threads. Max one per sample. [%d]\n"
    " -i, --title          STR   Assign a title to run.\n"
    " -c, --continue             Do not stop when a sample triggers a program error.\n"
    " -v, --verbose              Print progress messages during runtime.\n"
    "     --version              Print the version and exit.\n"
    " -h, --help                 Print this help message.\n"
    , QUAQC_VERSION, QUAQC_YEAR
    , DEFAULT_MITO, DEFAULT_PLTD
    , DEFAULT_MAPQ, DEFAULT_MIN_QLEN, DEFAULT_MIN_FLEN, DEFAULT_MAX_QLEN, DEFAULT_MAX_FLEN
    , DEFAULT_MAX_DEPTH
    , DEFAULT_TSS_SIZE, DEFAULT_TSS_QLEN
    , DEFAULT_OUT_EXT, DEFAULT_BAM_EXT
    , DEFAULT_THREADS
  );
}

// Utilities ------------------------------------------------------------ 

// This allows the use of scientific notation, which atoi cannot handle.
static int parse_num(const char *str) {
  double num = strtod(str, NULL);
  if (errno == ERANGE || num > (double)INT_MAX) {
    quit("Parameter value '%s' is too large to be parsed as an integer.", str);
  }
  return (int)num;
}

static void check_inputs(char **files, const int n) {
  char **rpaths = alloc(sizeof(char *) * n);
  for (int i = 0; i < n; i++) {
    rpaths[i] = realpath(files[i], NULL);
    if (rpaths[i] == NULL) {
      quit("Failed to open file '%s': %s", files[i], strerror(errno));
    }
  }
  int d;
  void *h = str_hash_init(rpaths, 1, &d);
  for (int i = 1; i < n; i++) {
    if (str_hash_exists(h, rpaths[i])) {
      quit("Found duplicate input: '%s'", rpaths[i]);
    } else {
      str_hash_add(h, rpaths[i]);
    }
  }
  str_hash_destroy(h);
  for (int i = 0; i < n; i++) {
    free(rpaths[i]);
  }
  free(rpaths);
}

static bool check_dir(const char *dir) {
  struct stat statbuf;
  return stat(dir, &statbuf) != 0 ? false : !!S_ISDIR(statbuf.st_mode);
}

void *calloc_or_die(size_t size, const char *func_name) {
#ifdef TRACK_ALLOCS
  static size_t alloc_total = 0;
  if (alloc_total == 0) {
    printf("Function\tTotal\tNow\n");
    printf("------------------------------------------------------------\n");
  }
  alloc_total += size;
  printf("%s\t%f MB\t%zu B\n", func_name, (double)alloc_total / 1024.0 / 1024.0, size);
#endif
  void *result = calloc(1, size);
  if (result == NULL) {
    fprintf(stderr, "[E::%s] Out of memory (requested %lu B).\n",
      func_name, size);
    exit(EXIT_FAILURE);
  }
  return result;
}

long get_mem(void) {
  struct rusage mem;
  getrusage(RUSAGE_SELF, &mem);
  long bytes = mem.ru_maxrss;
#ifdef __linux__
  bytes *= 1024; // ru_maxrss is in KBytes on linux, Bytes on macOS
#endif
  return bytes;
}

static void print_mem(void) {
  long bytes = get_mem();
  fprintf(stderr, "Peak memory usage: %'.2f MB\n",
    ((double) bytes / 1024.0) / 1024.0);
}

void print_time(const time_t s) {
  if (s > 7200) {
    fprintf(stderr, "%'.2f hours.\n", ((double) s / 60.0) / 60.0);
  } else if (s > 120) {
    fprintf(stderr, "%'.2f minutes.\n", (double) s / 60.0);
  } else {
    fprintf(stderr, "%'ld second%s.\n", s, pluralize(s));
  }
}

static int str_split(char *str, char **res, const char sep, const int max_size) {
  const char delim[2] = { sep, '\0' };
  int field_n = 0;
  if (strlen(str) == 0) return 0;
  res[field_n++] = strtok(str, delim);
  while (field_n < max_size && res[field_n - 1] != NULL) {
    res[field_n++] = strtok(NULL, delim);
  }
  if (field_n == max_size && res[field_n - 1] != NULL) return -1;
  return field_n - 1;
}

static void *str_split_and_hash(char *str, int *n, int *dups, const char sep) {
  const char delim[2] = { sep, '\0' };
  *n = 0;
  if (strlen(str) == 0) return NULL;
  int d, d_n = 0;
  char *rg;
  void *h = str_hash_init(NULL, 0, &d);
  rg = strtok(str, delim);
  for (int i = 0; i < MAX_HASH_SIZE && rg != NULL; i++) {
    d = str_hash_add(h, rg);
    d_n += d == 0;
    *n += 1;
    rg = strtok(NULL, delim);   
  }
  *dups = d_n;
  if (rg != NULL) {
    *n = -1;
  } else {
    *n -= *dups;
  }
  return h;
}

static void *read_file_and_hash(const char *fn, int *n, int *dups) {
  *n = 0;
  int d, d_n = 0;
  void *h = str_hash_init(NULL, 0, &d);
  gzFile fp = gzopen(fn, "r");
  if (fp == 0) {
    *n = -2;
    return NULL;
  }
  kstring_t str = { 0, 0, NULL };
  kstream_t *ks = ks_init(fp);
  int ks_len, dret, l_n = 0;
  while ((ks_len = ks_getuntil(ks, KS_SEP_LINE, &str, &dret)) >= 0) {
    if (ks_len > 0 && str.s[0] != '\0') {
      l_n++;
      if (l_n > MAX_HASH_SIZE) {
        *n = -1;
        str_hash_free_and_destroy(h);
        return NULL;
      }
      fprintf(stderr, "RG:%s\n", str.s);
      d = str_hash_add(h, strdup(str.s));
      d_n += d == 0;
    }
  }
  ks_destroy(ks);
  free(str.s);
  if (gzclose(fp) != Z_OK) {
    fp = NULL;
    *n = -3;
    str_hash_free_and_destroy(h);
    return NULL;
  }
  *n = l_n;
  if (l_n == 0) {
    str_hash_free_and_destroy(h);
    return NULL;
  } else {
    *dups = d_n;
    *n -= *dups;
    return h;
  }
}

// params_t ----------------------------------------------------------------------

static params_t *params; // quaqc_thread_handler() needs to be able to see this globally

static params_t *init_params(int argc, char **argv) {
  char **argv_clone = alloc((argc + 1) * sizeof(char *));
  for (int i = 0; i < argc; i++) {
    argv_clone[i] = strdup(argv[i]);
  }
  params_t *params  = alloc(sizeof(params_t));
  params->argv      = argv_clone;
  params->argc      = argc;
  params->mapq      = DEFAULT_MAPQ;
  params->qlen_max  = DEFAULT_MAX_QLEN;
  params->flen_max  = DEFAULT_MAX_FLEN;
  params->qlen_min  = DEFAULT_MIN_QLEN;
  params->flen_min  = DEFAULT_MIN_FLEN;
  params->depth_max = DEFAULT_MAX_DEPTH;
  params->tss_size  = DEFAULT_TSS_SIZE;
  params->tss_qlen  = DEFAULT_TSS_QLEN;
  params->qerr      = true;
  params->threads   = DEFAULT_THREADS;
  return params;
}

static void destroy_params(params_t *params) {
  for (int i = 0; i < params->argc + 1; i++) {
    free(params->argv[i]);
  }
  free(params->argv);
  free(params->peaks);
  free(params->tss);
  free(params->blist);
  free(params->tlist);
  str_hash_destroy(params->mito);           // v
  str_hash_destroy(params->pltd);           // -> Keys owned by argv
  str_hash_destroy(params->tseqs);          // ^
  str_hash_free_and_destroy(params->trg);   // --> Owns the memory
  free(params);
}

// globals_t ----------------------------------------------------------------------

static globals_t *init_globals(const params_t *params) {
  globals_t *globals = alloc(sizeof(globals_t));
  check_for_large_i32(params->qhist_max, "--max-qhist");
  globals->read_sizes = alloc(sizeof(int32_t) * (params->qhist_max + 1));
  check_for_large_i32(params->fhist_max, "--max-fhist");
  globals->frag_sizes = alloc(sizeof(int32_t) * (params->fhist_max + 1));
  check_for_large_i32(params->tss_size, "--tss-size");
  globals->tss = alloc(sizeof(int32_t) * params->tss_size);
  check_for_large_i32(params->depth_max, "--max-depth");
  globals->depths = alloc(sizeof(int32_t) * (params->depth_max + 1));
  return globals;
}

static void reset_globals(globals_t *globals, const params_t *params) {
  memset(globals->mapqs, 0, sizeof(int32_t) * 256);
  memset(globals->read_gc, 0, sizeof(int32_t) * 101);
  memset(globals->read_sizes, 0, sizeof(int32_t) * (params->qhist_max + 1));
  memset(globals->frag_sizes, 0, sizeof(int32_t) * (params->fhist_max + 1));
  memset(globals->tss, 0, sizeof(int32_t) * params->tss_size);
  memset(globals->depths, 0, sizeof(int32_t) * (params->depth_max + 1));
}

static void destroy_globals(globals_t *globals) {
  free(globals->depths);
  free(globals->read_sizes);
  free(globals->frag_sizes);
  free(globals->tss);
  free(globals);
}

// results_t ----------------------------------------------------------------------

static results_t **results;  // quaqc_thread_handler() needs to be able to see this globally

static results_t *init_results(const params_t *params) {
  results_t *results = alloc(sizeof(results_t));
  results->nucl_shared = init_globals(params);
  return results;
}

static void reset_results(results_t *results, const params_t *params) {
  stats_t *tmp;
  while (results->seqs != NULL) {
    tmp = results->seqs;
    results->seqs = results->seqs->next;
    free(tmp);
  }
  free(results->nucl);
  free(results->pltd);
  free(results->mito);
  globals_t *globals = results->nucl_shared;
  reset_globals(globals, params);
  memset(results, 0, sizeof(results_t));
  results->nucl_shared = globals;
}

static void destroy_results(results_t *results) {
  if (results == NULL) return;
  stats_t *tmp;
  while (results->seqs != NULL) {
    tmp = results->seqs;
    results->seqs = results->seqs->next;
    free(tmp);
  }
  destroy_globals(results->nucl_shared);
  free(results->nucl);
  free(results->pltd);
  free(results->mito);
  free(results);
}

// main program loop ----------------------------------------------------------------------

static void *quaqc_thread_handler(void *tind) {
  for (int i = 0; i < sample_count; i++) {
    if (*((int *) tind) == threads_ind[i]) {
      char *sample = params->argv[argv_ind[i]];

      htsFile *bam = hts_open(sample, "r");
      if (bam == NULL) {
        error(params->qerr, "Failed to open file '%s'.", sample);
        goto loop_end;
      }

      msg("Starting file: %s [thread#%d]\n", sample, *((int *) tind) + 1);

      bam->idx = sam_index_load(bam, sample);
      if (bam->idx == NULL) {
        msg("Building missing index ...");
        if (sam_index_build(sample, 0) < 0) {
          error(params->qerr, "Failed to build index for file '%s'.", sample);
          warn("Giving up on file '%s'.", sample);
          goto loop_end;
        }
        bam->idx = sam_index_load(bam, sample);
        msg("done.\n");
      }

      if (!bam) {
        error(params->qerr, "Bad input for file '%s'.", sample);
        warn("Giving up on file '%s'.", sample);
        goto loop_end;
      }

      quaqc_run(bam, results[threads_ind[i]], params);

      if (!results[threads_ind[i]]->success) {
        error(params->qerr, "Failed to run QC on file '%s'.", sample);
        goto loop_end;
      }

      if (!params->no_out && print_results(sample, results[threads_ind[i]], params)) {
        error(params->qerr, "Failed to output QC results on file '%s'.", sample);
        goto loop_end;
      }

      pthread_mutex_lock(&sample_success_lock);
      sample_success++;
      pthread_mutex_unlock(&sample_success_lock);

loop_end:

      if (params->json != NULL) {
        pthread_mutex_lock(&json_lock);
        if (results[threads_ind[i]]->success) {
          if (append_json_result(sample, results[threads_ind[i]], params)) {
            error(params->qerr, "Failed to write to JSON for sample '%s'", bam->fn);
            // TODO: Should I give up on creating the JSON and set params->json = NULL? Or quit even if -c is set?
          }
        } else {
          if (append_json_fail(sample)) {
            error(params->qerr, "Failed to write to JSON for sample '%s'", bam->fn);
          }
        }
        pthread_mutex_unlock(&json_lock);
      }

#if 0
      for (int j = 0; j < params->tss_size; j++) {
        printf("%d\t%d\n", j - params->tss_size / 2, results[threads_ind[i]]->nucl_shared->tss[j]);
      }
#endif
#if 0
      for (int j = 0; j <= results[threads_ind[i]]->nucl->max_depth; j++) {
        printf("%d\t%d\n", j, results[threads_ind[i]]->nucl_shared->depths[j]);
      }
#endif
#if 0
      for (int j = 0; j <= params->fhist_max; j++) {
        printf("%d\t%d\n", j, results[threads_ind[i]]->nucl_shared->frag_sizes[j]);
      }
#endif

      reset_results(results[threads_ind[i]], params);
      if (bam != NULL) hts_close(bam);

    } // END: if (*((int *) tind) == threads_ind[i])
  } // END: for (int i = 0; i < sample_count; i++) 
  return NULL;
}

// Program entry ----------------------------------------------------------------------

static int quaqc_main(int argc, char *argv[]) {

  /* fprintf(stderr, "%zu\n", sizeof(stats_t)); */

  setlocale(LC_NUMERIC, "en_US.UTF-8");  // For thousandths sep

  char default_pltd[sizeof(DEFAULT_PLTD)] = DEFAULT_PLTD;
  char default_mito[sizeof(DEFAULT_MITO)] = DEFAULT_MITO;
  char *pltd[1024], *mito[1024];
  int pltd_n = -1, mito_n = -1, tseqs_n, trg_n;
  params = init_params(argc, argv);
  char *peak_bed = NULL;

  int opt, mapq_tmp, d;
  while ((opt = getopt_long(argc, argv, opts_short, opts_long, NULL)) != -1) {
    switch (opt) {
      case MITOCHONDRIA:
        if (params->mito != NULL) {
          quit("--mitochondria has already been set.");
        }
        params->mito = str_split_and_hash(optarg, &mito_n, &d, ',');
        if (mito_n == 0) {
          quit("Failed to parse any mitochondria (--mitochondria).");
        } else if (mito_n == -1) {
          quit("Too many mitochondria specified, max = %'d (--mitochondria).", MAX_HASH_SIZE);
        }
        if (d > 0) warn("Found duplicate names in --mitochondria.");
        break;
      case PLASTIDS:
        if (params->pltd != NULL) {
          quit("--plastids has already been set.");
        }
        params->pltd = str_split_and_hash(optarg, &pltd_n, &d, ',');
        if (pltd_n == 0) {
          quit("Failed to parse any plastids (--plastids).");
        } else if (pltd_n == -1) {
          quit("Too many plastids specified, max = %'d (--plastids).", MAX_HASH_SIZE);
        }
        if (d > 0) warn("Found duplicate names in --plastids.");
        break;
      case PEAKS:
        if (params->peaks != NULL) {
          quit("--peaks has already been set.");
        }
        peak_bed = optarg;
        params->peaks = bed_read(optarg);
        if (params->peaks == NULL) {
          quit("Failed to read --peaks.");
        }
        bed_unify(params->peaks);
        params->peaks_n = bed_n(params->peaks);
        break;
      case TSS:
        if (params->tss != NULL) {
          quit("--tss has already been set.");
        }
        params->tss = bed_read(optarg);
        if (params->tss == NULL) {
          quit("Failed to read --tss.");
        }
        params->tss_n = bed_n(params->tss);
        break;
      case TARGET_NAMES:
        if (params->tseqs != NULL) {
          quit("--target-names has already been set.");
        }
        params->tseqs = str_split_and_hash(optarg, &tseqs_n, &d, ',');
        if (tseqs_n == 0) {
          quit("Failed to parse any target names (--target-names).");
        } else if (tseqs_n == -1) {
          quit("Too many target names specified, max = %'d (--target-names).", MAX_HASH_SIZE);
        }
        if (d > 0) warn("Found duplicate names in --target-names.");
        break;
      case TARGET_LIST:
        if (params->tlist != NULL) {
          quit("--target-list has already been set.");
        }
        params->tlist = bed_read(optarg);
        if (params->tlist == NULL) {
          quit("Failed to read --target-list.");
        }
        bed_unify(params->tlist);
        params->tlist_n = bed_n(params->tlist);
        break;
      case BLACKLIST:
        if (params->blist != NULL) {
          quit("--blacklist has already been set.");
        }
        params->blist = bed_read(optarg);
        if (params->blist == NULL) {
          quit("Failed to read --blacklist.");
        }
        bed_unify(params->blist);
        params->blist_n = bed_n(params->blist);
        break;
      case RG_NAMES:
        if (params->trg != NULL) {
          quit("Target read groups have already been set (--rg-names).");
        }
        params->trg = str_split_and_hash(optarg, &trg_n, &d, ',');
        if (trg_n == 0) {
          quit("Failed to parse any read groups (--rg-names).");
        } else if (trg_n == -1) {
          quit("Too many read groups specified, max = %'d (--rg-names).", MAX_HASH_SIZE);
        } else {
          params->trg_n = trg_n;
        }
        if (d > 0) warn("Found duplicate names in --rg-names.");
        break;
      case RG_LIST:
        if (params->trg != NULL) {
          quit("Target read groups have already been set (--rg-list).");
        }
        params->trg = read_file_and_hash(optarg, &trg_n, &d);
        if (trg_n == 0) {
          quit("Failed to parse any read groups (--rg-list).");
        } else if (trg_n == -1) {
          quit("Too many read groups specified, max = %'d (--rg-list).", MAX_HASH_SIZE);
        } else if (trg_n == -2) {
          quit("Failed to open file (--rg-list).");
        } else if (trg_n == -3) {
          warn("Failed to close file (--rg-list).");
        }
        if (d > 0) warn("Found duplicate names in --rg-list.");
        break;
      case USE_SECONDARY:
        if (params->use_2nd) {
          quit("--use-secondary has already been set.");
        }
        params->use_2nd = true;
        break;
      case USE_CHIMERIC:
        if (params->use_chi) {
          quit("--use-chimeric has already been set.");
        }
        params->use_chi = true;
        break;
      case USE_NOMATE:
        if (params->use_nomate) {
          quit("--use-nomate has already been set.");
        }
        params->use_nomate = true;
        break;
      case USE_DUPS:
        if (params->use_dups) {
          quit("--use-dups has already been set.");
        }
        params->use_dups = true;
        break;
      case MAPQ:
        if (params->mapq != DEFAULT_MAPQ) {
          quit("--mapq has already been set.");
        }
        mapq_tmp = parse_num(optarg);
        if (mapq_tmp < 0 || mapq_tmp > 255) {
          quit("--mapq must be in the range 0-255.");
        }
        params->mapq = (uint8_t) mapq_tmp;
        break;
      case MIN_QLEN:
        if (params->qlen_min != DEFAULT_MIN_QLEN) {
          quit("--min-qlen has already been set.");
        }
        params->qlen_min = parse_num(optarg);
        if (params->qlen_min < 0) {
          quit("--min-qlen cannot be <0.");
        }
        break;
      case MIN_FLEN:
        if (params->flen_min != DEFAULT_MIN_FLEN) {
          quit("--min-flen has already been set.");
        }
        params->flen_min = parse_num(optarg);
        if (params->flen_min < 0) {
          quit("--min-flen cannot be <0.");
        }
        break;
      case MAX_QLEN:
        if (params->qlen_max != DEFAULT_MAX_QLEN) {
          quit("--max-qlen has already been set.");
        }
        params->qlen_max = parse_num(optarg);
        if (params->qlen_max < 1) {
          quit("--max-qlen cannot be <1.");
        }
        break;
      case MAX_FLEN:
        if (params->flen_max != DEFAULT_MAX_FLEN) {
          quit("--max-flen has already been set.");
        }
        params->flen_max = parse_num(optarg);
        if (params->flen_max < 1) {
          quit("--max-flen cannot be <1.");
        }
        break;
      case MAX_DEPTH:
        if (params->depth_max != DEFAULT_MAX_DEPTH) {
          quit("--max-depth has already been set.");
        }
        params->depth_max = parse_num(optarg);
        if (params->depth_max < 1) {
          quit("--max-depth cannot be <1.");
        }
        break;
      case MAX_QHIST:
        if (params->qhist_max != 0) {
          quit("--max-qhist has already been set.");
        }
        params->qhist_max = parse_num(optarg);
        if (params->qhist_max < 1) {
          quit("--max-qhist cannot be <1.");
        }
        break;
      case MAX_FHIST:
        if (params->fhist_max != 0) {
          quit("--max-fhist has already been set.");
        }
        params->fhist_max = parse_num(optarg);
        if (params->fhist_max < 1) {
          quit("--max-fhist cannot be <1.");
        }
        break;
      case TSS_SIZE:
        if (params->tss_size != DEFAULT_TSS_SIZE) {
          quit("--tss-size has already been set.");
        }
        params->tss_size = parse_num(optarg);
        if (params->tss_size < 1) {
          quit("--tss-size cannot be <1.");
        }
        break;
      case TSS_QLEN:
        if (params->tss_qlen != DEFAULT_TSS_QLEN) {
          quit("--tss-qlen has already been set.");
        }
        params->tss_qlen = parse_num(optarg);
        if (params->tss_qlen < 0) {
          quit("--tss-qlen cannot be a negative value.");
        }
        break;
      case TN5_SHIFT:
        if (params->tn5_shift) {
          quit("--tss-tn5 has already been set.");
        }
        params->tn5_shift = true;
        break;
      case USE_ALL:
        if (params->use_all) {
          quit("--use-all has already been set.");
        }
        params->use_all = true;
        break;
      case LENIENT:
        if (params->lenient) {
          quit("--lenient has already been set.");
        }
        params->lenient = true;
        break;
      case NFR:
        if (params->nfr) {
          quit("--nfr has already been set.");
        }
        params->nfr = true;
        break;
      case NBR:
        if (params->nbr) {
          quit("--nbr has already been set.");
        }
        params->nbr = true;
        break;
      case LOW_MEM:
        if (params->low_mem) {
          quit("--low-mem has already been set.");
        }
        params->low_mem = true;
        break;
      case FOOTPRINT:
        if (params->footprint) {
          quit("--footprint has already been set.");
        }
        params->footprint = true;
        break;
      case CHIP:
        if (params->chip) {
          quit("--chip has already been set.");
        }
        params->chip = true;
        break;
      case NO_SE:
        if (params->no_se) {
          quit("--no-se has already been set.");
        }
        params->no_se = true;
        break;
      case USE_DOVETAILS:
        if (params->use_dovetail) {
          quit("--use-dovetail has already been set.");
        }
        params->use_dovetail = true;
        break;
      case OUTPUT_DIR:
        if (params->out_dir != NULL) {
          quit("--output-dir has already been set.");
        }
        params->out_dir = optarg;
        if (!check_dir(params->out_dir)) {
          quit("Incorrect directory provided to --output-dir.");
        }
        break;
      case OUTPUT_EXT:
        if (params->out_ext != NULL) {
          quit("--output-ext has already been set.");
        }
        params->out_ext = optarg;
        if (strlen(params->out_ext) < 1) {
          quit("--output-ext cannot be an empty string.");
        }
        break;
      case NO_OUTPUT:
        if (params->no_out) {
          quit("--no-output has already been set.");
        }
        params->no_out = true;
        break;
      case KEEP_FILTERED:
        if (params->save) {
          quit("--keep-filtered has already been set.");
        }
        params->save = true;
        break;
      case KEEP_DIR:
        if (params->keep_dir != NULL) {
          quit("--keep-dir has already been set.");
        }
        params->keep_dir = optarg;
        if (!check_dir(params->keep_dir)) {
          quit("Incorrect directory provided to --keep-dir.");
        }
        break;
      case KEEP_EXT:
        if (params->keep_ext != NULL) {
          quit("--keep-ext has already been set.");
        }
        params->keep_ext = optarg;
        if (strlen(params->keep_ext) < 1) {
          quit("--keep-ext cannot be an empty string.");
        }
        break;
      case JSON:
        if (params->json != NULL) {
          quit("--json has already been set.");
        }
        params->json = optarg;
        if (strlen(params->json) < 1) {
          quit("--json cannot be an empty string.");
        }
        break;
      case OMIT_GC:
        if (params->omit_gc) {
          quit("--omit-gc has already been set.");
        }
        params->omit_gc = true;
        break;
      case OMIT_DEPTH:
        if (params->omit_depth) {
          quit("--omit-depth has already been set.");
        }
        params->omit_depth = true;
        break;
      case FAST:
        if (params->fast) {
          quit("--fast has already been set.");
        }
        params->omit_gc = true;    // About ~7-8% of the runtime with PE150.
        params->omit_depth = true; // About ~6-7% of the runtime with PE150.
        params->fast = true;       // About ~15-16% of the runtime with PE150.
        break;
      case THREADS:
        if (params->threads != 1) {
          quit("--threads has already been set.");
        }
        params->threads = parse_num(optarg);
        if (params->threads < 1) {
          quit("--threads must be a positive integer.");
        }
        break;
      case TITLE:
        if (params->title != NULL) {
          quit("--title has already been set.");
        }
        params->title = optarg;
        if (strlen(params->title) < 1) {
          quit("--title cannot be an empty string.");
        }
        break;
      case CONTINUE:
        if (!params->qerr) {
          quit("--continue has already been set.");
        }
        params->qerr = false;
        break;
      case VERBOSE:
        if (params->v && !params->vv) {
          params->vv = true;
        } else if (!params->v) {
          params->v = true;
        } else {
          quit("--verbose can only be set twice at most.");
        }
        break;
      case VERSION:
        printf("quaqc v%s\n", QUAQC_VERSION);
        exit(EXIT_SUCCESS);
      case HELP:
        help();
        exit(EXIT_SUCCESS);
      default:
        fputs("Encountered fatal error, exiting. Run quaqc -h for usage.\n", stderr);
        exit(EXIT_FAILURE);
    }
  }

  if (params->fhist_max == 0) {
    params->fhist_max = params->flen_max;
  }
  if (params->qhist_max == 0) {
    params->qhist_max = params->qlen_max;
  }

  if (params->lenient) {
    params->use_dups = true;
    params->use_nomate = true;
    params->mapq = 10;
    params->use_dovetail = true;
  }

  if (params->nbr && params->nfr) {
    quit("--nfr and --nbr cannot be set simultaneously.");
  }
  if (params->nbr && params->footprint) {
  // TODO: Potentially change this if I remove --tss-qlen=0 from --nbr.
    quit("--nbr and --footprint cannot be set simultaneously.");
  }
  if (params->low_mem && params->footprint) {
    quit("--low-mem and --footprint cannot be set simultaneously.");
  }
  if (params->nfr && params->chip) {
    quit("--nfr and --chip cannot be used simultaneously.");
  }
  if (params->nbr && params->chip) {
    quit("--nbr and --chip cannot be used simultaneously.");
  }
  if (params->footprint && params->chip) {
    quit("--footprint and --chip cannot be used simultaneously.");
  }
  if (params->low_mem && params->chip) {
    quit("--low-mem and --chip cannot be used simultaneously.");
  }

  if (params->nfr) {
    params->flen_max = 120;
    params->tn5_shift = true;
  }
  if (params->nbr) {
    params->flen_min = 150;
    params->flen_max = 1000;
    params->tss_qlen = 0;
  }
  if (params->footprint) {
    params->tss_size = 501;
    params->tss_qlen = 1;
    params->tn5_shift = true;
  }
  if (params->low_mem) {
    params->tss_size = 1;
    params->qhist_max = 1;
    params->fhist_max = 1;
    params->depth_max = 1;
  }
  if (params->chip) {
    if (params->tss != NULL) {
      quit("--chip and --tss cannot be used simultaneously.");
    }
    if (peak_bed == NULL) {
      quit("--peaks must be set when using --chip.");
    }
    params->tss = bed_read(peak_bed);
    if (params->tss == NULL) {
      quit("Failed to read --peaks.");
    }
    params->tss_n = bed_n(params->tss);
    params->tss_size = 5001;
    params->tss_qlen = 0;
  }

  if (params->tn5_shift && params->tss_qlen == 0) {
    msg("Warning: --tss-tn5 is ignored when --tss-qlen=0.\n");
  }

  if (params->tlist != NULL && params->blist != NULL) {
    quit("--target-list and --blacklist cannot both be used.");
  }

  if (params->omit_depth) params->depth_max = 0;

  if (params->tss != NULL) {
    int up = params->tss_size, down = params->tss_size;
    up = up / 2; down = down / 2 + down % 2 - 1;
    bed_resize(params->tss, (hts_pos_t) up, (hts_pos_t) down);
  } else {
    params->tss_size = 0;
  }

  if (mito_n == -1) {
    mito_n = str_split(default_mito, mito, ',', 1024);
    int d;
    params->mito = str_hash_init(mito, mito_n, &d);
  }
  if (pltd_n == -1) {
    pltd_n = str_split(default_pltd, pltd, ',', 1024);
    int d;
    params->pltd = str_hash_init(pltd, pltd_n, &d);
  }
  if (params->keep_ext == NULL) params->keep_ext = DEFAULT_BAM_EXT;
  if (params->out_ext == NULL) params->out_ext = DEFAULT_OUT_EXT;

  time_t time_start = time(NULL);

  if (optind == argc) quit("Missing input BAMs.");
  params->flag_n = optind - 1;

  if (params->json != NULL && init_json(params)) {
    quit("Failed to create JSON file.");
  }

  sample_count = argc - optind;
  if (params->threads > sample_count) params->threads = sample_count;
  threads = alloc(sizeof(pthread_t) * params->threads);
  threads_ind = alloc(sizeof(int) * sample_count);
  argv_ind = alloc(sizeof(int) * sample_count);
  for (int i = 0; i < sample_count; i++) {
    threads_ind[i] = (float) i / sample_count * params->threads;
    argv_ind[i] = i + optind;
  }

  check_inputs(argv + optind, sample_count);

  results = alloc(sizeof(results) * params->threads);
  for (int i = 0; i < params->threads; i++) {
    results[i] = init_results(params);
  }

  int *tind = alloc(sizeof(int *) * params->threads);
  for (int i = 0; i < params->threads; i++) {
    tind[i] = i;
    if (pthread_create(&threads[i], NULL, quaqc_thread_handler, &(tind[i]))) {
      quit("Failed to create thread#%d", i + 1);
    }
  }
  for (int i = 0; i < params->threads; i++) {
    if (pthread_join(threads[i], NULL)) {
      quit("Failed to terminate thread#%d", i + 1);
    }
  }
  free(tind);

  if (params->json != NULL && finish_json()) {
    quit("Failed to finish generation of JSON file.");
  }

  if (params->v) {
    time_t time_end = time(NULL);
    fputc('\n', stderr);
      fprintf(stderr, "Done. Processed %'d sample%s (%'d successfully) in ", sample_count, pluralize(sample_count), sample_success);
    print_time(difftime(time_end, time_start));
    if (params->json != NULL) {
      msg("All QC results saved in JSON output: ");
      if (strcmp(params->json, "-") == 0) {
        msg("/dev/stdout\n");
      } else {
        msg("%s\n", params->json);
      }
    }
    print_mem();
  }

  for (int i = 0; i < params->threads; i++) {
    destroy_results(results[i]);
  }
  free(results);
  destroy_params(params);
  free(threads);
  free(threads_ind);
  free(argv_ind);

  return EXIT_SUCCESS;

}

int main(int argc, char *argv[]) {
  return quaqc_main(argc, argv);
}

