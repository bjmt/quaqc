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
#include "quaqc.h"
#include "depth.h"

#ifndef kroundup32
#define kroundup32(x) (--(x), (x)|=(x)>>1, (x)|=(x)>>2, (x)|=(x)>>4, (x)|=(x)>>8, (x)|=(x)>>16, ++(x))
#endif

#define INIT_HIST_SIZE            64
#define BEDGRAPH_MAX_READ_SIZE   300

typedef struct depths_t {
  hts_pos_t beg, end;
  int size;
  int32_t *hist;
} depths_t;

#ifdef DEBUG0
static void print_depths(void *depths) {
  depths_t *d = (depths_t *) depths;
  for (int i = 0; i < d->size; i++) {
    fprintf(stderr, "[%d: %*d]", i, 4, d->hist[i]);
  }
  fputc('\n', stderr);
}
#endif
#ifdef DEBUG0
static void print_hist(int32_t *hist, int max) {
  for (int i = 0; i < max; i++) {
    fprintf(stderr, "[%d: %*d]", i, 4, hist[i]);
  }
  fputc('\n', stderr);
}
#endif

void *init_depths(void) {
  depths_t *depths = alloc(sizeof(depths_t));
  depths->size = INIT_HIST_SIZE;
  depths->hist = alloc(depths->size * sizeof(int32_t));
  depths->end = -1;
  return (void *) depths;
}

void *init_bedGraph(const params_t *params) {
  depths_t *depths = alloc(sizeof(depths_t));
  depths->size = params->bg_qlen > 0 ? params->bg_qlen : BEDGRAPH_MAX_READ_SIZE;
  depths->size *= 2;
  depths->size += BEDGRAPH_MAX_READ_SIZE + 3;
  kroundup32(depths->size);
  depths->hist = alloc(depths->size * sizeof(int32_t));
  depths->end = -1;
  return (void *) depths;
}

void destroy_depths(void *depths) {
  if (depths != NULL) {
    depths_t *d = (depths_t *) depths;
    if (d->size && d->hist != NULL) free(d->hist);
    free(d);
  }
}

static inline void grow_depths(depths_t *d, const int size) {
  const int prev = d->size;
  d->size = size;
  kroundup32(d->size);
  int *hist = alloc(d->size * sizeof(int32_t));
  if (hist == NULL) quit("Out of memory.");
  for (int i = 0; i < prev; i++) {
    hist[(d->beg + i) & (d->size - 1)] = d->hist[(d->beg + i) & (prev - 1)];
  }
  free(d->hist);
  d->hist = hist;
}

void purge_and_reset_depths(void *depths, int32_t *hist, const int depths_max) {
  depths_t *d = (depths_t *) depths;
#ifdef DEBUG0
  print_depths(depths);
#endif
  for (int j, i = d->beg; i < d->end; i++) {
    j = i & (d->size - 1);
    hist[min(d->hist[j], depths_max)]++;
    d->hist[j] = 0;
  }
  d->beg = 0;
  d->end = -1;
#ifdef DEBUG0
  print_hist(hist, depths_max);
#endif
}

void purge_and_reset_bedGraph(gzFile bgfile, void *bedGraph, const char *chr) {
  depths_t *b = (depths_t *) bedGraph;
  for (int i = b->beg; i < b->end; i++) {
    if (b->hist[i & (b->size - 1)] != 0) {
      gzprintf(bgfile, "%s\t%d", chr, i);
      while (i + 1 < b->end && b->hist[i & (b->size - 1)] == b->hist[(i + 1) & (b->size - 1)]) i++;
      gzprintf(bgfile, "\t%d\t%d\n", i + 1, b->hist[i & (b->size - 1)]);
    }
  }
  for (int i = 0; i < b->size; i++) {
    b->hist[i] = 0;
  }
  b->beg = 0;
  b->end = -1;
}

void add_read_to_bedGraph(gzFile bgfile, void *bedGraph, const hts_pos_t qbeg0, const hts_pos_t qbeg, const hts_pos_t qend, const char *chr) {
  depths_t *b = (depths_t *) bedGraph;
  if (b->end == -1) {
    b->beg = qbeg0;
    b->end = qend;
  }
  if (qend - qbeg0 > b->size) grow_depths(b, qend - qbeg0);
  if (qbeg0 > b->beg) {
    for (int i = b->beg; i < qbeg0; i++) {
      if (b->hist[i & (b->size - 1)] != 0) {
        gzprintf(bgfile, "%s\t%d", chr, i);
        while (i + 1 < qbeg0 && b->hist[i & (b->size - 1)] == b->hist[(i + 1) & (b->size - 1)]) i++;
        gzprintf(bgfile, "\t%d\t%d\n", i + 1, b->hist[i & (b->size - 1)]);
      }
    }
    for (int i = b->beg; i < qbeg0; i++) {
      b->hist[i & (b->size - 1)] = 0;
    }
    b->beg = qbeg0;
  }
  if (qend > b->end) b->end = qend;
  for (int i = qbeg; i < qend; i++) {
    b->hist[i & (b->size - 1)]++;
  }
}

void add_read_to_depths(const bam1_t *aln, const hts_pos_t qend, void *depths, int32_t *hist, const int depths_max) {
  const hts_pos_t beg = aln->core.pos, qlen = 1 + qend - aln->core.pos;
  const uint32_t *cigar = bam_get_cigar(aln);
  depths_t *d = (depths_t *) depths;
  if (beg > d->beg) {
    // I think this still works when the gap between reads is greater than the size of the
    // buffer, since it is zeroed as it travels along the ring buffer and loops to just zeros.
    for (int j, i = d->beg; i < beg; i++) {
      j = i & (d->size - 1);
      hist[min(d->hist[j], depths_max)]++;
      d->hist[j] = 0; 
    }
    d->beg = beg;
#ifdef DEBUG0
    print_hist(hist, depths_max);
#endif
  }
  if (qlen > d->size) grow_depths(d, qlen);
  if (qend > d->end) d->end = qend;
  static const int query[16] = {
  //M I D N  S H P =  X B ? ?  ? ? ? ?
    1,0,0,0, 0,0,0,1, 1,0,0,0, 0,0,0,0
  };
  for (int k = 0, r = 0, b; k < aln->core.n_cigar; k++) {
    b = bam_cigar_oplen(cigar[k]);
    if (query[bam_cigar_op(cigar[k])]) {
      for (int i = r; i < r + b; i++) {
        d->hist[(d->beg + i) & (d->size - 1)]++;
      }
    }
    r += b;
  }
#ifdef DEBUG0
  print_depths(depths);
#endif
}

