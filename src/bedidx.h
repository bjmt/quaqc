/*  bedidx.h -- BED file indexing header file.

    Copyright (C) 2017 Genome Research Ltd.
    Copyright (C) 2024 Benjamin Jean-Marie Tremblay

    Author: Valeriu Ohan <vo2@sanger.ac.uk>

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
DEALINGS IN THE SOFTWARE.  */

/*   Modifications to bedidx.h:
 *
 * - Make bed_print() available
 * - Delete functions unused by quaqc.c
 * - Added new functions
 *
 *
 *     - Benjamin Jean-Marie Tremblay
 */

#ifndef BEDIDX_H
#define BEDIDX_H

#include "htslib/hts.h"

#define LIDX_SHIFT 13
#define ALL 0
#define FILTERED 1

void *bed_read(const char *fn);
void bed_destroy(void *_h);
void bed_unify(void *reg_hash);
void bed_resize(void *_h, const hts_pos_t up, const hts_pos_t down);
int bed_overlap_within(const void *_h, const char *chr, hts_pos_t beg, hts_pos_t end);
int bed_overlap(const void *_h, const char *chr, hts_pos_t beg, hts_pos_t end);
int bed_overlap_offset(void *_h, const char *seq, const hts_pos_t beg, const hts_pos_t end);
hts_pos_t bed_total(void *reg_hash, const char *chr, const hts_pos_t size);
int bed_n_in_chr(void *_h, const char *chr);
int bed_n(void *_h);
int bed_remove_overlaps(void *reg_hash);

#if 0
void bed_print(void *reg_hash);
#endif

#endif
