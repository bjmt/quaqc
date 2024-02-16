/*  bedidx.c -- BED file indexing.

    Copyright (C) 2011 Broad Institute.
    Copyright (C) 2014, 2017-2019 Genome Research Ltd.
    Copyright (C) 2024 Benjamin Jean-Marie Tremblay

    Author: Heng Li <lh3@sanger.ac.uk>

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

/*   Modifications to bedidx.c:
 *
 * - Don't include config.h
 * - Remove end arg from  bed_minoff()
 * - Delete functions unused by quaqc.c
 * - Added new functions
 * - bed_read() now reads the strand column as well
 *
 *
 *     - Benjamin Jean-Marie Tremblay
 */

/* #include <config.h> */

#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <stdio.h>
#include <errno.h>
#include <limits.h>
#include "zlib.h"
#include "bedidx.h"

#include "htslib/ksort.h"

#include "htslib/kseq.h"
KSTREAM_INIT(gzFile, gzread, 8192)

static inline int lt_pair_pos(hts_pair_pos_t a, hts_pair_pos_t b) {
    if (a.beg == b.beg) return a.end < b.end;
    return a.beg < b.beg;
}
KSORT_INIT_STATIC(hts_pair_pos_t, hts_pair_pos_t, lt_pair_pos)

/*! @typedef
 * @abstract bed_reglist_t - value type of the BED hash table
 * This structure encodes the list of intervals (ranges) for the regions provided via BED file or
 * command line arguments.
 * @field *a           pointer to the array of intervals.
 * @field n            actual number of elements contained by a
 * @field m            number of allocated elements to a (n <= m)
 * @field *idx         index array for computing the minimum offset
 */
typedef struct {
    int n, m;
    hts_pair_pos_t *a;
    char *s;
    int *idx;
    int filter;
} bed_reglist_t;

#include "htslib/khash.h"
KHASH_MAP_INIT_STR(reg, bed_reglist_t)

typedef kh_reg_t reghash_t;

#if 0
// Debug function
void bed_print(void *reg_hash) {
    reghash_t *h = (reghash_t *)reg_hash;
    bed_reglist_t *p;
    khint_t k;
    int i;
    const char *reg;

    if (!h) {
        printf("[E::bed_print] Hash table is empty!\n");
        return;
    }
    for (k = kh_begin(h); k < kh_end(h); k++) {
        if (kh_exist(h,k)) {
            reg = kh_key(h,k);
            printf("Region: '%s'\n", reg);
            if ((p = &kh_val(h,k)) != NULL && p->n > 0) {
                printf("Filter: %d\n", p->filter);
                for (i=0; i<p->n; i++) {
                    printf("\tinterval[%d]: %"PRIhts_pos"-%"PRIhts_pos"(%c)\n",
                           i,p->a[i].beg,p->a[i].end,p->s[i]);
                }
            } else {
                printf("Region '%s' has no intervals!\n", reg);
            }
        }
    }
}
#endif

static int *bed_index_core(int n, hts_pair_pos_t *a)
{
    int i, j, l, *idx, *new_idx;
    l = 0; idx = 0;
    for (i = 0; i < n; ++i) {
        hts_pos_t beg, end;
        beg = a[i].beg >> LIDX_SHIFT; end = a[i].end >> LIDX_SHIFT;
        if (l < end + 1) {
            int old_l = l;
            l = end + 1;
            kroundup32(l);
            new_idx = realloc(idx, l * sizeof(*idx));
            if (!new_idx) {
                free(idx);
                return NULL;
            }
            idx = new_idx;

            for (j = old_l; j < l; ++j)
                idx[j] = -1;
        }

        for (j = beg; j < end+1; ++j)
            if (idx[j] < 0)
                idx[j] = i;
    }
    return idx;
}

static void bed_index(void *_h)
{
    reghash_t *h = (reghash_t*)_h;
    khint_t k;
    for (k = 0; k < kh_end(h); ++k) {
        if (kh_exist(h, k)) {
            bed_reglist_t *p = &kh_val(h, k);
            if (p->idx) free(p->idx);
            ks_introsort(hts_pair_pos_t, p->n, p->a);
            p->idx = bed_index_core(p->n, p->a);
        }
    }
}

static void bed_index_without_sort(void *_h)
{
    reghash_t *h = (reghash_t*)_h;
    khint_t k;
    for (k = 0; k < kh_end(h); ++k) {
        if (kh_exist(h, k)) {
            bed_reglist_t *p = &kh_val(h, k);
            if (p->idx) free(p->idx);
            p->idx = bed_index_core(p->n, p->a);
        }
    }
}

static inline void coord_resize(hts_pos_t *s1, hts_pos_t *s2, const hts_pos_t up, const hts_pos_t down, const char s) {
#if 0
    fprintf(stderr, "bedidx.c: Before: %"PRId64"\t%"PRId64"\t%c\n", *s1, *s2, s);
#endif
  if (s == '-') {
    *s1 = down > *s1 ? 0: *s1 - down;
    *s2 += up;
  } else {
    *s1 = up > *s1 ? 0 : *s1 - up;
    *s2 += down;
  }
#if 0
    fprintf(stderr, "bedidx.c: After:  %"PRId64"\t%"PRId64"\t%c\t(-%"PRId64", +%"PRId64")\n", *s1, *s2, s, up, down);
#endif
}

// Resize ranges from the center.
void bed_resize(void *_h, const hts_pos_t up, const hts_pos_t down) {
    reghash_t *h = (reghash_t *)_h;
    bed_reglist_t *p;
    for (khint_t k = kh_begin(h); k < kh_end(h); k++) {
        if (kh_exist(h,k)) {
            if ((p = &kh_val(h,k)) != NULL && p->n > 0) {
                for (int i=0; i<p->n; i++) {
                    coord_resize(&(p->a[i].beg), &(p->a[i].end), up, down, p->s[i]);
                }
            }
        }
    }
    bed_index_without_sort(_h);
}

static int bed_minoff(const bed_reglist_t *p, hts_pos_t beg) {
    int i, min_off=0;

    if (p && p->idx) {
        min_off = (beg>>LIDX_SHIFT >= p->n)? p->idx[p->n-1] : p->idx[beg>>LIDX_SHIFT];
        if (min_off < 0) { // TODO: this block can be improved, but speed should not matter too much here
            hts_pos_t n = beg>>LIDX_SHIFT;
            if (n > p->n)
                n = p->n;
            for (i = n - 1; i >= 0; --i)
                if (p->idx[i] >= 0)
                    break;
            min_off = i >= 0? p->idx[i] : 0;
        }
    }

    return min_off;
}

static int bed_overlap_core(const bed_reglist_t *p, hts_pos_t beg, hts_pos_t end)
{
    int i, min_off;
    if (p->n == 0) return 0;
    min_off = bed_minoff(p, beg);

    for (i = min_off; i < p->n; ++i) {
        if (p->a[i].beg >= end) break; // out of range; no need to proceed
        if (p->a[i].end > beg && p->a[i].beg < end)
            return 1; // find the overlap; return
    }
    return 0;
}

int bed_n(void *_h) {
    int count = 0;
    reghash_t *h = (reghash_t *)_h;
    bed_reglist_t *p;
    for (khint_t k = kh_begin(h); k < kh_end(h); k++) {
        p = &kh_val(h, k);
        if (kh_exist(h, k)) {
            count += p->n;
        }
    }
    return count;
}

int bed_n_in_chr(void *_h, const char *chr) {
    const reghash_t *h = (const reghash_t*)_h;
    if (!h) return 0;
    khint_t k = kh_get(reg, h, chr);
    if (k == kh_end(h)) return 0;
    bed_reglist_t *p = &kh_val(h, k);
    return p->n;
}

int bed_overlap(const void *_h, const char *chr, hts_pos_t beg, hts_pos_t end)
{
    const reghash_t *h = (const reghash_t*)_h;
    khint_t k;
    if (!h) return 0;
    k = kh_get(reg, h, chr);
    if (k == kh_end(h)) return 0;
    return bed_overlap_core(&kh_val(h, k), beg, end);
}

// What about multiple overlapping ranges? Just return offset for first one?
int bed_overlap_offset(void *_h, const char *seq, const hts_pos_t beg, const hts_pos_t end) {
    const reghash_t *h = (const reghash_t*)_h;
    if (!h) return INT_MIN;
    khint_t k = kh_get(reg, h, seq);
    if (k == kh_end(h)) return INT_MIN;
    const bed_reglist_t *p = &kh_val(h, k);
    if (p->n == 0) return INT_MIN;
    int min_off = bed_minoff(p, beg);
    for (int i = min_off; i < p->n; ++i) {
        if (p->a[i].beg >= end) return INT_MIN;
        if (p->a[i].end > beg && p->a[i].beg <= end) {
            return p->s[i] == '-' ? p->a[i].end - end : beg - p->a[i].beg;
        }
    }
    return INT_MIN;
}

/** @brief Trim a sorted interval list, inside a region hash table,
 *   by removing completely contained intervals and merging adjacent or
 *   overlapping intervals.
 *  @param reg_hash    the region hash table with interval lists as values
 *
 *  Note: this function trashes the strand info.
 */

void bed_unify(void *reg_hash) {

    int i, j, new_n;
    reghash_t *h;
    bed_reglist_t *p;

    if (!reg_hash)
        return;

    h = (reghash_t *)reg_hash;

    for (i = kh_begin(h); i < kh_end(h); i++) {
        if (!kh_exist(h,i) || !(p = &kh_val(h,i)) || !(p->n))
            continue;

        for (new_n = 0, j = 1; j < p->n; j++) {
            if (p->a[new_n].end < p->a[j].beg) {
                p->a[++new_n] = p->a[j];
            } else {
                if (p->a[new_n].end < p->a[j].end)
                    p->a[new_n].end = p->a[j].end;
            }
        }

        p->n = ++new_n;

        for (int k = 0; k < p->n; k++) {
          p->s[k] = '.';
        }
    }
}

// Calculate the sum of ranges covering a particular sequence.
// BED coordinates are half open, so I think just (end - beg) is ok for sizes.
hts_pos_t bed_total(void *reg_hash, const char *chr, const hts_pos_t size) {
    reghash_t *h = (reghash_t *)reg_hash;
    if (!h) return 0;
    khint_t k = kh_get(reg, h, chr);
    if (k == kh_end(h)) return 0;
    const bed_reglist_t *p = &kh_val(h, k);
    hts_pos_t s = 0;
    for (int i = 0; i < p->n; i++) {
        if (p->a[i].beg > size) break;
        if (p->a[i].end > size) {
            s += size - p->a[i].beg;
            break;
        }
        s += p->a[i].end - p->a[i].beg;
    }
    return s;
}

/* "BED" file reader, which actually reads two different formats.

   BED files contain between three and nine fields per line, of which
   only the first three (reference, start, end) are of interest to us.
   BED counts positions from base 0, and the end is the base after the
   region of interest.  While not properly documented in the specification,
   it is also possible to have 'browser' and 'track' lines in BED files that
   do not follow the standard format and should be ignored.  Examination
   of the BED file reading code in
   http://genome-source.cse.ucsc.edu/gitweb/?p=kent.git shows that BED
   files can also have comment lines starting with '#', leading whitespace
   is stripped, and that fields are separated by one or more consecutive
   whitespace characters.

   The alternative format was originally for reading positions in VCF
   format.  This expects two columns, which indicate the reference and
   a position.  The position corresponds to a single base, and unlike
   BED counts from 1.

   Which format is in use is determined based on whether one or two
   numbers can be decoded on the line.  As this choice is made line-by-line
   in this implementation, it is possible (but probably a bad idea) to mix
   both formats in the same file.  If trying to read a VCF file by this
   method, it would be important to ensure that the third column (ID) does
   not contain any entries that start with a digit, to avoid the line
   erroneously being parsed as a BED file entry.

   The BED specification is at http://www.genome.ucsc.edu/FAQ/FAQformat.html
   The VCF specification is at https://github.com/samtools/hts-specs
 */

void *bed_read(const char *fn)
{
    reghash_t *h = kh_init(reg);
    gzFile fp;
    kstream_t *ks = NULL;
    int dret;
    unsigned int line = 0, save_errno;
    kstring_t str = { 0, 0, NULL };

    if (NULL == h) return NULL;
    // read the list
    fp = strcmp(fn, "-")? gzopen(fn, "r") : gzdopen(fileno(stdin), "r");
    if (fp == 0) return 0;
    ks = ks_init(fp);
    if (NULL == ks) goto fail;  // In case ks_init ever gets error checking...
    int ks_len;
    while ((ks_len = ks_getuntil(ks, KS_SEP_LINE, &str, &dret)) >= 0) { // read a line
        char *ref = str.s, *ref_end;
        char strand;
        uint64_t beg = 0, end = 0;
        int num = 0;
        khint_t k;
        bed_reglist_t *p;

        if (ks_len == 0)
            continue; // skip blank lines

        line++;
        while (*ref && isspace(*ref)) ref++;
        if ('\0' == *ref) continue;  // Skip blank lines
        if ('#'  == *ref) continue;  // Skip BED file comments
        ref_end = ref;   // look for the end of the reference name
        while (*ref_end && !isspace(*ref_end)) ref_end++;
        if ('\0' != *ref_end) {
            *ref_end = '\0';  // terminate ref and look for start, end
            num = sscanf(ref_end + 1, "%"SCNu64" %"SCNu64" %*s %*s %c",
                         &beg, &end, &strand);
        }
        if (1 == num) {  // VCF-style format
            end = beg--; // Counts from 1 instead of 0 for BED files
        }
        if (num < 1 || end < beg) {
            // These two are special lines that can occur in BED files.
            // Check for them here instead of earlier in case someone really
            // has called their reference "browser" or "track".
            if (0 == strcmp(ref, "browser")) continue;
            if (0 == strcmp(ref, "track")) continue;
            if (num < 1) {
                fprintf(stderr,
                        "[E::bed_read] Parse error reading \"%s\" at line %u\n",
                        fn, line);
            } else {
                fprintf(stderr,
                        "[E::bed_read] Parse error reading \"%s\" at line %u : "
                        "end (%"PRIu64") must not be less "
                        "than start (%"PRIu64")\n",
                        fn, line, end, beg);
            }
            errno = 0; // Prevent caller from printing misleading error messages
            goto fail;
        }
        if (num < 3) strand = '.';

        // Put reg in the hash table if not already there
        k = kh_get(reg, h, ref);
        if (k == kh_end(h)) { // absent from the hash table
            int ret;
            char *s = strdup(ref);
            if (NULL == s) goto fail;
            k = kh_put(reg, h, s, &ret);
            if (-1 == ret) {
                free(s);
                goto fail;
            }
            memset(&kh_val(h, k), 0, sizeof(bed_reglist_t));
        }
        p = &kh_val(h, k);

        // Add begin,end to the list
        if (p->n == p->m) {
            p->m = p->m ? p->m<<1 : 4;
            hts_pair_pos_t *new_a = realloc(p->a, p->m * sizeof(p->a[0]));
            if (NULL == new_a) goto fail;
            char *new_s = realloc(p->s, p->m * sizeof(char));
            if (NULL == new_s) goto fail;
            p->a = new_a;
            p->s = new_s;
        }
        p->s[p->n] = strand;
        p->a[p->n].beg = beg;
        p->a[p->n++].end = end;
    }
    // FIXME: Need to check for errors in ks_getuntil.  At the moment it
    // doesn't look like it can return one.  Possibly use gzgets instead?

    if (gzclose(fp) != Z_OK) {
        fp = NULL;
        goto fail;
    }
    ks_destroy(ks);
    free(str.s);
    bed_index(h);
    /* bed_unify(h); */
    return h;
 fail:
    save_errno = errno;
    if (ks) ks_destroy(ks);
    if (fp) gzclose(fp);
    free(str.s);
    bed_destroy(h);
    errno = save_errno;
    return NULL;
}

void bed_destroy(void *_h)
{
    reghash_t *h;
    khint_t k;

    if (!_h)
        return;

    h = (reghash_t*)_h;
    for (k = 0; k < kh_end(h); ++k) {
        if (kh_exist(h, k)) {
            free(kh_val(h, k).a);
            free(kh_val(h, k).idx);
            free((char*)kh_key(h, k));
        }
    }
    kh_destroy(reg, h);
}

