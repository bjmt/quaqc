#ifndef DEPTH_H
#define DEPTH_H

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

#include <stdint.h>
#include "htslib/hts.h"
#include "htslib/sam.h"
#include "quaqc.h"
#include "zlib.h"

void *init_depths(void);
void add_read_to_depths(const bam1_t *aln, const hts_pos_t qlen, void *depths, int32_t *hist, const int depths_max);
void purge_and_reset_depths(void *depths, int32_t *hist, const int depths_max);
void destroy_depths(void *depths);
void *init_bedGraph(const params_t *params);
void purge_and_reset_bedGraph(gzFile bgfile, void *bedGraph, const char *chr);
void add_read_to_bedGraph(gzFile bgfile, void *bedGraph, const hts_pos_t qbeg, const hts_pos_t qend, const char *chr);

#endif
