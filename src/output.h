#ifndef OUTPUT_H
#define OUTPUT_H

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

#include "quaqc.h"
#include "zlib.h"

samFile *init_filtered_bam(bam_hdr_t *hdr, const char *fn, const params_t *params);

int print_results(char *fn, results_t *results, const params_t *params);

int init_json(const params_t *params);

int finish_json(void);

int append_json_fail(char *fn);

int append_json_result(char *fn, results_t *results, const params_t *params);

gzFile init_bedGraph_f(const char *fn, const params_t *params);

#endif
