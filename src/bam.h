#ifndef BAM_H
#define BAM_H

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

#include "htslib/hts.h"
#include "htslib/sam.h"
#include "quaqc.h"

void quaqc_run(htsFile *bam, results_t *results, const params_t *params);

#endif
