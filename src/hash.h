#ifndef HASH_H
#define HASH_H

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

#include <stdbool.h>

void *str_hash_init(char *str_arr[], const int str_n, int *dups);
bool str_hash_exists(void *h, const char *str);
void str_hash_add(void *h, char *str_one);
void str_hash_destroy(void *h);
int str_hash_size(void *h);

#endif
