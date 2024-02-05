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
#include "htslib/khash.h"
KHASH_SET_INIT_STR(str)

void *str_hash_init(char *str_arr[], const int str_n, int *dups) {;
  khash_t(str) *h = kh_init(str);
  khint_t k;
  int absent;
  *dups = 0;
  for (int i = 0; i < str_n; i++) {
    k = kh_put(str, h, str_arr[i], &absent);
    *dups += absent == 0;
  }
  return (void *) h;
}

void str_hash_add(void *h, char *str_one) {
  if (h != NULL) {
    int ret;
    kh_put(str, (khash_t(str) *) h, str_one, &ret);
  }
}

bool str_hash_exists(void *h, const char *str) {
  return h != NULL && kh_get(str, (khash_t(str) *) h, str) != kh_end((khash_t(str) *) h);
}

void str_hash_destroy(void *h) {
  if (h != NULL) kh_destroy(str, (khash_t(str) *) h);
}

int str_hash_size(void *h) {
  return h == NULL ? 0 : (int) kh_size((khash_t(str) *) h);
}
