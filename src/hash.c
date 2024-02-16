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

/* unsigned int str_hash_start(void *h) { */
/*   return (unsigned int) kh_begin(h); */
/* } */

unsigned int str_hash_end(void *h) {
  return (unsigned int) kh_end((khash_t(str) *) h);
}

const char *str_hash_key(void *h, const unsigned int i) {
  return kh_key((khash_t(str) *) h, i);
}

bool str_ind_exists(void *h, const unsigned int i) {
  return kh_exist((khash_t(str) *) h, i);
}

void *str_hash_init(char *str_arr[], const int str_n, int *dups) {;
  khash_t(str) *h = kh_init(str);
  int absent;
  *dups = 0;
  for (int i = 0; i < str_n; i++) {
    kh_put(str, h, str_arr[i], &absent);
    *dups += absent == 0;
  }
  return (void *) h;
}

int str_hash_add(void *h, char *str_one) {
  int ret = 1;
  if (h != NULL) {
    kh_put(str, (khash_t(str) *) h, str_one, &ret);
  }
  return ret;
}

bool str_hash_exists(void *h, const char *str) {
  return h != NULL && kh_get(str, (khash_t(str) *) h, str) != kh_end((khash_t(str) *) h);
}

void str_hash_free_and_destroy(void *h) {
  if (h != NULL) {
    khash_t(str) *h_ = (khash_t(str) *) h;
    for (khint_t i = 0; i < kh_end(h_); i++) {
      if (kh_exist(h_, i)) free((char *) kh_key(h_, i));
    }
    kh_destroy(str, h_);
  }
}

void str_hash_destroy(void *h) {
  if (h != NULL) kh_destroy(str, (khash_t(str) *) h);
}

int str_hash_size(void *h) {
  return h == NULL ? 0 : (int) kh_size((khash_t(str) *) h);
}
