#   quaqc: Quick Unassuming Atac-seq Quality Control
#   Copyright (C) 2023  Benjamin Jean-Marie Tremblay
#
#   This program is free software: you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation, either version 3 of the License, or
#   (at your option) any later version.
#
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with this program.  If not, see <https://www.gnu.org/licenses/>.

# User-adjustable variables
CC      ?=cc
CFLAGS  +=-std=gnu99
LDLIBS  +=-lm -lpthread
ZDIR    ?=./libs/zlib
HTSDIR  ?=./libs/htslib
PREFIX  ?=/usr/local
BINDIR  ?=bin
MANDIR  ?=share/man

debug: CFLAGS+=-g3 -Og -Wall -Wextra -Wdouble-promotion -Wno-sign-compare \
	-fsanitize=address,undefined -fno-omit-frame-pointer -DDEBUG
debug: quaqc

profile: CFLAGS+=-fprofile-instr-generate -fcoverage-mapping
profile: release

release: CFLAGS+=-DNDEBUG -O3 -flto
release: LDFLAGS+=-dead_strip
release: quaqc

release-full: libz libhts release

ZOPS=
ZLIB=

HTSOPS=
HTSLIB=

# TODO: What happens if libhts is built with libdeflate? Do I need to add something to LDLIBS?
# HTSOPS+=--without-libdeflate  # To disable

ifneq ($(native),)
	CFLAGS+=-march=native
	HTSOPS+=CFLAGS=-march=native
endif

ifeq ($(with_lzma),)
	HTSOPS+=--disable-lzma
else
	LDLIBS+=-llzma
endif

ifeq ($(with_bz2),)
	HTSOPS+=--disable-bz2
else
	LDLIBS+=-lbz2
endif

ifeq ($(with_curl),)
	HTSOPS+=--disable-libcurl
else
	LDLIBS+=-lcurl
## Need to test if this is necessary on linux
# ifeq ($(shell uname -s),Darwin)
#         LDLIBS+=-lssl -lcrypto
# endif
endif

ifeq ($(z_dyn),)
	HTSOPS+=LDFLAGS=-L$(ZDIR) CPPFLAGS=-I$(ZDIR)
	ZLIB=$(ZDIR)/libz.a
else
	LDLIBS+=-lz
endif

ifeq ($(hts_dyn),)
	HTSLIB=$(HTSDIR)/libhts.a
else
	LDLIBS+=-lhts
endif

libz/libz.a:
	(cd $(ZDIR) && ./configure --prefix=./ --static)
	$(MAKE) -C $(ZDIR)

libz: libz/libz.a

libhts/libhts.a: 
	(cd $(HTSDIR) && ./configure $(HTSOPS))
	$(MAKE) -C $(HTSDIR) lib-static

libhts: libhts/libhts.a

clean/libz:
	$(MAKE) -C $(ZDIR) clean

clean/libhts:
	$(MAKE) -C $(HTSDIR) clean

clean/quaqc:
	-rm -f ./src/*.o
	-rm -f ./quaqc

clean: clean/libz clean/libhts clean/quaqc

mostlyclean/quaqc:
	-rm -f ./src/*.o

mostlyclean: clean/libz clean/libhts mostlyclean/quaqc

src/%.o: src/%.c
	$(CC) $(CFLAGS) -c $^ -o $@

objects := $(patsubst %.c,%.o,$(wildcard src/*.c))

quaqc: $(objects)
	$(CC) $(CFLAGS) $(LDFLAGS) $(LDLIBS) $(objects) -o $@ $(HTSLIB) $(ZLIB)

install: quaqc
	install -p ./quaqc $(PREFIX)/$(BINDIR)

uninstall:
	-rm -f $(PREFIX)/$(BINDIR)/quaqc

