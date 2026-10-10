# Builds the crackle test harness against the grackle library built in src/clib
# (i.e. the working tree, not the installed library). Build grackle first with `make` in src/clib.

GRACKLE_CLIB ?= ../..
CC ?= gcc
CFLAGS ?= -O2 -Wall -DOMIT_LEGACY_INTERNAL_GRACKLE_FUNC

LIBDIR := $(abspath $(GRACKLE_CLIB)/.libs)

crackle_harness: crackle_harness.c $(LIBDIR)/libgrackle.so
	$(CC) $(CFLAGS) -o $@ $< -I$(GRACKLE_CLIB) -L$(LIBDIR) -lgrackle -Wl,-rpath,$(LIBDIR) -lm

clean:
	rm -f crackle_harness

.PHONY: clean
