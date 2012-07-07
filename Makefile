LIBS_PATH=$(HOME)/c/libs

STRING_BUF_PATH := $(LIBS_PATH)/string_buffer

ifdef DEBUG
	CFLAGS := -DDEBUG=1 --debug
else
	CFLAGS := -O3
endif

ifeq ($(UNAME), Darwin)
	CFLAGS := $(CFLAGS) -fnested-functions
endif

CFLAGS := $(CFLAGS) -Wall -Wextra -I$(STRING_BUF_PATH) -L$(STRING_BUF_PATH)
LIBFLAGS := -lstrbuf -lz

all:
	gcc $(CFLAGS) -o cortex.o -c cortex.c
	ar -csru libcortex.a cortex.o
	gcc $(CFLAGS) $(LIBFLAGS) -o cortex_test cortex_test.c cortex.o

clean:
	if test -e cortex.o; then rm cortex.o; fi
	if test -e libcortex.a; then rm libcortex.a; fi
	if test -e cortex_test; then rm cortex_test; fi
	if test -e cortex_test.dSYM; then rm -r cortex_test.dSYM; fi
	if test -e cortex_test.greg; then rm cortex_test.greg; fi
