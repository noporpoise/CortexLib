LIBS_PATH=../

ifdef DEBUG
	FLAGS=-DDEBUG=1 --debug
endif

STRING_BUF_PATH := $(LIBS_PATH)/string_buffer

all:
	gcc -o cortex_test $(FLAGS) -Wall -lz -I $(STRING_BUF_PATH) \
	  cortex_test.c cortex.c $(STRING_BUF_PATH)/string_buffer.c

clean:
	if test -e cortex_test; then rm cortex_test; fi
	if test -e cortex_test.dSYM; then rm -r cortex_test.dSYM; fi
	if test -e cortex_test.greg; then rm cortex_test.greg; fi
