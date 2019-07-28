.PHONY: clean all
struve.o: struve.c struve.h Makefile
	$(CC) -c struve.c -std=c99 -Wall -Wextra -pedantic -fPIC -O3
libstruve.so: struve.o
	$(CC) -shared -o $@ $^
libstruve.a: struve.o
	$(AR) rvs $@ $^
all: libstruve.a libstruve.so
clean:
	-$(RM) struve.o libstruve.so libstruve.a
