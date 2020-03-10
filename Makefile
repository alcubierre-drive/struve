.PHONY: clean all
all: libstruve.a libstruve.so
test: test.c
	$(CC) $^ -o $@ -std=c99 -Wall -Wextra -pedantic -O3 -lm -lstruve
struve.o: struve.c struve.h Makefile
	$(CC) -c $< -o $@ -std=c99 -Wall -Wextra -pedantic -fPIC -O3
libstruve.so: struve.o
	$(CC) -shared -o $@ $^
libstruve.a: struve.o
	$(AR) rvs $@ $^
clean:
	-$(RM) -f struve.o libstruve.so libstruve.a test
