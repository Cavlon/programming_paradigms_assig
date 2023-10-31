CC=gcc
DEPS=structs.h

all:  gs.o parse.c
	$(CC) -o parse parse.c gs.o

%.o: %.c $(DEPS)
	$(CC) -c -o $@ $<

clean:
	rm -rf out