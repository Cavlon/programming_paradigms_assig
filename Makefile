CC=gcc

all:  program.c
	$(CC) -o out program.c

%.o: %.c $(DEPS)
	$(CC) -c -o $@ $<

clean:
	rm -rf out