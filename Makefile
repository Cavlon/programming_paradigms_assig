CC=g++
DEPS=linalg.h parse.h svp.h

all: linalg.o parse.o svp.o main.cpp
	$(CC) -o runme main.cpp linalg.o parse.o svp.o

%.o: %.cpp $(DEPS)
	$(CC) -c -o $@ $<

clean:
	rm -rf runme.exe test.exe time.exe linalg.o parse.o svp.o

test: linalg.o parse.o svp.o test.cpp
	$(CC) -o test test.cpp linalg.o parse.o svp.o

time: linalg.o parse.o svp.o time.cpp
	$(CC) -o time time.cpp linalg.o parse.o svp.o

profile: linalg.o parse.o svp.o main.cpp
	$(CC) -pg -o runme main.cpp linalg.o parse.o svp.o