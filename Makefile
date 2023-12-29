CC=g++
DEPS=linalg.h parse.h svp.h

all: linalg.o parse.o svp.o main.cpp
	$(CC) -o runme.exe main.cpp linalg.o parse.o svp.o

%.o: %.cpp $(DEPS)
	$(CC) -c -o $@ $<

clean:
	rm -rf runme.exe test.exe time.exe profile.exe linalg.o parse.o svp.o result.txt

test: linalg.o parse.o svp.o test.cpp
	$(CC) -o test.exe test.cpp linalg.o parse.o svp.o

time: linalg.o parse.o svp.o time.cpp
	$(CC) -o time.exe time.cpp linalg.o parse.o svp.o

profile: linalg.o parse.o profile.cpp
	$(CC) -pg -o profile.exe profile.cpp linalg.o parse.o