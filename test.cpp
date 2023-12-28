#include "svp.h"
#include <cassert>
#include <cstring>

using namespace std;

Matrix StringToMatrix(string input, int args);
void Test(string input, int args);

int main(int argc, char** argv){

    cout << "Testing" << endl;

    Test("[1 0 0] [0 1 0] [0 0 1]", 10);    // Basic Execution Check

    Test("[1 1 1] [-1 0 2] [3 5 6]", 10);   // Simple Lattice Chack

    Test("[3 0 0] [2 0 0] [1 1 1]", 10);    // Linear Dependence Check

    // 10-Dimensional Lattice Check
    Test("[1 0 0 0 0 0 0 0 0 0] [0 1 0 0 0 0 0 0 0 0] [0 0 1 0 0 0 0 0 0 0] [0 0 0 1 0 0 0 0 0 0] [0 0 0 0 1 0 0 0 0 0] [0 0 0 0 0 1 0 0 0 0] [0 0 0 0 0 0 1 0 0 0] [0 0 0 0 0 0 0 1 0 0] [0 0 0 0 0 0 0 0 1 0] [0 0 0 0 0 0 0 0 0 1]", 101);
    
    return 0;
}
void Test(string input, int args){
    Matrix m = StringToMatrix(input, args);
    double res = SVP(m);
    assert(res == 1);
    cout << input << " Successful\n" << endl;
}

Matrix StringToMatrix(string input, int args){
    char** argvals = new char*[args];
    size_t start = 0, end = 0;

    size_t index = 1;
    while ((end = input.find(' ', start)) != string::npos) {
        string argument = input.substr(start, end - start);
        if (!argument.empty()){
            argvals[index] = new char[argument.length()];
            strcpy(argvals[index], argument.c_str());
        }
        start = end + 1;
        ++index;
    }

    string argument = input.substr(start, input.size() - start);
    if (!argument.empty()){
        argvals[index] = new char[argument.length()];
        strcpy(argvals[index], argument.c_str());
    }

    return Parse(args, argvals);
}