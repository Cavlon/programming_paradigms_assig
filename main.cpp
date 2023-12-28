#include "svp.h"
#include <cstring>

using namespace std;

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

int main(int argc, char** argv){

    Matrix m = StringToMatrix("[1 0 0] [0 1 0] [0 0 1]", 10);
    double res = SVP(m);

    return 0;
}