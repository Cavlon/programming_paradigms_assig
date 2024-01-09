#include <cassert>
#include <cstring>
#include <fstream>
#include "svp.h"
#include "parse.h"

using namespace std;

Matrix StringToMatrix(string input, int args);
void Test(string input, int args);

int main(int argc, char** argv) {
    cout << "Testing" << endl;

    Test("[1 0 0] [0 1 0] [0 0 1]", 10);    // Basic Execution Check

    Test("[1 1 1] [-1 0 2] [3 5 6]", 10);   // Simple Lattice Chack

    return 0;
}

// Run the program and check for correctness
void Test(string input, int args) {
    Matrix m = StringToMatrix(input, args);

    SVP(&m);

    ifstream inFile("result.txt");

    if (!inFile.is_open()) {
        throw runtime_error("Error opening the result file");
    }

    double res;
    inFile >> res;

    inFile.close();

    assert(res == 1);
    cout << input << " Successful\n" << endl;
}

// Simulate terminal input using a string
Matrix StringToMatrix(string input, int args) {
    // Simulated argument array
    char** argvals = new char*[args];
    // Indices in the string to split arguments
    size_t start = 0, end = 0;

    size_t index = 1;
    // Finds the next whitespace
    while ((end = input.find(' ', start)) != string::npos) {
        // Splits the string into an argument
        string argument = input.substr(start, end - start);
        // Cast the string into a char array
        if (!argument.empty()) {
            size_t len = argument.length();
            argvals[index] = new char[len];
            // Copy the argument into the argument array
            snprintf(argvals[index], len+1, "%s", argument.c_str());
        }
        start = end + 1;
        ++index;
    }

    // There won't be a space at the end of the string
    string argument = input.substr(start, input.size() - start);
    if (!argument.empty()) {
        size_t len = argument.length();
        argvals[index] = new char[len];
        snprintf(argvals[index], len+1, "%s", argument.c_str());
    }

    return Parse(args, argvals);
}
