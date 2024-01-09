#include "svp.h"
#include "parse.h"
#include <chrono>

using namespace std;

int main(int argc, char** argv) {
    Matrix m = Parse(argc, argv);
    SVP(&m);

    // size_t samples = 100000;
    // chrono::microseconds time (0);

    // for (size_t i = 0; i < samples; ++i){
    //     auto start = chrono::high_resolution_clock::now();
    //     Matrix m = Parse(argc, argv);
    //     SVP(&m);
    //     auto end = chrono::high_resolution_clock::now();
    //     time += chrono::duration_cast<chrono::microseconds>(end - start);
    // }

    // cout << '\n' << "Time taken: " << time.count()/samples << " microseconds" << endl;
    return 0;
}
