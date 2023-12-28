#include "svp.h"
#include <chrono>

using namespace std;

int main(int argc, char** argv){

    Matrix m = Parse(argc, argv);

    size_t samples = 1;
    chrono::microseconds time (0);

    for (size_t i = 0; i < samples; ++i){
        auto start = chrono::high_resolution_clock::now();
        SVP(m);
        auto end = chrono::high_resolution_clock::now();
        time += chrono::duration_cast<chrono::microseconds>(end - start);
    }

    cout << '\n' << "Time taken: " << time.count()/samples << " microseconds" << endl;
    
    return 0;
}