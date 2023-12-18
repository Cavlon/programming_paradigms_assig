#include <iostream>
#include <string>
#include <cmath>
#include <limits>
#include "linalg.h"
#include "gs.h"
#include "parse.h"

using namespace std;

void LLL(Matrix& m){
    int dim = m.getDim();
    Matrix G = GS(m);

    int k = 1;
    while (k < dim){
        for (int j = k-1; j>=0; j--){
            float mu = GetMu(m(k), G(j));

            // cout << "k=" << k << " " << "j=" << j << " " << "mu=" << mu << endl;

            if (fabsf(mu) > 0.5){
                vector<float> subvec = m(j) * roundf(mu);
                m(k) = m(k) - subvec;
                G = GS(m);
            }
        }

        vector<float> kcheck = G(k);
        vector<float> kchecklow = G(k-1);
        float mu = GetMu(m(k), G(k-1));
        if ((kcheck * kcheck) > ((0.75 - (mu * mu)) * (kchecklow * kchecklow))){
            k++;
        } else {
            m(k).swap(m(k-1));
            G = GS(m);
            k = (k-1 > 1) ? k-1 : 1;
        }
    }
}

int main(int argc, char** argv){
    Matrix m = Parse(argc, argv);
    cout << "Input Matrix:\n";
    Print(m);

    LLL(m);

    cout << "Reduced Matrix:\n";
    Print(m);

    float shortest = std::numeric_limits<float>::max();
    int shrtInd = 0;
    for (int i=0; i<m.getDim(); i++){
        vector<float> basisVec = m(i);
        float sqrnorm = basisVec * basisVec;
        if (sqrnorm < shortest){
            shortest = sqrnorm;
            shrtInd = i;
        }
    }
    cout << "Shortest Vector is:\n";
    Print(m(shrtInd));
    cout << "\nShortest Length: " << sqrtf(shortest);
    return 0;
}