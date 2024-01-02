#include "svp.h"

using namespace std;

int main(int argc, char** argv) {
    Matrix m = Parse(argc, argv);
    SVP(&m);
    return 0;
}
