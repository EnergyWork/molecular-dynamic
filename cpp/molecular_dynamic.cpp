#include "Headers.hpp"
#include "MDSystemClass.hpp"

using namespace std;

#define CLOCKS_PER_SEC 1000
#define STEPS 10000
#define N_ATOMS 25
#define CUBE_SIZE 3
#define DIM 3
#define SPEED 2

int main()
{   
    cout << "Start" << endl;
    clock_t start = clock();
    MDSystem *mdsys = new MDSystem(N_ATOMS, CUBE_SIZE, DIM, SPEED);
    mdsys->init_system(true);
    for (size_t i = 0; i < STEPS; i++) {
        mdsys->calc_forces();
        mdsys->integrate();
        mdsys->print_to_file();
    }
    clock_t end = clock();
    double secs = (double)(end - start) / CLOCKS_PER_SEC;
    cout << "Done! " << secs << " sec" << endl;
}