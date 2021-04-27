#include <iostream>
#include <vector>
#include <string>
#include <time.h>
#include "MDSystemClass.hpp"

using namespace std;

#define CLOCKS_PER_SEC 1000
#define STEPS 1000

int main()
{
    clock_t start = clock();
    MDSystem *mdsys = new MDSystem(10, 3, 3, 2); //n_atoms=10, cube_size=3, dim=3, speed=2
    mdsys->init_system(false);
    for (size_t i = 0; i < STEPS; i++) {
        mdsys->calc_forces();
        mdsys->integrate();
        mdsys->print_to_file();
    }
    clock_t end = clock();
    double secs = (double)(end - start) / CLOCKS_PER_SEC;
    cout << "Done! " << secs << " sec" << endl;
}