#include <iostream>
#include <vector>
#include <string>
#include <time.h>
#include "MDSystemClass.hpp"

using namespace std;

#define CLOCKS_PER_SEC 1000
#define STEPS 10

void print_vector(vector<double> vec)
{
    for (auto el: vec)
        cout << el << " ";
    cout << endl;
}

void test(MDSystem *sys)
{
    auto r1 = vector<double>{3, 4, 5};
    auto r2 = vector<double>{6, 2, 4};
    cout << sys->lenght(r1) << endl;
    cout << sys->force_LD(0.78) << endl;
    sys->normalize(r1);
    print_vector(r1);
}

int main()
{
    cout << "Start" << endl;
    clock_t start = clock();
    MDSystem *mdsys = new MDSystem(10, 3, 3, 2); //n_atoms=10, cube_size=3, dim=3, speed=2
    mdsys->init_system(false);
    //test(mdsys);
    for (size_t i = 0; i < STEPS; i++) {
        mdsys->calc_forces();
        mdsys->integrate();
        mdsys->print_to_file();
    }
    clock_t end = clock();
    double secs = (double)(end - start) / CLOCKS_PER_SEC;
    cout << "Done! " << secs << " sec" << endl;
}