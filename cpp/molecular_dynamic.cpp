#include "Headers.hpp"
#include "MDSystemClass.hpp"

using namespace std;

#define CLOCKS_PER_SEC 1000
#define STEPS 5000

void print_vector(vector<double> vec)
{
    for (auto el: vec)
        cout << el << " ";
    cout << endl;
}

void test(MDSystem *sys)
{
    auto r1 = vector<double>{0.75, -1.5, 3};
    //auto r2 = vector<double>{3, 2, 4};
    cout << 0.75 / sys->lenght(r1) << endl;
    print_vector(sys->normalize(r1));
}

int main()
{   
    cout << "Start" << endl;
    clock_t start = clock();
    MDSystem *mdsys = new MDSystem(5, 3, 3, 2); // example: n_atoms=10, cube_size=3, dim=3, speed=2
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