#include <iostream>
#include <vector>
#include <string>
using namespace std;

class MDSystem
{
private:
    uint32_t N, SIZE, DIM, SPEED;
    bool large_motion, b_nim_error;
    double L_FREE_MOTION, dt, eps, sigma, rcut, rmin, fulltime;
    vector<vector<double>> r, dr, f, m, v;

    void init_vars();
    void clear_file();
    void new_frame();
    double force_LD(vector<double> vec);
    vector<double> verle_R(vector<double>, vector<double>, vector<double>, double);
    vector<double> verle_V(vector<double>, double);
    vector<double> NIM(vector<double>, vector<double>, double);
    double NIM_fix(double, double);
    void PBC(vector<double>& coords, double s);
    double correct_coord(double coord, double left_bound, double right_bound);
    double lenght(vector<double> vec);
    double normalise(vector<double> coords);
    double get_random_number(int min, int max);
public:
    MDSystem(uint32_t, uint32_t, uint32_t, uint32_t);
    ~MDSystem();

    void init_system(bool);
    void calc_forces();
    void integrate();
    void print_to_file();
};