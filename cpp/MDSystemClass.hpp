#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <iomanip>
using namespace std;

class MDSystem
{
public:
    uint32_t N, DIM;
    double SIZE, SPEED;
    bool large_motion, b_nim_error;
    double L_FREE_MOTION, dt, eps, sigma, rcut, rmin, fulltime;
    vector<vector<double>> r, dr, f, v;
    vector<double> m;
    ofstream out, dbg;

    void init_vars();
    void clean_file();
    void new_frame();
    double force_LD(double vec_lenght);
    vector<double> verle_R(vector<double> r, vector<double> dr, vector<double> f, double m, double dt);
    vector<double> verle_V(vector<double> dr, double dt);
    vector<double> NIM(vector<double> r1, vector<double> r2, double s);
    double NIM_fix(double coord, double s);
    void PBC(vector<double> &coords, double s);
    double correct_coord(double coord, double left_bound, double right_bound);
    double lenght(vector<double> vec);
    void normalize(vector<double> &coords);
    double get_random_number(int min, int max);
    void add_force(vector<double> &force, vector<double> forces, double ff);
    vector<double> sub(vector<double> r1, vector<double> r2);

public:
    MDSystem(uint32_t n_atoms, double cube_size, uint32_t dim, double speed);
    void init_system(bool zero_v);
    void calc_forces();
    void integrate();
    void print_to_file();
};