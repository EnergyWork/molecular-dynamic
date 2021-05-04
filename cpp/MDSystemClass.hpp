#pragma once
#include "Headers.hpp"
#include "AtomClass.hpp"

class MDSystem
{
public: // replace with private
    uint32_t N, DIM;
    double SIZE, SPEED;
    bool large_motion, b_nim_error;
    double L_FREE_MOTION, dt, eps, sigma, rcut, rmin, fulltime;
    ofstream out, dbg;
    Atoms atoms;

    void init_vars();
    void clean_file();
    double force_LD(double vec_lenght);
    double NIM_fix(double coord, double s);
    double correct_coord(double coord, double left_bound, double right_bound);
    double lenght(Vector vect);
    double get_random_number(int min, int max);
    Vector verle_R(Vector r, Vector dr, Vector f, double m, double dt); // -> Vector verle_R(Atom a, double dt);
    Vector verle_V(Vector dr, double dt);
    Vector NIM(Vector r1, Vector r2, double s);
    Vector PBC(Vector coords, double s); // -> remove into AtomClass
    Vector normalize(Vector vect); // -> remove into Vector
    Vector add_force(Vector force, Vector dr, double ff);  // -> remove into AtomClass
    Vector sub(Vector r1, Vector r2); // -> remove into Vector

public:
    MDSystem(uint32_t n_atoms, double cube_size, uint32_t dim, double speed);
    void init_system(bool zero_v);
    void calc_forces();
    void integrate();
    void print_to_file();
};