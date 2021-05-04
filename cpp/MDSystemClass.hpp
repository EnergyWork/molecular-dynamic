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
    Vector3d
    void clean_file();
    double force_LD(double vec_lenght);
    double NIM_fix(double coord, double s);
    //double correct_coord(double coord, double left_bound, double right_bound);
    //double lenght(Vector3d vect);
    double get_random_number(int min, int max);
    Vector3d verle_R(Atom a, double dt);
    Vector3d verle_V(Atom atom, double dt);
    Vector3d NIM(Vector3d r1, Vector3d r2, double s);
    Vector3d generate_v(bool zero_v);
    pair<Vector3d, Vector3d> generate_r_dr();
    //Vector3d PBC(Vector3d coords, double s); // -> remove into AtomClass
    //Vector3d normalize(Vector3d vect); // -> remove into Vector3d
    //Vector3d add_force(Vector3d force, Vector3d dr, double ff);  // -> remove into AtomClass
    //Vector3d sub(Vector3d r1, Vector3d r2); // -> remove into Vector3d

public:
    MDSystem(uint32_t n_atoms, double cube_size, uint32_t dim, double speed);
    void init_system(bool zero_v);
    void calc_forces();
    void integrate();
    void print_to_file();
};