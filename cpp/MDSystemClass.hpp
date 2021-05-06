#pragma once
#include "Headers.hpp"
#include "AtomClass.hpp"

class MDSystem
{
public: // replace with private
    uint32_t N, DIM, N_THREADS;
    double SIZE, SPEED;
    bool large_motion, b_nim_error;
    double L_FREE_MOTION, dt, eps, sigma, rcut, rmin, fulltime;
    ofstream out, dbg;
    vector<Atom> atoms;
    vector<thread> threads;
    mutex mtx;

    void init_vars();
    void clean_file();
    double force_LD(double vec_lenght);
    double NIM_fix(double coord, double s);
    double get_random_number(int min, int max);
    Vector3d verle_R(Atom a, double dt);
    Vector3d verle_V(Atom atom, double dt);
    Vector3d NIM(Vector3d r1, Vector3d r2, double s);
    Vector3d generate_v(bool zero_v);
    pair<Vector3d, Vector3d> generate_r_dr(size_t x, size_t y, size_t z, double dh, Vector3d v);
    vector<pair<uint32_t, uint32_t>> spread(uint32_t n_atoms, uint32_t n_ths);
    void simulate(size_t from, size_t to);
    void print_to_file();
public:
    MDSystem(uint32_t n_atoms, double cube_size, uint32_t dim, double speed, uint32_t n_threads);
    void init_system(bool zero_v);
    void solve(size_t steps);
};