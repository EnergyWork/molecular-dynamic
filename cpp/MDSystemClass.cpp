#include "MDSystemClass.hpp"
#include <iostream>
#include <vector>
#include <string>
#include <time.h>
#include <math.h>
#include <algorithm>
#include <numeric>
using namespace std;

MDSystem::MDSystem(uint32_t n_atoms, double cube_size, uint32_t dim, double speed)
{
   N = n_atoms;
   SIZE = cube_size;
   DIM = dim;
   SPEED = speed;
   init_vars();
}

void MDSystem::clean_file()
{
   out.open(".\\cpp\\coords.xmol", ios::trunc);
   out.close();
   out.open(".\\cpp\\coords.xmol", ios::app);

   dbg.open(".\\cpp\\debug.txt", ios::trunc);
   dbg.close();
   dbg.open(".\\cpp\\debug.txt", ios::app);
}

void MDSystem::print_to_file(bool debg = false)
{
    out << N << endl << "***** time = * " << fulltime << " *****\n";
    for (size_t i = 0; i < N; i++) {
        out << "Ar ";
        out << fixed << setprecision(4) << r[i][0] << " " << r[i][1] << " " << r[i][2] << " ";
        out << v[i][0] << " " << v[i][1] << " " << v[i][2] << endl;
    }
}

double MDSystem::lenght(vector<double> vec)
{   
    double sum = 0;
    for (double c: vec)
        sum += pow(c, 2);
    return sqrt(sum);
}

double MDSystem::get_random_number(int min, int max)
{
    return (double)min + (rand() / (RAND_MAX / ((double)max - (double)min)));
}

void MDSystem::init_vars()
{
    large_motion = false;
    b_nim_error = false;
    L_FREE_MOTION = pow(pow(SIZE, 3.), 1./3.) / (2. * N);
    dt = 0.001;
    eps = 1.;
    sigma = 1.;
    rcut = 2.5 * sigma;
    rmin = 0.00001;
    fulltime = dt;
}

void MDSystem::init_system(bool zero_v)
{
    srand(42); // srand(time(NULL));
    clean_file();
    vector<double> v_vec;
    for (size_t i = 0; i < N; i++) {
        if (zero_v) {
            v_vec = {0., 0., 0.};
        } else {
            v_vec = {
               get_random_number(-SPEED, SPEED),
               get_random_number(-SPEED, SPEED), 
               get_random_number(-SPEED, SPEED)
            };
        }
        v.push_back(v_vec);
        m.push_back(1.0);
    }
    double k = ceil(pow(N, 1. / 3.));
    double dh = SIZE / k;
    size_t counter = 0;
    for (size_t x = 0; x < k; x++) {
        for (size_t y = 0; y < k; y++) {
            for (size_t z = 0; z < k; z++){
                if (counter < N) {
                    vector<double> r = { (x + 1. / 2.) * dh, (y + 1. / 2.) * dh, (z + 1. / 2.) * dh };
                    vector<double> dr = { v[counter][0] * 2. * dt, v[counter][1] * 2. * dt, v[counter][2] * 2. * dt };
                    this->r.push_back(r);
                    this->dr.push_back(dr);
                    if (lenght(dr) > L_FREE_MOTION) {
                        cout << "init system error" << endl;
                    }
                    counter++;
                } else {
                    break;
                }
            }
        }
    }
    print_to_file();
}

double MDSystem::force_LD(double len)
{
    if (len > rcut)
        return 0.;
    if (len < rmin)
        return force_LD(rmin);
    double x = sigma / len;
    return -48. * eps * (pow(x, 13.) - 0.5 * pow(x, 7.));
}

vector<double> MDSystem::verle_R(vector<double> r, vector<double> dr, vector<double> f, double m, double dt)
{
    vector<double> tmpr;
    for (size_t i = 0; i < DIM; i++) {
        tmpr.push_back(r[i] + (dr[i] + (f[i] / (2. * m)) * pow(dt, 2.)));
    }
    return tmpr;
}

vector<double> MDSystem::verle_V(vector<double> dr, double dt)
{
    vector<double> tmpv;
    for (double el: dr)
        tmpv.push_back(el / (2. * dt));
    return tmpv;
}

vector<double> MDSystem::add_force(vector<double> force, vector<double> dr, double ff)
{
    vector<double> forces_sum;
    for (size_t i = 0; i < DIM; i++)
        forces_sum.push_back(force[i] + dr[i] * ff);
    return forces_sum;
}

vector<double> MDSystem::normalize(vector<double> dr)
{
    vector<double> tmp;
    double dr_len = lenght(dr);
    for (double el: dr)
        tmp.push_back(el / dr_len);
    return tmp;
}

vector<double> MDSystem::NIM(vector<double> r1, vector<double> r2, double s)
{
    vector<double> tmp_crds;
    for (size_t i = 0; i < DIM; i++) {
        double n = -(r1[i] - r2[i]);
        double fixed_coord = NIM_fix(n, s);
        tmp_crds.push_back(fixed_coord);
    }
    vector<double> dst = vector<double>{s, s, s};
    if (lenght(tmp_crds) > lenght(dst)) {
        if (!(b_nim_error)) {
            b_nim_error = true;
            cout << "nim error" << endl;
        }
    }
    return tmp_crds;
}

double MDSystem::NIM_fix(double coord, double s)
{
    if (coord >= s / 2.) {
        coord = s - coord;
    } else if (coord <= -s / 2.) {
        coord = coord + s;
    }
    return coord;
}

void MDSystem::calc_forces()
{
    vector<double> rij, _dr;
    double ff = 0, _rij = 0;
    for (size_t i = 0; i < N; i++) {
        f.push_back(vector<double>(3));
        for (size_t j = 0; j < N; j++) {
            if (i != j) {
                rij = NIM(r[i], r[j], SIZE);
                //dbg << fixed << setprecision(5) << rij[0] << "  " << rij[1] << "  " << rij[2] << endl;
                _rij = lenght(rij);
                ff = force_LD(_rij);
                _dr = rij;
                _dr = normalize(_dr);
                f[i] = add_force(f[i], _dr, ff);
            }
        }
    }
}

vector<double> MDSystem::sub(vector<double> r1, vector<double> r2)
{
    vector<double> tmp;
    for (size_t i = 0; i < DIM; i++)
        tmp.push_back(r1[i] - r2[i]);
    return tmp;
}

vector<double> MDSystem::PBC(vector<double> r, double s)
{   
    vector<double> tmp;
    for (double coord: r)
        tmp.push_back(correct_coord(coord, 0., s));
    return tmp;
}

double MDSystem::correct_coord(double coord, double left_bound, double right_bound)
{
    double len = right_bound - left_bound;
    double d;
    if (coord >= right_bound) {
        d = floor(coord / len);
        coord = coord - d * len;
    } else if (coord < left_bound) {
        d = floor(abs(coord) / len);
        coord = right_bound - (abs(coord) - d * len);
    }
    return coord;
}

void MDSystem::integrate()
{
    for (size_t i = 0; i < N; i++) {
        vector<double> r_tmp = r[i];
        r[i] = verle_R(r[i], dr[i], f[i], m[i], dt);
        dbg << fixed << setprecision(9) << r[i][0] << "  " << r[i][1] << "  " << r[i][2] << endl;
        dr[i] = sub(r[i], r_tmp);
        if (lenght(dr[i]) > L_FREE_MOTION) {
            if (!(large_motion)) {
                large_motion = true;
                cout << "too large motion detected" << endl;
            }
        }
        v[i] = verle_V(dr[i], dt);
        r[i] = PBC(r[i], SIZE);
        fulltime += dt;
    }
}