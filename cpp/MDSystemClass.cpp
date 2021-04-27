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
    this->N = n_atoms;
    this->SIZE = cube_size;
    this->DIM = dim;
    this->SPEED = speed;
    this->init_vars();
}

void MDSystem::clean_file()
{
    this->out.open(".\\cpp\\coords.xmol", ios::trunc);
    this->out.close();
    this->out.open(".\\cpp\\coords.xmol", ios::app);
}

void MDSystem::print_to_file()
{
    out << this->N << endl
        << "***** time = *" << this->fulltime << "*****\n";
    for (size_t i = 0; i < this->N; i++) {
        out << "Ar ";
        out << this->r[i][0] << " " << this->r[i][1] << " " << this->r[i][2] << " ";
        out << this->v[i][0] << " " << this->v[i][1] << " " << this->v[i][2] << endl;
    }
}

double MDSystem::lenght(vector<double> vec)
{   
    for (double &c: vec)
        c = pow(c, 2);
    return accumulate(vec.begin(), vec.end(), 0);
}

double MDSystem::get_random_number(int min, int max)
{
    return (double)min + (rand() / (RAND_MAX / ((double)max - (double)min)));
}

void MDSystem::init_vars()
{
    this->large_motion = false;
    this->b_nim_error = false;
    this->L_FREE_MOTION = pow(pow(this->SIZE, 3), 1/3.) / (2. * this->N);
    this->dt = 0.001;
    this->eps = 1;
    this->sigma = 1;
    this->rcut = 2.5 * this->sigma;
    this->rmin = 0.00001;
    this->fulltime = this->dt;
}

void MDSystem::init_system(bool zero_v)
{
    srand(time(NULL));
    this->clean_file();
    vector<double> v_vec;
    for (size_t i = 0; i < this->N; i++) {
        if (zero_v){
            v_vec = {0., 0., 0.};
        } else {
            v_vec = {
                this->get_random_number(-this->SPEED, this->SPEED),
                this->get_random_number(-this->SPEED, this->SPEED), 
                this->get_random_number(-this->SPEED, this->SPEED)
            };
        }
        this->v.push_back(v_vec);
    }
    double k = ceil(pow(this->N, 1.0 / 3));
    double dh = this->SIZE / k;
    for (size_t i = 0; i < this->N; i++) {
        this->m.push_back(1.0);
    }
    size_t counter = 0;
    for (size_t x = 0; x < k; x++) {
        for (size_t y = 0; y < k; y++) {
            for (size_t z = 0; z < k; z++){
                if (counter < this->N) {
                    vector<double> r = {(x + 1.0 / 2) * dh, (y + 1.0 / 2) * dh, (z + 1.0 / 2) * dh};
                    vector<double> dr = {this->v[counter][0] * 2 * this->dt, this->v[counter][1] * 2 *this->dt, this->v[counter][2] * 2 * this->dt};
                    this->r.push_back(r);
                    this->dr.push_back(dr);
                    if (this->lenght(dr) > this->L_FREE_MOTION) {
                        cout << "init system error" << endl;
                    }
                    counter++;
                }
            }
        }
    }
    this->print_to_file();
}

double MDSystem::force_LD(double len)
{
    if (len > this->rcut)
        return 0;
    if (len < this->rmin)
        return this->force_LD(this->rmin);
    double x = this->sigma / len;
    return -48. * this->eps * (pow(x, 13.) - 0.5 * pow(x, 7.));
}

vector<double> MDSystem::verle_R(vector<double> r, vector<double> dr, vector<double> f, double m, double dt)
{
    vector<double> tmpr;
    for (size_t i = 0; i < this->DIM; i++) {
        tmpr.push_back(r[i] + (dr[i] + (f[i] / (2 * m)) * pow(dt, 2)));
    }
    return tmpr;
}

vector<double> MDSystem::verle_V(vector<double> dr, double dt)
{
    vector<double> tmpv;
    for (auto el: dr)
        tmpv.push_back(el / (2 * dt));
    return tmpv;
}

void MDSystem::add_force(vector<double> &force, vector<double> dr, double ff)
{
    for (auto &el: dr)
        el *= ff;
    for (size_t i = 0; i < force.size(); i++)
        force[i] += dr[i];
}

void MDSystem::normalize(vector<double> &dr)
{
    for (auto &el: dr)
        el /= this->lenght(dr);
}

vector<double> MDSystem::NIM(vector<double> r1, vector<double> r2, double s)
{
    vector<double> tmp_crds, dist;
    for (size_t i = 0; i < this->DIM; i++)
        tmp_crds.push_back(this->NIM_fix(-(r1[i] - r2[i]), s));
    dist = vector<double>{s, s, s};
    if (this->lenght(tmp_crds) > this->lenght(dist))
        if (!(this->b_nim_error)) {
            this->b_nim_error = true;
            cout << "nim error" << endl;
        }
    return tmp_crds;
}

double MDSystem::NIM_fix(double coord, double s)
{
    if (coord >= s / 2.)
        coord -= s;
    else if (coord <= -s / 2.)
        coord += s;
    return coord;
}

void MDSystem::calc_forces()
{
    double ff, _rij;
    for (size_t i = 0; i < this->N; i++) {
        this->f.push_back(vector<double>(3));
        for (size_t j = 0; j < this->N; j++){
            if (i != j) {
                vector<double> rij = this->NIM(r[i], r[j], this->SIZE);
                _rij = this->lenght(rij);
                ff = this->force_LD(_rij);
                vector<double> dr = rij;
                this->normalize(dr);
                this->add_force(this->f[i], dr, ff);
            }
        }
    }
}

vector<double> MDSystem::sub(vector<double> r1, vector<double> r2)
{
    vector<double> tmp;
    for (size_t i = 0; i < this->DIM; i++)
        tmp.push_back(r1[i] - r2[i]);
    return tmp;
}

void MDSystem::PBC(vector<double> &r, double s)
{
    for (auto &coord: r)
        coord = this->correct_coord(coord, 0, s);
}

double MDSystem::correct_coord(double coord, double left_bound, double right_bound)
{
    double len = right_bound - left_bound;
    double d;
    if (coord >= right_bound) {
        d = coord - left_bound;
        coord -= len * floor(d / len);
    } else if (coord < left_bound) {
        d = left_bound - coord;
        coord = right_bound - len * (d / len - floor(d / len));
    }
    return coord;
}

void MDSystem::integrate()
{
    vector<double> r_tmp;
    for (size_t i = 0; i < this->N; i++) {
        r_tmp = this->r[i];
        this->r[i] = this->verle_R(this->r[i], this->dr[i], this->f[i], this->m[i], this->dt);
        this->dr[i] = sub(this->r[i], r_tmp);
        if (this->lenght(dr[i]) > this->L_FREE_MOTION) {
            if (!(this->large_motion)) {
                this->large_motion = true;
                cout << "too large motion detected" << endl;
            }
        }
        this->v[i] = this->verle_V(this->dr[i], this->dt);
        this->fulltime += this->dt;
        this->PBC(this->r[i], this->SIZE);
    }
}