#include "MDSystemClass.hpp"
#include <iostream>
#include <vector>
#include <string>
#include <time.h>
#include <math.h>
using namespace std;

MDSystem::MDSystem(uint32_t n_atoms, uint32_t cube_size, uint32_t dim, uint32_t speed)
{
    this->N = n_atoms;
    this->SIZE = cube_size;
    this->DIM = dim;
    this->SPEED = speed;
    this->init_vars();
}

MDSystem::~MDSystem()
{
}

void MDSystem::clear_file()
{

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
    this->clear_file();
    for (size_t i = 0; i < this->N; i++) {
        this->v.push_back(vector<double>{
            this->get_random_number(-this->SPEED, this->SPEED), 
            this->get_random_number(-this->SPEED, this->SPEED), 
            this->get_random_number(-this->SPEED, this->SPEED)
        });
    }
} 