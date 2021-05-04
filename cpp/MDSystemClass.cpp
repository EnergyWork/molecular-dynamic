#include "MDSystemClass.hpp"

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
    rmin = 0.000001;
    fulltime = dt;
}
void MDSystem::print_to_file()
{
    out << N << endl << "***** time = * " << fulltime << " *****\n";
    for (size_t i = 0; i < N; i++) {
        out << "Ar ";
        out << fixed << setprecision(4) << atoms[i].r[0] << " " << atoms[i].r[1] << " " << atoms[i].r[2] << " ";
        out << atoms[i].v[0] << " " << atoms[i].v[1] << " " << atoms[i].v[2] << endl;
    }
}
double MDSystem::get_random_number(int min, int max)
{
    return (double)min + (rand() / (RAND_MAX / ((double)max - (double)min)));
}
Vector3d MDSystem::generate_v(bool zero_v)
{
    if (zero_v) {
        return Vector3d();
    } else {
        vector<double> tmp;
        tmp = {
            get_random_number(-SPEED, SPEED),
            get_random_number(-SPEED, SPEED), 
            get_random_number(-SPEED, SPEED)
        };
        return Vector3d(tmp);
    }
}
pair<Vector3d, Vector3d> MDSystem::generate_r_dr(size_t x, size_t y, size_t z, double dh)
{
    Vector3d r = Vector3d(vector<double>({ (x + 1. / 2.) * dh, (y + 1. / 2.) * dh, (z + 1. / 2.) * dh }));
    Vector3d dr = Vector3d(vector<double>({ v[counter][0] * 2. * dt, v[counter][1] * 2. * dt, v[counter][2] * 2. * dt }));
    return make_pair(r, dr);
}
void MDSystem::init_system(bool zero_v = false)
{
    srand(time(NULL));
    clean_file();
    double k = ceil(pow(N, 1. / 3.));
    double dh = SIZE / k;
    size_t counter = 0;
    for (size_t x = 0; x < k; x++) {
        for (size_t y = 0; y < k; y++) {
            for (size_t z = 0; z < k; z++) {
                if (counter < N) {
                    Atom tmp;
                    tmp.v = generate_v(zero_v);
                    tmp.m = 1.0;
                    tmp.f = Vector3d();
                    pair<Vector3d, Vector3d> gen_r_dr = generate_r_dr(x, y, z, dh);
                    tmp.r = gen_r_dr.first;
                    tmp.dr = gen_r_dr.second;
                    if (tmp.dr.length() > L_FREE_MOTION) {
                        if (!large_motion) {
                            large_motion = true;
                            cout << "init system error" << endl;
                        }
                    }
                    atoms.push_back(tmp);
                    counter++;
                } else break;
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
Vector3d MDSystem::verle_R(Atom atom, double dt)
{
    return atom.r + (atom.dr + (atom.f / (2. * atom.m)) *  pow(dt, 2.));
}
Vector3d MDSystem::verle_V(Atom atom, double dt)
{
    return atom.dr / (2. * dt);
}
Vector3d MDSystem::NIM(Vector3d r1, Vector3d r2, double s)
{
    Vector3d tmp_crds;
    for (size_t i = 0; i < DIM; i++) {
        double n = -(r1[i] - r2[i]);
        double fixed_coord = NIM_fix(n, s);
        tmp_crds.push_back(fixed_coord);
    }
    Vector3d dst = Vector3d{s, s, s};
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
        coord = coord - s;
    } else if (coord <= -s / 2.) {
        coord = coord + s;
    }
    return coord;
}
void MDSystem::calc_forces()
{
    for (size_t i = 0; i < N; i++) {
        atoms[i].f.clear();
        for (size_t j = 0; j < N; j++) {
            if (i != j) {
                Vector3d rij = NIM(atoms[i].r, atoms[j].r, SIZE);
                double _rij = rij.lenght();
                double ff = force_LD(_rij);
                Vector3d _dr = rij.normalize();
                atoms[i].f += (rij * ff);
            }
        }
    }
}
void MDSystem::integrate()
{
    for (size_t i = 0; i < N; i++) {
        Atom tmp = atoms[i];
        atoms[i].r = verle_R(atoms[i], dt);
        atoms[i].dr = atoms[i].r - tmp;
        if (atoms[i].dr.lenght() > L_FREE_MOTION) {
            if (!(large_motion)) {
                large_motion = true;
                cout << "too large motion detected" << endl;
            }
        }
        atoms[i].v = verle_V(atoms[i], dt);
        atoms[i].pbc(SIZE); 
        fulltime += dt;
    }
}