#include "MDSystemClass.hpp"

MDSystem::MDSystem(uint32_t n_atoms, double cube_size, uint32_t dim, double speed, uint32_t n_threads) :
    N(n_atoms), SIZE(cube_size), DIM(dim), SPEED(speed), N_THREADS(n_threads)
{
    threads = vector<thread>(N_THREADS);
    init_vars();
}
void MDSystem::clean_file()
{
    out.open(".\\cpp\\coords.xmol", ios::trunc);
    out.close();
    out.open(".\\cpp\\coords.xmol", ios::app);

    // dbg.open(".\\cpp\\debug.txt", ios::trunc);
    // dbg.close();
    // dbg.open(".\\cpp\\debug.txt", ios::app);
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
pair<Vector3d, Vector3d> MDSystem::generate_r_dr(size_t x, size_t y, size_t z, double dh, Vector3d v)
{
    Vector3d r = Vector3d(vector<double>({ (x + 1. / 2.) * dh, (y + 1. / 2.) * dh, (z + 1. / 2.) * dh }));
    Vector3d dr = Vector3d(vector<double>({ v[0] * (2. * dt), v[1] * (2. * dt), v[2] * (2. * dt) }));
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
                    pair<Vector3d, Vector3d> gen_r_dr = generate_r_dr(x, y, z, dh, tmp.v);
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
        tmp_crds[i] = fixed_coord;
    }
    Vector3d dst = Vector3d(vector<double>({s, s, s}));
    if (tmp_crds.length() > dst.length()) {
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
        coord = coord - s; // too many questions
    } else if (coord <= -s / 2.) {
        coord = coord + s;
    }
    return coord;
}
vector<pair<uint32_t, uint32_t>> MDSystem::spread(uint32_t n_atoms, uint32_t n_ths)
{
    vector<pair<uint32_t, uint32_t>> tmp;
    vector<int> t = vector<int>(n_ths);
    for (size_t i = 0, j = 0; i < n_atoms; i++, j++) {
        t[j]++;
        if ((j+1) % (n_ths) == 0) 
            j = -1;
    }
    uint32_t f = 0, st;
    for (size_t i = 0; i < n_ths; i++) {
        st = t[i];
        tmp.push_back(make_pair(f, f + st));
        f += st;
    }
    return tmp;
}
void MDSystem::simulate(uint32_t from, uint32_t to) 
{
    for (size_t i = from; i < to; i++) {
        // вычисляем силу для i-го атома
        atoms[i].f = Vector3d(vector<double>(DIM));
        for (size_t j = 0; j < N; j++) {
            if (i != j) {
                Vector3d rij;
                mtx.lock();
                rij = NIM(atoms[i].r, atoms[j].r, SIZE);
                mtx.unlock();
                double _rij = rij.length();
                double ff = force_LD(_rij);
                rij.normalize();
                atoms[i].f += (rij * ff);
            }
        }
        // вычисляем его новые координаты и скорость
        Atom tmp = atoms[i];
        atoms[i].r = verle_R(atoms[i], dt);
        atoms[i].dr = (atoms[i].r - tmp.r);
        mtx.lock();
        {
        if (atoms[i].dr.length() > L_FREE_MOTION) {
            if (!(large_motion)) {
                large_motion = true;
                cout << "too large motion detected" << endl;
            }
        }
        }
        mtx.unlock();
        atoms[i].v = verle_V(atoms[i], dt);
        atoms[i].pbc(SIZE); 
        mtx.lock();
        fulltime += dt;
        mtx.unlock();
    }
}
void MDSystem::solve(size_t steps)
{
    vector<pair<uint32_t, uint32_t>> intervals = spread(N, N_THREADS);
    for (size_t step = 0; step < steps; step++) {
        for (size_t i = 0; i < N_THREADS; i++)
            threads[i] = thread(&MDSystem::simulate, this, intervals[i].first, intervals[i].second);
        for (auto &th : threads) 
            th.join();
        print_to_file();
    }
}