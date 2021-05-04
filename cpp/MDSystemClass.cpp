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

//rewrite
void MDSystem::print_to_file()
{
    out << N << endl << "***** time = * " << fulltime << " *****\n";
    for (size_t i = 0; i < N; i++) {
        out << "Ar ";
        out << fixed << setprecision(4) << r[i][0] << " " << r[i][1] << " " << r[i][2] << " ";
        out << v[i][0] << " " << v[i][1] << " " << v[i][2] << endl;
    }
}

//delete, removed into Vector3d structutre
double MDSystem::lenght(Vector3d vect)
{   
    double sum = 0;
    for (double c: vect)
        sum += pow(c, 2);
    return sqrt(sum);
}

double MDSystem::get_random_number(int min, int max)
{
    return (double)min + (rand() / (RAND_MAX / ((double)max - (double)min)));
}

//rewrite
void MDSystem::init_system(bool zero_v = false)
{
    srand(time(NULL));
    clean_file();
    Vector3d v_vec;
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
        f.push_back(Vector3d(DIM));
    }
    double k = ceil(pow(N, 1. / 3.));
    double dh = SIZE / k;
    size_t counter = 0;
    for (size_t x = 0; x < k; x++) {
        for (size_t y = 0; y < k; y++) {
            for (size_t z = 0; z < k; z++){
                if (counter < N) {
                    Vector3d r = { (x + 1. / 2.) * dh, (y + 1. / 2.) * dh, (z + 1. / 2.) * dh };
                    Vector3d dr = { v[counter][0] * 2. * dt, v[counter][1] * 2. * dt, v[counter][2] * 2. * dt };
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
//rewrite
Vector3d MDSystem::verle_R(Vector3d r, Vector3d dr, Vector3d f, double m, double dt)
{
    Vector3d tmpr;
    for (size_t i = 0; i < DIM; i++) {
        tmpr.push_back(r[i] + (dr[i] + (f[i] / (2. * m)) * pow(dt, 2.)));
    }
    return tmpr;
}
//rewrite
Vector3d MDSystem::verle_V(Vector3d dr, double dt)
{
    Vector3d tmpv;
    for (double el: dr)
        tmpv.push_back(el / (2. * dt));
    return tmpv;
}

//delete
Vector3d MDSystem::add_force(Vector3d force, Vector3d dr, double ff)
{
    Vector3d forces_sum;
    for (size_t i = 0; i < DIM; i++)
        forces_sum.push_back(force[i] + dr[i] * ff);
    return forces_sum;
}

//delete, removed into Vector3d
Vector3d MDSystem::normalize(Vector3d dr)
{
    Vector3d tmp;
    double dr_len = lenght(dr);
    for (double el: dr)
        tmp.push_back(el / dr_len);
    return tmp;
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

//deleet, removed into Vector3d structure
Vector3d MDSystem::sub(Vector3d r1, Vector3d r2)
{
    Vector3d tmp;
    for (size_t i = 0; i < DIM; i++)
        tmp.push_back(r1[i] - r2[i]);
    return tmp;
}

//delete, removed into AtomClass
Vector3d MDSystem::PBC(Vector3d r, double s)
{   
    Vector3d tmp;
    for (double coord: r)
        tmp.push_back(correct_coord(coord, 0., s));
    return tmp;
}

//delete, removed into AtomClass
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

void MDSystem::calc_forces()
{
    for (size_t i = 0; i < N; i++) {
        f[i] = Vector3d(DIM); // delete or -> Atoms[i].f = Vector3d(DIM);
        for (size_t j = 0; j < N; j++) {
            if (i != j) {
                Vector3d rij = NIM(r[i], r[j], SIZE); // -> NIM(Atoms[i].r, Atoms[j].r, SIZE);
                double _rij = lenght(rij); // -> rij.lenght();
                double ff = force_LD(_rij);
                Vector3d _dr = normalize(rij); // -> rij.normalize()
                f[i] = add_force(f[i], _dr, ff); // -> Atom[i].f += rij * ff
            }
        }
    }
}

void MDSystem::integrate()
{
    for (size_t i = 0; i < N; i++) {
        Vector3d r_tmp = r[i]; // -> Atom tmp = Atoms[i];
        r[i] = verle_R(r[i], dr[i], f[i], m[i], dt); // -> Atoms[i].r = verle_R(Atoms[i], dt);
        dr[i] = sub(r[i], r_tmp); // -> Atoms[i].dr = Atoms[i].r - tmp;
        if (lenght(dr[i]) > L_FREE_MOTION) { // -> Atoms[i].dr.lenght();
            if (!(large_motion)) {
                large_motion = true;
                cout << "too large motion detected" << endl;
            }
        }
        v[i] = verle_V(dr[i], dt); // -> Atoms[i].v = verle_V(Atoms[i].dr, dt);
        r[i] = PBC(r[i], SIZE); // -> Atoms[i].pbc(); 
        fulltime += dt;
    }
}