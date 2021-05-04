#include "Headers.hpp"

class Atom 
{
private:
    double correct_coord(double coord, double left_bound, double right_bound)
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

public:
    Vector3d r, dr, f, v;
    double m;

    Atom(){ }
    Atom(Vector3d r, Vector3d dr, Vector3d f, Vector3d v, double m)
    {
        this->r = r;
        this->dr = dr;
        this->f = f;
        this->v = v;
        this->m = m;
    }
    void pbc(double s)
    {
        for (size_t i = 0; i < r.size(); i++) {
            r[i] = correct_coord(r[i], 0., s);
        }
    }
    Atom& operator = (const Atom& atom) 
    {
        this->r = atom.r;
        this->dr = atom.dr;
        this->f = atom.f;
        this->v = atom.v;
        this->m = atom.m;
        return *this;
    }
};