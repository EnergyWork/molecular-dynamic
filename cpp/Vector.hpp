#pragma once
#include <vector>
#include <math.h>
using namespace std;

class Vector3d
{
private:
    size_t DIM;
    vector<double> data;
public:
    Vector3d() {
        DIM = 3;
        data = vector<double>(DIM);
    }
    Vector3d(vector<double> vec) {
        DIM = 3;
        data = vec;
    }
    double& operator [] (const int i) const {
        return data[i];
    }
    Vector3d operator+(const Vector3d &v) const {
        Vector3d tmp;
        for (size_t i = 0; i < DIM; i ++)
            tmp.push(data[i] + v[i]);
        return tmp;
    }
    Vector3d operator - (const Vector3d &v) const {
        Vector3d tmp;
        for (size_t i = 0; i < DIM; i ++)
            tmp.push(data[i] - v[i]);
        return tmp;
    }
    Vector3d operator * (const double &n) const {
        Vector3d tmp;
        for (size_t i = 0; i < DIM; i ++)
            tmp.push(data[i] * n);
        return tmp;
    }
    Vector3d operator / (const double &n) const {
        Vector3d tmp;
        for (size_t i = 0; i < DIM; i ++)
            tmp.push(data[i] / n);
        return tmp;
    }
    Vector3d& operator = (const Vector3d &v) {
        data.clear();
        data = v;
        return *this;
    }
    Vector3d& operator += (const Vector3d &v) {
        for (size_t i = 0; i < DIM; i ++)
            data[i] += v[i];
        return *this;
    }
    Vector3d& operator -= (const Vector3d &v) {
        for (size_t i = 0; i < DIM; i ++)
            data[i] -= v[i];
        return *this;
    }
    void push(const double n) {
        data.push_back(n);
    }
    void clear() {
        data.clear();
    }
    size_t size() {
        return data.size();
    }
    double length() {
        double sum = 0;
        for (double c: data)
            sum += pow(c, 2);
        return sqrt(sum);
    }
    void normalize() {
        double dr_len = this->length();
        for (double& el: data)
            el = el / dr_len;
    }
};