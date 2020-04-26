#ifndef REAL_H
#define REAL_H
#include <cmath>
namespace TEST {
typedef double real;
class real3 {
public:
  real x[3];

  real3() {
    x[0] = 0;
    x[1] = 0;
    x[2] = 0;
  }

  real3(double *x_) {
    x[0] = x_[0];
    x[1] = x_[1];
    x[2] = x_[2];
  }
  real3(double x1, double x2, double x3) {
    x[0] = x1;
    x[1] = x2;
    x[2] = x3;
  }

  inline real3 operator*(real scale) const {
    return real3(x[0] * scale, x[1] * scale, x[2] * scale);
  }
  inline real operator*(const real3 &y) const {
    return (x[0] * y.x[0] + x[1] * y.x[1] + x[2] * y.x[2]);
  }
  friend inline real3 operator*(real scale, const real3 &x) {
    return real3(x.x[0] * scale, x.x[1] * scale, x.x[2] * scale);
  }
  template <typename T> inline real3 operator/(T scale) const {
    return real3(x[0] / scale, x[1] / scale, x[2] / scale);
  }
  inline real3 operator+(const real3 &y) const {
    return real3(x[0] + y.x[0], x[1] + y.x[1], x[2] + y.x[2]);
  }
  inline real3 operator-(const real3 &y) const {
    return real3(x[0] - y.x[0], x[1] - y.x[1], x[2] - y.x[2]);
  }
  inline real3 &operator-=(const real3 &r) {
    x[0] -= r.x[0];
    x[1] -= r.x[1];
    x[2] -= r.x[2];
  }
  inline real3 &operator+=(const real3 &r) {
    x[0] += r.x[0];
    x[1] += r.x[1];
    x[2] += r.x[2];
  }
  inline const real &operator()(int i) const { return x[i]; }
  inline const real &operator[](int i) const { return x[i]; }
  inline real &operator()(int i) { return x[i]; }
  inline real &operator[](int i) { return x[i]; }
  inline real3 &operator*=(double scale) {
    x[0] *= scale;
    x[1] *= scale;
    x[2] *= scale;
    return *this;
  }
  friend inline real3 cross(const real3 &l, const real3 &r) {
    return real3(l.x[1] * r.x[2] - l.x[2] * r.x[1],
                 l.x[2] * r.x[0] - l.x[0] * r.x[2],
                 l.x[0] * r.x[1] - l.x[1] * r.x[0]);
  }
  inline real norm() const {
    return sqrt(x[0] * x[0] + x[1] * x[1] + x[2] * x[2]);
  }
};
} // namespace MFast
#endif