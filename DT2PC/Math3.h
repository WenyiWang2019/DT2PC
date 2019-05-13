#ifndef MATH3_H
#define MATH3_H

#include <assert.h>
#include <math.h>
#include <algorithm> 
#include <string.h>

/// Vector dim 3
template <typename T>
class Vector3 {
 public:
  T &operator[](size_t i) {
    assert(i < 3);
    return data[i];
  }
  const T &operator[](size_t i) const {
    assert(i < 3);
    return data[i];
  }

  size_t getElementCount() const { return 3; }
  T &r() { return data[0]; }
  T &g() { return data[1]; }
  T &b() { return data[2]; }
  const T &r() const { return data[0]; }
  const T &g() const { return data[1]; }
  const T &b() const { return data[2]; }
  T &Y() { return data[0]; }
  T &U() { return data[1]; }
  T &V() { return data[2]; }
  const T &Y() const { return data[0]; }
  const T &U() const { return data[1]; }
  const T &V() const { return data[2]; }
  T &x() { return data[0]; }
  T &y() { return data[1]; }
  T &z() { return data[2]; }
  const T &x() const { return data[0]; }
  const T &y() const { return data[1]; }
  const T &z() const { return data[2]; }
  void normalize() {
    const T norm2 = data[0] * data[0] + data[1] * data[1] + data[2] * data[2];
    if (norm2 != 0.0) {
      T invNorm = static_cast<T>(1.0 / sqrt(norm2));
      (*this) *= invNorm;
    }
  }
  T getNorm() const { return static_cast<T>(sqrt(getNorm2())); }
  T getNorm2() const { return (*this) * (*this); }
  Vector3 &operator=(const Vector3 &rhs) {
    memcpy(data, rhs.data, sizeof(data));
    return *this;
  }
  Vector3 &operator+=(const Vector3 &rhs) {
    data[0] += rhs.data[0];
    data[1] += rhs.data[1];
    data[2] += rhs.data[2];
    return *this;
  }
  Vector3 &operator-=(const Vector3 &rhs) {
    data[0] -= rhs.data[0];
    data[1] -= rhs.data[1];
    data[2] -= rhs.data[2];
    return *this;
  }
  Vector3 &operator-=(const T a) {
    data[0] -= a;
    data[1] -= a;
    data[2] -= a;
    return *this;
  }
  Vector3 &operator+=(const T a) {
    data[0] += a;
    data[1] += a;
    data[2] += a;
    return *this;
  }
  Vector3 &operator/=(const T a) {
    assert(a != 0);
    data[0] /= a;
    data[1] /= a;
    data[2] /= a;
    return *this;
  }
  Vector3 &operator*=(const T a) {
    data[0] *= a;
    data[1] *= a;
    data[2] *= a;
    return *this;
  }
  Vector3 &operator=(const T a) {
    data[0] = a;
    data[1] = a;
    data[2] = a;
    return *this;
  }
  Vector3 &operator=(const T *const rhs) {
    data[0] = rhs[0];
    data[1] = rhs[1];
    data[2] = rhs[2];
    return *this;
  }
  T operator*(const Vector3 &rhs) const {
    return (data[0] * rhs.data[0] + data[1] * rhs.data[1] + data[2] * rhs.data[2]);
  }
  Vector3 operator^(const Vector3 &rhs) const {
    return Vector3<T>(data[1] * rhs.data[2] - data[2] * rhs.data[1],
                         data[2] * rhs.data[0] - data[0] * rhs.data[2],
                         data[0] * rhs.data[1] - data[1] * rhs.data[0]);
  }
  Vector3 operator-() const { return Vector3<T>(-data[0], -data[1], -data[2]); }
  friend Vector3 operator+(const Vector3 &lhs, const Vector3 &rhs) {
    return Vector3<T>(lhs.data[0] + rhs.data[0], lhs.data[1] + rhs.data[1],
                         lhs.data[2] + rhs.data[2]);
  }
  friend Vector3 operator+(const T lhs, const Vector3 &rhs) {
    return Vector3<T>(lhs + rhs.data[0], lhs + rhs.data[1], lhs + rhs.data[2]);
  }
  friend Vector3 operator+(const Vector3 &lhs, const T rhs) {
    return Vector3<T>(lhs.data[0] + rhs, lhs.data[1] + rhs, lhs.data[2] + rhs);
  }
  friend Vector3 operator-(const Vector3 &lhs, const Vector3 &rhs) {
    return Vector3<T>(lhs.data[0] - rhs.data[0], lhs.data[1] - rhs.data[1],
                         lhs.data[2] - rhs.data[2]);
  }
  friend Vector3 operator-(const T lhs, const Vector3 &rhs) {
    return Vector3<T>(lhs - rhs.data[0], lhs - rhs.data[1], lhs - rhs.data[2]);
  }
  friend Vector3 operator-(const Vector3 &lhs, const T rhs) {
    return Vector3<T>(lhs.data[0] - rhs, lhs.data[1] - rhs, lhs.data[2] - rhs);
  }
  friend Vector3 operator*(const T lhs, const Vector3 &rhs) {
    return Vector3<T>(lhs * rhs.data[0], lhs * rhs.data[1], lhs * rhs.data[2]);
  }
  friend Vector3 operator*(const Vector3 &lhs, const T rhs) {
    return Vector3<T>(lhs.data[0] * rhs, lhs.data[1] * rhs, lhs.data[2] * rhs);
  }
  friend Vector3 operator/(const Vector3 &lhs, const T rhs) {
    assert(rhs != 0);
    return Vector3<T>(lhs.data[0] / rhs, lhs.data[1] / rhs, lhs.data[2] / rhs);
  }
  bool operator<(const Vector3 &rhs) const {
    if (data[0] == rhs.data[0]) {
      if (data[1] == rhs.data[1]) {
        return (data[2] < rhs.data[2]);
      }
      return (data[1] < rhs.data[1]);
    }
    return (data[0] < rhs.data[0]);
  }

  bool operator>(const Vector3 &rhs) const {
    if (data[0] == rhs.data[0]) {
      if (data[1] == rhs.data[1]) {
        return (data[2] > rhs.data[2]);
      }
      return (data[1] > rhs.data[1]);
    }
    return (data[0] > rhs.data[0]);
  }
  bool operator==(const Vector3 &rhs) const {
	  return (int(data[0]) == int(rhs.data[0]) && int(data[1]) == int(rhs.data[1]) && int(data[2]) == int(rhs.data[2]));
  }
  bool operator!=(const Vector3 &rhs) const {
    return (data[0] != rhs.data[0] || data[1] != rhs.data[1] || data[2] != rhs.data[2]);
  }
  friend std::ostream &operator<<(std::ostream &os, const Vector3 &vec) {
    os << vec[0] << " " << vec[1] << " " << vec[2] << std::endl;
    return os;
  }
  friend std::istream &operator>>(std::istream &is, Vector3 &vec) {
    is >> vec[0] >> vec[1] >> vec[2];
    return is;
  }

  void SetXYZ(const T x, const T y, const T z)
  {
	data[0] = x;
	data[1] = y;
	data[2] = z;
  }

  Vector3(const T a) { data[0] = data[1] = data[2] = a; }
  Vector3(const T x, const T y, const T z) {
    data[0] = x;
    data[1] = y;
    data[2] = z;
  }
  Vector3(const Vector3 &vec) {
    data[0] = vec.data[0];
    data[1] = vec.data[1];
    data[2] = vec.data[2];
  }
  Vector3(const T *vec) { memcpy(data, vec, sizeof(data)); }
  Vector3() = default;
  ~Vector3(void) = default;

 private:
  T data[3];
};

template <typename T>
struct Box3 {
  Vector3<T> min;
  Vector3<T> max;
  bool contains(const Vector3<T> point) const {
    return !(point.x() < min.x() || point.x() > max.x() || point.y() < min.y() ||
             point.y() > max.y() || point.z() < min.z() || point.z() > max.z());
  }

  Box3 merge(const Box3 &box) {
    min.x() = std::min(min.x(), box.min.x());
    min.y() = std::min(min.y(), box.min.y());
    min.z() = std::min(min.z(), box.min.z());
    max.x() = std::max(max.x(), box.max.x());
    max.y() = std::max(max.y(), box.max.y());
    max.z() = std::max(max.z(), box.max.z());
    return box;
  }

  bool intersects(const Box3 &box) {
    return max.x() >= box.min.x() && min.x() <= box.max.x() && max.y() >= box.min.y() &&
           min.y() <= box.max.y() && max.z() >= box.min.z() && min.z() <= box.max.z();
  }

  friend std::ostream &operator<<(std::ostream &os, const Box3 &box) {
    os << box.min[0] << " " << box.min[1] << " " << box.min[2] << " " << box.max[0] << " "
       << box.max[1] << " " << box.max[2] << std::endl;
    return os;
  }
  friend std::istream &operator>>(std::istream &is, Box3 &box) {
    is >> box.min[0] >> box.min[1] >> box.min[2] >> box.max[0] >> box.max[1] >> box.max[2];
    return is;
  }
};

//!    3x3 Matrix
template <typename T>
class Matrix3 {
 public:
  T *operator[](const size_t rowIndex) {
    assert(rowIndex < 3);
    return data[rowIndex];
  }
  const T *operator[](const size_t rowIndex) const {
    assert(rowIndex < 3);
    return data[rowIndex];
  }
  size_t getColumnCount() const { return 3; }
  size_t getRowCount() const { return 3; }
  Matrix3 &operator=(const Matrix3 &rhs) {
    memcpy(data, rhs.data, sizeof(data));
    return *this;
  }
  void operator+=(const Matrix3 &rhs) {
    for (int i = 0; i < 3; ++i) {
      for (int j = 0; j < 3; ++j) {
        this->data[i][j] += rhs.data[i][j];
      }
    }
  }
  void operator-=(const Matrix3 &rhs) {
    for (int i = 0; i < 3; ++i) {
      for (int j = 0; j < 3; ++j) {
        this->data[i][j] -= rhs.data[i][j];
      }
    }
  }
  void operator-=(const T a) {
    for (int i = 0; i < 3; ++i) {
      for (int j = 0; j < 3; ++j) {
        this->data[i][j] -= a;
      }
    }
  }
  void operator+=(const T a) {
    for (int i = 0; i < 3; ++i) {
      for (int j = 0; j < 3; ++j) {
        this->data[i][j] += a;
      }
    }
  }
  void operator/=(const T a) {
    assert(a != 0);
    for (int i = 0; i < 3; ++i) {
      for (int j = 0; j < 3; ++j) {
        this->data[i][j] /= a;
      }
    }
  }
  void operator*=(const T a) {
    for (int i = 0; i < 3; ++i) {
      for (int j = 0; j < 3; ++j) {
        this->data[i][j] *= a;
      }
    }
  }
  Vector3<T> operator*(const Vector3<T> &rhs) const {
    Vector3<T> res;
    for (int i = 0; i < 3; ++i) {
      res[i] = 0;
      for (int j = 0; j < 3; ++j) {
        res[i] += this->data[i][j] * rhs[j];
      }
    }
    return res;
  }
  Matrix3 operator*(const Matrix3 &rhs) const {
    Matrix3<T> res;
    for (int i = 0; i < 3; ++i) {
      for (int j = 0; j < 3; ++j) {
        res.data[i][j] = 0;
        for (int k = 0; k < 3; ++k) {
          res.data[i][j] += this->data[i][k] * rhs.data[k][j];
        }
      }
    }
    return res;
  }
  Matrix3 operator+(const Matrix3 &rhs) const {
    Matrix3<T> res;
    for (int i = 0; i < 3; ++i) {
      for (int j = 0; j < 3; ++j) {
        res.data[i][j] = this->data[i][j] + rhs.data[i][j];
      }
    }
    return res;
  }
  Matrix3 operator-(const Matrix3 &rhs) const {
    Matrix3<T> res;
    for (int i = 0; i < 3; ++i) {
      for (int j = 0; j < 3; ++j) {
        res.data[i][j] = this->data[i][j] - rhs.data[i][j];
      }
    }
    return res;
  }
  Matrix3 operator-() const {
    Matrix3<T> res;
    for (int i = 0; i < 3; ++i) {
      for (int j = 0; j < 3; ++j) {
        res.data[i][j] = -this->data[i][j];
      }
    }
    return res;
  }
  Matrix3 operator*(T rhs) const {
    Matrix3<T> res;
    for (int i = 0; i < 3; ++i) {
      for (int j = 0; j < 3; ++j) {
        res.data[i][j] = this->data[i][j] * rhs;
      }
    }
    return res;
  }
  Matrix3 operator/(T rhs) const {
    assert(rhs != 0);
    Matrix3<T> res;
    for (int i = 0; i < 3; ++i) {
      for (int j = 0; j < 3; ++j) {
        res.data[i][j] = this->data[i][j] / rhs;
      }
    }
    return res;
  }
  Matrix3 transpose() const {
    Matrix3<T> res;
    for (int i = 0; i < 3; ++i) {
      for (int j = 0; j < 3; ++j) {
        res.data[i][j] = this->data[j][i];
      }
    }
    return res;
  }
  Matrix3() = default;
  Matrix3(const T a) {
    for (int i = 0; i < 3; ++i) {
      for (int j = 0; j < 3; ++j) {
        data[i][j] = a;
      }
    }
  }
  Matrix3(const Matrix3 &rhs) { memcpy(data, rhs.data, sizeof(data)); }
  ~Matrix3(void) = default;

  friend inline Matrix3<T> operator*(T lhs, const Matrix3<T> &rhs) {
    Matrix3<T> res;
    for (int i = 0; i < 3; ++i) {
      for (int j = 0; j < 3; ++j) {
        res.data[i][j] = lhs * rhs.data[i][j];
      }
    }
    return res;
  }
  static void makeIdentity(Matrix3<T> &mat) {
    memset(mat.data, 0, sizeof(mat.data));
    for (int i = 0; i < 3; ++i) {
      mat[i][i] = 1;
    }
  }

  static void makeScale(const T sx, const T sy, const T sz, Matrix3<T> &mat) {
    makeIdentity(mat);
    mat[0][0] = sx;
    mat[1][1] = sy;
    mat[2][2] = sz;
  }
  static void makeUniformScale(const T s, Matrix3<T> &mat) { makeScale(s, s, s, mat); }
  static void makeRotation(const T angle, const T ax, const T ay, const T az, Matrix3<T> &mat) {
    T c = cos(angle);
    T l_c = 1 - c;
    T s = sin(angle);
    mat[0][0] = ax * ax + (1 - ax * ax) * c;
    mat[0][1] = ax * ay * l_c - az * s;
    mat[0][2] = ax * az * l_c + ay * s;
    mat[1][0] = ax * ay * l_c + az * s;
    mat[1][1] = ay * ay + (1 - ay * ay) * c;
    mat[1][2] = ay * az * l_c - ax * s;
    mat[2][0] = ax * az * l_c - ay * s;
    mat[2][1] = ay * az * l_c + ax * s;
    mat[2][2] = az * az + (1 - az * az) * c;
  }

 private:
  T data[3][3];
};

// Slightly modified version of http://www.melax.com/diag.html?attredirects=0
// A must be a symmetric matrix.
// returns Q and D such that
// Diagonal matrix D = QT * A * Q;  and  A = Q*D*QT
void Diagonalize(const Matrix3<double> &A, Matrix3<double> &Q, Matrix3<double> &D) {
  const int maxsteps = 24;  // certainly wont need that many.
  int k0, k1, k2;
  double o[3], m[3];
  double q[4] = {0.0, 0.0, 0.0, 1.0};
  double jr[4];
  double sqw, sqx, sqy, sqz;
  double tmp1, tmp2, mq;
  Matrix3<double> AQ;
  double thet, sgn, t, c;
  for (int i = 0; i < maxsteps; ++i) {
    // quat to matrix
    sqx = q[0] * q[0];
    sqy = q[1] * q[1];
    sqz = q[2] * q[2];
    sqw = q[3] * q[3];
    Q[0][0] = (sqx - sqy - sqz + sqw);
    Q[1][1] = (-sqx + sqy - sqz + sqw);
    Q[2][2] = (-sqx - sqy + sqz + sqw);
    tmp1 = q[0] * q[1];
    tmp2 = q[2] * q[3];
    Q[1][0] = 2.0 * (tmp1 + tmp2);
    Q[0][1] = 2.0 * (tmp1 - tmp2);
    tmp1 = q[0] * q[2];
    tmp2 = q[1] * q[3];
    Q[2][0] = 2.0 * (tmp1 - tmp2);
    Q[0][2] = 2.0 * (tmp1 + tmp2);
    tmp1 = q[1] * q[2];
    tmp2 = q[0] * q[3];
    Q[2][1] = 2.0 * (tmp1 + tmp2);
    Q[1][2] = 2.0 * (tmp1 - tmp2);

    // AQ = A * Q;
    AQ[0][0] = Q[0][0] * A[0][0] + Q[1][0] * A[0][1] + Q[2][0] * A[0][2];
    AQ[0][1] = Q[0][1] * A[0][0] + Q[1][1] * A[0][1] + Q[2][1] * A[0][2];
    AQ[0][2] = Q[0][2] * A[0][0] + Q[1][2] * A[0][1] + Q[2][2] * A[0][2];
    AQ[1][0] = Q[0][0] * A[0][1] + Q[1][0] * A[1][1] + Q[2][0] * A[1][2];
    AQ[1][1] = Q[0][1] * A[0][1] + Q[1][1] * A[1][1] + Q[2][1] * A[1][2];
    AQ[1][2] = Q[0][2] * A[0][1] + Q[1][2] * A[1][1] + Q[2][2] * A[1][2];
    AQ[2][0] = Q[0][0] * A[0][2] + Q[1][0] * A[1][2] + Q[2][0] * A[2][2];
    AQ[2][1] = Q[0][1] * A[0][2] + Q[1][1] * A[1][2] + Q[2][1] * A[2][2];
    AQ[2][2] = Q[0][2] * A[0][2] + Q[1][2] * A[1][2] + Q[2][2] * A[2][2];

    // D  = Q.transpose() * AQ;
    D[0][0] = AQ[0][0] * Q[0][0] + AQ[1][0] * Q[1][0] + AQ[2][0] * Q[2][0];
    D[0][1] = AQ[0][0] * Q[0][1] + AQ[1][0] * Q[1][1] + AQ[2][0] * Q[2][1];
    D[0][2] = AQ[0][0] * Q[0][2] + AQ[1][0] * Q[1][2] + AQ[2][0] * Q[2][2];
    D[1][0] = AQ[0][1] * Q[0][0] + AQ[1][1] * Q[1][0] + AQ[2][1] * Q[2][0];
    D[1][1] = AQ[0][1] * Q[0][1] + AQ[1][1] * Q[1][1] + AQ[2][1] * Q[2][1];
    D[1][2] = AQ[0][1] * Q[0][2] + AQ[1][1] * Q[1][2] + AQ[2][1] * Q[2][2];
    D[2][0] = AQ[0][2] * Q[0][0] + AQ[1][2] * Q[1][0] + AQ[2][2] * Q[2][0];
    D[2][1] = AQ[0][2] * Q[0][1] + AQ[1][2] * Q[1][1] + AQ[2][2] * Q[2][1];
    D[2][2] = AQ[0][2] * Q[0][2] + AQ[1][2] * Q[1][2] + AQ[2][2] * Q[2][2];

    o[0] = D[1][2];
    o[1] = D[0][2];
    o[2] = D[0][1];
    m[0] = fabs(o[0]);
    m[1] = fabs(o[1]);
    m[2] = fabs(o[2]);

    k0 = (m[0] > m[1] && m[0] > m[2])
             ? 0
             : (m[1] > m[2]) ? 1 : 2;  // index of largest element of offdiag
    k1 = (k0 + 1) % 3;
    k2 = (k0 + 2) % 3;
    if (o[k0] == 0.0) {
      break;  // diagonal already
    }
    thet = (D[k2][k2] - D[k1][k1]) / (2.0 * o[k0]);
    sgn = (thet > 0.0) ? 1.0 : -1.0;
    thet *= sgn;  // make it positive
    t = sgn /
        (thet + ((thet < 1.E6) ? sqrt(thet * thet + 1.0) : thet));  // sign(T)/(|T|+sqrt(T^2+1))
    c = 1.0 / sqrt(t * t + 1.0);                                    //  c= 1/(t^2+1) , t=s/c
    if (c == 1.0) {
      break;  // no room for improvement - reached machine precision.
    }
    jr[0] = jr[1] = jr[2] = jr[3] = 0.0;
    jr[k0] = sgn * sqrt((1.0 - c) / 2.0);  // using 1/2 angle identity sin(a/2) = sqrt((1-cos(a))/2)
    jr[k0] *= -1.0;  // since our quat-to-matrix convention was for v*M instead of M*v
    jr[3] = sqrt(1.0 - jr[k0] * jr[k0]);
    if (jr[3] == 1.0) {
      break;  // reached limits of floating point precision
    }
    q[0] = (q[3] * jr[0] + q[0] * jr[3] + q[1] * jr[2] - q[2] * jr[1]);
    q[1] = (q[3] * jr[1] - q[0] * jr[2] + q[1] * jr[3] + q[2] * jr[0]);
    q[2] = (q[3] * jr[2] + q[0] * jr[1] - q[1] * jr[0] + q[2] * jr[3]);
    q[3] = (q[3] * jr[3] - q[0] * jr[0] - q[1] * jr[1] - q[2] * jr[2]);
    mq = sqrt(q[0] * q[0] + q[1] * q[1] + q[2] * q[2] + q[3] * q[3]);
    q[0] /= mq;
    q[1] /= mq;
    q[2] /= mq;
    q[3] /= mq;
  }
}
typedef Vector3<double> Vector3D;
typedef Vector3<double> Point3D;
typedef Box3<double> Box3D;
typedef Vector3<uint32_t> Color3B;
typedef Matrix3<double> Matrix3D;

template <typename T>
T Clip(const T &n, const T &lower, const T &upper) {
  return std::max(lower, std::min(n, upper));
}
template <typename T>
bool ApproximatelyEqual(T a, T b, T epsilon = std::numeric_limits<double>::epsilon()) {
  return fabs(a - b) <= ((fabs(a) < fabs(b) ? fabs(b) : fabs(a)) * epsilon);
}

#endif /* MATH3_H */
