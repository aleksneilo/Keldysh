#pragma once
#include <algorithm>
#include <cmath>
#include <complex>
#include <stdexcept>
#include <vector>

namespace sns {
using Complex = std::complex<double>;
// Dense column-major matrix; Floquet indices are translated by +NF.
struct Matrix {
    int rows = 0, cols = 0;
    std::vector<Complex> data;
    Matrix() = default;
    Matrix(int r, int c) : rows(r), cols(c), data(checked_size(r,c)) {}
    static size_t checked_size(int r, int c) {
        if (r < 0 || c < 0) throw std::invalid_argument("negative matrix size");
        return size_t(r)*size_t(c);
    }
    Complex& operator()(int r, int c) { return data.at(size_t(c)*rows+r); }
    Complex operator()(int r, int c) const { return data.at(size_t(c)*rows+r); }
    static Matrix identity(int n) { Matrix a(n,n); for(int j=0;j<n;++j) a(j,j)=1.; return a; }
};
inline void same_shape(const Matrix& a, const Matrix& b) {
    if(a.rows!=b.rows || a.cols!=b.cols) throw std::invalid_argument("matrix shape mismatch");
}
inline Matrix operator+(Matrix a,const Matrix& b) { same_shape(a,b); for(size_t j=0;j<a.data.size();++j) a.data[j]+=b.data[j]; return a; }
inline Matrix operator-(Matrix a,const Matrix& b) { same_shape(a,b); for(size_t j=0;j<a.data.size();++j) a.data[j]-=b.data[j]; return a; }
inline Matrix operator*(Matrix a,Complex s) { for(auto& v:a.data) v*=s; return a; }
inline Matrix operator*(Complex s,Matrix a) { return a*s; }
inline Matrix operator*(const Matrix& a,const Matrix& b) {
    if(a.cols!=b.rows) throw std::invalid_argument("matrix product mismatch");
    Matrix c(a.rows,b.cols);
    for(int j=0;j<b.cols;++j) for(int k=0;k<a.cols;++k) {
        Complex v=b(k,j); if(v==Complex{}) continue;
        for(int i=0;i<a.rows;++i) c(i,j)+=a(i,k)*v;
    } return c;
}
inline Matrix adjoint(const Matrix& a) { Matrix b(a.cols,a.rows); for(int j=0;j<a.cols;++j) for(int i=0;i<a.rows;++i) b(j,i)=std::conj(a(i,j)); return b; }
inline double norm(const Matrix& a) { double v=0; for(auto z:a.data) v=std::hypot(v,std::abs(z)); return v; }
inline bool finite(const Matrix& a) { for(auto z:a.data) if(!std::isfinite(z.real()) || !std::isfinite(z.imag())) return false; return true; }
// Partial-pivot Gaussian elimination, multiple RHS. No explicit matrix inversion.
inline Matrix solve_dense(Matrix a,Matrix b) {
    if(a.rows!=a.cols || a.rows!=b.rows) throw std::invalid_argument("linear solve shape mismatch");
    const int n=a.rows; double scale=0; for(auto z:a.data) scale=std::max(scale,std::abs(z));
    for(int k=0;k<n;++k) {
        int p=k; for(int i=k+1;i<n;++i) if(std::abs(a(i,k))>std::abs(a(p,k))) p=i;
        if(std::abs(a(p,k))<=1e-14*scale || !finite(a)) throw std::runtime_error("singular/nonfinite linear system");
        for(int j=k;j<n;++j) std::swap(a(k,j),a(p,j));
        for(int j=0;j<b.cols;++j) std::swap(b(k,j),b(p,j));
        for(int i=k+1;i<n;++i) { Complex q=a(i,k)/a(k,k); a(i,k)=0.;
            for(int j=k+1;j<n;++j) a(i,j)-=q*a(k,j);
            for(int j=0;j<b.cols;++j) b(i,j)-=q*b(k,j);
        }
    }
    for(int k=n-1;k>=0;--k) for(int j=0;j<b.cols;++j) {
        for(int i=k+1;i<n;++i) b(k,j)-=a(k,i)*b(i,j);
        b(k,j)/=a(k,k);
    }
    if(!finite(b)) throw std::runtime_error("nonfinite linear solution");
    return b;
}
inline Matrix block(const Matrix& a,int r,int c,int n) { Matrix b(n,n); for(int j=0;j<n;++j) for(int i=0;i<n;++i) b(i,j)=a(r+i,c+j); return b; }
inline void set_block(Matrix& a,int r,int c,const Matrix& b) { for(int j=0;j<b.cols;++j) for(int i=0;i<b.rows;++i) a(r+i,c+j)=b(i,j); }
inline Matrix nambu(const Matrix& a,const Matrix& b,const Matrix& c,const Matrix& d) { int k=a.rows; Matrix v(2*k,2*k); set_block(v,0,0,a);set_block(v,0,k,b);set_block(v,k,0,c);set_block(v,k,k,d);return v; }
} // namespace sns
