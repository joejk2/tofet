///////////////////////////////////////////////////////////////////////
//  This file is part of ToFeT.
//  
//  ToFeT is free software: you can redistribute it and/or modify
//  it under the terms of the GNU Lesser General Public License as published by
//  the Free Software Foundation, either version 3 of the License, or
//  (at your option) any later version.
//  
//  ToFeT is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU Lesser General Public License for more details.
//  
//  You should have received a copy of the GNU Lesser General Public License
//  along with ToFeT.  If not, see <http://www.gnu.org/licenses/>.
///////////////////////////////////////////////////////////////////////
/* 
 * Author: Victor Ruehle
 * Created on July 11, 2007, 3:24 PM
*/

#ifndef _mat_H
#define	_mat_H

#include "types.h"
#include "ostream"
#include "vec.h"

class matrix
{
public:
        
    matrix() {};
    matrix(const real &v) { *this=v; }
    matrix(const matrix &m) { *this=m; }
    

    matrix &operator=(const real &v);
    matrix &operator=(const matrix &v);
    //vec &operator+=(const vec &v);
    //vec &operator-=(const vec &v);
    matrix &operator*=(const double &d){
        for(size_t i=0; i<9; ++i) _m[i] *=d;
        return *this;
    }
    matrix &operator/=(const double &d){
        for(size_t i=0; i<9; ++i) _m[i] /=d;
        return *this;
    }  
    matrix &operator-=(const matrix &v){
        for(size_t i=0; i<9; ++i) _m[i] -= v._m[i];
        return *this;
    }
    matrix &operator+=(const matrix &v){
        for(size_t i=0; i<9; ++i) _m[i] += v._m[i];
        return *this;
    }  
    
    void UnitMatrix();
    
    void set(const byte_t &i, const byte_t &j, const real &v) { _m[i*3+j] = v; }
    const real &get(const byte_t &i, const byte_t &j) const { return _m[i*3+j]; }
    
    real *operator[](size_t i) { return &_m[i*3]; }       
    
    struct eigensystem_t {
        real eigenvalues[3];
        vec eigenvecs[3];
    };
    void SolveEigensystem(eigensystem_t &out);
    
    matrix &Transpose(){
        std::swap( _m[1], _m[3]);
        std::swap( _m[2], _m[6]);
        std::swap( _m[5], _m[7]);
        return *this;
    }
    ///matrix multiplication
    matrix  operator * (const matrix & a){
        matrix r;
        r._m[0] = _m[0] * a._m[0] + _m[1] * a._m[3] + _m[2] * a._m[6];
        r._m[1] = _m[0] * a._m[1] + _m[1] * a._m[4] + _m[2] * a._m[7];
        r._m[2] = _m[0] * a._m[2] + _m[1] * a._m[5] + _m[2] * a._m[8];
        
        r._m[3] = _m[3] * a._m[0] + _m[4] * a._m[3] + _m[5] * a._m[6];
        r._m[4] = _m[3] * a._m[1] + _m[4] * a._m[4] + _m[5] * a._m[7];
        r._m[5] = _m[3] * a._m[2] + _m[4] * a._m[5] + _m[5] * a._m[8];
        
        r._m[6] = _m[6] * a._m[0] + _m[7] * a._m[3] + _m[8] * a._m[6];
        r._m[7] = _m[6] * a._m[1] + _m[7] * a._m[4] + _m[8] * a._m[7];
        r._m[8] = _m[6] * a._m[2] + _m[7] * a._m[5] + _m[8] * a._m[8];
        
        return r;
    }
       
    vec operator * ( const vec & a){
       return vec( _m[0] * a.getX() + _m[1] * a.getY() + _m[2] * a.getZ(), 
                _m[3] * a.getX() + _m[4] * a.getY() + _m[5] * a.getZ(),
                _m[6] * a.getX() + _m[7] * a.getY() + _m[8] * a.getZ() ); 
    }
    
private:
    real _m[9];
};

inline matrix &matrix::operator=(const real &v)
{
    for(size_t i=0; i<9; ++i)
        _m[i] = v;
    return *this;
}

inline matrix &matrix::operator=(const matrix &m)
{
    for(size_t i=0; i<9; ++i)
        _m[i] = m._m[i];
    return *this;
}

inline void matrix::UnitMatrix()
{
    for(size_t i=0; i<9; ++i)
        _m[i] = 0.;//(*this) = 0.;
    _m[0] = _m[4] = _m[8] = 1.0;
}

inline std::ostream &operator<<(std::ostream &out, matrix& m)
{
      out << '|' << m[0][0] << ',' << m[0][1] << ',' << m[0][2] << '|' << std::endl;
      out << '|' << m[1][0] << ',' << m[1][1] << ',' << m[1][2] << '|' << std::endl;
      out << '|' << m[2][0] << ',' << m[2][1] << ',' << m[2][2] << '|' << std::endl;
      return out;
}

inline matrix operator|( const vec & a, const vec & b){
    matrix res;
    res.set(0,0, a.getX() * b.getX());
    res.set(0,1, a.getX() * b.getY());
    res.set(0,2, a.getX() * b.getZ());
    res.set(1,0, a.getY() * b.getX());
    res.set(1,1, a.getY() * b.getY());
    res.set(1,2, a.getY() * b.getZ());
    res.set(2,0, a.getZ() * b.getX());
    res.set(2,1, a.getZ() * b.getY());
    res.set(2,2, a.getZ() * b.getZ());
    return res;
}


inline matrix operator*(const matrix & r, const double &d){
    return ( matrix(r) *= d);
}
inline matrix operator/(const matrix & r,const double &d){
    return ( matrix(r) /= d);
}  
inline matrix operator+(const matrix & r, const matrix & v){
    return ( matrix(r) += v);
}
inline matrix operator-(const matrix & r, const matrix & v){
    return ( matrix(r) -= v);
}
    

/* provided the matrix a diagonalizes it and returns the eigenvalues 
   lambda0 = a[0][0], lambda1 = a[1][1], lambda2= a[2][2], ...
   as well as the corresponding eigenvectors v[][0], v[][1], v[][2] 
*/
int cjcbi(matrix &a, matrix &v, double eps=1e-10, int jt=100);



#endif	/* _mat_H */

