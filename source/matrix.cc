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
// 
// File:   matrix.cc
// Author: ruehle
//
// Created on July 23, 2007, 5:39 PM
//

#include "matrix.h"

int cjcbi(matrix &a, matrix &v, double eps, int jt)
{ 
    int n=3;
    int i,j,p,q,l;
    double fm,cn,sn,omega,x,y,d;
    l=1;
    
    v.UnitMatrix();
    
    while (true) {
        fm=0.0;
        for (i=1; i<=n-1; i++)
            for (j=0; j<=i-1; j++) {
                d=fabs(a[i][j]);
                if ((i!=j)&&(d>fm))
                    { fm=d; p=i; q=j;}
            }
        if (fm<eps)  return(1);
        if (l>jt)  return(-1);
        l=l+1;
        x=-a[p][q]; y=(a[q][q]-a[p][p])/2.0;
        omega=x/sqrt(x*x+y*y);
        if (y<0.0) omega=-omega;
        sn=1.0+sqrt(1.0-omega*omega);
        sn=omega/sqrt(2.0*sn);
        cn=sqrt(1.0-sn*sn);
        fm=a[p][p];
        a[p][p]=fm*cn*cn+a[q][q]*sn*sn+a[p][q]*omega;
        a[q][q]=fm*sn*sn+a[q][q]*cn*cn-a[p][q]*omega;
        a[p][q]=0.0; a[q][p]=0.0;
        for (j=0; j<=n-1; j++)
        if ((j!=p)&&(j!=q))
          { //u=p*n+j; w=q*n+j;
            fm=a[p][j];
            a[p][j]=fm*cn+a[q][j]*sn;
            a[q][j]=-fm*sn+a[q][j]*cn;
          }
        for (i=0; i<=n-1; i++)
          if ((i!=p)&&(i!=q))
            { //u=i*n+p; w=i*n+q;
              fm=a[i][p];
              a[i][p]=fm*cn+a[i][q]*sn;
              a[i][q]=-fm*sn+a[i][q]*cn;
            }
        for (i=0; i<=n-1; i++)
          { //u=i*n+p; w=i*n+q;
            fm=v[i][p];
            v[i][p]=fm*cn+v[i][q]*sn;
            v[i][q]=-fm*sn+v[i][q]*cn;
          }
      }
    return(1);
}
#include <iostream>
using namespace std;
void matrix::SolveEigensystem(eigensystem_t &out)
{
    matrix m(*this);
    matrix v;
    cjcbi(m, v);

    for(int i=0; i<3; ++i) {
        out.eigenvalues[i] = m[i][i];
        out.eigenvecs[i] = vec(v[0][i], v[1][i], v[2][i]);
        out.eigenvecs[i].normalize();
    }

    //cout << v << endl;
    // sort by eigenvalues
    if(out.eigenvalues[0] > out.eigenvalues[1]) {
        std::swap(out.eigenvalues[0], out.eigenvalues[1]);
        std::swap(out.eigenvecs[0], out.eigenvecs[1]);
    }
    if(out.eigenvalues[1] > out.eigenvalues[2]) {
        std::swap(out.eigenvalues[1], out.eigenvalues[2]);
        std::swap(out.eigenvecs[1], out.eigenvecs[2]);
    }
    if(out.eigenvalues[0] > out.eigenvalues[1]) {
        std::swap(out.eigenvalues[0], out.eigenvalues[1]);
        std::swap(out.eigenvecs[0], out.eigenvecs[1]);
    }    
}
    
