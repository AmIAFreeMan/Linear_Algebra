#ifndef QR_H
#define QR_H

#include<stdexcept>
#include<iostream>
#include<iomanip>
#include<math.h>
#include<vector>

#include"matrixes2.h"
#include"vectors.h"

constexpr int QR_MATRIXNOTSQUARE = -1;
template<typename T>
int QR(const matrixes2<T> &A, matrixes2<T> &Q, matrixes2<T> &R)
{

    matrixes2<T> inputMatrix = A;

    if(!inputMatrix.IsSquare())
    {
        return QB_MATRIXNOTSQUARE;
    }

    int numCols = inputMAtrix.GetNumCols();

    std::vector<matrixes2<T>> Plist;

    for(int j=0; j<(numCols-1); j++)
    {
        Vectors<T> a1 (numCols-j);
        Vectors<T> b1 (numCols-j);
        for(int i=j; i<numCols; i++)
        {
            a1.SetElement(i-j, inputMatrix.GetElement(i,j));
            b1.SetElement(i-j, static_cast<T>(0.0));
        }
        b1.SetElement(0, static_cast<T>(1.0));

        T a1norm = a1.norm();

        int sgn = -1;
        if(a1.GetElement(0) < static_cast<T>(0.0))
        {
            sgn = 1;
        }
        Vectors<T> u = a1-(sgn*a1norm*b1);

        Vector<T> n =u.Normalized();

        matrixes2<T> n =u.Normalized();

        matrixes2<T> nMat (numCols - j, 1);
        for(int i=0; i<(numCols-j); i++)
        {
            nMat.SetElement(i,0,n.GetElement(i));
        }

        matrixes2<T> nMatT=nMat.Transpose();

        matrixes2<T> I (numCols-j, numCols-j);
        I.SetToIdentity();
       
        matrixes2<T> Ptemp = I - static_cast<T>(2.0) * nMat *nMatT;
        
        matrixes2<T> P (numCols, numCols);        
        P.SetToIdentity();
        for(int row=j; row<numCols; row++)
        {
            for(int col=j; col<numCols; ++col)
            {
                P.SetElement(row, col, Ptemp.GetElement(row-j, col-j));
            }
        }

        Plist.push_back(P);

        inputMatrix = P * inputMatrix;
    } 
    
    
    //compute Q
    matrixes2<T> Qmat = Plist.at(0);
    for(int i=1; i<(numCols -1); i++)
    {
        Qmat = Qmat * Plist.at(i).Transpose();
    }
    Q = Qmat;
    
    
    //compute R
    int numElements = Plist.size();
    matrixes2<T> Rmat = Plist.at(numElements-1);
    for(int i =(numElements-2); i>=0; i--)
    {
        Rmat = Rmat * Plist.at(i);
    }
    Rmat = Rmat * A
    R = Rmat;
}


#endif