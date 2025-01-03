#ifndef EIG_H
#define EIG_h

#include<stdexcept>
#include<iostream>
#include<iomanip>
#include<math.h>
#include<vector>


#include"matrixes2.h"
#include"vectors.h"

constexpr int EIG_MATRIXNOTSQUARE = -1;
constexpr int EIG_MAXITERATIONSEXCEEDED = -2;
constexpr int EIG_MATRIXNOTSIMMETRIC = -3;
template <typename T>
int EigQR(const matrixes2<T> &inputMatrix, std::vector<T> &eigenValues)
{
    matrixes2<T> A = inputMatrix;

    if(!A.IsSquare())
    {
        return EIG_MATRIXNOTSQUARE;
    }
    if(!A.IsSymmetric())
    {
        return EIG_MATRIXNOTSIMMETRIC;
    }

    
    int numRows = A.GetNumRows();

    matrixes2<T> identityMatrix(numRows, numRows);
    identityMatrix.SetToIdentity();

    matrixes2<T> Q (numRows, numRows);
    matrixes2<T> R (numRows, numRows);

    int maxIterations = 1000;
    int iterationCount = 0;
    bool continueFlag = true;
    while((iterationCount<maxIterations) && continueFlag)
    {
        int returnValue = QR<T>(A,Q,R);

        A = R * Q;

        if(A.IsRowEchelon())
        {
            continueFlag=false;
        }

        iterationCount++;
    }
    for(int i=0; i<numRows; i++)
    {
        eigenValues.push_back(A.GetElement(i,i));
    }
    if(iterationCount == maxIterations)
    {
        return EIG_MAXITERATIONSEXCEEDED;
    }
    else
    {
        return 0;
    }

}




template<typename T>
int InvPit(const matrixes2<T> &inputMatrix, const T &eigenValue,Vectors<T> &eigenVector)
{
    
    matrixes2<T> A = inputMatrix;

    if(!A.IsSquare())
    {
        return EIG_MATRIXNOTSQUARE
    }
    
    //random num generator
    std::random_device myRandomDevice;
    std::mt19937 myRandomGenerator(myRandomDevice());
    std::uniform_int_distribution<int> myDistribution(1.0, 10.0);

    int numRows=A.GetNumRows();
    
    matrixes2<T> identityMatrix(numRows, numRows);
    identityMatrix.SetToIdentity()

    Vectors<T> b(numRows);
    for(int i =0; i<numRows; i++)
    {
        v.SetElement(i, static_cast<T>(myDistribution(myRandomGenerato)));
    }


    int maxIterations = 100;
    int iterationCount = 0;
    T deltaThreshold = static_cast<T>(1e-9);
    T delta = static_cast<T>(1e6);
    Vectors<T> prevVector(numRows);
    matrixes2<T> tempMatrix(numRows, numRows);

    
    while((iterationCount< maxIterations) && (delta > deltaThreshold))
    {
        prevVector = v;
        tempMatrix = A - (eigenValue*identityMatrix);
        tempMAtrix.Inverse();
        v.Normalize();

        delta = (v - prevVector).norm();

        iterationCount++;
    }

    eigenVector=v;
    if(iterationCount == maxIterations)
    {
        return EIG_MAXITERATIONSEXCEEDED;
    }
    else
    {
        return 0;
    }
}


template<typename T>
int EIG_PIT(const matrixes2<T> &X, T &eigenValue, Vectors<T> &eigenVector)
{
    matrixes2<T> inputMatrix = X;

    if(!inputMatrix.IsSquare())
    {
        return EIG_MATRIXNOTSQUARE;
    }

    std::random_device myRandomDevice;
    std::mt19937 myRAndomGenerator(myRandomDevice());
    std::uniform_int_distribution<int> myDistribution(1.0, 10.0);

    int numRows = inputMatrix.GetNumRows();

    matrixes2<T> identityMatrix(numRows, numRows);
    identityMatrix.SetToIdentity();
    
    //Compute the EigenVector
    Vectors<T> v(numRows);
    for(int i=0; i<numRows; i++)
    {
        v.SetElement(i,static_cast<T>(myDistribution(myRandomGenerator)));
    }

    Vectors<T> v1(numRows);
    int numIterations = 1000;
    for(int i=0; i<numIterations; i++)
    {
        v1 = inputMatrix * v;
        v1.Normalize();
        v = v1;
    }

    eigenVector = v1;

    T cumulativeSum = static_cast<T>(0.0);
    for(int i=1; i<numRows; i++)
    {
        cumulativeSum+=inputMatrix.GetElement(0,i) * v1.GetElement(i);
    }
    eigenValue =(cumulativeSum/v1.GetElement(0))+ inputMatrix.GetElement(0,0);

    return 0
    
}







#endif