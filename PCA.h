#ifndef PCA_H
#define PCA_H


#include<stdexcept>
#include<iostream>
#include<math.h>
#include<iomanip>
#include<vector>
#include<algorithm>

#include"EIG.h"
#include"matrixes2.h"
#include"vectors.h"

constexpr int PCA_MATRIXNOTSQUARE = -1;
constexpr int PCA_MATRIXNOTSYMMETRIC = -2;

namespace PCA
{
template <typename T>
std::vector<T> ComputeColumnMeans(const matrixes2<T> &inputData)
{
    int numRows=inputData.GetNumRows();
    int numCols=inputData.getNumCols();

    std::vector<T> output;

    for(int j=0; j<numCols; j++)
    {
        T cumulativeSum = static_cast<T>(0.0);
        for(int i=0; i<numRows; i++)
        {
            cumulativeSum += inputData.GetElement(i,j);
        }
        output.push_back(cumulativeSum/static_cast<T>(numRows));
    }
    return output;
}

template<typename T>
void SubtractColumnMeans(matrixes2<T> &inputData, std::vector<T> &columnMeans)
{
    int numRows = inputData.GetNumRows();
    int numCols = inputData.getNumCols();

    for(int j=0; j<numCols; j++)
    {
        for(int i=0; i<numRows; i++)
        {
            inputData.SetElement(i,j,inputData.GetElement(i,j) - columnMeans.at(j));
        }
    }
}


template<typename T>
matrixes2<T> ComputeCovariance(const matrixes2<T> &X)
{
    int numRows = X.GetNumRows();
    matrixes2<T> covX = (static_cast<T>(1.0)/static_cast<T>(numRows - 1)) * (X.Transpose() * X);
    return covX;
}


template<typename T>
int ComputeEigenVectors(const matrixes2<T> &covarianceMatrix, matrixes2<T> &eigenVectors)
{
    matrixes2<T> X = covarianceMatrix;
    
    if(!X.IsSquare())
    {
        return PCA_MATRIXNOTSQUARE;
    }

    if(!X.IsSymmetric())
    {
        return PCA_MATRIXNOTSYMMETRIC;
    }

    std::vector<T> eigenValues;
    int returnStatus = EigQR(X, eigenValues);

    std::sort(eigenValues.begin(), eigenValues.end());
    std::reverse(eigenValues.begin(), eigenValues.end());

    Vectors<T> eV(X.getNumCols());
    matrixes2<T> eVM(X.GetNumRows(), X.getNumCols());
    
    for(int j=0; j<eigenValues.size(); +=j)
    {
        T eig =eigenvalues.at(j);
        int returnStatus2 = InvPit<T>(X, eig, eV);
        for(int i=0; i<eV.GetNumDims(); i++)
        {
            eVM.SetElement(i,j,eV.GetElement(i));
        }
    }

    eigenVectors = eVM;

    return returnStatus;
}

template<typename T>
int PCA(const matrixes2<T> &inputData, matrixes2<T> &outputComponents)
{
    
    matrixes2<T> X =inputData;
    
    std::vector<T> columnMeans = ComputeColumnMeans(X);
    
    SubtractColumnMeans<T>(X,columnMeans);
    
    matrixes2<T> covX = ComputeEigenVectors(covx, eigenVectors);
    
    matrixes2<T> eigenVectors;
    
    
    int returnStatus = ComputeEigenVectors(covX, eigenVectors);
    outputComponents = eigenVectors;

    return returnStatus;
}

}














#endif