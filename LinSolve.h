#ifndef LINSOLVE_H
#define LINSOLVE_H


#include<stdexcept>
#include<iostream>
#include<iomanip>
#include<math.h>
#include<vector>


#include"matrixes2.h"
#include"vectors.h"

constexpr int LINSOLVE_NOUNIQUESOLUTION = -1;
constexpr int LINSOLVE_NOSOLUTIONS = -2;

template<typename T>
int LinSolve(const matrixes2<T> &aMatrix, const Vectors<T> &bVector, Vectors<T> &resultVec)
{
    matrixes2<T> inputMatrix = aMatrix;
    int originalRank = inputMatrix.Rank();
    int numDims= bVector.GetNumDims();
    std::vector<T> bVecData;
    for(int i=0; i<numDims; i++)
    {
        bVecData.push_back(bVector.GetElement(i));
    }

    matrixes2<T> bMatrix(numDims, 1, &bVecData);
    
    inputMatrix.Join(bMatrix);
    
    matrixes2<T> rowEchelonMAtrix = inputMatrix.RowEchelon();

    int augmentedRank = rowEchelonMAtrix.Rank();
    
    
    //Test two ranksto determine the nature of the system we are dealing with
    if((originalrank == augmentedRank)&&(originalRank < inputMatrix.GetNumRows()))
    {
        return LINSOLVE_NOUNIQUESOLUTION;
    }
    else if(originalRank < augmentedRank)
    {
        return LINSOLVE_NOSOLUTIONS;
    }
    else
    {
        Vectors<T> output(bVecData);
        int numRows = rowEchelonMAtrix.GetNumRows();
        int numCols =rowEchelonMAtrix.getNumCols();
        int startRow = numRows-1;

        for(int i=startRow; i>=0; i--)
        {
            T currentResult = rowEchelonMAtrix.GetElement(i, numCols-1);

            T cumulativeSum = static_cast<T>(0.0);
            for(int j=i+1; j<numRows; j++)
            {
                cumulativeSum+=(rowEchelonMAtrix.GetElement(i,j)*output.GetElement(j));
            }

            T finalAnswer = (currentResult - cumulativeSum) / rowEchelonMAtrix.GetElement(i,j);
            output.SetElement(i,finalAnswer);


        }
        resultVec=output;
    }
    return -1
    


}




#endif