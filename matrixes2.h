#ifndef MATRIXES2_H
#define MATRIXES2_H

#include <stdexcept>
#include <iostream>
#include <iomanip>
#include <math.h>
#include <vector>

#include "vectors.h"


template <class T>
class matrixes2
{
  public:
      
      
    // Define constructors
    matrixes2();
    matrixes2(int nRows, int nCols);
    matrixes2(int nRows, int nCols, const T *inputData);
    matrixes2(const matrixes2<T>& inputMatrix);
    matrixes2(int nRows, int nCols, const std::vector<T> *inputData);
      
      
    //destructur
    ~matrixes2();

      
    //configuration methods
    bool Resize(int numRows, int numCols);
    void SetToIdentity();

      
      
    //elements acces methods
    T GetElement(int row, int col);
    bool SetElement(int row, int col, T elementValue);
    int GetNumRows();
    int getNumCols();

      

    bool Inverse();

    matrixes2<T> RowEchelon();

    matrixes2<T> Transpose();

    T Determinant();

      

      
    //overload == operator
    bool operator== (const matrixes2<T>& rhs);
    bool Compare (const matrixes2<T>& matrix1, double tollerance);
      
    //overload the  assign operator
    matrixes2<T> operator= (const matrixes2<T> &rhs)
      
    //overload +, - and * operators (friends)
    template <class U> friend matrixes2 operator+ (const matrixes2<U>& lhs, const matrixes2<U>& rhs);
    template <class U> friend matrixes2 operator+ (const U& lhs, const matrixes2<U>& rhs);
    template <class U> friend matrixes2 operator+ (const matrixes2<U>& lhs, const U& rhs);

    template <class U> friend matrixes2 operator- (const matrixes2<U>& lhs, const matrixes2<U>& rhs);
    template <class U> friend matrixes2 operator- (const U& lhs, const matrixes2<U>& rhs);
    template <class U> friend matrixes2 operator- (const matrixes2<U>& lhs, const U& rhs);

    template <class U> friend matrixes2 operator* (const matrixes2<U>& lhs, const matrixes2<U>& rhs);
    template <class U> friend matrixes2 operator* (const U& lhs, const matrixes2<U>& rhs);
    template <class U> friend matrixes2 operator* (const matrixes2<U>& lhs, const U& rhs);

    template <class U> friend Vectors<U> operator* (const matrixes2<U>& lhs, const Vectors<U>& rhs);

    bool Join(const matrixes2<T>& matrix2);
    bool Separete(matrixes2<T> *matrix1, matrixes2<T> *matrix2, int colNum);
      
    bool IsRowEchelon();
    bool IsSquare();
    bool IsSymmetric();

    bool IsNonZero();

    static int Rank(const matrixes2<T> &inputMatrix);
      
    
    
  private:
    int Sub2Ind(int row, int col);
    bool CloseEnough(T f1, T f2);
    void SwapRow(int i, int j);
    void MultAdd(int i, int j, T multfactor);
    void multRow(int i, T multfactor);
    int FindRowWhitMaxElement(int colNumber, int startingRow);
    matrixes2<T> FindSubMatrix(int rowNum, int colNum);
    void PrintMatrix();
    
    
  private:
    T *m_matrixData;
    int m_nRows, m_nCols, m_nElements;
    
};



/* **************************************************************************************************************
CONSTRUCTOR / DESTRUCTOR FUNCTIONS
*****************************************************************************************************************/



// The default constructor .
template <class T>
matrixes2<T>::matrixes2()
{
  m_nRows = 1;
  m_nCols = 1;
  m_nElements = 1;
  m_matrixData = new T[m_nElements];
  m_matrixData[0] = 0.0;
}





//Construct empty matrix (all elements 0)
template <class T>
matrixes2<T>::matrixes2(int nRows, int nCols)
{
  m_nRows = nRows;
  m_nCols = nCols;
  m_nElemnts = m_nRows * m_nCols;
  m_matrixData = new T[m_nElements];
  for(int i = 0; i<m_nElements; i++)
  {
    m_matrixData[i] = 0.0;
  }
}



//Construct from const linear array
template <class T>
matrixes2<T>::matrixes2(int nRows, int nCols, const T *inputData)
{
  m_nRows = nRows;
  m_nCols = nCols;
  m_nElemnts = m_nRows * m_nCols;
  m_matrixData = new T[m_nElements];
  for(int i = 0; i<m_nElements; i++)
  {
    m_matrixData[i] = inputData[i];
  }
}




//Copy constructor
template<class T>
matrixes2<T>::matrixes2(const matrixes2<T>& inputMatrix)
{
  m_nRows = inputMatrix.m_nRows;
  m_nCols = inputMatrix.m_nCols;
  m_nElements = inputMatrix.m_nElements;
  m_matrixData = new T[m_nElements];
  for(int i = 0; i<m_nElements; i++)
  {
    m_matrixData[i] = inputData[i];
  }
}



// counstructor from std::vector
template<class T>
matrixes2<T>::matrixes2(int nRows, int nCols, const std::vector<T> *inputData)
{
  m_nRows = nRows;
  m_nCols = nCols;
  m_nElements = m_nRows * m_nCols;
  m_nmatrixData = new T[m_nElements];
  for(int i=0; i<m_nElemnts; i++)
  {
    m_nmatrixData[i] = inputData->at(i);
  }
}




//Destructoring
template <class T>
matrixes2<T>::~matrixes2()
{
  if(m_matrixData)
  {
    delete[] m_matrixData;
  }
  m_matrixData = nullptr;
}




// converts existing matrix in identity matrix
template<class T>
void matrixes2<T>::SetToIdentity()
{
  if(!isSquare())
  {
    throw std::invalid_argument("Cannot form an identity matrix that is not square");
  }
  for(int row=0; row<m_nRows; row++)
  {
    for(int col=0; col<m_nCols; ++col)
    {
      if(col == row)
      {
        m_matrixData[Sub2Ind(row,col)] = 1.0;
      }
      else
      {
        m_matrixData[Sub2Ind(row col)] = 0.0;
      }
    }
  }
}









/*********************************************************************************************************
CONFIGURATION FUNCTIONS
**********************************************************************************************************/

template <class T>
bool matrixes2<T>::Resize(int numRows, int numCols)
{
  m_nRows = numRows;
  m_nCols = numCols;
  m_nElements = (m_nRows * m_nCols);
  delete[] m_matrixData;
  m_matrixData = new T[m_nElements];
  if(m_matrixData != nullptr)
  {
    for(int i=0; i<m_nElements; i++)
    {
      m_matrixData[i] = 0.0;
    }
    return true;
  }
  else
  {
    return false;
  }
}



/**********************************************************************************************************
ELEMENT FUNCTIONS 
***********************************************************************************************************/

template <class T>
T matrixes2<T>::GetElement(int row, int col)
{
  int linearIndex = Sub2Ind(row, col);
  if(linearIndex >= 0)
  {
    return m_matrixData[linearIndex];
  }
  else
  {
    return 0.0;
  }
}


template <class T>
bool matrixes2<T>::SetElement(int row, int col, T elementValue)
{
  int linearIndex = Sub2Ind(row, col);
  if(linearIndex >= 0)
  {
    m_matrixData[linearIndex] = elementValue;
    return true;
  }
  else
  {
    return false;
  }
}

template <class T>
int matrixes2<T>::GetNumRows()
{
  return m_nRows;
}

template<class T>
int matrixes2<T>::getNumCols()
{
  return m_nCols;
}

template<class T>
bool matrixes2<T>::Compare(const matrixes2<T>& matrix1, double tollerance)
{
  // chech if matrixes have the same dimensions
  int numRows1 = matrix1.m_nRows;
  int numCols1 = matrix1.m_nCols;
  if((numRows1 != m_nRows) || (numCols1 != m_nCols))
  {
    return false;
  }

  // loop over all elements and compute the sum ofsquared differences
  double comulativeSum = 0.0;
  for(int i = 0; i<m_nElements; i++)
  {
    T elements1 = matrix1.m_matrixData[i];
    T elements2 = m_matrixData[i];
    cumulativeSum += ((element1 - element2) * (elements1 - elements2));  
  }
  double finalValue = sqrt(cumulativeSum / ((numRows1 * numCols1)-1));
  if(finalValue < tollerance)
  {
    return true;
  }
  else
  {
    return false;
  }
}


/*************************************************************************************
OVERLOAD OPERATOR FUNCTIONS
************************************************************************************/





/****************************************************************************************
THE + OPERATOR
*****************************************************************************************/

// matrix + matrix
template <class T>
matrixes2<T> operator+ (const matrixes2<T>& lhs, const matrixes2<T>& rhs)
{
  int numRows =lhs.m_nRows;
  int numCols = lhs.m_nCols;
  int numElements = numRows * numCols;
  T *tempResult = new T[numElements];
  for(int i = 0; i<numElements; i++)
  {
    tempResult[i] = lhs.m_matrixData[i] + rhs.m_matrixData[i];
  }
  matrixes2<T> result(numRows, numCols, tempResult);
  delete[] tempResult;
  return result;
}



//scaler + matrix
template <class T>
matrixes2<T> operator+ (const matrixes2<T>& lhs, const matrixes2<T>& rhs)
{
  int numRows =lhs.m_nRows;
  int numCols = lhs.m_nCols;
  int numElements = numRows * numCols;
  T *tempResult = new T[numElements];
  for(int i = 0; i<numElements; i++)
  {
    tempResult[i] = lhs + rhs.m_matrixData[i];
  }
  matrixes2<T> result(numRows, numCols, tempResult);
  delete[] tempResult;
  return result;
}



//matrix + scaler
template <class T>
matrixes2<T> operator+ (const matrixes2<T>& lhs, const matrixes2<T>& rhs)
{
  int numRows =lhs.m_nRows;
  int numCols = lhs.m_nCols;
  int numElements = numRows * numCols;
  T *tempResult = new T[numElements];
  for(int i = 0; i<numElements; i++)
  {
    tempResult[i] = lhs.m_matrixData[i] + rhs;
  }
  matrixes2<T> result(numRows, numCols, tempResult);
  delete[] tempResult;
  return result;
}




/*************************************************************************************
THE - OPERATOR
*************************************************************************************/

// matrix - matrix
template <class T>
matrixes2<T> operator- (const matrixes2<T>& lhs, const matrixes2<T>& rhs)
{
  int numRows =lhs.m_nRows;
  int numCols = lhs.m_nCols;
  int numElements = numRows * numCols;
  T *tempResult = new T[numElements];
  for(int i = 0; i<numElements; i++)
  {
    tempResult[i] = lhs.m_matrixData[i] - rhs.m_matrixData[i];
  }
  matrixes2<T> result(numRows, numCols, tempResult);
  delete[] tempResult;
  return result;
}

//scaler - matrix
template <class T>
matrixes2<T> operator- (const T& lhs, const matrixes2<T>& rhs)
{
  int numRows =lhs.m_nRows;
  int numCols = lhs.m_nCols;
  int numElements = numRows * numCols;
  T *tempResult = new T[numElements];
  for(int i = 0; i<numElements; i++)
  {
    tempResult[i] = rhs - lhs.m_matrixData[i];
  }
  matrixes2<T> result(numRows, numCols, tempResult);
  delete[] tempResult;
  return result;
}

//matrix - scaler
template <class T>
matrixes2<T> operator+ (const matrixes2<T>& lhs, const matrixes2<T>& rhs)
{
  int numRows =lhs.m_nRows;
  int numCols = lhs.m_nCols;
  int numElements = numRows * numCols;
  T *tempResult = new T[numElements];
  for(int i = 0; i<numElements; i++)
  {
    tempResult[i] = lhs.m_matrixData[i] - rhs;
  }
  matrixes2<T> result(numRows, numCols, tempResult);
  delete[] tempResult;
  return result;
}


/*******************************************************************************************************
THE * OPERATOR
*******************************************************************************************************/
//matrix * vector
template <class T>
matrixes2<T> operator* (matrixes2<T>& lhs, const Vectors<T>& rhs)
{
  if(lhs.m_nCols != rhs.GetNumDims())
  {
    throw std::invalid_argument("Number of columns in matrix must be equal to the number of rows in vector.");
  }

  Vectors<T> result(lhs.m_nRows);

  for(int row=0; row<lhs.m_nRows; ++row)
  {
    T cumulativeSum += static_cast<T>(0.0);
    for(int col=0; col<lhs.m_nCols; ++col)
    {
      cumulativeSum += (lhs.GetElement(row, col) * rhs.GetElement(col));
    }
    result.SetElement(row, cumulativeSum);
  }
  return result;
}






//scaler * matrix
template <class T>
matrixes2<T> operator* (const T& lhs, const matrixes2<T>& rhs)
{
  int numRows =lhs.m_nRows;
  int numCols = lhs.m_nCols;
  int numElements = numRows * numCols;
  T *tempResult = new T[numElements];
  for(int i = 0; i<numElements; i++)
  {
    tempResult[i] = rhs * lhs.m_matrixData[i];
  }
  matrixes2<T> result(numRows, numCols, tempResult);
  delete[] tempResult;
  return result;
}

//matrix * scaler
template <class T>
matrixes2<T> operator* (const T& lhs, const matrixes2<T>& rhs)
{
  int numRows =lhs.m_nRows;
  int numCols = lhs.m_nCols;
  int numElements = numRows * numCols;
  T *tempResult = new T[numElements];
  for(int i = 0; i<numElements; i++)
  {
    tempResult[i] = rhs.m_matrixData[i] * lhs;
  }
  matrixes2<T> result(numRows, numCols, tempResult);
  delete[] tempResult;
  return result;
}


//matrix * matrix
template <class T>
matrixes2<T> operator* (const matrixes2<T>& lhs, const matrixes2<T>& rhs)
{
  int l_numRows =lhs.m_nRows;
  int l_numCols = lhs.m_nCols;
  
  int r_numRows =rhs.m_nRows;
  int r_numCols = rhs.m_nCols;
  
  //standard matrix multiplication condition
  //size of output == size of RHS
  if(l_numCols == r_numRols)
  {
    T *tempResult = new T[lhs.m_nRows * rhs.m_nCols];
    //Loop through each row of LHS
    for(int lhsRow = 0; lhsRow<l_numRows; lhsRow++)
      {
      //Loop through each row of RHS
      for(int rhsCol=0; rhsCol<r_numCols; rhsCol++)
        {
          T elementResult = 0.0;
        //Loop through each of THIS LHS row
        for (int lhsCol =0; lhsCol<l_numCols; lhsCol++)
        {
          //compute LHS linear index
          int lhsLinearIndex = (lhsRow * l_numCols) + lhsCol;
          
          //compute the RHS linear index (based on LHS col)
          //rhs row number == to this column number
          int rhsLinearIndex = (lhsCol * r_numCols) + rhsCol;
          
          //perform the calculation on those elements
          
          elementResult += (lhs.matrixData[lhsLinearIndex] * rhs.m_matrixData[rhsLinearIndex]);
        }
        //result
        int resultLinearIndex = (lhsRow * r_numCols) + rhsCol;
        tempResult[resultLinearIndex] = elementResult;
      
      }
    }
    matrixes2<T> result(l_numRows, r_numCols, tempResult);
    delete[] tempResult;
    return result;
  }
  else
  {
  matrixes2<T> result(1,1);
  return result;
  }
} 

/***************************************************************************************************************
THE == OPERATOR
****************************************************************************************************************/


template <class T>
bool matrixes2<T>::operator== (const matrixes2<T>& rhs)
{
  if((this->m_nRows != rhs.m_nRows) && (this->m_nCols != rhs.m_nCols))
  {
    return false;
  }

  bool flag = true;
  for(int i=0; i<this->m_nElements; i++)
  {
    if(!CloseEnough(this->m_matrixData[i], rhs.m_matrixData[i]))
    {
      flag = false;
    }
  }
  return flag;
}

/***************************************************************************************************************
THE ASSIGNMENT (=) OPERATOR
****************************************************************************************************************/
template<class T>
matrixes2<T> matrixes2<T>::operator= (const matrixes2<T> &rhs)
{
  if(this-> &rhs)
  {
    m_nRows=rhs.m_nRows;
    m_nCols=rhs.m_nCols;
    m_nElements=rhs.m_nElements;

    if(m_matrixData)
    {
      delete[] m_matrixData;
    }

    m_matrixData=new T[m_nElements];
    for(int i=0; i<m_nElements; i++)
    {
      m_matrixData[i] = rhs.m_matrixData[i];
    }

    return *this;
  }
}




/************************************************************************************************************
SEPARETE THE MATRIXE INTO TWO PARTS, AROUND THE COLUMN NUMBER PROVIDED
(Note that the output is returned into the matrixes2<T> pointers in the input argument lists) 
*************************************************************************************************************/
template<class T>
bool matrixes2<T>::Separete(matrixes2<T> *matrix1, matrixes2<T> *matrix2, int colNum)
{
  //compute the size of new matrixes
  int numRows = m_nRows;
  int numCols1 = m_nCols;
  int numCols2 = m_nCols - colNum;


  //resize two matrices to proper dimensions
  matrix1->Resize(numRows, numCols1);
  matrix2->Resize(numRows, numCols1);

  //Loop over original matrix and store into appropriate elemnts of the two
  //output matrices
  for(int row = 0; row<m_nRows; row++)
  {
    for(int col =0; col<m_nCols; col++)
    {
      if(col<colNum)
      {
        matrix1->SetElement(row,col,this->GetElement(row,col));
      }
      else
      {
        matrix2->SetElement(row, col-colNum, this->GetElement(row, col));
      }
    }
  }
  return true;
}

/***************************************************************************************************************
JOIN TWO MATRIXES TOGETHER
****************************************************************************************************************/

template<class T>
bool matrixes2<T>::Join(const matrixes2<T>& matrix2)
{
  int numRows1 = m_nRows;
  int numRows2 = matrix2.m_nRows;
  int numCols1 = m_nCols;
  int numCols2 = matrix2.m_nCols;


  //if matrices have different numbers of rows, then this operation makes no sense
  if(numRows1 != numRows2)
  {
    throw std::invalid_argument("Attempt to join matrices with different numbers of rows is invalid.");
  }

  //Allocate memory for the result
  //Note that only the number of columns increases
  T*newmatrixData = new T[numRows1*(numCols1+numCols2)];
  
  //Copy the two matrices into the new one
  int linearIndex, resultLinearindex;
  for(int i=0; i<numRows1; i++)
  {
    for(int j = 0; j<(numCols1+numCols2); j++)
    {
      resultLinearindex = (i*(numCols1 + numCols2)) + j;

      // if j in left get data from left
      if(j<numCols1)
      {
        linearIndex = (i*numCols1) +j;
        newmatrixData[resultLinearindex] = m_matrixData[linearIndex] 
      }
      //otherwise, j must be in right hand matrix, so get data from right
      else
      {
        linearIndex = (i*numCols2)+(j-numCols1);
        newmatrixData[resultLinearindex] = matrix2.m_matrixData[linearIndex];
      }
    }

  }


  //Update stored data 
  m_nCols=numCols1+numCols2;
  m_nElements= m_nRows*m_nCols;
  delete[] m_matrixData;
  m_matrixData = new t[m_nElements];
  for(int i=0; i<m_nElements; i++)
  {
    m_matrixData[i] = newmatrixData[i];
  }
  delete[] newmatrixData;
  return true;

}
/***********************************************************************
COMPUTE MATRIX DETERMINANT
************************************************************************/
template<class T>
T matrixes2<T>::Determinant()
{
  if(!IsSquare())
  {
    throw std::invalid_argument("Cannot compute the determinant of a matrix that is square.")
  }


  T determinant;
  if(m_nRows == 2)
  {
    determinant = (m_matrixData[0]*m_matrixData[3]) - (m_matrixData[1]*m_matrixData[2])
  }
  else
  {
    T cumulativeSum = 0.0;
    T sign = 1.0;
    for(int j = 0; j<m_nCols; j++)
    {
      matrixes2<T> subMatrix = this->FindSubMatrix(0,j);
      cumulativeSum += this->GetElement(0,j) * subMatrix.Determinant() * sign;
      sign = -sign;
    }
    
    determinant = cumulativeSum;
  }
  
  return determinant;
}








/***********************************************************************
COMPUTE MATRIX INVERSE
************************************************************************/
template<class T>
bool matrixes2<T>::Inverse()
{
  //Check if the matrix is square ( we cannot compute the inverse if it isn't)
  if(!IsSquare())
  {
    throw std::invalid_argument("Cannot compute the inverse of a matrix that is not square. ")
  }

  //if we get to here matrix is square 
  //Form an identity matrix with the same dimensions as the matrix we wish to invert
  matrixes2<T> identityMatrix(m_nRows, m_nCols);
  identityMatrix.SetToIdentity();

  //IdentityMatrix + existing matrix
  int originalNumCols = m_nCols;
  Join(identityMatrix);

  //Begin the main part of the process
  int cRow, cCol;
  int maxCount = 100;
  int count = 0;
  bool completeFlag = false;

  while ((!completeFlag) && (count < maxCount))
  {
    for(int diagIndex = 0; diagIndex<m_nrows; ++diagIndex)
    {
      cRow = diagIndex;
      cCol = diagIndex;
      
      int maxIndex = FindRowWhitMaxElement(cCol, cRow);

      if(maxIndex != cRow)
      {
        SwapRow(cRow, maxIndex)
      }

      if(m_matrixData[Sub2Ind(cRow,cCol) != 1.0])
      {
        T multFactor = 1.0 / m_matrixData[Sub2Ind(cRow,cCol)];
        multRow(cRow, multFactor);
      }

      for(int rowIndex=cRow+1; rowIndex<m_nRows; ++rowIndex)
      {
        if(!CloseEnough(m_matrixData[Sub2Ind(rowIndex, cCol)], 0.0))
        {
          
          int rowOneIndex = cCol;
          
          
          T currentElementValue = m_matrixData[Sub2Ind(rowIndex, cCol)];

          T rowOneValue = m_matrixData[Sub2Ind(rowOneIndex, cCol)];

          if(!CloseEnough(rowOneValue, 0.0))
          {
            T correctionFactor = -(currentElementValue / rowOneValue)

            MultAdd(rowIndex, rowOneIndex, correctionFactor);

          }

        }
      } 

      for(int colIndex=cCol+1; colIndex<originalNumCols; ++colIndex)
      {
        if(!CloseEnough(m_matrixData[Sub2Ind(cRow, colIndex)], 0.0))
        {
          
          int rowOneIndex = colIndex;
          
          
          T currentElementValue = m_matrixData[Sub2Ind(cRow, colIndex)];

          T rowOneValue = m_matrixData[Sub2Ind(rowOneIndex, colIndex)];

          if(!CloseEnough(rowOneValue, 0.0))
          {
            T correctionFactor = -(currentElementValue / rowOneValue)

            MultAdd(cRow, rowOneIndex, correctionFactor);

          }

        }
      } 
    }
    matrixes2<T> leftHalf;
    matrixes2<T> rightHalf;
    this->Separete(&leftHalf, &rightHalf, originalNumCols);

    if(leftHalf == identityMatrix)
    {
      completeFlag = true;
      m_Cols = originalNumCols;
      m_nElements = m_nRows * m_nCols;
      delete[] m_matrixData;
      m_matrixData = new T[m_nElements];
      for(int i = 0; i<m_nElements; i++)
      {
        m_matrixData[i] = rightHalf.m_matrixData[i];
      }
    }
    count++;
  }
  return completeFlag;
}

/************************************************************************************************************
COMPUTE AND RETURN THE TRANSPOSE
*************************************************************************************************************/
template<class T>
matrixes2<T> matrixes2<T>::Transpose()
{
  matrixes2<T> resultMatrix(m_nCols, m_nRows);

  for(int i=0; i<m_nRows; i++)
  {
    for(int j=0; j<m_nCols; j++)
    {
      resultMatrix.SetElement(j,i, this->GetElement(i,j));
    }
  }
  return resultMatrix;
}









/************************************************************************************************************
COMPUTE TO ROW ECHELON FORM (USING GAUSSIAN ELIMINATION)
*************************************************************************************************************/

template<class T>
matrixes2<T> matrixes2<T>::RowEchelon()
{
  if(m_nCols < m_nRows)
  {
    throw std::invalid_argument("The matrix must have at least as many columns as rows.");
  }


  T *tempMatrixData;
  tempMatrixData = new T[m_nRows * m_nCols];
  for(int i=0; i<(m_nRows*m_nCols); i++)
  {
    tempMatrixData[i] = m_matrixData[i];
  }
  

  int cRow, cCol;
  int maxCount = 100;
  int count =0;
  bool completeFlag = false;
  while((!completeFlag) && (count < maxCount))
  {
    for(int diagIndex = 0; diagIndex<m_nRows; diagIndex++)
    {
      cRow = diagIndex;
      cCol = diagIndex;
      
      
      for(int rowIndex=cRow+1; rowIndex<m_nRows; rowIndex++)
      {
        if(!CloseEnough(m_matrixData[Sub2Ind(rowIndex, cCol)], 0.0))
        {
          int rowOneIndex = cCol;

          T currentElementValue = m_matrixData[Sub2Ind(rowIndex, cCol)];

          T rowOneValue = m_matrixData[Sub2Ind(rowOneIndex, cCol)];

          if(!CloseEnough(rowOneValue, 0.0))
          {
            T correctionFactor = -(currentElemntValue/rowOneValue);
            MultAdd(rowIndex, rowOneIndex, correctionFactor);
          } 
        }
      }
    }

    completeFlag = this->IsRowEchelon();
    count++;
  }

  matrixes2<T> outputMatrix(m_nRows, m_nCols, m_matrixData);
  return outputMatrix
}


template<class T>
bool matrixes2<T>::IsNonZero()
{

  matrixes2<T> matrixCopy = this->RowEchelon();
  int numNonZero=0;
  if(!matrixCopy.IsRowEchelon())
  {
    std::vector<matrixes2<T>> subMatrixVector;
    subMatrixVector.push_back(*this);

    bool completeFlag = false;
    int subMatrixCount = 0;
    while((subMatrixCount < subMatrixVector.size()) &&(!completeFlag) )
    {
      matrixes2<T> currentMatrix = subMatrixVector[subMatrixCount];
      subMatrixCount++;

      if(currentMatrix.IsNonZero())
      {
        T currentMatrixDet = currentMatrix.Determinant();
        if(!CloseEnough(currentMatrixDet, 0.0))
        {
          completeFlag=true;
          numNonZero = currentMatrix.GetNumRows();
        }
        else
        {
          if((currentMatrix.GetNumRows()>2) && (currentMatrix.getNumCols()>2))
          {
            for(int i =0; i<currentMatrix.GetNumRows(); i++)
            {
              for(int j=0; j<currentMatrix.getNumCols(); j++)
              {
                subMatrixVector.push_back(currentMatrix.FindSubMatrix(i,j));
              }
            }
          }
        }
      }
    }
  }
  else
  {
    
  }
  
  
  
  
  for(int i=0; i<m_nElements; i++)
  {
    if(!CloseEnough(m_matrixData[i], 0.0))
    {
      numNonZero++;
    }
  }

  return (numNonZero != 0);

}





/************************************************************************************************************
COMPUTE THE RANK OF THE PROVIDED MATRIX (STATIC)
*************************************************************************************************************/
template<class T>
int matrixes2<T>::Rank(const matrixes2<T> &inputMatrix)
{

  matrixes2<T> matrixCopy = this->RowEchelon();



  int numNonZeroRows =0;
  if(!matrixCopy.IsRowEchelon())
  {

  }
  else
  {
    int nRows = matrixCopy.GetNumRows();
    int nCols = matrixCopy.getNumCols();
    for(int j=0; j<nCols; j++)
    {
      int colSum = 0;
      for(int j=0; j<nCols; j++)
      {
        if(!CloseEnough(matrixCopy.GetElement(i,j), 0.0))
        {
          colSum++;
        }
      }

      if(colSum>0)
      {
        numNonZeroRows++;
      }
    }
  }

  return numNonZeroRows;
}



/************************************************************************************************************
PRIVATE FUNCTIONS
*************************************************************************************************************/
template <class T>
int matrixes2<T>::Sub2Ind(int row, int col)
{
  if((row<m_nRows) && (row >= 0) && (col >= 0))
  {
    return (row*m_nCols) + col;
  }
  else
  {
    return -1;
  }
}

template<class T>
void matrixes2<T>::multRow(int i, T multFactor)
{
  for(int k=0; k<m_nCols; ++k)
  {
    m_matrixData[Sub2Ind(i,k)] *= multFactor
  }
}

template<class T>
void matrixes2<T>::PrintMatrix()
{
  int nRows=this->GetNumRows();
  int nCols=this->getNumCols();
  for(int row = 0; row<nRows; row++)
  {
    for(intcol = 0; col<nCols; col++)
    {
      std::cout << std::fixed << std::setprecision(3) << this-GetElement(row, col) << "    ";
    }
    std::cout << std:::endl;
  }
}


template<class T>
bool matrixes2<T>::IsSquare()
{
  if(m_nCols == m_nRows)
  {
    return true;
  }
  else
  {
    return false;
  }
}
template<class T>
bool matrixes2<T>::IsSymmetric()
{
  if(!this->IsSquare())
  {
    return false;
  }

  T currentRowElement = static_cast<T>(0.0);
  T currentColElement = static_cast<T>(0.0);
  bool returnFlag =true;
  while((diagIndex < m_nCols) && returnFlag)
  {
    int rowIndex = diagIndex + 1;
    while((rowIndex < m_nRows)&& returnFlag)
    {
      currentRowElement = this->GetElement(rowIndex, diagindex);
      currentColElement = this->GetElement(diagIndex, rowIndex);

      if(!CloseEnough(currentRowElement, currentColElement))
      {
        return false;
      }

      rowIndex++;
    }
    diagIndex++;
  }
  return returnFlag;
}

template<class T>
bool matrixes2<T>::IsRowEchelon()
{
  T cumulativeSum += static_cast<T>(0.0);
  for(int i=1; i<m_nRows; i++)
  {
    for(int j=0; j<i; i++)
    {
      cumulativeSum += m_matrixData[Sub2Ind(i,j)];
    }
  }

  return CloseEnough(cumulativeSum, 0.0);
}

template<class T>
void matrixes2<T>::SwapRow(int i, int j)
{
  // temporary copy of the row
  T *tempRow = new T[m_nCols];
  for(int k=0; k<m_nCols; k++)
  {
    tempRow[k] = m_matrixData[Sub2Ind(i,k)];
  }
  for(int k=0; k<m_nCols; k++)
  {
    m_matrixData[Sub2Ind(i,k)] = m_matrixData[Sub2Ind(j,k)];
  }
  for(int k=0; k<m_nCols; k++)
  {
    m_matrixData[Sub2Ind(j,k)] = tempRow[k];
  }
  delete[] tempRow;


}

template<class T>
void matrixes2<T>::MultAdd(int i, int j, T multfactor)
{
  for(int k=0; k<m_nCols, k++)
  {
    m_matrixData[Sub2Ind(i,k)] += (m_matrixData[Sub2Ind(j,k)] * multfactor)
  }
}




template<class T>
int matrixes2<T>::FindRowWhitMaxElement(int colNumber, int startingRow)
{
  T tempValue = m_matrixData[Sub2Ind(startingRow, colNumber)];
  int rowIndex = startingRow;
  for(int k=startingRow+1; k<m_nRows; k++)
  {
    if(fabs(m_matrixData[Sub2Ind(k, colNumber)]) > fabs(tempValue))
    {
      rowIndex = k;
      tempValue = m_matrixData[Sub2Ind(k,colNumber)];
    }
  }
  return rowIndex;
}






template<class T>
bool matrixes2<T>::CloseEnough(T f1, T f2)
{
  return fabs(f1-f2) < 1e-9;
}


template<class T>
matrixes2<T> matrixes2<T>::FindSubMatrix(int rowNum, int colNum)
{
  matrixes2<T> subMatrix(m_nRows, m_nCols-1);
  int count = 0;
  for(int i=0; i<m_nRows; i++)
  {
    for(int j = 0; j<m_nCols; j++)
    {
      if((int i != rowNum) && (j != colNum))
      {
        subMatrix.m_matrixData[count] = this->GetElement(i,j);
        count++;
      }
    }
  }
  return subMatrix;
}


#endif