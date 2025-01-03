#ifndef VECTORS_H
#define VECTORS_H

#include <stdexcept>
#include <iostream>
#include <iomanip>
#include <math.h>
#include <vector>

template<class T>
class Vectors
{
    public:
       Vectors();

       Vectors(std::vector<T> inputData);

       Vectors(int numDims);

       ~Vectors();

       int GetNumDims() const;

       T GetElement(int index) const;
       void SetElement(int index, T vlaue);

       T norm();
       Vectors<T> Normalized();
       void Normalize();

       Vectors<T> operator+ (const Vectors<T> &rhs) const;
       Vectors<T> operator- (const Vectors<T> &rhs) const;
       Vectors<T> operator* (const T &rhs) const;

       template <class U> friend Vectors<U> operator* (const U &lhs, const Vectors<U> &rhs);

       static T dot(const Vectors<T> &a, const Vectors<T> &b);
       static Vectors<T> cross(const Vectors<T> &a, const Vectors<T> &b);

    private:
       std::vector<T> m_vectorData;
       int m_nDims;

};

template<class T>
Vectors<T>::Vectors(int numDims)
{
    m_nDims = numDims;
    m_vectorData = std::vector<T>(numDims, static_cast<T>(0.0));
}


template<class T>
Vectors<T>::Vectors()
{
    m_nDims = 0;
    m_vectorData = std::vector<T>();
}

template<class T>
Vectors<T>::Vectors(std::vector<T> inputData)
{
    m_nDims = inputData.size();
    m_vectorData = inputData;
}

template<class T>
Vectors<T>::~Vectors(){}

/***************************************************************************************************************
FUNCTIONS TO RETURN PARAMETRS
****************************************************************************************************************/

template<class T>
int Vectors<T>::GetNumDims() const
{
    return m_nDims;
}

/***************************************************************************************************************
FUNCTIONS TO HANDEL ELEMENTS OF THE VECTOR
****************************************************************************************************************/


template<class T>
T Vectors<T>::GetElement(int index) const
{
    return m_vectorData[index];
}


template<class T>
void Vectors<T>::SetElement(int index, T value)
{
    m_vectorData[index] = value;
}

/***************************************************************************************************************
FUNCTIONS TO PERFORM COMPUTATIONS ON THE VECTOR
****************************************************************************************************************/

template<class T>
T Vectors<T>::norm()
{
    T cumulativeSum = static_cast<T>(0.0);
    for(int i =0; i<m_nDims; i++)
    {
        cumulativeSum+=(m_vectorData.at(i) * m_vectorData.at(i));

        return sqrt(cumulativeSum);
    }
}

template<class T>
Vectors<T> Vectors<T>::Normalized()
{
    T vecNorm = this->norm();

    Vectors<T> result(m_vectorData);
    return result*(static_cast<T>(1.0) / vecNorm);
}

template<class T>
void Vectors<T>::Normalize()
{
    T vecNorm=this->norm();
    for(int i = 0; i<m_nDims; i++)
    {
        T temp = m_vectorData.at(i) * static_cast<T>(1.0) ? vecNorm;
        m_vectorData.at(i) = temp;
    }
}








/***************************************************************************************************************
OVERLOAD OPERATORS
****************************************************************************************************************/
template<class T>
Vectors<T> Vectors<T>::operator+ (const Vectors<T> &rhs) const 
{
    if(m_nDims != rhs.n_mDims)
    {
        throw std::invalid_argument("Vector dimensions do not match");
    }

    std::vector<T> resultData;
    for(int i=0; i<m_nDims; i++)
    {
        resultData.push_back(m_vectorData.at(i) + rhs.m_vectorData.at(i));
    }

    Vectors<T> result(resultData);
    return result;
}

template<class T>
Vectors<T> Vectors<T>::operator- (const Vectors<T> &rhs) const 
{
    if(m_nDims != rhs.n_mDims)
    {
        throw std::invalid_argument("Vector dimensions do not match");
    }

    std::vector<T> resultData;
    for(int i=0; i<m_nDims; i++)
    {
        resultData.push_back(m_vectorData.at(i) - rhs.m_vectorData.at(i));
    }

    Vectors<T> result(resultData);
    return result;
}


template<class T>
Vectors<T> Vectors<T>::operator* (const T &rhs) const 
{

    std::vector<T> resultData;
    for(int i=0; i<m_nDims; i++)
    {
        resultData.push_back(m_vectorData.at(i) * rhs);
    }

    Vectors<T> result(resultData);
    return result;
}

/***************************************************************************************************************
FRIEND FUNCTIONS
****************************************************************************************************************/

template<class T>
Vectors<T> operator* (const T &lhs, const Vectors<T> &rhs)
{
    std::vector<T> resultData;
    for(int i=0; i<m_nDims; i++)
    {
        resultData.push_back(lhs*rhs.m_vectorData.at(i));
    }

    Vectors<T> result(resultData);
    return result;
}


/***************************************************************************************************************
STATIC FUNCTIONS
****************************************************************************************************************/
template<class T>
T Vectors<T>::dot(const Vectors<T> &a, const Vectors<T> &b)
{
    if(a.m_nDims != b.n_mDims)
    {
        throw std::invalid_argument("Vector dimensions do not match for the dot-product to be computed.");
    }

    T cumulativeSum = static_cast<T>(0.0);
    for(int i=0; i<m_nDims; i++)
    {
        cumulativeSum += a.m_vectorData.at(i) * b.m_vectorData.at(i);
    }

    
    return cumulativeSum;
}


template<class T>
Vectors<T> Vectors<T>::cross(const Vectors<T> &a, const Vectors<T> &b)
{
    if(a.m_nDims != b.m_nDims)
    {
        throw std::invalid_argument("Vector dimensions must match for the cross-product to be complited.");
    }
    
    if(a.m_nDims != 3)
    {
        throw std::invalid_argument("The cross-product can only be computed for three-dimensional vectors.");
    }

    std::vector<T> resultData;
    resultData.push_back(a.m_vectorData.at(1) * b.m_vectorData.at(2) - (a.m_vectorData.at(2) * b.m_vectorData.at(1)));
    resultData.push_back(-((a.m_vectorData.at(0) * b.m_vectorData.at(2)) - (a.m_vectorData.at(2) * b.m_vectorData.at(0))));
    resultData.push_back((a.m_vectorData.at(0) * b.m_vectorData.at(1)) - (a.m_vectorData.at(1) * b.m_vectorData.at(0)));

    Vectors<T> result(resultData);
    return result;

}

#endif