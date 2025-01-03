#ifndef LSQ_H
#define LSQ_H


#include<stdexcept>
#include<iostream>
#include<iomanip>
#include<math.h>
#include<vector>

#include"vectors.h"
#include"matrixes2.h"

constexpr int LSQ_NOINVERSE = -1;

template<typename T>
int LSQ(matrixes2<T> &Xin, const Vectors<T> &yin, Vectors<T> &result)
{
    matrixes2<T> X =Xin;
    Vectors<T> y = yin;

    matrixes2<T> XT = X.Transpose();
    matrixes2<T> XXT = XT * X;

    if(!XXT.Inverse())
    {
        return LSQ_NOINVERSE;
    }
    
    matrixes2<T> XXTXT = XXT * XT;

    result = XXTXT * y;

    return 1;
    
}



#endif