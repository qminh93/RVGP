#include "Result.h"

Result::Result()
{
    this->rSize = 0;
}

Result::Result(int rSize)
{
    this->rSize = rSize;
    allTerms = field <mat> (rSize,4);
}

Result::~Result()
{
    allTerms.clear();
}

void Result::save(const char* datafile)
{
    allTerms.save(datafile,arma_binary);
}

void Result::load(const char* datafile)
{
    allTerms.load(datafile,auto_detect);
}

void Result::set(int t, mat M, mat S, Prediction* P)
{
    if (t < rSize)
    {
        allTerms(t,0) = M;
        allTerms(t,1) = S;
        allTerms(t,2) = P->pred;
        mat temp      = ones<mat>(1,2);
        temp(0,0)    *= P->rse;
        temp(0,1)    *= P->mnlp;
        allTerms(t,3) = temp;
    }
}

mat Result::meanQ(int i)
{
    return allTerms(i,0);
}

mat Result::sigmaQ(int i)
{
    return allTerms(i,1);
}

mat Result::pred(int i)
{
    return allTerms(i,2);
}

double Result::rse(int i)
{
    return allTerms(i,3)(0,0);
}

double Result::mnlp(int i)
{
    return allTerms(i,3)(0,1);
}


