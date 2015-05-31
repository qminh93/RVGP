#ifndef RESULT_H
#define RESULT_H

#include "DataClass.h"


class Result : public DataClass
{
    public:
        int rSize;
        field <mat> allTerms;

        Result();
        Result(int rSize);
        ~Result();

        void    save(const char* datafile);
        void    load(const char* datafile);
        void    set(int i, mat M, mat S, Prediction* P);

        mat     meanQ(int i);
        mat     sigmaQ(int i);
        mat     pred(int i);
        double  rse(int i);
        double  mnlp(int i);
    protected:
    private:
};

#endif // RESULT_H
