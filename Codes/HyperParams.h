#ifndef HYPERPARAMS_H
#define HYPERPARAMS_H

#include "DataClass.h"

class HyperParams : public DataClass
{
    public:
        HyperParams();
        ~HyperParams();

        void    save(const char* datafile);
        void    load(const char* datafile);
        void    loadcsv(const char* datafile);

        void    setmean(double mean);
        void    setnoise(double noise);
        void    setkparams(vd &k_params);

        double  mean();
        double  noise();
        double  signal();
        double  ls(int i);
        mat     kparams();
    protected:
    private:
        field <mat> params;
};

#endif // HYPERPARAMS_H
