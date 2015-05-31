#ifndef KERNEL_H
#define KERNEL_H

#include "HyperParams.h"

class Kernel
{
    public:
        Kernel();
        Kernel(HyperParams* hyper);
        ~Kernel();

        mat kmat(mat &A, mat &B);
        double k(vec &xi, vec &xj);
        double k(rowvec &xi, rowvec &xj);
        double dk(int t, vec &xi, vec &xj);
    protected:
        HyperParams* hyper;
    private:
};

#endif // KERNEL_H
