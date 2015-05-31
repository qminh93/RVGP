#include "Kernel.h"

Kernel::Kernel()
{
    hyper = new HyperParams();
}

Kernel::Kernel(HyperParams* hyper)
{
    this->hyper = hyper;
}

Kernel::~Kernel()
{
    delete hyper;
}

mat Kernel::kmat(mat &A, mat &B)
{
    mat result(A.n_rows,B.n_rows);
    NFOR(i,j,A.n_rows,B.n_rows)
    {
        rowvec Ai   = A.row(i),
               Bj   = B.row(j);
        result(i,j) = k(Ai,Bj);
    }
    return result;
}

double Kernel::k(rowvec &xi, rowvec &xj)
{
    double  result = 0.0;

	SFOR(i,(int)xi.n_elem)
	{
		result += pow(((xi[i] - xj[i]) / hyper->ls(i)) , 2.0);
	}
	result = pow(hyper->signal(), 2.0) * exp(-0.5 * result);

	return result;
}

double Kernel::k(vec &xi, vec &xj)
{
    double  result = 0.0;

	SFOR(i,(int)xi.n_elem)
		result += pow(((xi[i] - xj[i]) / hyper->ls(i)) , 2.0);
	result = pow(hyper->signal(), 2.0) * exp(-0.5 * result);

	return result;
}

double Kernel::dk(int t,vec &xi, vec &xj)
{
    double  tmp     = 0.0,
            result  = 0.0;

	SFOR(i,(int)xi.n_elem)
		tmp += pow(((xi[i] - xj[i]) / hyper->ls(i)),2.0);
    tmp = exp(-0.5 * tmp);

	if (t > 0)
	{
		result = pow((xi[t-1] - xj[t-1]) * hyper->signal() , 2.0) * tmp / pow(hyper->ls(t-1), 3.0);
	}
	else
    {
        result = 2.0 * hyper->signal() * tmp;
    }

	return result;
}
