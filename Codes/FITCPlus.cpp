#include "FITCPlus.h"

FITCPlus::FITCPlus(OrganizedData* data, SGPSetting* setting):SGPPlus(data,setting)
{
    this->SGP_name = "FITC";
}

FITCPlus::~FITCPlus()
{
    // do nothing
}

Pair <mat>* FITCPlus::subgradient(int i)
{
    mat xi, yi, ki, Li0, Li1;

    xi  = mat(data->getxb(i));
    yi  = mat(data->getyb(i) - data->mean() * ones <mat> (xi.n_rows,1));

    ki  = ker->kmat(xm,xi);
    //Li0 = pow(data->noise(),-2.0) * ki * trans(ki);
    //Li1 = pow(data->noise(),-2.0) * ki * yi;

    Li0 = zeros<mat>(data->nSupport, data->nSupport);
    Li1 = zeros<mat>(data->nSupport, 1);

    double total1 = 0.0;
    double total2 = 0.0;

    mat kmj, qjj, kmj_t;
    for (int j = 0; j < (int)xi.n_rows; j++)
    {
         clock_t tB1 = clock();
         rowvec xj   = xi.row(j);
         kmj         = ki.submat(0,j,data->nSupport - 1,j);
         kmj_t       = kmj.t();
         qjj         = kmj_t * Kmm_inv * kmj;
         double sj   = 1.0 / (ker->k(xj,xj) - qjj(0,0) + pow(data->noise(),2.0));
         clock_t tB2 = clock();
         Li0         = Li0 + sj * kmj * kmj_t;
         Li1         = Li1 + sj * kmj * yi(j,0);
         clock_t tEnd = clock();
         total1 += lapse(tB1,tEnd);
         total2 += lapse(tB2,tEnd);
         //cout << "j = " << j << " : " << lapse(tB1,tEnd) << " " << lapse(tB2,tEnd) << endl;
    }

    //cout << "Total : " << total1 << " " << total2 << endl;

    Pair <mat>* result = new Pair <mat>(Li0, Li1);

    // Prevent memory leaks
    Li0.clear();
    Li1.clear();
    xi.clear();
    yi.clear();
    ki.clear();

    return result;
}


Pair<mat>* FITCPlus::subgradient(mat &sample)
{
    mat xi, yi, ki, Li0, Li1;

    xi  = sample.submat(0,0,sample.n_rows - 1,sample.n_cols - 2);
    yi  = sample.submat(0,sample.n_cols - 1,sample.n_rows - 1,sample.n_cols - 1) - data->mean();
    ki  = ker->kmat(xm,xi);

    Li0 = zeros<mat>(data->nSupport, data->nSupport);
    Li1 = zeros<mat>(data->nSupport, 1);

    for (int j = 0; j < (int)xi.n_rows; j++)
    {
        rowvec xj   = xi.row(j);
        mat kmj     = (mat)ki.col(j),
        qjj         = trans(kmj) * Kmm_inv * kmj;
        double sj   = 1.0 / (ker->k(xj,xj) - qjj(0,0) + pow(data->noise(),2.0));
        Li0         = Li0 + sj * kmj * trans(kmj);
        Li1         = Li1 + sj * kmj * yi(j,0);
    }

    Pair <mat>* result = new Pair <mat>(Li0,Li1);

    // Prevent memory leaks
    Li0.clear();
    Li1.clear();
    xi.clear();
    yi.clear();
    ki.clear();
    return result;
}
