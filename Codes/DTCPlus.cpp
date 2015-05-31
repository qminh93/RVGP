#include "DTCPlus.h"

DTCPlus::DTCPlus(OrganizedData* data, SGPSetting* setting):SGPPlus(data,setting)
{
    this->SGP_name = "DTC";
}

DTCPlus::~DTCPlus()
{
    //dtor
}

Pair <mat>* DTCPlus::subgradient(int i)
{
    mat xi, yi, ki, kit, Li0, Li1;

    xi  = mat(data->getxb(i));
    yi  = mat(data->getyb(i) - data->mean() * ones <mat> (xi.n_rows,1));

    ki  = ker->kmat(xm,xi);
    kit = ki.t();
    Li0 = pow(data->noise(),-2.0) * ki * kit;
    Li1 = pow(data->noise(),-2.0) * ki * yi;

    Pair <mat>* result = new Pair <mat>(Li0,Li1);

    // Prevent memory leaks
    Li0.clear();
    Li1.clear();
    xi.clear();
    yi.clear();
    ki.clear();

    return result;
}

Pair <mat>* DTCPlus::subgradient(mat &sample)
{
    mat xi, yi, ki, kit, Li0, Li1;

    xi  = sample.submat(0,0,sample.n_rows - 1,sample.n_cols - 2);
    yi  = sample.submat(0,sample.n_cols - 1,sample.n_rows - 1,sample.n_cols - 1) - data->mean();
    ki  = ker->kmat(xm,xi);
    kit = ki.t();
    Li0 = pow(data->noise(),-2.0) * ki * kit;
    Li1 = pow(data->noise(),-2.0) * ki * yi;

    Pair <mat>* result = new Pair <mat>(Li0,Li1);

    // Prevent memory leaks
    Li0.clear();
    Li1.clear();
    xi.clear();
    yi.clear();
    ki.clear();
    return result;
}
