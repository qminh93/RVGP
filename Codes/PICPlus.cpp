#include "PICPlus.h"

PICPlus::PICPlus(OrganizedData* data,SGPSetting* setting):SGPPlus()
{
    this->SGP_name   = "PIC";
    this->data       = data;
    this->setting    = setting;
    this->usePrecomp = true;
    initialize();
    this->pre->precompute(data);
}

PICPlus::PICPlus(OrganizedData* data,SGPSetting* setting, PrecomputedData* precomp):SGPPlus()
{
    this->SGP_name   = "PIC";
    this->data       = data;
    this->setting    = setting;
    this->usePrecomp = true;
    initialize();
    this->pre        = precomp;
}

PICPlus::~PICPlus()
{
    //dtor
}

Prediction* PICPlus::predict(mat &qfm, mat &sigma)
{
    mat pred(data->nTest,1), predvar(data->nTest,1);
    int cur = 0;

    SFOR(i,data->nBlock) if (data->getxt(i).n_rows)
    {
        int bTest = data->getxt(i).n_rows;
        mat mean  = data->mean() * ones <mat> (bTest,1),
            xt    = data->getxt(i),
            predi = pre->getM(i) * qfm + pre->getL(i),
            // MNLP
            temp  = pre->getM(i) * sigma * trans(pre->getM(i)),
            vari  = pre->getV(i) + temp.diag();

        SFOR(j,(int)predi.n_rows)
        {
            pred(cur,0) = predi(j,0);

            // MNLP
            predvar(cur++,0) = vari(j,0);
        }
    }

    double _rse         = rse(pred,yt);
    double _mnlp        = mnlp(pred,yt,predvar);
    Prediction* result  = new Prediction(pred,_rse,_mnlp);

    pred.clear();
    cout << "rse: " << _rse << endl;
    cout << "mnlp: " << _mnlp << endl;
    return result;
}

