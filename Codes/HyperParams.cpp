#include "HyperParams.h"

HyperParams::HyperParams()
{
    params = field <mat> (1,3);
}
HyperParams::~HyperParams()
{
    params.clear();
}

void HyperParams::save(const char* datafile)
{
    params.save(datafile,arma_binary);
}

void HyperParams::loadcsv(const char* datafile)
{
    ifstream fin(datafile);
    string   token;
    vd       _kparams;

    getline(fin,token,',');
    setmean(log(atof(token.c_str())));
    getline(fin,token,',');
    setnoise(atof(token.c_str()));

    while (getline(fin,token,','))
        _kparams.push_back(atof(token.c_str()));

    setkparams(_kparams);

    fin.close();
}

void HyperParams::load(const char* datafile)
{
    params.load(datafile,auto_detect);
}

void HyperParams::setmean(double _mean)
{
    params(0,0) = mat(1,1).fill(_mean);
}

void HyperParams::setnoise(double _noise)
{
    params(0,1) = mat(1,1).fill(_noise);
}

void HyperParams::setkparams(vd &_kparams)
{
    params(0,2) = mat(1,(int)_kparams.size());
    for (int i = 0; i < (int)_kparams.size(); i++)
        params(0,2)(0,i) = _kparams[i];
}

double HyperParams::mean()
{
    return exp(params(0,0)(0,0));
}

double HyperParams::noise()
{
    return exp(params(0,1)(0,0));
}

double HyperParams::signal()
{
    return exp(kparams()(0,0));
}

double HyperParams::ls(int i)
{
    return exp(kparams()(0,i+1));
}

mat HyperParams::kparams()
{
    return params(0,2);
}
