#include "PrecomputedData.h"
#include "Kernel.h"

PrecomputedData::PrecomputedData()
{
    allTerms.clear();
}

PrecomputedData::~PrecomputedData()
{
    allTerms.clear();
}

void PrecomputedData::save(const char* datafile)
{
    allTerms.save(datafile,arma_binary);
}

void PrecomputedData::save(vector <string> dataset)
{
    save(dataset[4].c_str());
}

void PrecomputedData::load(const char* datafile)
{
    allTerms.load(datafile,auto_detect);
}

void PrecomputedData::load(vector <string> dataset)
{
    load(dataset[4].c_str());
}

mat PrecomputedData::getRbb(int i)
{
    return allTerms(i,0);
}

mat PrecomputedData::getM(int i)
{
    return allTerms(i,1);
}

mat PrecomputedData::getL(int i)
{
    return allTerms(i,2);
}

mat PrecomputedData::getV(int i)
{
    return allTerms(i,3);
}

void PrecomputedData::precompute(OrganizedData* data)
{
    int nBlock = data->nBlock;
    mat xm = data->getxm();
    Kernel* ker = new Kernel(data->hyper);
    mat Kmm_inv = inv(ker->kmat(xm,xm));
    allTerms = field<mat>(nBlock,4);

    SFOR(i,data->nBlock)
    {
        cout << "Precomputing " << i << "-th Block!!" << endl;
        mat xb    = data->getxb(i),
            xt    = data->getxt(i);

        int bTest = xt.n_rows,
            bSize = xb.n_rows;

        mat noise = pow(data->noise(),2.0) * eye<mat>(bSize,bSize),
            mean  = data->mean() * ones<mat>(bSize,1),
            Kbb   = ker->kmat(xb,xb),
            Kbm   = ker->kmat(xb,xm),
            Pb    = Kbm * Kmm_inv,
            Qbb   = Pb * trans(Kbm),
            Jbb   = Kbb - Qbb;

        allTerms(i,0) = inv(Jbb + noise);

        if (bTest)
        {
            mat Ktb = ker->kmat(xt,xb),
                Ktm = ker->kmat(xt,xm),
                Ktt = ker->kmat(xt,xt),
                Pt  = Ktm * Kmm_inv,
                Qtb = Pt * trans(Kbm),
                Jtb = Ktb - Qtb,
                Lt = Jtb * allTerms(i,0);

            allTerms(i,1) = Pt - Jtb * allTerms(i,0) * Pb;
            allTerms(i,2) = Lt * (data->getyb(i) - mean);

            mat temp = Ktt - allTerms(i,1) * trans(Ktm) - Lt * trans(Ktb);
            // V[ft]
            allTerms(i,3) = temp.diag();
        }
    }
}
