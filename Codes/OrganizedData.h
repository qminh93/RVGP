#ifndef ORGANIZEDDATA_H
#define ORGANIZEDDATA_H

#include "DataClass.h"
#include "HyperParams.h"
#include "RawData.h"

class OrganizedData : public DataClass
{
    public:
        int nTrain, nTest, nSupport, nDim, nBlock;
        HyperParams* hyper;
        field <mat> train, test, support;

        OrganizedData();
        ~OrganizedData();

        void save(vs dataset);
        void load(vs dataset);
        void loadHyp(string hypfile, string mode);
        void process(RawData* raw, int nBlock, double pTest, int nSupport);

        mat getxb(int i);
        mat getyb(int i);
        mat getxt(int i);
        mat getyt(int i);
        mat getxm();
        mat getym();

        mat kparams();
        double noise();
        double mean();
        double signal();
        double ls(int i);

    protected:
    private:
};

#endif // ORGANIZEDDATA_H
