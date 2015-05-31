#ifndef PRECOMPUTEDDATA_H
#define PRECOMPUTEDDATA_H

#include "DataClass.h"
#include "OrganizedData.h"

class PrecomputedData : public DataClass
{
    public:
        PrecomputedData();
        ~PrecomputedData();

        mat getRbb(int i);
        mat getM(int i);
        mat getL(int i);
        mat getV(int i);

        void precompute(OrganizedData* data);

        void save(const char* datafile);
        void save(vector <string> dataset);
        void load(const char* datafile);
        void load(vector <string> dataset);
    protected:
        field <mat> allTerms;
    private:
};

#endif // PRECOMPUTEDDATA_H
