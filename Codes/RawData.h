#ifndef RAWDATA_H
#define RAWDATA_H

#include "DataClass.h"

class RawData : public DataClass
{
    public:
        mat X;
        int nDim, nData;

        RawData();
        RawData(mat& X);
        ~RawData();

        void save(const char* datafile,const char* mode);
        void load(const char* datafile,const char* mode);
    protected:
    private:
};

#endif // RAWDATA_H
