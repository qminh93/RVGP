#ifndef PICPLUS_H
#define PICPLUS_H

#include "SGPPlus.h"


class PICPlus : public SGPPlus
{
    public:
        PICPlus(OrganizedData* data, SGPSetting* setting);
        PICPlus(OrganizedData* data, SGPSetting* setting, PrecomputedData* pre);
        ~PICPlus();

        Prediction* predict(mat &qfm, mat &sigma);
    protected:
    private:
};

#endif // PICPLUS_H
