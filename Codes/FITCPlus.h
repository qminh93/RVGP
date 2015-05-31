#ifndef FITCPLUS_H
#define FITCPLUS_H

#include "SGPPlus.h"


class FITCPlus : public SGPPlus
{
    public:
        FITCPlus(OrganizedData* data, SGPSetting* setting);
        ~FITCPlus();

        Pair <mat>* subgradient(mat &sample);
        Pair <mat>* subgradient(int i);

    protected:
    private:
};

#endif // FITCPLUS_H
