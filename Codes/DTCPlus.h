#ifndef DTCPLUS_H
#define DTCPLUS_H

#include "SGPPlus.h"


class DTCPlus : public SGPPlus
{
    public:
        DTCPlus(OrganizedData* data,SGPSetting* setting);
        ~DTCPlus();

        Pair <mat>* subgradient(mat &sample);
        Pair <mat>* subgradient(int i);
    protected:
    private:
};

#endif // DTCPLUS_H
