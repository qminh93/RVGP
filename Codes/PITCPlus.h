#ifndef PITCPLUS_H
#define PITCPLUS_H

#include "SGPPlus.h"


class PITCPlus : public SGPPlus
{
    public:
        PITCPlus(OrganizedData* data,SGPSetting* setting);
        ~PITCPlus();
    protected:
    private:
};

#endif // PITCPLUS_H
