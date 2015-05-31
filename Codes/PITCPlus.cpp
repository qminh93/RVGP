#include "PITCPlus.h"

PITCPlus::PITCPlus(OrganizedData* data, SGPSetting* setting):SGPPlus(data,setting)
{
    this->SGP_name = "PITC";
    //ctor
}

PITCPlus::~PITCPlus()
{
    //dtor
}
