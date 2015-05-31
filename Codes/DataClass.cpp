#include "DataClass.h"

DataClass::DataClass()
{
    //ctor
}

DataClass::~DataClass()
{
    //dtor
}

void DataClass::save(const char* datafile)
{
    save(datafile,"bin");
}

void DataClass::save(const char* datafile, const char* mode)
{

}

void DataClass::save(vector <string> dataset)
{

}

void DataClass::load(const char* datafile)
{
    load(datafile, "bin");
}

void DataClass::load(const char* datafile,const char* mode)
{

}

void DataClass::load(vector <string> dataset)
{

}
