#ifndef DATACLASS_H
#define DATACLASS_H

#include "libhead.h"

class DataClass
{
    public:
        DataClass();
        virtual ~DataClass();

        virtual void save(const char* datafile);
        virtual void save(const char* datafile,const char* mode);
        virtual void save(vector <string> dataset);

        virtual void load(const char* datafile);
        virtual void load(const char* datafile,const char* mode);
        virtual void load(vector <string> dataset);
    protected:
    private:
};

#endif // DATACLASS_H
