#ifndef SGPPLUS_H
#define SGPPLUS_H
#include "libhead.h"
#include "Kernel.h"
#include "OrganizedData.h"
#include "PrecomputedData.h"
#include "Result.h"
#define TAU 9
#define LZ 0.1
#define LPOW 0.75
#define CONVERGE_THRESHOLD 30
#define CONVERGE_EPS 0.001

class SGPPlus
{
    public:
        OrganizedData*      data;
        PrecomputedData*    pre;
        Kernel*             ker;
        SGPSetting*         setting;

        Pair <mat>*         E;
        mat                 Kmm_inv,xm,yt;
        bool                usePrecomp;
        int                 popSize;
        Result*             exactRes;
        Result*             approxRes;
        string              SGP_name;

        SGPPlus();
        SGPPlus(OrganizedData* data, SGPSetting* setting);
        SGPPlus(OrganizedData* data, SGPSetting* setting, PrecomputedData* pre);
        ~SGPPlus();

        void                exact();
        void                approx();

        virtual void        initialize();
        double              rse(mat &X, mat &Y);
        double              mnlp(mat &X, mat &Y, mat &V);
        Pair <mat>*         computeSM(Pair <mat>* L);


        void                update(int t);
        Pair <mat>*         gradient(vi &RS);
        virtual Pair <mat>* subgradient(int index);
        virtual Pair <mat>* subgradient(mat &sample);
        virtual Prediction* predict(mat &qfm, mat &sigma);

        void                saveExact();
        void                saveApprox();
    protected:
    private:

};

#endif // SGPPLUS_H
