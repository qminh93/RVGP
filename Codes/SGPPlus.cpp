#include "SGPPlus.h"

SGPPlus::SGPPlus()
{

}

SGPPlus::SGPPlus(OrganizedData* data,SGPSetting* setting)
{
    this->data       = data;
    this->setting    = setting;
    this->usePrecomp = false;

    cout << "Running SGP with Hyper Parameters: " << endl;

    cout << "Mean : " << data->mean() << endl;

    cout << "Noise: " << data->noise() << endl;

    cout << "Signal: " << data->signal() << endl;

    for (int i = 0; i < data->nDim; i++)
        cout << i << "-th length scale: " << data->ls(i) << endl;

    initialize();
}

SGPPlus::SGPPlus(OrganizedData* data, SGPSetting* setting, PrecomputedData* pre)
{
    SGPPlus(data,setting);

    this->pre        = pre;
    this->usePrecomp = true;
}

SGPPlus::~SGPPlus()
{
    //Destructor
}

void SGPPlus::initialize()
{
    xm         = data->getxm();
    ker        = new Kernel(data->hyper);
    Kmm_inv    = inv(ker->kmat(xm,xm));
    pre        = new PrecomputedData();
    E          = new Pair <mat>();
    exactRes   = new Result(0);
    approxRes  = new Result(setting->nPred + 1);
    popSize    = (setting->blockSampling ? data->nBlock : data->nTrain);

    yt         = mat(data->nTest,1);

    int blk = 0, cur = 0;
    while (blk < data->nBlock)
    {
        for (int i = 0; i < (int)data->getyt(blk).n_rows; i++)
        {
            yt(cur++,0) = data->getyt(blk)(i,0) - data->mean();
        }
        blk++;
    }
}

// N x 1 mat
double SGPPlus::rse(mat &X,mat &Y)
{
    mat diff = X - Y;
    return sqrt(dot(diff,diff) / data->nTest);
}

double SGPPlus::mnlp(mat &X, mat &Y, mat &V)
{
    double MNLP = 0;

    SFOR(i,X.n_rows)
        MNLP += pow(Y(i,0) - X(i,0),2.0) / abs(V(i,0)) + log(2 * PI * abs(V(i,0)));

    MNLP = MNLP / (2.0 * data->nTest);

    return MNLP;
}

Pair <mat>* SGPPlus::computeSM(Pair <mat>* L)
{
    mat S = -0.5 * inv(L->first),
        M = S * L->second;

    return new Pair <mat> (S,M);
}

void SGPPlus::exact()
{
    cout << SGP_name << " exact:" << endl;

    clock_t tStart  = clock();
    vi RS           = randsample(popSize,popSize);
    Pair <mat>* L   = gradient(RS);
    Pair <mat>* SM  = computeSM(L);

    exactRes->set(0,SM->first,SM->second,predict(SM->second,SM->first));
    clock_t tEnd    = clock();

    exactRes->save(setting->logs[2].c_str());

    if (setting->measureTime)
    {
        FILE* tlog = fopen(setting->logs[0].c_str(),"w");
        fprintf(tlog,"%.3f\n",lapse(tStart,tEnd));
        fclose(tlog);
    }
}

void SGPPlus::approx()
{
    cout << SGP_name <<  " approximate:" << endl;
    vd  ite_rt,cum_rt;
    double ave_pred = 0;

    int maxIte      = setting->nPred * setting->interval,
        unchange    = 0,
        t           = 0,
		nIter = 0;

    E->first        = setting->gamma * eye <mat> (data->nSupport,data->nSupport);
    E->second       = setting->alpha * randu <mat> (data->nSupport,1) + (setting->beta - 0.5 * setting->alpha) * ones <mat> (data->nSupport,1);

    Pair<mat>* SM   = computeSM(E);
    approxRes->set(nIter++,SM->first,SM->second,predict(SM->second,SM->first));

    while (unchange < CONVERGE_THRESHOLD && t++ < maxIte)
    {
        cout << "Iteration " << t << endl;

        mat prev    = E->second;

        clock_t tStart = clock();
        update(t);
        clock_t tEnd   = clock();
        ite_rt.push_back(lapse(tStart,tEnd));
        if (cum_rt.size() == 0)
            cum_rt.push_back(ite_rt.back());
        else
            cum_rt.push_back(cum_rt.back() + ite_rt.back());

        if (t % setting->interval == 0)
        {
            SM = computeSM(E);
            clock_t pStart = clock();
            approxRes->set(nIter++,SM->first,SM->second,predict(SM->second,SM->first));
            clock_t pEnd   = clock();
            ave_pred += lapse(pStart,pEnd);
        }

        unchange = (rse(E->second,prev) * data->nTest < CONVERGE_EPS ? unchange + 1 : 0);
    }

    approxRes->save(setting->logs[3].c_str());

    cout << "Average prediction time = " << (ave_pred / setting->nPred) << endl;
    cout << "Result summary:" << endl;
    for (int i = 0; i < nIter; i++)
    {
        cout << "At iteration " << i * setting->interval << ", mrse = " << approxRes->rse(i) << ", mnlp = " << approxRes->mnlp(i) << endl;
    }
    if (setting->measureTime)
    {
        FILE* tlog = fopen(setting->logs[1].c_str(),"w");
        for (int i = 0; i < (int)cum_rt.size(); i++)
            fprintf(tlog,"%.3f %.3f\n",ite_rt[i],cum_rt[i]);
        fclose(tlog);
    }
}

void SGPPlus::update(int t)
{
    vi RS            = randsample(popSize,setting->sSize);
    Pair <mat>* L    = gradient(RS);

    double learnRate = LZ / pow(1 + LZ * TAU * t, LPOW);

    E->first        += learnRate * (- E->first - 0.5 * L->first - TAU * E->first);
    E->second       += learnRate * (- E->second + L->second - TAU * E->second);

    L->first.clear();
    L->second.clear();
    delete(L);
}

Pair <mat>* SGPPlus::gradient(vi &RS)
{
    mat L0 = zeros <mat> (data->nSupport,data->nSupport),
        L1 = zeros <mat> (data->nSupport,1);

    double scale = (double) popSize / (double) setting->sSize;

    if (setting->blockSampling)
    {
        for (int i = 0; i < (int)RS.size(); i++)
        {
            cout << i << "-th Block/Point: " << RS[i] << endl;
            Pair <mat>* Li   = subgradient(RS[i]);
            L0              += Li->first;
            L1              += Li->second;

            Li->first.clear();
            Li->second.clear();
            delete(Li);
        }

        L0 = Kmm_inv + scale * Kmm_inv * L0 * Kmm_inv;
        L1 = scale * Kmm_inv * L1;
    }

    else
    {
        mat sample = mat((int)RS.size(),data->nDim + 1);
        for (int i = 0; i < (int)RS.size(); i++)
        {
            int blk = 0, j = RS[i];
            while (j >= (int)data->getxb(blk).n_rows)
                j  -= (int)data->getxb(blk++).n_rows;

            for (int k = 0; k < data->nDim; k++)
                sample(i,k) = data->getxb(blk)(j,k);
            sample(i,data->nDim) = data->getyb(blk)(j,0);
        }

        Pair <mat>* Li  = subgradient(sample);
        L0 = Kmm_inv + scale * Kmm_inv * Li->first * Kmm_inv;
        L1 = scale * Kmm_inv * Li->second;

        Li->first.clear();
        Li->second.clear();
        delete(Li);
    }

    Pair <mat>* result = new Pair <mat>(L0,L1);
    L0.clear();
    L1.clear();

    return result;
}

Pair <mat>* SGPPlus::subgradient(int i)
{
    int bSize   = data->getxb(i).n_rows;
    mat mean    = data->mean() * ones <mat> (bSize,1),
        xbi     = data->getxb(i),
        ybi     = data->getyb(i) - mean;

    mat ki      = ker->kmat(xm,xbi),
        kit     = ki.t(),
        kii     = ker->kmat(xbi,xbi),
        qii     = kit * Kmm_inv * ki,
        noise   = pow(data->noise(),2.0) * eye <mat> (bSize,bSize),
        l       = inv(kii - qii + noise);

    mat Li0     = ki * l * kit,
        Li1     = ki * l * ybi;

    Pair <mat>* result = new Pair<mat> (Li0,Li1);

    ki.clear();
    kii.clear();
    qii.clear();
    l.clear();
    noise.clear();
    mean.clear();
    xbi.clear();
    ybi.clear();
    Li0.clear();
    Li1.clear();

    return result;
}

Pair <mat>* SGPPlus::subgradient(mat &sample)
{
    int bSize   = sample.n_rows;
    mat mean    = data->mean() * ones <mat> (bSize,1),
        xbi     = sample.submat(0,0,sample.n_rows - 1,sample.n_cols - 2),
        ybi     = sample.submat(0,sample.n_cols - 1,sample.n_rows - 1,sample.n_cols - 1) - data->mean();

    mat ki      = ker->kmat(xm,xbi),
        kit     = ki.t(),
        kii     = ker->kmat(xbi,xbi),
        qii     = kit * Kmm_inv * ki,
        noise   = pow(data->noise(),2.0) * eye <mat> (bSize,bSize),
        l       = inv(kii - qii + noise);

    mat Li0     = ki * l * kit,
        Li1     = ki * l * ybi;

    Pair <mat>* result = new Pair <mat> (Li0,Li1);

    // Prevent memory leaks
    ki.clear();
    kii.clear();
    qii.clear();
    noise.clear();
    l.clear();
    mean.clear();
    xbi.clear();
    ybi.clear();
    Li0.clear();
    Li1.clear();

    return result;
}

Prediction* SGPPlus::predict(mat &qfm,mat &sigma)
{
    mat pred(data->nTest,1), predvar(data->nTest,1);
    int cur = 0;

    SFOR(i,data->nBlock) if (data->getxt(i).n_rows)
    {
        int bTest = data->getxt(i).n_rows;
        mat mean  = data->mean() * ones <mat> (bTest,1),
            xt    = mat(data->getxt(i)),
            Ktm   = ker->kmat(xt,xm),
            Ktt   = ker->kmat(xt,xt),
            Ptm   = Ktm * Kmm_inv,
            predi = Ptm * qfm,
            temp  = Ktt - Ptm * (eye<mat>(data->nSupport,data->nSupport) + sigma * Kmm_inv) * trans(Ktm),
            vari  = temp.diag();

        SFOR(j,(int)predi.n_rows)
        {
            pred(cur,0) = predi(j,0);
            predvar(cur++,0) = vari(j,0);
        }
    }

    double _rse         = rse(pred,yt),
           _mnlp        = mnlp(pred,yt,predvar);
    Prediction* result  = new Prediction(pred,_rse,_mnlp);

    pred.clear();
    cout << "rse: " << _rse << endl;
    cout << "mnlp: " << _mnlp << endl;

    return result;
}





