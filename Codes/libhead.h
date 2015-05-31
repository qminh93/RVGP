#ifndef LIBHEAD_H
#define LIBHEAD_H

#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <omp.h>
#include <armadillo>
#include <algorithm>
#include <fstream>
#include <vector>
#include <queue>
#include <deque>
#include <set>
#include <map>
#include <cstring>
#include <ctime>
#include <cmath>
#include <iomanip>
#include <stack>
#include <bitset>
#include <functional>
#include <numeric>

using namespace std;
using namespace arma;

#define SEED_DEFAULT 20051987
#define SEED(x) srand(x)
#define IRAND(a, b) (rand() % ((b) - (a) + 1) + (a)) // randomly generate an integer from [a, b]
#define DRAND(a, b) (((b) - (a)) * ((double) rand() / (RAND_MAX)) + (a)) // randomly generate a real number from [a, b]
#define EPS 0.00000000000000001
#define PI 3.141592653589793238462
#define EN 2.71828182845904509080
#define INFTY 1000000000000000.0
#define N(m, mean, var) ((1 / sqrt(2 * PI * (var))) * exp(-0.5 * (pow((m) - (mean), 2.0) / (var)))) // the pdf function for N(mean, var = sigma^2)
#define ENT(var) (0.5 * log(2.0 * PI * E * (var))) // the entropy of N(., var)
#define LN(m, mean, var) (-0.5 * log(2 * PI * (var)) - 0.5 * (pow((m) - (mean), 2.0)) / (var) ) // the log pdf function for N(mean, var = sigma^2)
#define SFOR(i,n) for (int i = 0; i < (int)n; i++)
#define NFOR(i,j,n,m) SFOR(i,n) SFOR(j,m)
#define vi  vector < int >
#define vd  vector < double >
#define vm  vector < mat >
#define vs  vector < string >
#define vvi vector < vi >
#define vvd vector < vd >
#define vvm vector < vm >


struct Partition
{
	int nBlock; // the number of partitions
	vvd C; // store the estimated centroids
	vvi member; // lists of data indices belonging to each cluster
	vi  nAssign; // nAssign[i] -- the cluster which the ith data point belongs to

	Partition(int nBlock, vvd &C, vvi &member, vi &nAssign)
	{
		this->nBlock  = nBlock;
		this->C       = C;
		this->member  = member;
		this->nAssign = nAssign;
	}

	~Partition()
	{
		for (int i = 0; i < (int) C.size(); i++)
			C[i].clear();
		for (int i = 0; i < (int) member.size(); i++)
			member[i].clear();
		C.clear(); member.clear(); nAssign.clear();
	}
};

template <class T> struct Pair
{
    T first;
    T second;

    Pair() {}
    Pair(T first, T second) { this->first = first; this->second = second; }
    Pair(Pair <T>* another) { first = another->first; second = another->second; }
    ~Pair() { first.clear(); second.clear(); }
};

struct SGPSetting
{
    bool  blockSampling, // false: point sampling;
          measureTime;
    int   sSize,         // number of blocks/points per sample
          nPred,         // number of prediction to be made
          interval,      // interval between predictions
          seed;          // random seed
    vs    logs;          // 0: timeLog, 1: exactLog, 2: approxLog

    double alpha,beta,gamma;


    SGPSetting()
    {
        // Do nothing
    }

    SGPSetting(bool blockSampling, bool measureTime, int sSize, int nPred, int interval, vs &logs, double alpha, double beta, double gamma,int seed)
    {
        this->blockSampling  = blockSampling;
        this->measureTime    = measureTime;
        this->sSize          = sSize;
        this->nPred          = nPred;
        this->interval       = interval;
        this->logs           = logs;
        this->alpha          = alpha;
        this->beta           = beta;
        this->gamma          = gamma;
        this->seed           = seed;
    }

    ~SGPSetting()
    {
        logs.clear();
    }
};

struct Prediction
{
    mat pred;
    double rse;
    double mnlp;

    Prediction()
    {
        this->rse  = 0.0;
        this->mnlp = 0.0;
    }

    Prediction(mat &pred, double rse, double mnlp)
    {
        this->pred = pred;
        this->rse  = rse;
        this->mnlp = mnlp;
    }

    ~Prediction()
    {
        pred.clear();
    }
};

string num2str(double x);
string num2str(int x);
double lapse (clock_t tStart, clock_t tEnd);
template <class ForwardIterator, class T> void iota(ForwardIterator first, ForwardIterator last, T value);
vi  randsample(int popSize, int sampleSize);
vd  r2v (rowvec &R);
vd  c2v (colvec &C);
mat v2m (vvd &A);
mat v2mat (vvd &A);
void csv2bin_blkdata(string src,string dest);
void csv2bin_support(string src,string dest);
#endif /* LIBHEAD_H_ */
