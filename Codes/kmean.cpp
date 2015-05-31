/*
 * kmean.cpp
 *
 *  Created on: 1 Dec, 2014
 *      Author: nghiaht
 */

#include "kmean.h"

KMean::KMean(RawData *data)
{
	this->nIter = KMEAN_ITERATION_DEFAULT;
	this->nThread = omp_get_num_procs();
	this->data = data;
	this->nBlock = data->nData;
	this->t1 = this->t2 = time(NULL);
}

KMean::KMean(int nIter, RawData *data)
{
	this->nIter = nIter;
	this->nThread = omp_get_num_procs();
	this->nBlock = nBlock;
	this->data = data;
	this->nBlock = this->data->nData;
	this->t1 = this->t2 = time(NULL);
}

KMean::~KMean()
{
	for (int i = 0; i < (int) C.size(); i++)
		C[i].clear();
	for (int i = 0; i < (int) member.size(); i++)
		member[i].clear();
	C.clear(); member.clear(); nAssign.clear();
}

void KMean::allocate()
{
	for (int t = 0; t < this->nBlock; t++)
		this->member[t].clear(); // free memory

	int chunk = this->data->nData / this->nThread;

	#pragma omp parallel for schedule(dynamic, chunk)
	for (int i = 0; i < this->data->nData; i++)
	{
		this->nAssign[i] = -1;
		double closest = INFTY;

		for (int j = 0; j < this->nBlock; j++)
		{
			double dist = Dist(j, i);
			if (dist < closest)
			{
				closest = dist;
				this->nAssign[i] = j;
			}
		}
	}

	for (int i = 0; i < this->data->nData; i++)
		this->member[this->nAssign[i]].push_back(i);
}

void KMean::reestimate()
{
	int chunk = this->nBlock / this->nThread;
	if (chunk == 0) chunk++;

	#pragma omp parallel for schedule(dynamic, chunk)
	SFOR(t,nBlock)
    {
		if ((int) this->member[t].size() == 0)
		{
			int pos = IRAND(0, this->data->nData - 1);
			C[t].clear(); // free memory
			rowvec R = data->X.row(pos);
			C[t] = r2v(R);
		}
		else
		{
			C[t].clear(); // free memory
			C[t] = vector <double> (data->nDim, 0.0);
			NFOR(i,j,(int)member[t].size(),data->nDim)
                C[t][j] += data->X(member[t][i],j);
			SFOR(i,data->nDim)
                C[t][i] /= (double) member[t].size();
		}
    }
}

void KMean::initialize()
{
	SEED(SEED_DEFAULT);
	vector <int> mark(data->nData, 0);

	nAssign.clear(); // free memory
	nAssign = vector <int> (data->nData, 0);

	SFOR(i,(int)C.size()) C[i].clear(); C.clear(); // free memory
	SFOR(i,(int)member.size()) member[i].clear(); member.clear(); // free memory

	member = vector < vector <int> > (nBlock);
	int pos;

	SFOR(i,nBlock)
	{
		pos = IRAND(0, data->nData - 1);
		while (mark[pos] > 0)
			pos = IRAND(0, data->nData - 1);
		mark[pos] = 1;
		rowvec R = data->X.row(pos);
		C.push_back(r2v(R));
	}

	mark.clear(); // free memory
}

Partition* KMean::cluster(int nBlock)
{
	this->nBlock = nBlock; this->t1 = time(NULL);

	cout << "Initializing Clusters" << endl;
	initialize();

	for (int t = 0; t < this->nIter; t++)
	{
	    cout << "Clustering Iteration " << t + 1 << endl;
		allocate();
		reestimate();
	}

	this->t2 = time(NULL);
	cout << "Done! Clustering Time = " << (double) (t2 -t1) << endl;

	Partition *result = new Partition(this->nBlock, this->C, this->member, this->nAssign);

	return result;
}

double KMean::Dist(int i, int j)
{
	double dist = 0.0;
	SFOR(t,data->nDim)
		dist += pow(C[i][t] - data->X(j,t), 2.0);
	return sqrt(dist);
}


