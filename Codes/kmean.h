/*
 * kmean.h
 *
 *  Created on: 1 Dec, 2014
 *      Author: nghiaht
 */

#ifndef KMEAN_H_
#define KMEAN_H_

#include "libhead.h"
#include "RawData.h"

#define KMEAN_ITERATION_DEFAULT 300

class KMean
{
	private :

	vector < vector <double> > C; // store the estimated centroids
	vector < vector <int> > member; // lists of data indices belonging to each cluster
	vector <int> nAssign; // nAssign[i] -- the cluster which the ith data point belongs to

	int nIter, nThread, nBlock;
	RawData *data;

	time_t t1, t2;

	void initialize();
	void allocate();
	void reestimate();
	double Dist(int i, int j);

	public:

	KMean(RawData *data);
	KMean(int nIter, RawData *data);
	~KMean();

	Partition* cluster(int nBlock);
};

#endif /* KMEAN_H_ */
