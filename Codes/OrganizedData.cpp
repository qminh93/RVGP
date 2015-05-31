#include "OrganizedData.h"
#include "kmean.h"

OrganizedData::OrganizedData()
{
    train.clear();
    test.clear();
    support.clear();
    hyper = new HyperParams();

    nTrain = nTest = nSupport = nDim = nBlock = 0;
}

OrganizedData::~OrganizedData()
{
    train.clear();
    test.clear();
    support.clear();
}

void OrganizedData::save(vs dataset)
{
    train.save   (dataset[0].c_str(),arma_binary);
    test.save    (dataset[1].c_str(),arma_binary);
    support.save (dataset[2].c_str(),arma_binary);
    hyper->save  (dataset[3].c_str());
}

void OrganizedData::load(vs dataset)
{
    train.load   (dataset[0].c_str(),auto_detect);
    test.load    (dataset[1].c_str(),auto_detect);
    support.load (dataset[2].c_str(),auto_detect);
    hyper->load  (dataset[3].c_str());

    nBlock      = train.n_rows;
    nSupport    = getxm().n_rows;

    for (int i = 0; i < nBlock; i++)
    {
        nDim    = getxb(i).n_cols;
        nTrain += getxb(i).n_rows;
        nTest  += getxt(i).n_rows;
    }
}

void OrganizedData::loadHyp(string hypfile, string mode)
{
    if (mode == "csv")
        hyper->loadcsv(hypfile.c_str());
    else
        hyper->load(hypfile.c_str());
}

void OrganizedData::process(RawData* raw, int nBlock, double pTest, int nSupport)
{
    int nData = raw->nData, nDim = raw->nDim - 1;

    this->nSupport = nSupport;
    this->nBlock   = nBlock;
    this->nDim     = nDim;

    train   = field<mat>(nBlock,2);
    test    = field<mat>(nBlock,2);
    support = field<mat>(1,2);

    mat xm(nSupport,nDim),
        ym(nSupport,1);

    vec mark(nData); mark.fill(0);

    printf("Randomly selecting %d supporting point ...\n", nSupport);

    for (int i = 0; i < nSupport; i++)
	{
		int pos = IRAND(0, nData - 1);
		while (mark[pos] > 0)
			pos = IRAND(0, nData - 1);
		mark[pos] = 1;
		for (int j = 0; j < nDim; j++)
			xm(i, j) = raw->X(pos,j);
		ym(i,0) = raw->X(pos,nDim);
	}

	support(0,0) = xm; xm.clear();
	support(0,1) = ym; ym.clear();

    cout << "Partitioning the remaining data into " << nBlock << " cluster using K-Mean ..." << endl;

    vvd _remain;

    for (int i = 0; i < nData; i++) if (!mark(i))
    {
        rowvec R = raw->X.row(i);
        _remain.push_back(r2v(R));
    }

    mat remaining = v2m(_remain);

    mark.clear();

    RawData* remain = new RawData(remaining);

    KMean* partitioner = new KMean(remain);

    Partition* clusters = partitioner->cluster(nBlock);

    cout << "Packaging training/testing data points into their respective cluster" << endl;

    for (int i = 0; i < nBlock; i++)
    {
        cout << "Processing block " << i + 1 << endl;

        int bSize   = (int) clusters->member[i].size(),
            tSize   = (int) floor(bSize * pTest),
            pos     = 0,
            counter = 0;

        mark = vec(bSize); mark.fill(0);

        if (bSize > tSize)  // if we can afford to draw tSize test points from this block without depleting it ...
        {
            mat xt(tSize,nDim),
                yt(tSize,1);

            for (int j = 0; j < tSize; j++)
            {
                pos = IRAND(0, bSize - 1);
				while (mark[pos] > 0)
					pos = IRAND(0, bSize - 1);
				mark[pos] = 1; pos = clusters->member[i][pos];

				for (int t = 0; t < nDim; t++)
					xt(j, t) = remain->X(pos,t);
				yt(j,0) = remain->X(pos,nDim);
            }

            bSize  -= tSize;
            nTest  += tSize;

            test(i,0) = xt; xt.clear();
            test(i,1) = yt; yt.clear();
        }

        nTrain += bSize;

        mat xb(bSize,nDim),
            yb(bSize,1);

        //cout << remain->X.n_rows << endl;

        for (int j = 0; j < (int)mark.n_elem; j++) if (mark[j] < 1)
        {
            for (int t = 0; t < nDim; t++) {
                xb(counter,t) = remain->X(clusters->member[i][j],t);
            }
            yb(counter++,0) = remain->X(clusters->member[i][j],nDim);
        }

        train(i,0) = xb; xb.clear();
        train(i,1) = yb; yb.clear();

        mark.clear();

        printf("Done ! nData[%d] = %d, nTrain[%d] = %d, nTest[%d] = %d .\n", i, (int) clusters->member[i].size(), i, train(i,0).n_rows, i, (int) test(i,0).n_rows);
    }
}

mat OrganizedData::getxb(int i)
{
    return train(i,0);
}

mat OrganizedData::getyb(int i)
{
    return train(i,1);
}

mat OrganizedData::getxt(int i)
{
    return test(i,0);
}

mat OrganizedData::getyt(int i)
{
    return test(i,1);
}

mat OrganizedData::getxm()
{
    return support(0,0);
}

mat OrganizedData::getym()
{
    return support(0,1);
}

mat OrganizedData::kparams()
{
    return hyper->kparams();
}

double OrganizedData::mean()
{
    return hyper->mean();
}

double OrganizedData::noise()
{
    return hyper->noise();
}

double OrganizedData::signal()
{
    return hyper->signal();
}

double OrganizedData::ls(int i)
{
    return hyper->ls(i);
}
