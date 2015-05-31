#include "RawData.h"

RawData::RawData()
{
    nDim = 0; nData = 0;
    X.clear();
}

RawData::RawData(mat& X)
{
    nDim = X.n_cols; nData = X.n_rows;
    this->X = X;
}

RawData::~RawData()
{
    X.clear();
}

void RawData::save(const char* datafile, const char* mode)
{
    if (!strcmp(mode,"csv"))
    {
        FILE* fout = fopen(datafile,"w");

        for (int i = 0; i < (int)X.n_rows; i++) {
            fprintf(fout,"%f",X(i,0));

            for (int j = 1; j < (int)X.n_cols; j++)
                fprintf(fout,",%f",X(i,j));

            fprintf(fout,"\n");
        }

        fclose(fout);
    }

    else if (!strcmp(mode,"bin"))
    {
        X.save(datafile,arma_binary);
    }

    else if (!strcmp(mode,"ascii"))
    {
        X.save(datafile,arma_ascii);
    }
}

void RawData::load(const char* datafile, const char* mode)
{
    X.clear(); nDim = 0; nData = 0;
    string temp = "csv"; int lapse = 1000;

    vector < vector <double> > buffer; buffer.clear();

    if (!strcmp(mode,temp.c_str()))
    {
        ifstream fin(datafile);
        string line, token;

        cout << "Loading raw data ..." << endl;

        while (getline(fin,line))
        {
        	nData++;
        	if (nData % lapse == 0) cout << nData << " data points have been loaded ..." << endl;

            stringstream parser(line);

            vd temp; temp.clear();

            while (getline(parser,token,','))
                temp.push_back(atof(token.c_str()));

            buffer.push_back(temp);
        }

        X = mat((int) buffer.size(), (int) buffer[0].size());
        for (int i = 0; i < (int) X.n_rows; i++)
        	for (int j = 0; j < (int) X.n_cols; j++) X(i, j) = buffer[i][j];

        for (int i = 0; i < (int) X.n_rows; i++) buffer[i].clear(); buffer.clear();

        cout << "Total number of data points : " << nData << endl;
    }
    else
    {
    	cout << "Loading binary formatted data ..." << endl;
        X.load(datafile,auto_detect);
        cout << "Total number of data points : " << X.n_rows << endl;
    }

    nDim = X.n_cols; nData = X.n_rows;
}
