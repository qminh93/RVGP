#include "libhead.h"

string num2str(double x)
{
	return static_cast <ostringstream*> (&(ostringstream() << x))->str();
}

string num2str(int x)
{
	return static_cast <ostringstream*> (&(ostringstream() << x))->str();
}

template<class ForwardIterator, class T>
void iota(ForwardIterator first, ForwardIterator last, T value)
{
    while(first != last)
    {
        *first++ = value;
        ++value;
    }
}

vi randsample(int popSize, int sampleSize)
{
    vi temp(popSize);

    iota(temp.begin(),temp.end(),0);
    random_shuffle(temp.begin(),temp.end());

    vi result(temp.begin(),temp.begin() + sampleSize);

    return result;
}

double lapse (clock_t tStart, clock_t tEnd)
{
    return (1000.0 * (tEnd - tStart) / CLOCKS_PER_SEC);
}

vd r2v (rowvec &R)
{
    vd r;
    SFOR(i, R.n_elem) r.push_back(R(i));

    return r;
}

vd c2v (colvec &C)
{
    vd c;
    SFOR(i, C.n_elem) c.push_back(C(i));

    return c;
}

mat v2m (vvd &A)
{
    if (A.size() == 0) return mat(0,0);
    mat M((int)A.size(),(int)A[0].size());
    for (int i = 0; i < (int)A.size(); i++)
        for (int j = 0; j < (int)A[0].size(); j++)
            M(i,j) = A[i][j];
    return M;
}

void csv2bin_support(string src,string dest)
{
    ifstream fin(src.c_str());
    string line,token;
    vvd Sx,Sy;
    while (getline(fin,line))
    {
        stringstream ss(line);
        vd point,pointy;
        while (getline(ss,token,','))
    		point.push_back(atof(token.c_str()));

        pointy.push_back(point.back()); point.pop_back();

    	Sx.push_back(point); Sy.push_back(pointy);
    }

    field <mat> M (1,2);
    M(0,0) = v2m(Sx);
    M(0,1) = v2m(Sy);
    M.save(dest.c_str(),arma_binary);
    fin.close();
}

void csv2bin_blkdata(string src,string dest)
{
	ifstream fin(src.c_str());

	string line,token;
	vector <vvd> Sx,Sy;
	int curblk = 1;
	vvd blkx,blky;
	vd pointx,pointy;

	while (getline(fin,line))
	{
		stringstream ss(line);
		pointx.clear();
		pointy.clear();

		getline(ss,token,',');
		int blkno = atoi(token.c_str());
		while (curblk < blkno)
		{
			Sx.push_back(blkx); blkx.clear();
			Sy.push_back(blky); blky.clear();
			curblk++;
		}

		while (getline(ss,token,','))
			pointx.push_back(atof(token.c_str()));

        pointy.clear(); pointy.push_back(pointx.back()); pointx.pop_back();
		blkx.push_back(pointx);
		blky.push_back(pointy);
	}

	Sx.push_back(blkx);
	Sy.push_back(blky);

	field <mat> M((int)Sx.size(),2);

	for (int i = 0; i < (int)Sx.size(); i++)
	{
		M(i,0) = v2m(Sx[i]);
		M(i,1) = v2m(Sy[i]);
	}

	M.save(dest.c_str(),arma_binary);
    Sx.clear(); Sy.clear(); blkx.clear(); blky.clear(); pointx.clear(); pointy.clear();
	fin.close();
}
