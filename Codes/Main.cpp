#include    "libhead.h"
#include    "PITCPlus.h"
#include    "DTCPlus.h"
#include    "PICPlus.h"
#include    "FITCPlus.h"


// Command Flags
string  WILDCARD = "@wild",
        DATAFLAG = "@data",
        PITCFLAG = "@pitc",
        PICFLAG  = "@pic",
        DTCFLAG  = "@dtc",
        FITCFLAG = "@fitc",
        EXITFLAG = "@exit",
        CONVFLAG = "@conv";

// Global variables to make life easier
ifstream config;
OrganizedData* data;
PrecomputedData* pre;
bool usePrecompute;

// Declare functions
bool                same(char* s1, char* s2);
string              readln();
vs                  get_dataset();
SGPSetting*         get_setting(bool &exact, bool &approx);
PrecomputedData*    get_precomp();
OrganizedData*      get_data();
RawData*            get_raw();
OrganizedData*      process_raw(RawData* raw,int nBlock, double pTest,int nSupport);


// Main functions
void wild_phase()
{
    // Do what ever you want here!!
    // Call @wild follow by @exit in configuration file
    string sup_file = "./Dataset/emulate_support.csv",
           sup_bin = "./Dataset/emulate_support.bin";
    csv2bin_support(sup_file,sup_bin);
}

void conv_phase()
{
    HyperParams*   hyper  = new HyperParams();
    vs             csvdata,bindata;

    csvdata.push_back(readln()); // train.csv file
    csvdata.push_back(readln()); // test.csv file
    csvdata.push_back(readln()); // support.csv file
    csvdata.push_back(readln()); // hyper.csv file

    bindata.push_back(readln()); // train.bin file
    bindata.push_back(readln()); // test.bin file
    bindata.push_back(readln()); // support.bin file
    bindata.push_back(readln()); // hyper.bin file

    csv2bin_blkdata(csvdata[0],bindata[0]);
    csv2bin_blkdata(csvdata[1],bindata[1]);
    csv2bin_support(csvdata[2],bindata[2]);
    hyper->loadcsv(csvdata[3].c_str());
    hyper->save(bindata[3].c_str());
}

void data_phase()
{
    //SEED(SEED_DEFAULT);
    string cmd = readln();

    // Use raw or processed data
    if (cmd == "[raw = on]")
    {
        RawData* rd       = get_raw();
        int      nBlock   = atoi(readln().c_str()),
                 nSupport = atoi(readln().c_str());
        double   pTest    = atof(readln().c_str());
                 data     = process_raw(rd,nBlock,pTest,nSupport);
    }
    else if (cmd == "[raw = off]")
    {
        data = get_data();
    }
    else
    {
        cout << "INVALID SYNTAX" << endl;
        exit(EXIT_FAILURE);
    }

    // Use precomputed data or not
    cmd = readln();
    if (cmd == "[precomp = on]")
    {
        usePrecompute = true;
        pre = get_precomp();
    }
    else if (cmd == "[precomp = off]")
    {
        usePrecompute = false;
    }
    else if (cmd == "[precomp = save]")
    {
        usePrecompute = true;
        pre = new PrecomputedData();
        pre->precompute(data);
        pre->save(readln().c_str());
    }
    else
    {
        cout << "INVALID SYNTAX" << endl;
        exit(EXIT_FAILURE);
    }
}

void pitc_phase()
{
	SEED(SEED_DEFAULT);

    bool exact, approx;

    SGPSetting* setting = get_setting(exact,approx);

    if (setting->seed >= 0) SEED(setting->seed);

    SGPPlus* pitc = new PITCPlus(data,setting);

    if (approx) pitc->approx();
    if (exact)  pitc->exact();
}

void pic_phase()
{
	SEED(SEED_DEFAULT);

    bool exact, approx;

    SGPSetting* setting = get_setting(exact,approx);

    if (setting->seed >= 0) SEED(setting->seed);

    SGPPlus* pic;

    if (usePrecompute)
        pic = new PICPlus(data,setting,pre);
    else
        pic = new PICPlus(data,setting);

    if (approx) pic->approx();
    if (exact)  pic->exact();
}

void dtc_phase()
{
	SEED(SEED_DEFAULT);

    bool exact, approx;

    SGPSetting* setting = get_setting(exact,approx);

    if (setting->seed >= 0) SEED(setting->seed);

    SGPPlus* dtc = new DTCPlus(data,setting);

    if (approx) dtc->approx();
    if (exact)  dtc->exact();
}

void fitc_phase()
{
    SEED(SEED_DEFAULT);

    bool exact, approx;

    SGPSetting* setting = get_setting(exact,approx);

    if (setting->seed >= 0) SEED(setting->seed);

    SGPPlus* fitc = new FITCPlus(data,setting);

    if (approx) fitc->approx();
    if (exact)  fitc->exact();
}
int main(int argc, char* argv[]) {
    string line;

    if (argc < 2)
    {
        cout << "NOT ENOUGH ARGUMENTS" << endl;
        config.open("C:/Coding/RVGP/Config/test2.txt");
        conv_phase();
        return 0;
    }

    config.open(argv[1]);

    while (getline(config,line))
    {
        if (!strcmp(line.c_str(),WILDCARD.c_str()))
            wild_phase();
        if (!strcmp(line.c_str(),CONVFLAG.c_str()))
            conv_phase();
        if (!strcmp(line.c_str(),DATAFLAG.c_str()))
            data_phase();
        if (!strcmp(line.c_str(),PITCFLAG.c_str()))
            pitc_phase();
        if (!strcmp(line.c_str(),PICFLAG.c_str()))
            pic_phase();
        if (!strcmp(line.c_str(),DTCFLAG.c_str()))
            dtc_phase();
        if (!strcmp(line.c_str(),FITCFLAG.c_str()))
            fitc_phase();
        if (!strcmp(line.c_str(),EXITFLAG.c_str()))
        {
            config.close();
            return 0;
        }
    }

    config.close();
}

// Supporting functions

bool same(char* s1, char* s2)
{
    return (!strcmp(s1,s2));
}

string readln()
{
    string str;
    getline(config,str);
    return str;
}

vs get_dataset()
{
    vs          dataset;
    SFOR(i,4)   dataset.push_back(readln());
    return      dataset;
}

RawData* get_raw()
{
    string      datafile = readln();
    RawData*    raw = new RawData();
                raw->load(datafile.c_str(),"csv");
    return      raw;
}

OrganizedData* process_raw(RawData* raw,int nBlock, double pTest, int nSupport)
{
    string          hypfile = readln(),
                    mode    = "csv";
    vs              dataset = get_dataset();
    OrganizedData*  od = new OrganizedData();

    od->loadHyp(hypfile,mode);
    od->process(raw,nBlock,pTest,nSupport);
    od->save(dataset);

    return od;
}

OrganizedData* get_data()
{
    vs              dataset = get_dataset();
    OrganizedData*  od = new OrganizedData();
                    od->load(dataset);
    return          od;
}

PrecomputedData* get_precomp()
{
    PrecomputedData*    pd = new PrecomputedData();
                        pd->load(readln().c_str());
    return              pd;
}

SGPSetting* get_setting(bool &exact, bool &approx)
{
    bool        blockSampling, // false: point sampling;
                measureTime;
    int         sSize,         // number of blocks/points per sample
                nPred,         // number of prediction to be made
                interval,      // interval between predictions
                seed;          // random seed
    vs          logs;          // 0: timeLog, 1: exactLog, 2: approxLog
    double      alpha,         // Tuning parameters
                beta,
                gamma;


    exact = approx = true;

    string cmd = readln();
    if (cmd == "[mode = exact]")  approx = false;
    if (cmd == "[mode = approx]") exact  = false;

    // Block sampling or Point sampling
    // blockSampling = ((readln() == "block") ? true : false);

    // Record running time
    measureTime   = ((readln() == "[time = yes]") ? true : false);

    // sample size
    sSize         =  atoi(readln().c_str());

    // number of predictions to be made
    nPred         =  atoi(readln().c_str());

    // Interval between predictions
    interval      =  atoi(readln().c_str());

    // Random seed
    seed          =  atoi(readln().c_str());

    // Tuning parameters
    alpha         =  atof(readln().c_str());
    beta          =  atof(readln().c_str());
    gamma         =  atof(readln().c_str());

    // exact time log
    logs.push_back(readln());

    // approx time log
    logs.push_back(readln());

    // exact log
    logs.push_back(readln());

    // approx log
    logs.push_back(readln());

    return new SGPSetting(true,measureTime,sSize,nPred,interval,logs,alpha,beta,gamma,seed);
}
