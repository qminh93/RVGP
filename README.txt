CONTENTS OF THIS FILE
---------------------
 * Introduction
 * Requirements
 * Usage
 * Data format
 * Configuration template
 * FAQ
 * Maintainers

INTRODUCTION
---------------------
RVGP is a C++ package which implements various anytime sparse Gaussian process regression (SGPR) models based on 
low-rank covariance matrix approximations such as PIC, DTC, and FITC. The resulting anytime SGPRs includes PIC+, 
DTC+, and FITC+. For a full description of these models, please refer to our paper: 

A Unifying Framework of Anytime Sparse Gaussian Process Models with Stochastic Variational Inference for Big Data. 
Trong Nghia Hoang, Quang Minh Hoang and Kian Hsiang Low. 
In Proceedings of the 32nd International Conference on Machine Learning (ICML '15)

REQUIREMENTS
---------------------
1. System:
   - Compatible with Windows/OS X/Linux
     (Tested on Linux server)
   - Minimum amount of RAM depends on dataset's size. 
     (For our experiment with 2 million data points, at least 16GB of RAM is recommended)

2. Supporting softwares:
   - gcc 4.8.1 or above (with OpenMP enabled - refer to: https://gcc.gnu.org/onlinedocs/libgomp/Enabling-OpenMP.html)
   - Armadillo 4.500.0 or above (for installation instruction, refer to: http://arma.sourceforge.net/download.html)

USAGE (FOR LINUX SYSTEM ONLY)
---------------------
1. Compilation:
   g++ *.cpp -o <executable-filename> -O3 -fopenmp -larmadillo -I <armadillo-header-folder-directory>

2. Expected input:
   We accept three forms of input (see format specification in the INPUT FORMAT section below)
	1) Raw data in CSV format (unpartitioned, unclustered): <dataset>, <hyper-parameters>
	2) Organized data in CSV format: <training data>, <testing data>, <support set>, <hyper-parameters>
	3) Organized data in Armadillo Binary (BIN) format (same as the above, but in binary format).

3. Expected output:
   We produce three types of output:
	1) Real-time statistics (RMSE and MNLP at different iterations) printed on-screen.
	2) Total incurred time (for constructing the predictive model only) versus number of iterations returned in a .log text file.
	3) Predictive results and variances at different iterations returned in a .bin Armadillo binary file.

4. Running:
   The compiled program takes in one single argument which is the configuration file specifying the above expected input and output 
   (see syntax below in CONFIGURATION TEMPLATE section). The deploy command is as followed:
   <executable-filename> <configuration-file>

INPUT FORMAT
---------------------
1.  Raw data format (CSV):
    CSV file, each line contains a tuple of (x, y). The components of x and y are numerical values separated by comma (standard CSV format).

2a. Hyper-parameter format (CSV):
    CSV file, 1 line contains the following:
    <mean>, <log-noise>, <log-signal>, <log-lengthscale-1>, <log-lengthscale-2>, ... <log-lengthscale-d>

2b. Hyper-parameter format (BIN):
    Converted from 2a. using RVGP in-built functions.

3a. Organized data format (CSV):
    3a.1) Training data:
	  CSV file, each line contains a tuple of (block-no, x, y), where block-no is the index of the partition 
	  whom training datapoint (x,y) belongs to. Components of this tuple are numerical values separated by comma.

    3a.2) Testing data:
	  CSV file, each line contains a tuple of (block-no, x, y), where block-no is the index of the partition 
	  whom testing datapoint (x,y) belongs to. Components of this tuple are numerical values separated by comma.

    3a.3) Support set:
	  CSV file, each line contains a tuple of support point (x, y). The components of x and y are numerical values 
	  separated by comma (standard CSV format).
3b. Organized data format (BIN):
    Converted from 3a. using RVGP in-build functions.

CONFIGURATION TEMPLATE
---------------------
1.  Overall template:
    One text file, separated into phases, marked by phase flags. The general format is specified below:

    @PHASE_FLAG_1
    phase settings
    @PHASE_FLAG_2
    phase settings 
    ...
    @exit

2.  Convert phase:
    This phase takes in 4 CSV organized data files and convert them into Armadillo Binary format (.BIN) 
    Flag: @conv
    Phase settings:
	<training-data-csv>
	<testing-data-csv>
	<support-set-csv>
	<hyper-param-csv>
	<training-data-bin-destination>
	<testing-data-bin-destination>
	<support-set-bin-destination>
	<hyper-param-bin-destination>

3.  Data phase:
    This phase specifies the input method of the program:
    Flag: @data
    Phase settings:
	[raw = on/off]
	{raw-data-info-if-raw-on}
	<training-data-bin>
	<testing-data-bin>
	<support-set-bin>
	<hyper-param-bin>
  	[precomp = on/off/save]
	{precomp-data-if-on-or-save}

    3.1 Raw data info:
	- If [raw = on], the following info need to be provided
	  <raw-data-csv>
	  <number-of-blocks>
	  <number-of-supports>
	  <percentage-of-data-to-be-used-as-test>
	  <hyper-param-csv>
	  RVGP will automatically cluster the raw data, partition it into 
	  blocks and save the results in binary files provided below.
	- If [raw = off], no further info need to be provided.

    3.2 Precompute data info:
	- If [precomp = off], RVGP will proceed to carry out regression 
	  without precomputing certain terms.
	- If [precomp = on], the following info need to be provided
	  <precomputed-data-file>
	  RVGP will then load the precomputed terms from the provided file.
	- If [precomp = save], the following info need to be provided
	  <precomputed-data-file>
	  RVGP will then proceed to precompute certain terms and save the 
	  result into designated file (in Armadillo binary format) for future uses.
	
4.  SGP phase:
    This phase specifies the SGP model used for regression.
    Flag: @pic, @fitc, @dtc
    Phase settings:
    	[mode = exact/approx/both]
	[time = yes/no]
	{run-settings}
	{initialization-parameters}
	{output-file-destinations}

    4.1 Mode:
	- Exact:  Update model using all training data blocks.
	- Approx: Update model using a sampled number of training data blocks at every iteration.
	- Both:   Carry out regression for each mode once.

    4.2 Time:
	- Yes:    Measure cumulative time spent constructing predictive model at every predicting iteration
	- No:	  Do not measure time
    
    4.3 Run settings:
	<sample-size> : Number of blocks sampled per iteration
	<number-of-predictions> : Number of anytime predictions made
	<interval-between-predictions> : Interval (number of interations) between predictions.
	
	e.g: {sample size = 1, number of predictions = 200, interval between predictions = 5},
	This means RVGP will run for 5 * 200 = 1000 iteration. At each iteration RVGP samples 1 
	training block and update the predictive model. At every fifth iteration, RVGP performs 
	1 prediction.

    4.4 Initialization parameters:
	Used to configure the initial predictive model (pre-update).
	<alpha>
	<beta>
	<gamma>
	This is strictly for tuning purpose. As empirically inspected, we suggest setting
	{alpha = 5, beta = 0, gamma = 0.00001}.

    4.5 Output file destinations:
	<exact-time-log>  : .time text file recording cumulative time spent constructing predictive model at every predicting iteration for exact SGP
	<approx-time-log> : .time text file recording cumulative time spent constructing predictive model at every predicting iteration for approximate SGP
	<exact-result>	  : .log binary file recording predictive results and variances at every predicting iteration for exact SGP
	<approx-result>   : .log binary file recording predictive results and variances at every predicting iteration for approximate SGP

5.  Wild phase
    This phase is strictly for developer only. Edit the code in wildphase function in Main.cpp to your preference!
    Flag: @wild
    Phase settings: none!

6.  Exit tag:
    Marking end of config file
    Flag: @exit

7.  Example:
    A.  Input:  Raw data from CSV file to be partitioned into 1000 blocks, 512 support points to be 
		generated, test ratio = 5%, precomputation to be saved to file, PIC Exact with timing, 200 
		predictions at 5 iterations interval, sample 1 block per iteration
		Corresponding config file: See example1.txt

    B.  Input:  Organized data from CSV files to be converted to bin, recomputation file supplied, 
		FITC Approx with timing, 200 predictions at 5 iterations interval, sample 1 block per iteration
		Corresponding config file: See example2.txt

    C.  Input:  Organized data from BIN files, no precomputation required, DTC exact and approximate 
		without timing, PIC approx, 200 predictions at 5 iterations interval, sample 1 block per iteration
		Corresponding config file: See example3.txt
	
FAQ
---------------------

MAINTAINERS
---------------------
Quang Minh, Hoang (qminh93@gmail.com)
Trong Nghia, Hoang (trongnghia_hoang@yahoo.com)