/*
 * data.h
 *
 *  Created on: 25 Apr 2010
 *      Author: Michael
 */

#ifndef DATA_H_
#define DATA_H_

#define CHI_RANGE_MIN(x) 0
#define CHI_RANGE_MAX(x) ((x) + ((x) / 2.0))
#define CHI_NUM_TRIALS(x, step) (int) ((CHI_RANGE_MAX(x) - CHI_RANGE_MIN(x)) / step)

#define PI 3.141592654
#define GAUSSIAN(x, var, mean) (pow((2 * PI * (var)), -0.5) * exp(-(((x) - (mean)) * ((x) - (mean))) / (2 * (var))))

typedef struct {
	double *data;
	double sum, mean, stdDev, rms, max, sumSq, fwhm;
} Var;

typedef struct {
	Var *V, *A;
	int size;
} Data;


typedef struct {
	double stdDev, bestStdDev;
	double mean, bestMean;
	int numStdDevTrials;
	int numMeanTrials;
} ChiParam;

typedef void (fitStrategy(ChiParam *params, int peaks, Data *o, Data *e));
typedef struct {
	fitStrategy *fitStrat;
	ChiParam *params;

	int peaks;

	double step, bestChi;
} ChiFit;

typedef struct {
	Data *observed;
	Data *expected;
	double *residuals;

	ChiFit *fit;

	double chi;
	double residualSum;

	double obsFreq;
	double minFreq;

	char name[50];
	char inFile[50];
	char outFile[50];
	char statsFile[50];
} Project;


Project *ProjectNew(char *name, double minFreq);
void ProjectCalcStats(Project *p, int peaks, int precision);
void ProjectSave(Project *p);
int ProjectStatsToFile(Project *p);
void ProjectInfo(Project *p);

Data *DataNew();
char *DataFormat(Data *d, char *name);
int DataFromFile(Data *d, char *filename);
int DataToFile(Data *o, Data *e,  char *filename);

void DataGaussianOverlay(ChiParam *params, int peaks, Data *o, Data *e);
void DataGaussianFit(ChiParam *params, int peaks, Data *o, Data *e);

void DataCalcStat(Data *d);
double DataCalcResiduals(double *r, Data *o, Data *e); //returns residual sum

ChiFit *ChiFitNew(Data *o, fitStrategy f, int peaks, double step);
void ChiFitAnalyse(ChiFit *cbf, Data *o, Data *e);
char *ChiFitFormat(ChiFit *c);
double ChiFitCalcChi(Data *o, Data *e);

#endif /* DATA_H_ */
