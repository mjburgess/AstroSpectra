
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "project.h"

#define debug(format, args...) fprintf (stdout, format, args)

int main(int argc, char **argv) {
	int precision = 3;

	if(argc > 1) {
		sscanf(argv[1], "%d", &precision);
	}

	Project *p = ProjectNew("c18o_single", 109782.1734);
	ProjectCalcStats(p, 1, precision);
	ProjectSave(p);
	ProjectInfo(p);
//	ProjectClose(p);

	p = ProjectNew("dco_single", 144077.321);
	ProjectCalcStats(p, 1, precision);
	ProjectSave(p);
	ProjectInfo(p);
//	ProjectClose(p);

//	Project *p = ProjectNew("c18o_double", 109782.1734);
//	ProjectCalcStats(p, 2, precision);
//	ProjectSave(p);
//	ProjectInfo(p);
//	ProjectClose(p);

	return EXIT_SUCCESS;
}

//PROJECT.c
Project *ProjectNew(char *name, double minFreq) {
	Project *p = malloc(sizeof(Project));

	p->minFreq = minFreq;

	strcpy(p->inFile, name);
	strcpy(p->outFile, name);
	strcpy(p->statsFile, name);
	strcpy(p->name, name);

	strcat(p->inFile, "_input.dat");
	strcat(p->outFile, "_output.dat");
	strcat(p->statsFile, "_stats.dat");

	p->observed = DataNew();
	p->expected = DataNew();

	if(!DataFromFile(p->observed, p->inFile)) {
		debug("Could not read from Input File (%s)\n", p->inFile);
		perror("OS Says:");
		exit(EXIT_FAILURE);
	}

	DataCalcStat(p->observed);

	p->expected->A->data = malloc(p->observed->size * sizeof(double));
	p->expected->V->data  = p->observed->V->data;
	p->expected->size = p->observed->size;

	p->residuals = malloc(p->observed->size * sizeof(double));
	return p;
}

void ProjectCalcStats(Project *p, int peaks, int precision) {
	fitStrategy *f;

	if(peaks > 3) {
		f = &DataGaussianOverlay;
	} else {
		f = &DataGaussianFit;
	}

	p->fit = ChiFitNew(p->observed, f, peaks, 1.0 / (double) precision);

	ChiFitAnalyse(p->fit, p->observed, p->expected);
	f(p->fit->params, p->fit->peaks, p->observed, p->expected);
	DataCalcStat(p->expected);
	p->residualSum = DataCalcResiduals(p->residuals, p->observed, p->expected);

	p->obsFreq = p->minFreq * p->expected->A->max;
}

int ProjectStatsToFile(Project *p) {
	FILE *fp = fopen(p->statsFile, "w");

	if(fp == NULL) {
		return 0;
	}

	fprintf(fp, "<%s, Observed Freq. (%lf MHz)>\n\n", p->name, p->obsFreq);
	fputs(ChiFitFormat(p->fit), fp);
	fputs("\n\n", fp);
	fprintf(fp, "\t<Residual Sum: %lf>\n\n", p->residualSum);
	fputs(DataFormat(p->observed, "Observed"), fp);
	fputs("\n\n", fp);
	fputs(DataFormat(p->expected, "Expected"), fp);
	fputs("\n", fp);
	fclose(fp);

	return 1;
}

void ProjectSave(Project *p) {
	if(!DataToFile(p->observed, p->expected, p->outFile)) {
		debug("Could not write to Output File (%s)", p->outFile);
	}

	if(!ProjectStatsToFile(p)) {
		debug("Could not write to Stats File (%s)", p->statsFile);
	}
}

void ProjectInfo(Project *p) {
#define SEP "+------------------------------------------+"

	puts(SEP);
	printf("  <\t%s\t*\n", p->name);
	puts(SEP);
	puts(ChiFitFormat(p->fit));
	printf("\n\t<Residual Sum: %.2lf>\n\n", p->residualSum);
	puts(DataFormat(p->observed, "Observed"));
	puts(DataFormat(p->expected, "Expected"));
	puts(SEP);
	printf("  *  Observed Freq. (%.2lf MHz)  >  \n", p->obsFreq);
	puts(SEP);
}

void ProjectClose(Project *p) {
	free(p->expected->A->data);
	free(p->expected->V->data);

	free(p->observed->A->data);
	free(p->observed->V->data);

	free(p->observed);
	free(p->expected);
	free(p->fit);
}

//DATA.c
Data *DataNew() {
	Data *d = malloc(sizeof(Data));
	d->A = malloc(sizeof(Var));
	d->V  = malloc(sizeof(Var));

	d->size = d->A->sum = d->A->sumSq = d->A->mean = d->A->stdDev = d->A->rms = d->A->max = 0.0;
	d->V->sum = d->V->mean = d->V->sumSq = d->V->stdDev = d->V->rms = d->V->max = 0.0;
	return d;
}

char *DataFormat(Data *d, char *name) {
#define DT(x) "\t" x "\n"
#define DP(x) "\t\t" x "\t %10.2lf\n"
#define DI(x) "\t\t" x "\t %9d\n"

	char *fmt = malloc(1000 * sizeof(char));

	//snprintf buggy under gcc -std=gnu99
	sprintf(fmt, DT("Data (%s) >")
			DT("A (intensity, K):")
				DP("Max  ") DP("Sum  ") DP("SumSq") DI("Size ") DP("Mean ") DP("StDev") DP("RMSq ") DP("FWHM ")
			DT("V (velocity, m/s):")
				DP("Max  ") DP("Sum  ") DP("SumSq") DI("Size ") DP("Mean ") DP("StDev") DP("RMSq ") DP("FWHM "),

			name,
			d->A->max, d->A->sum, d->A->sumSq, d->size, d->A->mean, d->A->stdDev, d->A->rms, d->A->fwhm,
			d->V->max, d->V->sum, d->V->sumSq, d->size, d->V->mean, d->V->stdDev, d->V->rms, d->V->fwhm);


	return fmt;
}

int DataFromFile(Data *d, char *filename) {
	char buffer[100];
	int maxSize = 250;

	d->V->data  = malloc(maxSize * sizeof(double));
    d->A->data = malloc(maxSize-- * sizeof(double));

	FILE *fp = fopen(filename, "r");

	if(fp == NULL) {
		return 0;
	}

	while(fgets(buffer, 100, fp) != NULL) {
		if(d->size == maxSize) {
			//enlarge arrays if input exceeds allocated size
			d->V->data = realloc(d->V->data, 50 * sizeof(double));
			d->A->data = realloc(d->A->data, 50 * sizeof(double));
			maxSize += 50;
		}

		sscanf(buffer, "%lf%lf", &d->V->data[d->size], &d->A->data[d->size]);
		d->size++;
	}

	fclose(fp);

	return 1;
}

int DataToFile(Data *o, Data * e, char *filename) {
	FILE *fp = fopen(filename, "w");

	if(fp == NULL) {
			perror("OS Said:");
			puts("Correct Error & Press Return");
			getchar();
		}

	fp = fopen(filename, "w");

	if(fp == NULL) {
		perror("OS Said:");
		puts("Failed Again, Exiting...");
		exit(EXIT_FAILURE);
	}

	fprintf(fp, "\n\tV(O)\tI(O)\tI(E)\tI(E) - I(O)\n");
	for(int i = 0; i < o->size; i++) {
		fprintf(fp, "\t%lf\t%lf\t%lf\t%lf\n", o->V->data[i], o->A->data[i], e->A->data[i], o->A->data[i] - e->A->data[i]);
	}

	fclose(fp);

	return 1;
}

void DataCalcStat(Data *d) {
	double size = (double) d->size;
	double fwhm = sqrt(8 * log(2));

	for(int i = 0; i < d->size; i++) {
		d->A->sum += d->A->data[i];
		d->A->sumSq  += d->A->data[i] * d->A->data[i];
		if(d->A->max < d->A->data[i]) {
			d->A->max = d->A->data[i];
		}

		d->V->sum += d->V->data[i];
		d->V->sumSq  += d->V->data[i] * d->V->data[i];
		if(d->V->max < d->V->data[i]) {
			d->V->max = d->V->data[i];
		}
	}


	d->A->mean = d->A->sum / size;
	d->A->rms = sqrt(d->A->sumSq / size);
	d->A->stdDev = (d->A->sumSq / size) - (d->A->mean * d->A->mean);
	d->A->stdDev = sqrt(d->A->stdDev);
	d->A->fwhm = fwhm * d->A->stdDev;

	d->V->mean = d->V->sum / size;
	d->V->rms = sqrt(d->V->sumSq / size);
	d->V->stdDev = (d->V->sumSq / size) - (d->V->mean * d->V->mean);
	d->V->stdDev = sqrt(d->V->stdDev);
	d->V->fwhm = fwhm * d->V->stdDev;
}

void DataGaussianFit(ChiParam *params, int peaks, Data *o, Data *e) {
	double *gPoints = malloc(peaks * sizeof(double));
	double gMax;

	for(int i = 0; i < o->size; i++) {
		for(int j = 0; j < peaks; j++) {
			gPoints[j] = o->A->data[i] * GAUSSIAN(e->V->data[i], params[j].stdDev * params[j].stdDev, params[j].mean);
		}

		gMax = gPoints[0];
		for(int k = 0; k < peaks; k++) {
			if((o->A->data[i] - gPoints[k]) > (o->A->data[i] - gMax)) {
				gMax = gPoints[k];
			}
		}

		e->A->data[i] = gMax;
	}
}

void DataGaussianOverlay(ChiParam *params, int peaks, Data *o, Data *e) {
	for(int i = 0; i < e->size; i++) {
		e->A->data[i] = o->A->data[i] * GAUSSIAN(e->V->data[i], params[0].stdDev * params[0].stdDev, params[0].mean);
	}
}

//unused, calcs resid. when as writing output
double DataCalcResiduals(double *r, Data *o, Data *e) {
	double sum = 0.0;

	for(int i = 0; i < o->size; i++) {
		r[i] = o->A->data[i] - e->A->data[i];
		sum += r[i];
	}

	return sum;
}

//ChiFIT.c
ChiFit *ChiFitNew(Data *o, fitStrategy f, int peaks, double step) {
	ChiFit *cf = malloc(sizeof(ChiFit));
	cf->fitStrat = f;
	cf->peaks = peaks;

	cf->params = malloc(peaks * sizeof(ChiParam));

	cf->step = step;
	for(int i = 0; i < peaks; i++) {
		cf->params[i].mean = cf->params[i].stdDev = 0.0;
		cf->params[i].bestMean    = o->V->mean;
		cf->params[i].bestStdDev  = o->V->stdDev;
		cf->params[i].numMeanTrials  = CHI_NUM_TRIALS(o->V->mean, cf->step);
		cf->params[i].numStdDevTrials = CHI_NUM_TRIALS(o->V->stdDev, cf->step);
	}


	cf->bestChi     = 1.0;

	return cf;
}

char *ChiFitFormat(ChiFit *c) {
	char *fmt = malloc(1000);

#define CT(x) "\n\t" x " > \n"
#define CP(x) "\t\t" x "\t %5.2lf\n"
#define CI(x) "\t\t" x " { Means: %d, StdDevs: %d, Total: %d }\n"

	sprintf(fmt, CT("ChiFit") CI("numTrys") CP("Redu'd X"),
				c->params[0].numMeanTrials, c->params[0].numStdDevTrials, (int) pow(c->params[0].numMeanTrials * c->params[0].numStdDevTrials, (double) c->peaks),
				c->bestChi);

	for(int i = 0; i < c->peaks; i++) {
		sprintf(fmt + strlen(fmt), "\t\tPeak (%d)\n", i);
		sprintf(fmt + strlen(fmt), CP("Std Dev.") CP("Mean    "), c->params[i].bestStdDev, c->params[i].bestMean);
	}
	return fmt;
}

void ChiFitVary(ChiFit *cf, Data *o, Data *e, int dim) {
	double tchi;

	for(int i = 0; i < cf->params[dim].numStdDevTrials; i++, cf->params[dim].stdDev += cf->step) {
		cf->params[dim].mean = 0;
		for(int j = 0; j < cf->params[dim].numMeanTrials; j++, cf->params[dim].mean += cf->step) {
			if(dim == 0) {
				cf->fitStrat(cf->params, cf->peaks, o, e);
				tchi = ChiFitCalcChi(o, e);
				if(tchi < cf->bestChi) {
					cf->bestChi = tchi;
					for(int i = 0; i < cf->peaks; i++) {
						cf->params[i].bestMean   = cf->params[i].mean;
						cf->params[i].bestStdDev = cf->params[i].stdDev;
					}
				}
			} else {
				ChiFitVary(cf, o, e, dim - 1);
			}
		}
	}
}

void ChiFitAnalyse(ChiFit *cf, Data *o, Data *e) {
	ChiFitVary(cf, o, e, cf->peaks - 1);
	for(int i = 0; i < cf->peaks; i++) {
			cf->params[i].mean = cf->params[i].bestMean;
			cf->params[i].stdDev = cf->params[i].bestStdDev;
		}
}

//Reduced Chi
double ChiFitCalcChi(Data *o, Data *e) {
	double chi = 0.0;
	for(int i = 0; i < o->size; i++) {
		chi += pow(o->A->data[i] - e->A->data[i], 2) /  o->A->stdDev * o->A->stdDev;
	}

	return chi / (o->size - 2 - 1); //df = size - params - 1
}


