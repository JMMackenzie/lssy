/* Reads an input file of binary floats and quantizes into a given number of
   bins using one of three different strategies. Bin boundaries and bin
   averages are reported.

   Writes a binary output file with the bin upper bounds (the last value to go
   into that bin, <= in the test) and corresponding bin representative values
   as a set of pairs, followed by the bin freq counts as constructed here.

   Commandline arguments, must give all four, no defaults:

   nbins, number of bins to be formed [default 256]
	 bintype, one of 1|2|3 [default 2 = fixed width in range]
	 index.sidx, sorted list of floats with two size_t.s first
	 binsfile.bin, list of computed bins

   Example

   	quantize 256 3 index.sidx index.bins 

   And then use index.bin as a control file for encoder.c to use when
   reducing and representing floats. Also needs to be supplied to
   decoder.c to reconstructed a file of 32-bit binned floats.

  Written by Alistair Moffat (The University of Melbourne) as part
  of the paper "Lossy Compression Options for Dense Index Retention"
  at SIGIR-AP 2023.

*/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <assert.h>


#define BIN1_GEOM 1		// number of items in smallest geometric bin

#define ALL_COLS (-1)

#define EPS 1e-10		// doubles only, don't use this with floats

/* comparison function for sorting floats */
int
cmp(const void *x1, const void *x2) {
	float f1=*(float*)x1, f2=*(float*)x2;;
	if (f1<f2) return -1;
	if (f1>f2) return +1;
	return 0;
}

/* Compute entropy of given frequencies*/
double
entropy(size_t frqs[], size_t n) {
	double ent = 0.0;
	double sum=0.0;
	size_t i;
	for (i=0; i<n; i++) {
		sum += frqs[i];
	}
	for (i=0; i<n; i++) {
		if (frqs[i]) {
			ent += frqs[i] * log(sum/frqs[i]);
		}
	}
	ent = ent/log(2.0)/sum;
	return ent;
}


/* a simple linear quantization, equal numbers of domain values in each bin
 * "Fixed Domain" FD
*/
void
bins_fixed_domain(size_t C[], size_t num_bins, float *F, size_t nF) {
	size_t i, step;
	size_t sofar=0;
	step = nF / num_bins;
	for (i=0; i<(num_bins-1)/2; i++) {
		C[i] = C[num_bins-i-1] = step;
		sofar += 2*step;
	}
	if (num_bins%2 == 0) {
		C[num_bins/2-1] = (nF-sofar)/2;
		C[num_bins/2  ] = (nF-sofar) - C[num_bins/2-1];
	} else {
		C[num_bins/2  ] = (nF-sofar);
	}
	return;
}

/* a simple linear quantization, equal slices of the **range** in each bin
 * "Fixed Range" FR
*/
void
bins_fixed_range(size_t C[], size_t num_bins, float *F, size_t nF) {
	double minF, maxF;
	size_t i, iF;
	double interval;

	/* establish the range of values in F, and the range interval */
	minF = F[0]    - EPS;
	maxF = F[nF-1] + EPS;
	interval = (maxF - minF) / num_bins;

	/* now count how many values in F in each of those sub ranges */
	for (i=0, iF=0; i<num_bins; i++) {
		C[i] = 0;
		while (iF < nF && F[iF] < minF + (i+1)*interval) {
			iF++;
			C[i]++;
		}
	}
	return;
}

/* now a non-linear quantization, growing and then shrinking again as
   a carefully fitted geometric sequence 
   "Geometric Domain" GR   
*/
void
bins_geometric_domain(size_t C[], size_t num_bins, float *F, size_t nF) {

	/* first find the geometric parameter */
	double lo=1.00000001;
	double hi = 1000.00;
	double r, fmid;
	size_t loops=0;

	/* via bisection on the governing equation */
	while (1) {
		if (hi-lo < EPS) break;
		r = (lo+hi)/2;
		fmid = BIN1_GEOM * (pow(r, num_bins/2.0) - 1)/(r-1);
		loops += 1;
		if (fmid < nF/2.0) {
			lo = r;
		} else {
			hi = r;
		}
	}
	fprintf(stderr, "geom ratio   = %10.8f, %lu iterations required\n",
		r, loops);

	/* and now assign bin sizes using that geometric ratio */
	double this=BIN1_GEOM;
	size_t sofar=2*BIN1_GEOM;

	/* assign the two end points */
	C[0] = C[num_bins-1] = BIN1_GEOM;

	/* assign all the in-between points, outside towards the middle */
	for (size_t i=1; i<(num_bins-1)/2; i++) {
		this *= r;
		C[i] = C[num_bins-i-1] = (size_t)(this);
		sofar += 2*C[i];
	}
	if (num_bins%2 == 0) {
		C[num_bins/2-1] = (nF-sofar)/2;
		C[num_bins/2  ] = (nF-sofar) - C[num_bins/2-1];
	} else {
		C[num_bins/2  ] = (nF-sofar);
	}
	return;
}

/* fixed range, but with quarter at bottom and quarter at top of bins
   allocated to singletons, bit of a simple hack, but reduces compression
   cost by almost a bit, and hence allows num_bins todouble, a virtuous
   cycle
   "Central Fixed Range" CFR
*/

void
bins_fixed_skinny(size_t C[], size_t num_bins, float *F, size_t nF) {

	size_t i, singles;

	/* first have 1/4 singleton bins at each of beginning and end */
	singles = num_bins/4;
	for (i=0; i<singles; i++) {
		C[i] = 1;
		C[num_bins-i-1] = 1;
	}
	/* and then fill in the blanks in between the easy way! */
	bins_fixed_range(
		C+singles,
		num_bins - 2*singles,
		F+singles,
		nF - 2*singles
	);
	return;
}

/* and now a tabulation of the methods of interest */

#define NUM_METHODS 4		// index of last method enabled

const char *labels[] = {"",
	"FD",
	"FR",
	"GD",
	"CFR",
	""};

void ((*bin_funcs[])(size_t *, size_t, float *, size_t)) =
	{bins_fixed_domain,
	 bins_fixed_range,
	 bins_geometric_domain,
	 bins_fixed_skinny};

/* print out the bin boundaries and bin averages, text format to stdout
*/
void
print_bins(size_t *C, size_t num_bins, float *F, size_t nF) {
	size_t i=0, strt=0, empty=0;
	double binrep, error, maxerror=0.0, avgerror=0.0;

	/* lets just do a quick bin check, how many are empty? */
	strt = 0;
	for (i=0; i<num_bins-1; i++) {
		if (F[strt]==F[strt+C[i]]) {
			empty += 1;
		}
		strt += C[i];
	}
	if (empty) {
		fprintf(stderr, "empty bins   = %lu\n", empty);
	}
	i = 0;
	strt = 0;
	while (strt<nF) {
		printf("bin %3lu has %7lu vals: ", i, C[i]);
		if (C[i] > 0) {
			printf("%9.6f to %9.6f, ",
				F[strt], F[strt+C[i]-1]);
			/* compute bin representative as average of the
			   values actually in this bin */
			binrep = 0.0;
			for (size_t j=strt; j<strt+C[i]; j++) {
				binrep += F[j];
			}
			binrep /= C[i];
#if 0
			/* or could use bin medians rather bin means */
			if (C[i]%2==0) {
				binrep = (F[strt+(C[i]-1)/2] +
					  F[strt+C[i]/2])/2;
			} else {
				binrep = F[strt+C[i]/2];
			}
#endif
			printf("rep %9.6f, ", binrep);

#if 0
			/* measure average error per bin value */
			error = 0.0;
			for (size_t j=strt; j<strt+C[i]; j++) {
				error += fabs(F[j] - binrep);
			}
			error /= C[i];
			printf("avgerr %9.6f", error);
#else
			/* measure worst error per bin */
			error = binrep - F[strt];
			if (F[strt+C[i]-1] - binrep > error) {
				error = F[strt+C[i]-1] - binrep;
			}
			printf("maxerr %9.6f", error);
#endif
			if (error>maxerror) {
				maxerror = error;
			}
			for (size_t j=strt; j<strt+C[i]; j++) {
				avgerror += fabs(F[j] - binrep);
			}
		}
		printf("\n");
		strt += C[i];
		i++;
	}
	assert(strt==nF);

	fprintf(stderr, "maxerror     = %8.6f\n", maxerror);
	fprintf(stderr, "avgerror     = %8.6f\n", avgerror/nF);
	fprintf(stderr, "entropy      = %.2f bits per bin id\n",
		entropy(C, num_bins));
	fprintf(stderr, "\n");
	return;
}

/* write a file of bin boundaries and representative values, followed by
   the complete set of bin frequencies (to be used by encoder and decoder)
   as binary output to bins.bin
*/
void
write_bins(size_t C[], size_t num_bins, float F[], size_t nF, FILE *fb) {
	size_t i=0, strt=0;
	double binrep;
	float fbinrep;
	size_t value=2;

	assert(fb);

	/* the first two size values */
	fwrite(&value, sizeof(size_t), 1, fb);
	fwrite(&num_bins, sizeof(size_t), 1, fb);

	/* and now the table */
	for (strt=0, i=0; i<num_bins; i++) {
		if (C[i] > 0) {
			binrep = 0.0;
			for (size_t j=strt; j<strt+C[i]; j++) {
				binrep += F[j];
			}
			binrep /= C[i];
		} else {
			binrep = F[strt+C[i]-1];
		}
		fwrite(F+strt+C[i]-1, sizeof(float), 1, fb);
		fbinrep = binrep;
		fwrite(&fbinrep, sizeof(float), 1, fb);
		strt += C[i];
	}

	/* final checks */
	assert(strt==nF);

	/* second output component is the set of bin frequencies decided
	   on in connection with the input data */
	fwrite(C, sizeof(*C), num_bins, fb);
	
}

int
main(int argc, char *argv[]) {

	float *F=NULL;
	size_t nF;
	size_t *C=NULL;
	size_t num_bins;
	size_t bintype;
	size_t ncols;
	size_t nrows;

	FILE *fi, *fb;

	if (argc!=5) {
		fprintf(stderr, "Usage: %s nbins bintype sidx-file bins-file\n",
			argv[0]);
		exit(EXIT_FAILURE);
	}

	/* pick up and check the four parameters */
	num_bins = atoi(argv[1]);
	if (num_bins<4) {
		fprintf(stderr, "minimum nbins is 4\n");
		exit(EXIT_FAILURE);
	}
	bintype = atoi(argv[2]);
	if (bintype<1 || bintype>NUM_METHODS) {
		fprintf(stderr, "invalid binning method:\n");
		for (int k=1; k<=NUM_METHODS; k++) {
			fprintf(stderr, "  -- bintype=%d for %s\n", k, labels[k]);
		}
		exit(EXIT_FAILURE);
	}
	if ((fi = fopen(argv[3], "r")) == NULL) {
		fprintf(stderr, "unable to open %s\n", argv[3]);
		exit(EXIT_FAILURE);
	}
	if ((fb = fopen(argv[4], "w")) == NULL) {
		fprintf(stderr, "unable to open %s\n", argv[4]);
		exit(EXIT_FAILURE);
	}

	fprintf(stderr, "\nquantizing using %s (type %lu binning)\n",
		labels[bintype], bintype);
	fprintf(stderr, "forming %lu bins\n", num_bins);

	/* fetch metadata from input */
	if (fread(&ncols, sizeof(size_t), 1, fi) != 1) {
 		fprintf(stderr, "fread() failure\n");
		exit(EXIT_FAILURE);
  }
	if (fread(&nrows, sizeof(size_t), 1, fi) != 1) {
 		fprintf(stderr, "fread() failure\n");
		exit(EXIT_FAILURE);
  }
	nF = ncols*nrows;

	C = malloc(num_bins * sizeof(size_t));
	assert(C);
	F = malloc(nF*sizeof(float));
	assert(F);

	/* and then fetch the data */
	if (fread(F, sizeof(*F), nF, fi) != nF) {
		fprintf(stderr, "fread() failure\n");
		exit(EXIT_FAILURE);
	}


	float minmag=1e20;
	float maxmag=1e-20;
	size_t num_zero=0;
	size_t num_neg=0;
	size_t num_pos=0;

	for (int i=0; i<nF; i++) {
		if (fabs(F[i]) < minmag) {
			minmag = fabs(F[i]);
		}
		if (fabs(F[i]) > maxmag) {
			maxmag = fabs(F[i]);
		}
		if (F[i] < 0.0) {
			num_neg++;
		} else if (F[i]>0.0) {
			num_pos++;
		} else {
			num_zero++;
		}
	}

	/* put some stats to stderr */
	fprintf(stderr, "\n");
	fprintf(stderr, "data columns = %lu\n", ncols);
	fprintf(stderr, "data rows    = %lu\n", nrows);
	fprintf(stderr, "total vals   = %lu\n", nF);
	fprintf(stderr, "bin count    = %lu\n", num_bins);
	fprintf(stderr, "average bin  = %lu values\n", nF/num_bins);
	fprintf(stderr, "\n");

	fprintf(stderr, "smallest mag = %.7g\n", minmag);
	fprintf(stderr, "biggest mag  = %.7g\n", maxmag);
	fprintf(stderr, "number neg   = %lu\n", num_neg);
	fprintf(stderr, "number zero  = %lu\n", num_zero);
	fprintf(stderr, "number pos   = %lu\n", num_pos);
	fprintf(stderr, "\n");

#if 0
	/* index floats data is now assumed to be sorted upon arrival */
	qsort(F, nF, sizeof(float), cmp);
#endif
	/* but no harm done to check */
	for (int i=0; i<nF-1; i++) {
		assert(F[i] <= F[i+1]);
	}

#if 0
	for (int i=0; i<10; i++) {
		printf("%15.12f\n", F[i]);
	}
	for (int i=nF-10; i<nF; i++) {
		printf("%15.12f\n", F[i]);
	}
#endif

	/* and now get on and do the work via the selected matching
	   function */
	bin_funcs[bintype](C, num_bins, F, nF);
	print_bins(C, num_bins, F, nF);
	write_bins(C, num_bins, F, nF, fb);
	fclose(fi);
	fclose(fb);

	return 0;
}
