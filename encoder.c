/* Reads a file of bin boundary values and bin representative values,
   as created by quantize.c from a completely sorted input file of
   32-bit floats. The quantize.c process determines the number of bins
   and the bin boundaries, then with those decisions made, the bin 
   frequencies are then use to transform the original stream of raw floats
   into approximate values indicated by their corresponding bin numbers.

   That stream of bin numbers is then entropy coded.

   Written by Alistair Moffat (The University of Melbourne) as part
   of the paper "Lossy Compression Options for Dense Index Retention"
   at SIGIR-AP 2023.
*
*/

#define COMPRESS_FILE "compressed-approx-index.bin"

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <assert.h>

#include "helpers.c"

int
main(int argc, char *argv[]) {

	FILE *fb=NULL, *fi=NULL, *fo=NULL;

	if ((argc != 4) ||
		(fb=fopen(argv[1], "r")) == NULL ||
		(fi=fopen(argv[2], "r")) == NULL ||
		(fo=fopen(argv[3], "w")) == NULL) {
		fprintf(stderr, "Usage: %s bins-file index-file prox-file\n",
			argv[0]);
		exit(EXIT_FAILURE);
	}

	make_arrays_and_read_bin_data(fb);

        fprintf(stderr, "read descriptions for %lu bins, ", num_bins);
	fprintf(stderr, "covering %zu symbols\n", total);


	/* ok, have the bin data, now for the fun part, second file
	   is a sequence of float values, each must be searched for
	   and mapped to a bin number */

	float f;

	if (fread(head, sizeof(*head), HEADER, fi) != 1) {
    read_error();
  } 
	fwrite(head, sizeof(*head), HEADER, fo);

	size_t cnt=0;

	while (fread(&f, sizeof(f), 1, fi) == 1) {

		// printf("f = %10.7f, ", f);

		/* loop fetches and processes one number */
		int lo, hi, md;

		/* writing binary search, now that's brave */
		lo = 0; hi = num_bins-1;
		while (lo < hi) {
			md = lo + (hi-lo)/2;
			if (f <= U[md]) {
				hi = md;
			} else {
				lo = md+1;
			}
		}

		assert(lo==0 || U[lo-1]<f);
		assert(f <= U[lo]);

		cnt++;

		/* ok, so the bin number we want to code is "lo",
		   let's give it our best shot! */
		arith_encode(lo, c, num_bins, fo);
	}

	encoder_close(fo);
	fclose(fo);

	fprintf(stderr, "wrote %lu codes for floats to %s\n",
		cnt, COMPRESS_FILE);
	fprintf(stderr, "wrote %lu bytes of output ",
		bytes_out);
	fprintf(stderr, "including %d bytes of header\n", HEADER);
	fprintf(stderr, "corresponds to %.4f bits/float, ",
		8.0*bytes_out/cnt);
	fprintf(stderr, "or %.2f%% of raw float size\n",
		100*(8.0*bytes_out)/(32.0*cnt));
		
	return 0;
}
