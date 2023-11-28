/* Constant and variables common to the encoder and decoder.
 * Handled as globals in this file.
 *
 * Written by Alistair Moffat (The University of Melbourne) as part
 * of the paper "Lossy Compression Options for Dense Index Retention"
 * at SIGIR-AP 2023.
*/

#define HEADER 45       // bytes in index file to put straight through; FAISS has 45 byte headers

size_t num_bins=0;	    // the number of bins in this quantized model
float *U;               // the bin upper boundaries
float *S;               // the corresponding representative values
size_t *c;              // and the corresponding bin frequency counts

char head[HEADER+1];


/* constants that control the arithmetic coder */

#define BBYTES 7
        /* less than eight */
#define BBITS (BBYTES*8)
        /* multiple of 8, strictly less than 64 */
#define FULL ((1LL<<BBITS)-1)
#define FULLBYTE 255
#define PART ((1LL<<(BBITS-8)))
#define ZERO (0)
#define MINR ((1LL<<(BBITS-15)))

/* and then state variables for encoding */

uint64_t L=ZERO;
uint64_t R=FULL;
uint64_t low, high, total, scale;
uint8_t last_non_ff_byte=0, byte;
uint32_t num_ff_bytes=0;
int first=1;

size_t bytes_out=HEADER;

/* plus one additional state variable for decoding */

uint64_t D;


void read_error() {
    fprintf(stderr, "Did not read the expected number of bytes. Exiting, sorry!\n");
    exit(EXIT_FAILURE);
}

/* most of the setup and initializations are common to both
   encoder and decoder
*/
void
make_arrays_and_read_bin_data(FILE *fb) {

	int i;

	/* file fb is bin descriptions, has format:
		ncols:		size_t [should be 2]
		num_bins:	size_t
		(ubound, rep):	(float, float) [x numbins]
		bin_frqs	size_t [x numbins]
	*/

	if (fread(&num_bins, sizeof(size_t), 1, fb) != 1) {
    read_error();
  }
	assert(num_bins==2);
	if (fread(&num_bins, sizeof(size_t), 1, fb) != 1) {
    read_error();
  }
	U = malloc(num_bins*sizeof(*U));
	S = malloc(num_bins*sizeof(*S));
	c = malloc(num_bins*sizeof(*c));
	assert(U && S && c);

	for (i=0; i<num_bins; i++) {
    if (fread(U+i, sizeof(float), 1, fb) != 1) {
      read_error();
    }
		if (fread(S+i, sizeof(float), 1, fb) != 1) {
      read_error();
    }
	}
	for (i=0; i<num_bins; i++) {
		if (fread(c+i, sizeof(size_t), 1, fb) != 1) {
      read_error();
    }
	}
	fclose(fb);

	/* last setup step is to convert to cumfreqs, and assign total */
	for (i=1; i<num_bins; i++) {
		c[i] += c[i-1];
	}
	total = c[num_bins-1];
}

/* encode symbol 0<=s<n relative to comfreqs[0..n-1], send any output
   bytes that get generated to fp
*/
void
arith_encode(size_t s, size_t c[], size_t n, FILE *fp) {

	// printf("coding %lu, ", s);

	assert(R>total);

	/* allocated probability range for this symbol */
	if (s==0) {
		low = 0;
	} else {
		low = c[s-1];
	}
	high = c[s];
	// printf("low = %llu, high = %llu, ", low, high);
	
	/* the actual arithmetic coding step */
	scale = R/total;
	L += low*scale;
	if (high<total) {
		/* top symbol gets benefit of rounding gaps */
		R = (high-low)*scale;
	} else {
		R = R - low*scale;
	}

	/* now sort out the carry/renormalization process */
	if (L>FULL) {
		/* lower bound has overflowed, need first to push
		   a carry through the ff bytes and into the pending
		   non-ff byte */
		last_non_ff_byte += 1;
		L &= FULL;
		while (num_ff_bytes>0) {
			fputc(last_non_ff_byte, fp);
			num_ff_bytes--;
			last_non_ff_byte = ZERO;
			bytes_out++;
		}
	}

	/* more normal type of renorm step */
	while (R < PART)  {
		/* can output (or rather, save for later output)
		   a byte from the front of L */
		byte = (L>>(BBITS-8));
		if (byte!=FULLBYTE) {
			/* not ff, so can bring everything up to date */
			if (!first) {
				fputc(last_non_ff_byte, fp);
				bytes_out++;
			}
			while (num_ff_bytes) {
				fputc(FULLBYTE, fp);
				num_ff_bytes--;
				bytes_out++;
			}
			last_non_ff_byte = byte;
			first = 0;
		} else {
			/* ff bytes just get counted */
			num_ff_bytes++;
		}
		L <<= 8;
		L &= FULL;
		R <<= 8;
	}
}

/* finish off the output stream, then switch off the engine
*/
void
encoder_close(FILE *fp) {
	int i;
        if (!first) {
                fputc(last_non_ff_byte, fp);
		bytes_out++;
        }
        while (num_ff_bytes) {
                fputc(FULLBYTE, fp);
		bytes_out++;
                num_ff_bytes--;
        }

        /* then send the final bytes from L, to be sure to be sure */
        for (i=BBYTES-1; i>=0; i--) {
                fputc((L>>((8*i)))&FULLBYTE, fp);
		bytes_out++;
        }
}

/* when starting decoding, first thing required is to wind the handle
   and start the pump
*/
void
decoder_start(FILE *fp) {
	int i;
	D = 0;
        for (i=0; i<BBYTES; i++) {
                D <<= 8;
                D += fgetc(fp);
        }
}

/* decode symbol 0<=s<n relative to comfreqs[0..n-1], return the integer
   symbol number. All bytes are read from fp.
*/
size_t
arith_decode(size_t c[], size_t n, FILE *fp) {

	uint64_t target;
	uint64_t low, high, scale;
	size_t v=0;

	scale = R/total;
	assert(scale>0);
	target = D/scale;

	/* handle the rounding that might accrue at the top of the
	   range, and adjust downward if required */
	if (target>=total) target = total-1;

	// printf("target = %llu, ", target);

#if 0
	/* could use linear search in c[] */
	for (v=0; v<n && target>=c[v]; v++) {
	}
#else
	/* or can use binary search and maybe go a little faster */
	int lo=0, hi=n-1;
	/* elements c[lo..hi] inclusive being considered */
	while (lo<hi) {
		v = lo + (hi-lo)/2;
		if (c[v] <= target) {
			lo = v+1;
		} else {
			hi = v;
		}
	}
	v = lo;
#endif

	assert(v==0 || c[v-1]<=target);
	assert(v<n);
	assert(target<c[v]);

	// printf("decoded %lu\n", v);

	/* adjust, tracing the encoder, with D=V-L throughout */
	if (v==0) {
		low = 0;
	} else {
		low = c[v-1];
	}
	high = c[v];
	D -= low*scale;
	if (high<total) {
		R = (high-low)*scale;
	} else {
		R = R - low*scale;
	}
	assert(D<=R);

	while (R < PART) {
		/* range has shrunk, time to bring in another byte */
		R <<= 8;
		D <<= 8;
		D &= FULL;
		D += fgetc(fp);
	}
	assert(D<=R);

	return v;
}
