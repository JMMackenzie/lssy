# lssy cmprssn

This is the code repository for the paper **Lossy Compression Options for Dense Index Retention** by
Joel Mackenzie and Alistair Moffat (SIGIR-AP 2023).

The code is a very simple series of hacked together C/C++/Python tools - sorry about that.
I am very happy to help support any experiments you might want to run with our code, please make an issue if you need help.

## Citation
```
@inproceedings{lssy,
author = {Mackenzie, Joel and Moffat, Alistair},
title = {Lossy Compression Options for Dense Index Retention},
year = {2023},
url = {https://doi.org/10.1145/3624918.3625316},
doi = {10.1145/3624918.3625316},
booktitle = {Proceedings of the Annual International ACM SIGIR Conference on Research and Development in Information Retrieval in the Asia Pacific Region},
pages = {185–194},
series = {SIGIR-AP '23}
}
```


## Building the Code

The code, as mentioned, is very simple. The only dependencies are a modern C++ compiler with Intel TBB to
allow parallel sorting, a C compiler, and FAISS (etc) if you're using our Python scripts.

Simply run `make` after adjusting the makefile as suitable.

## Data

The data we used in our experiments is available at the following locations. If you cannot access it, please make an issue and we can share our copy.
Many of these indexes were part of [Pyserini's prebuilt indexes](https://github.com/castorini/pyserini/blob/master/pyserini/prebuilt_index_info.py)
or experimental artefacts from related work. We thank those authors for hosting the data.

- [Con-Arg](https://rgw.cs.uwaterloo.ca/pyserini/indexes/faiss.beir-v1.0.0-arguana.contriever.20230124.tar.gz) 
- [DPR-Wiki](https://www.dropbox.com/s/nq62qhodd237p9t/dindex-dpr-multi-pca128.tar.gz) 
- [ANCE-MARCO](https://rgw.cs.uwaterloo.ca/JIMMYLIN-bucket0/pyserini-indexes/dindex-msmarco-passage-ance-bf-20210224-060cef.tar.gz)


## Indexing

The general pipeline is as follows.


### Step 0: Setup

 We assume you have a flat FAISS index (`IndexFlatIP`) stored in a file called `my_flat.idx`.
See [FAISS](https://github.com/facebookresearch/faiss/wiki/Faiss-indexes) for more information.

### Step 1: Convert FAISS to a simple, sorted, float file.
Our quantization code assumes a specific file format that makes it easy to compute quantization boundaries.
We call these files `.sidx` files (simple sorted index).

The format is as follows:
```
<num_rows> <num_cols> <num_rows * num_columns>
``` 
The first two entries are unsigned 64-bit integers `size_t`; the remaining data is stored as 32-bit floats `float32_t`.

You can do this using the `faiss2simple` program which reads your flat index and outputs a corresponding sidx file.
```
./faiss2simple my_flat.idx my_resulting.sidx
```
NOTE: The `.sidx` file does not keep vectors together; it sorts each individual float in an increasing order.


### Step 2: Build your bins
The next step is to `quantize` the data into bins via the`sidx` file. You need to provide some arguments.
```
./quantize <number of bins> <bin type> <example.sidx> <your.bins>
```
- Number of bins: How many quantization bins to use? Default is 256.
- Bin type: This is the quantization scheme.
- `your.bins` is the output file with the bin ranges for quantization.

The bin type can be any of 1, 2, 3, or 4:
  -- bintype=1 for FD
  -- bintype=2 for FR
  -- bintype=3 for GD
  -- bintype=4 for CFR

### Step 3: Compress your index
Once you have the bins file, you are ready to encode your index; the program reads the bins file from
the quantizer and a FAISS index; it outputs the compressed index
```
./encoder <your.bins> <your-faiss-flat.idx> <your-faiss-flat.idx.compressed>
```

### Step 4: Decompress (lossy) index for querying
You can decode the index to get back to a lossy representation that can be queried.
```
./decoder <your.bins> <your-faiss-flat.idx.compressed> <your-lossy-faiss.idx>
```
That is, `<your-lossy-faiss.idx>` can be queried to generate a run file.

