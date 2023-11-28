import faiss
import numpy as np
import pickle
import argparse
import os
import pandas as pd
import random

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--faiss-index', type=str, help='path to FAISS index', required=True)
    parser.add_argument('--query-count', type=int, help='how many queries to generate', required=True)
    parser.add_argument('--output', type=str, help='output trec file', required=True)
    args = parser.parse_args()

    # (1): Load the index
    print ("Loading FAISS index from file")
    idx = faiss.read_index(args.faiss_index)
 
    # (2): Generate a range of random identifiers
    if args.query_count >= idx.ntotal:
        print ("Warning: We're going to sample every document as a query...")
    qids = sorted(random.sample(range(0, idx.ntotal), args.query_count))

    # (3): Generate the queries
    query_vectors = idx.reconstruct_batch(qids)
            
    # (4): Export them
    print ("Exporting the queries...")
    output = open(args.output + ".queries", 'wb')
    np.savez(output, qids=qids, query_vectors=query_vectors)

    # Get them back like this
    #d = np.load(args.output + ".npz")
    #print (d["qids"])
    #print (d["query_vectors"])
