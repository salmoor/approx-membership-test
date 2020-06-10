# Approximate Membership Test

The program performs approximate membership test of k-mers for given DNA sequences.

# Primer

The k-mer is a sequence of k consecutive letters. Within the scope of nucleotides 4-mer could be ACAG (Adenine-Cytosine-Adenine-Guanine).

# Motivation

Comparison of frequencies of k-mers of different sequences is computationally cheaper and can still give an insight of how much two sequences are related. This can also help determine if it is worth to perform sequence alignment.

# Installation

```bash
$ git clone https://github.com/salmoor/approx-membership-test.git
$ cd approx-membership-test
$ make
```

# How to Use

Program uses **Bloom Filter** to perform approximate membership test for DNA sequences. it calculates the overlaps between two [FASTA](https://en.wikipedia.org/wiki/FASTA_format) files. Program will take two separate FASTA files. First, *reference*, from which to query and second, *query* file, which will be used to query. Both FASTA files may include arbitrary number of DNA sequences of arbitrary length until they are in FASTA format.  The program implements 3 hash functions to encode k-mers in the Bloom Filer.

After the Installation step, you should have executable named *bloomFilter*.

You can run the program in the following way:

```bash
$ ./bloomFilter --ref reference.fasta --query query.fasta --kmer 4 --bloomsize 8
```

The input parameters are as follows:

- **−−ref** FASTA-formatted file that contains the sequences to index.
- **−−query** FASTA-formatted file that contains the sequences to search in reference.
- **−−kmer** k-mer length.
- **−−bloomsize** Size of the Bloom filter bit vector in bytes.

So, in the above example our reference file is *reference.fasta*, query file is *query.fasta*, we use k-mers of length *4,* i.e. *4-mers*. Size of the Bloom filter bit vector is *8* bytes (*32* bits).