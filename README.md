# CD-HIT - Representative Sequences

CD-HIT is a very widely used program for clustering and comparing protein or nucleotide sequences. CD-HIT was originally developed by [Dr. Weizhong Li](liwz@sdsc.edu) at [Dr. Adam Godzik's Lab](http://bioinformatics.burnham.org/) at the [Burnham Institute (now Sanford-Burnham Medical Research Institute)](http://www.sanfordburnham.org/)

CD-HIT is very fast and can handle extremely large databases. CD-HIT helps to significantly reduce the computational and manual efforts in many sequence analysis tasks and aids in understanding the data structure and correct the bias within a dataset.

- CD-HIT (CD-HIT-EST) clusters similar proteins (DNAs) into clusters that meet a user-defined similarity threshold.
- CD-HIT-2D (CD-HIT-EST-2D) compares 2 datasets and identifies the sequences in db2 that are similar to db1 above a threshold.
- CD-HIT-454 identifies natural and artificial duplicates from pyrosequencing reads.
- CD-HIT-OTU clusters rRNA tags into OTUs
- CD-HIT-DUP identifies duplicates from single or paired Illumina reads
- CD-HIT-LAP identifies overlapping reads 

## How to compile?
1. Compile with multi-threading support (default): 
```
$ make
```
2. Compile without multi-threading support (if you are on very old systems):
```
$ make openmp=no
```

### For cd-hit-auxtools
```
$ cd cd-hit-auxtools
$ make
```

### For psi-cd-hit
  Please download legacy BLAST (not BLAST+) and install the executables in your $PATH

## How to use

Basic usage:
```
$ cd-hit -i input_fasta.fasta -o output_fasta.fasta
```

With custom sequence identity threshold, number of threads and maximum memory (in MB):
```
$ cd-hit -i input_fasta.fasta -o output_fasta.fasta -T 4 -M 4000 -c 0.98
```
## Web Server

cd-hit is also available as web server, visit http://cd-hit.org for web server address.

## More Information

For more information, please visit http://cd-hit.org.

Most up-to-date documents are available at https://github.com/weizhongli/cdhit/wiki.