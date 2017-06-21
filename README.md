# CD-HIT 

## How to compile?
  1. Compile with multi-threading support (default): 
```bash
make
````
  2. Compile without multi-threading support (if you are on very old systems): 
```bash
make openmp=no
```

To compile cd-hit-auxtools:
```bash
cd cd-hit-auxtools
make
```

To run psi-cd-hit please first download legacy BLAST (not BLAST+) and install the executables in your $PATH

For more information, please visit http://cd-hit.org or read the [cdhit-users-guide.pdf](https://github.com/weizhongli/cdhit/blob/master/doc/cdhit-user-guide.pdf). 
Most up-to-date documents are available at http://weizhongli-lab.org/cd-hit/wiki/doku.php?id=cd-hit_user_guide.

cd-hit was originally hosted at Google Code, some of the old releases are still available from https://code.google.com/p/cdhit/.

cd-hit is also available as web server, visit http://cd-hit.org for web server address.
