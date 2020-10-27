# newspecmatch
`newspecmatch` is a small tool for the quantitativ comparison of IR spectra.
The comparison is based on different similarity measures:
- matchscore (MSC) based on the Cauchy-Schwarz inequality
- euclidean norm (EUC)
- Pearson correlation coefficient (PCC)
- Spearman rank correlation coefficient (SCC)

For usage and options see `newspecmatch --help`

## Installation
The source directory contains a `Makefile` that can be used with the `ifort` Intel Fortran compiler and `make`:
```bash
cd src/
make
```

Otherwise, the program can also be compiled using the `gfortran` GNU Fortran compiler:
```bash
cd src/
gfortran spectramod.f90 newspecmatch.f90 main.f90 -o newspecmatch
```

After having build the binary with either `ifort` or `gfortran`, simply move it to some place included in your `PATH` variable.


## Literature
- P. Pracht, D.F. Grant, S. Grimme, *Journal of Chemical Theory and Computation*, **2020**,
  DOI: [10.1021/acs.jctc.0c00877](https://pubs.acs.org/doi/abs/10.1021/acs.jctc.0c00877)
