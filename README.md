# GSLIB in modern Fortran
This is a working repository intended to build, prototype and test new functions for PyGSLIB in pure Fortran. 

The may effort at this time is gslib_modern_fortran. This is the GSLIB library written in modern Fortran, with compatibility with parallel processing with OpenMP. 

The objective is to create a clean GSLIB library appropriate for standalone executables, dynamic libraries, and Python-Fortran extensions (F2Py). 

TODO
-----
 - [ ] write simple interpolator for reference
 - [ ] implement interpolator with pointer functions, so it can allow stacking multiple interpolation functions.   
 - [ ] implement kdtree for neighborhood
 - [ ] implement efficient search for gridded data
 - [ ] add test with open OpenBLAS, and probably a helper function
 - [ ] find way to link OpenBLAS without absolute path


Note 
-----
Use LAPACK as a replacement for ksol.f90 and ktsol.f90 
to compile
- Install OpenBLAS precompiled. 
        ```conda install -c conda-forge openblas```
- Then, just link the library. Here an example 
        ``` gfortran -o .\test.exe .\test.f90 -L "OpenBLAS-0.3.15" "C:\Users\AMartinez\Miniconda3\pkgs\openblas-0.3.15-pthreads_h543f93c_0\Library\lib\openblas.lib" ```
