# gslib_new_code and gslib_modern_fortran
This is a working repository intended to build, prototype and test new functions for PyGSLIB in pure Fortran. 

The may effort at this time is gslib_modern_fortran. This is the GSLIB library written in modern Fortran, with compatibility with parallel processing with OpenMP. 

The objective is to create a clean GSLIB library appropriate for standalone executables, dynamic libraries, and Python-Fortran extensions (F2Py). 

TODO:
 [ ] write simple interpolator for reference
 [ ] implement interpolator with pointer functions, so it can allow stacking multiple interpolation functions.   
 [ ] implement kdtree for neighborhood
 [ ] implement efficient search for gridded data