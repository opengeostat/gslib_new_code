! Program to test neighborhood.f90

! to compile use command:
!    gfortran -fopenmp -o neighborhood test_neighborhood.f90 neighborhood.f90

! to run:
!    neighborhood > test_neighborhood_results.txt

program test_neighborhood


    use neighborhood
    implicit none

    real :: x(5) = [1,2,3,4,5]
    real :: y(5) = [1,2,3,4,5]
    real :: z(5) = [0,0,0,0,0]

    print *,  ball (x, y, z, x0 = 1., y0 = 1., z0 =1., r =2.)
end program