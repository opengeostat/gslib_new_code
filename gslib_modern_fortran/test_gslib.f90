! Program to test gslib.f90

! to compile use command:
!    gfortran -fopenmp -o gslib test_gslib.f90 gslib.f90

! to run:
!    gslib > test_gslib_results.txt

! TODO: check carefully the test for more possible errors and expected output
! issues, all subroutine using rotation matrix have risk of memory segmentation if ind > num of rotation matrix
! test functions

subroutine test_qsort_serial()
    use gslib
    implicit none
    integer, parameter :: na = 6
    real, dimension (na) :: A = [3.,5.,2.,1.,8.,2.0]
    integer, dimension (na) :: I = [1,2,3,4,5,6]

    call QSort(na,A,I,1,na)
    print *,  A
    print *,  I

end subroutine test_qsort_serial

subroutine test_qsort()
    use gslib
    implicit none
    integer, parameter :: na = 6
    real, dimension (na) :: A 
    integer, dimension (na) :: I 

    !$OMP PARALLEL private(A, I)
        A = [3.,5.,2.,1.,8.,2.0]
        I = [1,2,3,4,5,6]
        call QSort(na,A,I,1,na)
        print *,  A
        print *,  I
    !$OMP END PARALLEL

end subroutine test_qsort

subroutine test_ulcg()
    use gslib
    implicit none
    real, dimension (5) :: x
    integer, dimension (5) :: ind
    integer(i64) :: old_state
    integer :: i

    !$OMP PARALLEL private(old_state, x, i, ind)
        x = 0.
        old_state = 87687                        ! may set the seed within the parallel region
        do i=1 , 5
            old_state = lcg(old_state)
            x(i) = ulcg(old_state)
        end do
        call QSort(5,x,ind, 1, 5)
        print *,  x
    !$OMP END PARALLEL

end subroutine test_ulcg

subroutine test_ulcg_serial()
    use gslib
    implicit none
    real, dimension (5) :: x
    integer, dimension (5) :: ind
    integer(i64) :: old_state = 87687_i64
    integer :: i

    do i=1 , 5
        old_state = lcg(old_state)
        x(i) = ulcg(old_state)
    end do
    call QSort(5,x,ind, 1, 5)
    print *,  x

end subroutine test_ulcg_serial


subroutine test_ulcg_serial2()
    use gslib
    implicit none
    real, dimension (500) :: x
    integer(i64) :: old_state = 876878767
    integer :: i
    real :: mean, variance 

    do i=1 , 500
        old_state = lcg(old_state)
        x(i) = ulcg(old_state)
    end do
    mean = sum(x)/500
    variance = sum(x*x)/500-mean*mean

    print *, 'experimental mean, variance',  mean, variance
    print *, 'expected     mean, variance',  0.5, 1./12.

end subroutine test_ulcg_serial2

subroutine test_nscore_serial()
    use gslib
    implicit none
    integer, parameter :: nd = 10
    real, dimension(nd) :: vr, wt, prob, vrg
    integer, dimension(nd) :: ind
    integer :: ierror

    vr = [3.,2.,1.,4.,5.,6.,7.,8.,9.,10.]
    wt = [1.,1.,1.,1.,1.,1.,1.,1.,1.,1.]

    ! wt not defined
    call nscore(nd, vr,vrg = vrg , prob = prob, ind = ind, ierror = ierror)

    if (ierror == 0) then 
        print *, 'gauss value' , vrg
        print *, 'provavility' , prob
    else
        print *, 'error', ierror
    end if

end subroutine test_nscore_serial

subroutine test_nscore()
    use gslib
    implicit none
    integer, parameter :: nd = 10
    real, dimension(nd) :: vr, wt, prob, vrg
    integer, dimension(nd) :: ind
    integer :: ierror

    vr = [3.,2.,1.,4.,5.,6.,7.,8.,9.,10.]
    wt = [1.,1.,1.,1.,1.,1.,1.,1.,1.,1.]

    !$OMP PARALLEL private(vrg, prob, ind, ierror)
        ! wt not defined
        call nscore(nd, vr,vrg = vrg , prob = prob, ind = ind, ierror = ierror)
        if (ierror == 0) then 
            print *, 'gauss value' , vrg
            print *, 'provavility' , prob
            print *, ''
        else
            print *, 'error', ierror
        end if
   !$OMP END PARALLEL

end subroutine test_nscore


subroutine test_backtr_serial()
    
    use gslib
    implicit none
    ! inputs
    integer, parameter :: nt = 10
    integer :: ltail, utail, ierror
    real, dimension(nt) ::  vr, vrg, wt, prob
    integer, dimension(nt) ::  ind
    real :: ltpar, utpar, zmin, zmax, vrgs
    
    ! create transformation table
    vr = [3.,2.,1.,4.,5.,6.,7.,8.,9.,10.]
    wt = [1.,1.,1.,1.,1.,1.,1.,1.,1.,1.]

    call nscore(nt, vr, wt, vrg , prob, ind, ierror)

    if (ierror>0) then
        print *, "error "
        return
    end if

    !print sorted
    print *, "transformation table"
    print *, 'raw  ', vr(ind)   ! transformation table
    print *, 'gauss', vrg(ind)  ! transformation table

    !back transform
    vrgs = -1.1 ! normal value
    ltail = 1 
    ltpar = 1.
    utail = 1
    utpar = 1.
    zmin = 0
    zmax = 12
    print *, 'gaussian to back transform', vrgs
    print *, 'raw value back transformed', backtr(vrgs,nt,vr(ind),vrg(ind),zmin,zmax,ltail,ltpar, utail,utpar)

end subroutine test_backtr_serial


subroutine test_backtr()
    
    use gslib
    implicit none
    ! inputs
    integer, parameter :: nt = 10
    integer :: ltail, utail, ierror
    real, dimension(nt) ::  vr, vrg, wt, prob
    integer, dimension(nt) ::  ind
    real :: ltpar, utpar, zmin, zmax, vrgs


    !$OMP PARALLEL private(vr, vrg, wt, prob, ind, ierror, ltpar, utpar, zmin, zmax, vrgs, ltail, utail)

        ! create transformation table
        vr = [3.,2.,1.,4.,5.,6.,7.,8.,9.,10.]
        wt = [1.,1.,1.,1.,1.,1.,1.,1.,1.,1.]

        call nscore(nt, vr, wt, vrg , prob, ind, ierror)

        if (ierror>0) then
            print *, "error "
        else

            !print sorted
            print *, "transformation table"
            print *, 'raw  ', vr(ind)   ! transformation table
            print *, 'gauss', vrg(ind)  ! transformation table

            !back transform
            vrgs = -1.1 ! normal value
            ltail = 1 
            ltpar = 1.
            utail = 1
            utpar = 1.
            zmin = 0
            zmax = 12
            print *, 'gaussian to back transform', vrgs
            print *, 'raw value back transformed', backtr(vrgs,nt,vr(ind),vrg(ind),zmin,zmax,ltail,ltpar, utail,utpar)

        end if

    !$OMP END PARALLEL

end subroutine test_backtr


subroutine test_gcum_serial()
    
    use gslib
    implicit none
    ! inputs
    integer, parameter :: nt = 3
    real, dimension(nt) ::  vrg
    
    ! create array with gaussians
    vrg = [-1.64, -1.04, 10.64]

    print *, 'gaussian', vrg
    print *, 'cdf     ', gcum(vrg)

end subroutine test_gcum_serial

subroutine test_gcum()
    use gslib
    implicit none
    ! inputs
    integer, parameter :: nt = 3
    real, dimension(nt) ::  vrg
    
    !$OMP PARALLEL private(vrg)
        ! create array with gaussians
        vrg = [-1.64, -1.04, 10.64]

        print *, 'gaussian', vrg
        print *, 'cdf     ', gcum(vrg)

    !$OMP END PARALLEL

end subroutine test_gcum

subroutine test_powint()
    use gslib
    implicit none
    ! inputs
    real  , dimension (3):: rxlow, rxhigh, rylow, ryhigh, rxval, rpow
    real*8, dimension (3):: dxlow, dxhigh, dylow, dyhigh, dxval, dpow
    
    !$OMP PARALLEL private(rxlow, rxhigh, rylow, ryhigh, rxval, rpow, dxlow, dxhigh, dylow, dyhigh, dxval, dpow)
        ! create arrays and scalars
        rxlow  = [0.,1.,2.]
        rylow  = [0.,1.,2.] 
        rxhigh = [1.,2.,3.]
        ryhigh = [1.,2.,3.]  
        rxval =  [.5,1.5,2.5]
        rpow = 1.

        dxlow  = dble(rxlow)
        dylow  = dble(rylow) 
        dxhigh = dble(rxhigh)
        dyhigh = dble(ryhigh)  
        dxval =  dble(rxval)
        dpow = 1.
         
        print *, 'reals',  powint(rxlow,rxhigh,rylow,ryhigh,rxval,rpow)
        print *, 'doubl', dpowint(dxlow,dxhigh,dylow,dyhigh,dxval,dpow)

    !$OMP END PARALLEL

end subroutine test_powint

subroutine test_locate()
    use gslib
    implicit none
        
    ! input variables
    integer, parameter :: n = 5
    integer:: is, ie
    real:: x
    real, dimension(n) :: xx

    !$OMP PARALLEL private(xx, is, ie, x)
        ! create arrays and scalars
        xx  = [0.,1.,2.,3.,4.]
        is = 1
        ie = n
        x = 3.1
        print *, 'x = ', x
        print *, 'xx =', xx
        print *, 'x is between xx(j) and xx(j+1), where j =', locate(xx,n,is,ie,x)

    !$OMP END PARALLEL

end subroutine test_locate


subroutine test_dlocate()
    use gslib
    implicit none
        
    ! input variables
    integer, parameter :: n = 5
    integer:: is, ie
    real*8:: x
    real*8, dimension(n) :: xx

    !$OMP PARALLEL private(xx, is, ie, x)
        ! create arrays and scalars
        xx  = [0.,1.,2.,3.,4.]
        is = 1
        ie = n
        x = 3.1
        print *, 'x = ', x
        print *, 'xx =', xx
        print *, 'x is between xx(j) and xx(j+1), where j =', dlocate(xx,n,is,ie,x)

    !$OMP END PARALLEL

end subroutine test_dlocate


subroutine test_getindx() 
    use gslib
    implicit none
    real, dimension(3) :: loc
    real :: min,siz
    min = 0
    siz = 20
    loc = [3.8, 45.9, 300.5]

    print *, 'in a grid with xmin =', min, 'and cell size =', siz
    print *, 'the points with coordinate', loc
    print *, 'are in the grid index', getindx(min,siz,loc)

end subroutine test_getindx

subroutine test_ordrel
    use gslib
    implicit none

    ! inputs
    integer, parameter:: ncut = 6
    integer :: ivtype                            !input
    real, dimension(ncut) :: ccdf                !input
    real, dimension(ncut) :: ccdfo, aviol, xviol !oputput
    integer, dimension(ncut) :: nviol            !order violations

    ! example from https://geostatisticslessons.com/images/mikoverview/correction.png
    ccdf = [0.25 , 0.45, 0.30, 0.50, 0.75, 0.90]
    ivtype = 1                               

    call ordrel_cont(ncut,ccdf,ccdfo,nviol,aviol,xviol)

    print *, 'cdf raw                ', ccdf
    print *, 'cdf corrected          ', ccdfo
    print *, 'number of violations   ', nviol
    print *, 'average of violations  ', aviol
    print *, 'magnitude of violations', xviol

end subroutine test_ordrel


subroutine test_setrot()
    use gslib
    implicit none

    ! inputs
    integer, parameter :: nrotmat = 1
    integer :: ind
    real:: ang1,ang2,ang3,anis1,anis2

    !output
    real*8, dimension(nrotmat,3,3) :: rotmat

    ang1 = 45.
    ang2 = 35.
    ang3 = 70.
    anis1 = 2.
    anis2 = 3.

    ind = 1

    call setrot(ang1,ang2,ang3,anis1,anis2,ind,nrotmat,rotmat)
    print *, 'rotation matrix for ang1,ang2,ang3,anis1,anis2', ang1,ang2,ang3,anis1,anis2
    print *, real(rotmat(1,1,:))
    print *, real(rotmat(1,2,:))
    print *, real(rotmat(1,3,:))
    print *, 'rotation is a pibot point 0,0,0'

end subroutine test_setrot

subroutine test_sqdist()
    use gslib
    implicit none

    integer, parameter:: nrotmat = 1
    integer:: ind 
    real*8, dimension(nrotmat,3,3) :: rotmat
    real :: x1,y1,z1, x2,y2, z2
    real:: ang1,ang2,ang3,anis1,anis2

    ! get rotation matrix
    ang1 = 45.
    ang2 = 35.
    ang3 = 70.
    anis1 = 2.
    anis2 = 3.

    ind = 1

    call setrot(ang1,ang2,ang3,anis1,anis2,ind,nrotmat,rotmat)

    ! calc sqdist
    x1 = 0.
    y1 = 0.
    z1 = 0.
    x2 = 1.
    y2 = 1.
    z2 = 0.

    print *, 'coordinates x1,y1,z1,x2,y2,z2', x1,y1,z1,x2,y2,z2
    print *, 'parameter for rot matrix with ang1,ang2,ang3,anis1,anis2', ang1,ang2,ang3,anis1,anis2
    print *, 'isotropic dist',  real(sqdist(x1,y1,z1,x2,y2,z2,ind,nrotmat,rotmat))

end subroutine test_sqdist

subroutine test_sqdist_err()
    use gslib
    implicit none

    integer, parameter:: nrotmat = 1
    integer:: ind 
    real*8, dimension(nrotmat,3,3) :: rotmat
    real*8 :: distsq
    real :: x1,y1,z1, x2,y2, z2
    real:: ang1,ang2,ang3,anis1,anis2

    ! get rotation matrix
    ang1 = 45.
    ang2 = 35.
    ang3 = 70.
    anis1 = 2.
    anis2 = 3.

    ind = 2

    call setrot(ang1,ang2,ang3,anis1,anis2,ind,nrotmat,rotmat)
    
    ! calc sqdist
    x1 = 0.
    y1 = 0.
    z1 = 0.
    x2 = 1.
    y2 = 1.
    z2 = 0.

    distsq = sqdist(x1,y1,z1,x2,y2,z2,ind,nrotmat,rotmat)
    print *, 'rotation matrix',  rotmat
    print *, 'distance is nan?', ieee_is_nan(distsq)
    print *, 'isotropic dist',  distsq

end subroutine test_sqdist_err


subroutine test_beyond()
    use gslib
    implicit none
    real :: cdfval, zval

    cdfval = beyond_cdf(nccut = 5, ccut = [0.1,0.2, 0.3, 0.5, 0.8], &
                        ccdf = [0.1,0.2, 0.3, 0.5, 0.8], ncut = 1,cut = [0.5],cdf = [0.5], &
                        zmin = 0.,zmax = 10., ltail = CDF_TAIL_LINEAR, ltpar = 1., &
                        middle  = CDF_TAIL_LINEAR, mpar =1., &
                        utail = CDF_TAIL_LINEAR,utpar = 1.1,zval = 0.11) 

     zval = beyond_zval(ivtype = VARTYPE_CONTINOUS, nccut = 5, ccut = [0.1,0.2, 0.3, 0.5, 0.8], &
                        ccdf = [0.1,0.2, 0.3, 0.5, 0.8], ncut = 1,cut = [0.5],cdf = [0.5], &
                        zmin = 0.,zmax = 10., ltail = CDF_TAIL_LINEAR,ltpar = 1., &
                        middle  = CDF_TAIL_LINEAR, mpar =1., &
                        utail = CDF_TAIL_LINEAR,utpar = 1.1,cdfval = cdfval) 

    print *, 'cdfval for zval 0.11 is', cdfval
    print *, 'zval for this cdfval', zval
    
end subroutine test_beyond

subroutine test_cova3()
    use gslib
    implicit none

    integer, parameter:: nrotmat = 1
    integer:: ind 
    real*8, dimension(nrotmat,3,3) :: rotmat
    real :: x1,y1,z1, x2,y2, z2
    real:: ang1,ang2,ang3,anis1,anis2, c

    ! get rotation matrix
    ang1 = 45.
    ang2 = 35.
    ang3 = 70.
    anis1 = 2.
    anis2 = 3.

    ind = 1

    call setrot(ang1,ang2,ang3,anis1,anis2,ind,nrotmat,rotmat)

    ! calc sqdist
    x1 = 0.
    y1 = 0.
    z1 = 0.
    x2 = 1.
    y2 = 1.
    z2 = 0.

    c= cova3(x1,y1,z1,x2,y2,z2,ivarg =1 ,nst = 2,c0 = [0.25], &
             it =[COV_SPHERICAL,COV_SPHERICAL],cc=[0.25, 0.5], &
             aa = [1., 3.], irot = ind,nrotmat=nrotmat, rotmat = rotmat)

    print *, 'isotropic covariance', c

end subroutine test_cova3

program test_gslib


    use gslib
    implicit none

    print *, ''
    print *, 'test qsort'
    print *, 'expected result'
    print *, '   1.0       2.0       2.0      3.0      5.00       8.0'
    print *, '     4         6         3        1        2          5'
    call  test_qsort()

    print *, ''
    print *, 'test ulcg'
    print *, 'expected result'
    print *, "same sequence n core times."
    call  test_ulcg()

    print *, ''
    print *, 'test ulcg serial version'
    print *, 'expected result'
    print *, "7.47933537E-02  0.221471041      0.389764041      0.548619866      0.776747882"
    call  test_ulcg_serial()

    print *, ''
    print *, 'test ulcg serial version mean and variance'
    print *, 'expected result are the mean and variance of the uniform distribution with parameter a=0, b=1'
    call  test_ulcg_serial2()

    print *, ''
    print *, 'test nscore, gauinv, and qsort'
    print *, 'expected result: gaussian    {-0.67,  -1.036, ..., 1.04, 1.64} repated ncores'
    print *, 'expected result: probability {0.25,     0.15, ..., 0.85, 0.95} repated ncores'
    print *, 'expected result: gaussian and probability repated ncores'
    call  test_nscore()
    
    print *, ''
    print *, 'test nscore, gauinv, and qsort'
    print *, 'expected result: gaussian    {-0.67,  -1.036, ..., 1.04, 1.64} '
    print *, 'expected result: probability {0.25,     0.15, ..., 0.85, 0.95} '
    call  test_nscore_serial()

    print *, ''
    print *, 'test backtr'
    print *, 'expected result: TODO '
    call  test_backtr_serial()

    print *, ''
    print *, 'test backtr'
    print *, 'expected result: TODO '
    call  test_backtr()

    print *, ''
    print *, 'test gcum'
    print *, 'expected result:  gaussian  -1.63999999      -1.03999996       1.63999999 '
    print *, '                  cdf        5.05025983E-02  0.149170041      0.949497402 '
    call  test_gcum_serial()

    print *, ''
    print *, 'test gcum'
    print *, 'expected result:  gaussian  -1.63999999      -1.03999996       10.6400003 '
    print *, '                  cdf        5.05025983E-02  0.149170041       1.00000000 '
    call  test_gcum()

    print *, ''
    print *, 'test powint and dpowint'
    print *, 'expected result:   reals  0.50       1.50       2.50'
    print *, '                   doubl  0.50       1.50       2.50'
    call  test_powint()

    print *, ''
    print *, 'test powint and dpowint'
    print *, 'expected result:  >>>> x is between xx(j) and xx(j+1), where j = 4'
    call  test_locate()

    print *, ''
    print *, 'test powint and dpowint'
    print *, 'expected result:  >>>> x is between xx(j) and xx(j+1), where j = 4'
    call  test_dlocate()

    print *, ''
    print *, 'test getindx'
    print *, 'expected result:  are in the grid index  1   3   16'
    call  test_getindx()

    print *, ''
    print *, 'test ordrel'
    print *, 'expected result:   cdf raw         0.25      0.45      0.30     0.50     0.75     0.90'
    print *, '                   cdf corrected   0.25      0.375     0.375    0.50     0.75     0.90'
    print *, '                   number of violations      0    1         1         0     0      0 '
    print *, '                   average of violations     0    7.5E-02   7.5E-02   0.0   0.0    0.0 '
    print *, '                   magnitude of violations   0    7.5E-02   7.5E-02   0.0   0.0    0.0'
    call  test_ordrel()
    

    print *, ''
    print *, 'test ordrel'
    print *, 'expected result:    '
    print *, '                rotation matrix for ang1,ang2,ang3,anis1,anis2  45. 35. 70. 2.  3.'
    print *, '                    0.579      0.579    0.574'
    print *, '                   -0.311     -7.E-02   0.385'
    print *, '                    0.175     -0.268    9.3E-02'
    call  test_setrot()

    print *, ''
    print *, 'test ordrel'
    print *, 'expected result:    '
    print *, 'coordinates x1,y1,z1,x2,y2,z2   0. 0. 0. 1. 1. 0.'
    print *, 'parameter for rot matrix with ang1,ang2,ang3,anis1,anis2   45. 35. 70. 2. 3.'
    print *, 'isotropic dist   1.49582493'
    call  test_sqdist()

    print *, ''
    print *, 'test ordrel'
    print *, 'expected result:   array of nans, true, and nan '
    call  test_sqdist_err()

    print *, ''
    print *, 'test cova3'
    print *, 'expected result:   isotropic covariance  0.211179554'
    call  test_cova3()

    print *, ''
    print *, 'test beyond_cdf and beyond_zval '
    print *, 'expected result:   cdfval for zval 0.11 is  0.109999999'
    print *,  '                  zval for this cdfval  0.109999999'
    call  test_beyond()

end program