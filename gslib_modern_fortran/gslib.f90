module gslib

    implicit none

    real, parameter, private :: EPSLON=1.0e-20, UNEST=-1.0, PI=3.14159265,PMX=999., DEG2RAD=3.141592654/180.0

    contains

        !TODO: get sorting subroutine to replace sortem

        subroutine nscore(nd,vr,wt,tmp,vrg,ierror)
            !-----------------------------------------------------------------------
            !              Transform Univariate Data to Normal Scores
            !              ******************************************
            ! This subroutibe takes "nd" data "vr(i),i=1,...,nd" possibly weighted
            ! by "wt(i),i=,...,nd" and returns the normal scores transform N(0,1)
            ! as "vrg(i),i=1,...,nd".  The extra storage array "tmp" is required
            ! so that the data can be returned in the same order (just in case
            ! there are associated arrays like the coordinate location).
            !
            ! INPUT VARIABLES:
            !   nd               Number of data (no missing values)
            !   vr(nd)           Data values to be transformed
            !   tmin,tmax        data trimming limits             
            !   wt(nd)           Weight for each data (don't have to sum to 1.0)
            !                    this is an optional variable, if present will use it, 
            !                    otherwise will use equal weighted                   
            !   tmp(nd)          Temporary storage space for sorting
            !
            ! OUTPUT VARIABLES:
            !   vrg(nd)          normal scores
            !   ierror           error flag (0=error free,1=problem)
            !
            ! EXTERNAL REFERENCES:
            
            !   gauinv           Calculates the inverse of a Gaussian cdf
            !   sortem           sorts a number of arrays according to a key array
            !-----------------------------------------------------------------------

            ! inputs
            integer, intent(in) :: nd 
            real, intent(in), dimension(nd) :: vr 
            real, intent(in), dimension(nd), optional :: wt
            real, intent(inout), dimension(nd) :: tmp         ! this is a temp variable for sorting, an ugly solution but efficient

            !output 
            real, intent(out), dimension(nd) :: vrg 
            integer, intent(out) :: ierror


            ! internal variables
            real, dimension(nd) :: wtt
            real*8 :: pd
            real   :: twt, oldcp, cp
            integer :: i 

            ! Sort the data in ascending order and calculate total weight:
        
            ierror = 0
            twt    = 0.0
            do i=1,nd
                tmp(i) = real(i)
                if(present(wt)) then
                    twt = twt + 1.
                else
                    twt = twt + wt(i)
                end if
            end do
            if(nd < 1 .OR. twt < EPSLON) then
                ierror = 1
                return
            end if

            call sortem(1,nd,vr,2,wt,tmp,d,e,f,g,h)
        
            ! Compute the cumulative probabilities:
        
            oldcp = 0.0
            cp    = 0.0
            do i=1,nd
                cp     =  cp + wt(i) / twt
                wtt(i)  = (cp + oldcp)/ 2.0
                oldcp  =  cp
                call gauinv(dble(wt(i)),vrg(i),ierror)
                if (ierror>0) return
            end do
        
            ! Get the arrays back in original order:
        
            call sortem(1,nd,tmp,3,wtt,vr,vrg,e,f,g,h)
        
        end subroutine nscore

        subroutine gauinv(p,xp,ierror)
            !-----------------------------------------------------------------------
            ! Computes the inverse of the standard normal cumulative distribution
            ! function with a numerical approximation from : Statistical Computing,
            ! by W.J. Kennedy, Jr. and James E. Gentle, 1980, p. 95.
            !
            ! INPUT/OUTPUT:
            !   p    = double precision cumulative probability value: dble(psingle)
            !   xp   = G^-1 (p) in single precision
            !   ierr = 1 - then error situation (p out of range), 0 - OK
            !-----------------------------------------------------------------------
            ! inputs
            real*8, intent (in) :: p
            

            ! output
            real, intent (out) :: xp
            integer, intent (out) :: ierror

            ! internal variables
            real*8 :: y,pp
            real*8, parameter ::  lim=1.0e-10
            real*8, parameter ::  p0 = -0.322232431088, p1 = -1.0, p2=-0.342242088547
            real*8, parameter ::  p3 = -0.0204231210245, p4 = -0.0000453642210148
            real*8, parameter ::  q0 = 0.0993484626060, q1 = 0.588581570495, q2 = 0.531103462366
            real*8, parameter ::  q3 = 0.103537752850, q4 = 0.0038560700634
            
        
            ! Check for an error situation:
        
            ierror = 1
            if(p < lim) then
                xp = -1.0e10
                return
            end if
            if(p > (1.0-lim)) then
                xp =  1.0e10
                return
            end if
            ierror = 0
        
            ! Get k for an error situation:
        
            pp   = p
            if(p > 0.5) pp = 1 - pp
            xp   = 0.0
            if(p == 0.5) return
        
            ! Approximate the function:
        
            y  = dsqrt(dlog(1.0/(pp*pp)))
            
            xp = real( y + ((((y*p4+p3)*y+p2)*y+p1)*y+p0) / ((((y*q4+q3)*y+q2)*y+q1)*y+q0) )
            
            if(real(p) == real(pp)) xp = -xp
        
        end subroutine gauinv

        pure real function backtr(vrgs,nt,vr,vrg,zmin,zmax,ltail,ltpar, utail,utpar)
            !-----------------------------------------------------------------------
            !           Back Transform Univariate Data from Normal Scores
            !           *************************************************
            ! This subroutine backtransforms a standard normal deviate from a
            ! specified back transform table and option for the tails of the
            ! distribution.  Call once with "first" set to true then set to false
            ! unless one of the options for the tail changes.

            ! INPUT VARIABLES:
            !   vrgs             normal score value to be back transformed
            !   nt               number of values in the back transform tbale
            !   vr(nt)           original data values that were transformed
            !   vrg(nt)          the corresponding transformed values
            !   zmin,zmax        limits possibly used for linear or power model
            !   ltail            option to handle values less than vrg(1):
            !   ltpar            parameter required for option ltail
            !   utail            option to handle values greater than vrg(nt):
            !   utpar            parameter required for option utail
            !
            !-----------------------------------------------------------------------

            ! inputs
            integer, intent (in) :: nt, ltail, utail
            real, dimension(nt), intent (in) ::  vr, vrg
            real, intent (in) ::    ltpar, utpar, vrgs, zmin,zmax

            ! internal variables
            real :: lambda, cdflo, cdfbt, cpow, cdfhi
            integer :: j

            ! Value in the lower tail?    1=linear, 2=power, (3 and 4 are invalid)

            if(vrgs <= vrg(1)) then
                backtr = vr(1)
                cdflo  = gcum(vrg(1))
                cdfbt  = gcum(vrgs)
                if(ltail == 1) then
                    backtr = powint(0.0,cdflo,zmin,vr(1),cdfbt,1.0)
                else if(ltail == 2) then
                    cpow   = 1.0 / ltpar
                    backtr = powint(0.0,cdflo,zmin,vr(1),cdfbt,cpow)
                endif
            
            ! Value in the upper tail?     1=linear, 2=power, 4=hyperbolic:
            
            else if(vrgs >= vrg(nt)) then
                backtr = vr(nt)
                cdfhi  = gcum(vrg(nt))
                cdfbt  = gcum(vrgs)
                if(utail == 1) then
                    backtr = powint(cdfhi,1.0,vr(nt),zmax,cdfbt,1.0)
                else if(utail == 2) then
                    cpow   = 1.0 / utpar
                    backtr = powint(cdfhi,1.0,vr(nt),zmax,cdfbt,cpow)
                else if(utail == 4) then
                    lambda = (vr(nt)**utpar)*(1.0-gcum(vrg(nt)))
                    backtr = (lambda/(1.0-gcum(vrgs)))**(1.0/utpar)
                endif
            else
            
            ! Value within the transformation table:
            
                j = locate(vrg,nt,1,nt,vrgs)
                j = max(min((nt-1),j),1)
                backtr = powint(vrg(j),vrg(j+1),vr(j),vr(j+1),vrgs,1.0)

            endif

            ! there is a finction call getz, which is equal to this function except for this two lines
            !if(getz < zmin) getz = zmin
            !if(getz > zmax) getz = zmax
            ! getz is only used in trans.f90

        end function backtr

        pure elemental real function gcum(x)
            !-----------------------------------------------------------------------
            ! Evaluate the standard normal cdf given a normal deviate x.  gcum is
            ! the area under a unit normal curve to the left of x.  The results are
            ! accurate only to about 5 decimal places.
            !-----------------------------------------------------------------------

            !inputs
            real, intent(in) :: x

            !internal variables
            real :: z, t, e2

            z = x
            if(z < 0.) z = -z
            t    = 1./(1.+ 0.2316419*z)
            gcum = t*(0.31938153   + t*(-0.356563782 + t*(1.781477937 + &
            t*(-1.821255978 + t*1.330274429))))
            e2   = 0.
        
            !  6 standard deviations out gets treated as infinity:
        
            if(z <= 6.) e2 = exp(-z*z/2.)*0.3989422803
            gcum = 1.0- e2 * gcum
            if(x >= 0.) return
            gcum = 1.0 - gcum
            
        end function gcum


        pure elemental real function powint(xlow,xhigh,ylow,yhigh,xval,pow)
            !-----------------------------------------------------------------------
            ! Power interpolate the value of y between (xlow,ylow) and (xhigh,yhigh)
            !                 for a value of x and a power pow.
            !-----------------------------------------------------------------------

            ! input variables
            real, intent(in) :: xlow, xhigh, ylow, yhigh, xval, pow
        
            if((xhigh-xlow) < EPSLON) then
                powint = (yhigh+ylow)/2.0
            else
                powint = ylow + (yhigh-ylow)*(((xval-xlow)/(xhigh-xlow))**pow)
            end if
        
        end function powint


        pure elemental real*8 function dpowint(xlow,xhigh,ylow,yhigh,xval,pow)
            !-----------------------------------------------------------------------
            ! Power interpolate the value of y between (xlow,ylow) and (xhigh,yhigh)
            !                 for a value of x and a power pow.
            !-----------------------------------------------------------------------

            ! input variables
            real*8, intent(in) :: xlow, xhigh, ylow, yhigh, xval, pow
        
            if((xhigh-xlow) < EPSLON) then
                dpowint = (yhigh+ylow)/2.0
            else
                dpowint = ylow + (yhigh-ylow)*(((xval-xlow)/(xhigh-xlow))**pow)
            end if
        
        end function dpowint 

        pure integer function locate(xx,n,is,ie,x) result(j)
            !-----------------------------------------------------------------------
            ! Given an array "xx" of length "n", and given a value "x", this routine
            ! returns a value "j" such that "x" is between xx(j) and xx(j+1).  xx
            ! must be monotonic, either increasing or decreasing.  j=is-1 or j=ie is
            ! returned to indicate that x is out of range.
            !
            ! Bisection Concept From "Numerical Recipes", Press et. al. 1986  pp 90.
            !-----------------------------------------------------------------------
            
            ! input variables
            integer, intent (in):: n, is, ie
            real, intent (in) :: x
            real, intent (in), dimension(n) :: xx

            !internal variable
            integer :: jl, ju, jm, iss
            
            ! make internal copy of is
            iss = is

            ! Initialize lower and upper methods:
            
            if(iss <= 0) iss = 1
            jl = is-1
            ju = ie
            if(xx(n) <= x) then
                j = ie
                return
            end if
            
            ! If we are not done then compute a midpoint:
            
            10 if(ju-jl > 1) then
                jm = (ju+jl)/2
            
                ! Replace the lower or upper limit with the midpoint:
            
                if((xx(ie) > xx(iss)).eqv.(x > xx(jm))) then
                    jl = jm
                else
                    ju = jm
                endif
                go to 10
            endif
        
            ! Return with the array index:
        
            j = jl

            return
    
        end function locate


        pure integer function dlocate(xx,n,is,ie,x) result(j)
            !-----------------------------------------------------------------------
            ! Given an array "xx" of length "n", and given a value "x", this routine
            ! returns a value "j" such that "x" is between xx(j) and xx(j+1).  xx
            ! must be monotonic, either increasing or decreasing.  j=is-1 or j=ie is
            ! returned to indicate that x is out of range.
            !
            ! Bisection Concept From "Numerical Recipes", Press et. al. 1986  pp 90.
            !-----------------------------------------------------------------------
            
            ! input variables
            integer, intent (in):: n, is, ie
            real*8, intent (in) :: x
            real*8, intent (in), dimension(n) :: xx

            !internal variable
            integer :: jl, ju, jm, iss
            
            ! make internal copy of is
            iss = is

            ! Initialize lower and upper methods:
            
            if(iss <= 0) iss = 1
            jl = is-1
            ju = ie
            if(xx(n) <= x) then
                j = ie
                return
            end if
            
            ! If we are not done then compute a midpoint:
            
            10 if(ju-jl > 1) then
                jm = (ju+jl)/2
            
                ! Replace the lower or upper limit with the midpoint:
            
                if((xx(ie) > xx(iss)).eqv.(x > xx(jm))) then
                    jl = jm
                else
                    ju = jm
                endif
                go to 10
            endif
            
            ! Return with the array index:
        
            j = jl

            return
        
        end function dlocate

            
        pure subroutine beyond(ivtype,nccut,ccut,ccdf,ncut,cut,cdf,zmin,zmax, &
            ltail,ltpar,middle,mpar,utail,utpar,zval,cdfval,ierr)
            !-----------------------------------------------------------------------
            !                     Go Beyond a Discrete CDF
            !                     ************************
            ! This subroutine is a general purpose subroutine to interpolate within
            ! and extrapolate beyond discrete points on a conditional CDF.  If the
            ! Z value "zval" is specified then the corresponding CDF value "cdfval"
            ! will be computed, if the CDF value "cdfval" is specified the
            ! corresponding Z value "zval" will be computed.

            ! INPUT/OUTPUT VARIABLES:                
            !   ivtype           variable type (1=continuous, 0=categorical)
            !   nccut            number of cutoffs defining the conditional CDF
            !   ccut()           real array of the nccut cutoffs
            !   ccdf()           real array of the conditional cdf values
            !   ncut             number of cutoffs defining the global CDF
            !   cut()            real array of the ncut cutoffs
            !   cdf()            real array of the global cdf values
            !   zmin,zmax        minimum and maximum allowable data values
            !   ltail            option to handle values in lower tail
            !   ltpar            parameter required for option ltail
            !   middle           option to handle values in the middle
            !   mpar             parameter required for option middle
            !   utail            option to handle values in upper tail
            !   utpar            parameter required for option utail
            !   output variables
            !   zval             interesting cutoff (if -1 then it is calculated)
            !   cdfval           interesting CDF (if -1 then it is calculated)
            !   ierr             error code
            !       1            zval > UNEST .AND. cdfval > UNEST, 
            !                    or zval <= UNEST .AND. cdfval <= UNEST
            !       2            Unacceptable input option
            !-----------------------------------------------------------------------
            
            ! Input variables
            integer, intent (in) :: ivtype, nccut, ncut
            real, intent (in), dimension (ncut):: ccut, ccdf
            real, intent (in), dimension (1):: cut, cdf       !consider using scalar here
            real, intent (in):: zmin, zmax, ltpar, utpar, mpar
            integer, intent (in):: ltail, middle, utail
            
            !output variable
            integer, intent (out):: ierr
            real, intent (out):: zval, cdfval

            !  Internal variables
            integer ::   cclow,cchigh, i, ipart, idat, iupp, ilow
            real :: cum, powr, temp, lambda
        
            ! Check for both "zval" and "cdfval" defined or undefined:
        
            ierr  = 1
            if(zval > UNEST .AND. cdfval > UNEST) return
            if(zval <= UNEST .AND. cdfval <= UNEST) return
        
            ! Handle the case of a categorical variable:
        
            if(ivtype == 0) then
                cum = 0
                do i=1,nccut
                    cum = cum + ccdf(i)
                    if(cdfval <= cum) then
                        zval = ccut(i)
                        return
                    endif
                end do
                return
            end if
        
            ! Figure out what part of distribution: ipart = 0 - lower tail
            !                                       ipart = 1 - middle
            !                                       ipart = 2 - upper tail
            ierr  = 0
            ipart = 1
            if(zval > UNEST) then
                if(zval <= ccut(1))       ipart = 0
                if(zval >= ccut(nccut))   ipart = 2
            else
                if(cdfval <= ccdf(1))     ipart = 0
                if(cdfval >= ccdf(nccut)) ipart = 2
            endif
        
            ! ARE WE IN THE LOWER TAIL?
        
            if(ipart == 0) then
                if(ltail == 1) then
                
                    ! Straight Linear Interpolation:
                
                    powr = 1.0
                    if(zval > UNEST) then
                        cdfval = powint(zmin,ccut(1),0.0,ccdf(1), &
                        zval,powr)
                    else
                        zval = powint(0.0,ccdf(1),zmin,ccut(1), &
                        cdfval,powr)
                    endif
                else if(ltail == 2) then
                
                    ! Power Model interpolation to lower limit "zmin"?
                
                    if(zval > UNEST) then
                        cdfval = powint(zmin,ccut(1),0.0,ccdf(1), &
                        zval,ltpar)
                    else
                        powr = 1.0 / ltpar
                        zval = powint(0.0,ccdf(1),zmin,ccut(1), &
                        cdfval,powr)
                    endif
                
                    ! Linear interpolation between the rescaled global cdf?
                
                else if(ltail == 3) then
                    if(zval > UNEST) then
                    
                        ! Computing the cdf value. Locate the point and the class bound:
                    
                        idat = locate(cut,ncut,1,ncut,zval)
                        iupp = locate(cut,ncut,1,ncut,ccut(1))
                    
                        ! Straight linear interpolation if no data; otherwise, linear:
                    
                        if(idat <= 0 .OR. idat >= ncut .OR. &
                        iupp <= 0 .OR. iupp >= ncut) then
                            cdfval = powint(zmin,cut(1),0.0,cdf(1), &
                            zval,1.)
                        else
                            temp   = powint(cut(idat),cut(idat+1), &
                            cdf(idat),cdf(idat+1),zval,1.)
                            cdfval = temp*ccdf(1)/cdf(iupp)
                        endif
                    else
                    
                        ! Computing Z value: Are there any data out in the tail?
                    
                        iupp = locate(cut,ncut,1,ncut,ccut(1))
                    
                        ! Straight linear interpolation if no data; otherwise, local linear
                        ! interpolation:
                    
                        if(iupp <= 0 .OR. iupp >= ncut) then
                            zval = powint(0.0,cdf(1),zmin,cut(1), &
                            cdfval,1.)
                        else
                            temp = cdfval*cdf(iupp)/ccdf(1)
                            idat = locate(cdf,ncut,1,ncut,temp)
                            if(idat <= 0 .OR. idat >= ncut) then
                                zval = powint(0.0,cdf(1),zmin, &
                                cut(1),cdfval,1.)
                            else
                                zval = powint(cdf(idat),cdf(idat+1), &
                                cut(idat),cut(idat+1),temp,1.)
                            end if
                        endif
                    endif
                else
                
                    ! Error situation - unacceptable option:
                
                    ierr = 2
                    return
                endif
            endif
        
            ! FINISHED THE LOWER TAIL,  ARE WE IN THE MIDDLE?
        
            if(ipart == 1) then
            
                ! Establish the lower and upper limits:
            
                if(zval > UNEST) then
                    cclow = locate(ccut,nccut,1,nccut,zval)
                else
                    cclow = locate(ccdf,nccut,1,nccut,cdfval)
                endif
                cchigh = cclow + 1
                if(middle == 1) then
                
                    ! Straight Linear Interpolation:
                
                    powr = 1.0
                    if(zval > UNEST) then
                        cdfval = powint(ccut(cclow),ccut(cchigh), &
                        ccdf(cclow),ccdf(cchigh),zval,powr)
                    else
                        zval = powint(ccdf(cclow),ccdf(cchigh), &
                        ccut(cclow),ccut(cchigh),cdfval,powr)
                    endif
                
                    ! Power interpolation between class bounds?
                
                else if(middle == 2) then
                    if(zval > UNEST) then
                        cdfval = powint(ccut(cclow),ccut(cchigh), &
                        ccdf(cclow),ccdf(cchigh),zval,mpar)
                    else
                        powr = 1.0 / mpar
                        zval = powint(ccdf(cclow),ccdf(cchigh), &
                        ccut(cclow),ccut(cchigh),cdfval,powr)
                    endif
                
                    ! Linear interpolation between the rescaled global cdf?
                
                else if(middle == 3) then
                    ilow = locate(cut,ncut,1,ncut,ccut(cclow))
                    iupp = locate(cut,ncut,1,ncut,ccut(cchigh))
                    if(cut(ilow) < ccut(cclow))  ilow = ilow + 1
                    if(cut(iupp) > ccut(cchigh)) iupp = iupp - 1
                    if(zval > UNEST) then
                        idat = locate(cut,ncut,1,ncut,zval)
                    
                        ! Straight linear interpolation if no data; otherwise, local linear
                        ! interpolation:
                    
                        if(idat <= 0 .OR. idat >= ncut .OR. &
                        ilow <= 0 .OR. ilow >= ncut .OR. &
                        iupp <= 0 .OR. iupp >= ncut .OR. &
                        iupp <= ilow) then
                            cdfval=powint(ccut(cclow),ccut(cchigh), &
                            ccdf(cclow),ccdf(cchigh),zval,1.)
                        else
                            temp = powint(cut(idat),cut(idat+1), &
                            cdf(idat),cdf(idat+1),zval,1.)
                            cdfval=powint(cdf(ilow),cdf(iupp), &
                            ccdf(cclow),ccdf(cchigh),temp,1.)
                        endif
                    else
                    
                        ! Straight linear interpolation if no data; otherwise, local linear
                        ! interpolation:
                    
                        if(ilow <= 0 .OR. ilow >= ncut .OR. &
                        iupp <= 0 .OR. iupp >= ncut .OR. &
                        iupp <= ilow) then
                            zval=powint(ccdf(cclow),ccdf(cchigh), &
                            ccut(cclow),ccut(cchigh),cdfval,1.)
                        else
                            temp=powint(ccdf(cclow),ccdf(cchigh), &
                            cdf(ilow),cdf(iupp),cdfval,1.)
                            idat = locate(cdf,ncut,1,ncut,temp)
                            if(cut(idat) < ccut(cclow)) idat=idat+1
                            if(idat <= 0 .OR. idat >= ncut .OR. &
                            cut(idat+1) > ccut(cchigh)) then
                                zval = powint(ccdf(cclow), &
                                ccdf(cchigh),ccut(cclow), &
                                ccut(cchigh),cdfval,1.)
                            else
                                zval = powint(cdf(idat),cdf(idat+1), &
                                cut(idat),cut(idat+1),temp,1.)
                            end if
                            zval = powint(cdf(idat),cdf(idat+1), &
                            cut(idat),cut(idat+1),temp,1.)
                        endif
                    endif
                else
                
                    ! Error situation - unacceptable option:
                
                    ierr = 2
                    return
                endif
            endif
        
            ! FINISHED THE MIDDLE,  ARE WE IN THE UPPER TAIL?
        
            if(ipart == 2) then
                if(utail == 1) then
                    powr = 1.0
                    if(zval > UNEST) then
                        cdfval = powint(ccut(nccut),zmax,ccdf(nccut), &
                        &                                   1.0,zval,powr)
                    else
                        zval   = powint(ccdf(nccut),1.0,ccut(nccut), &
                        zmax,cdfval,powr)
                    endif
        
                else if(utail == 2) then
                
                    ! Power interpolation to upper limit "utpar"?
                
                    if(zval > UNEST) then
                        cdfval = powint(ccut(nccut),zmax,ccdf(nccut), &
                        &                                   1.0,zval,utpar)
                    else
                        powr = 1.0 / utpar
                        zval   = powint(ccdf(nccut),1.0,ccut(nccut), &
                        zmax,cdfval,powr)
                    endif
                
                    ! Linear interpolation between the rescaled global cdf?
                
                else if(utail == 3) then
                    if(zval > UNEST) then
                    
                        ! Approximately Locate the point and the class bound:
                    
                        idat = locate(cut,ncut,1,ncut,zval)
                        ilow = locate(cut,ncut,1,ncut,ccut(nccut))
                        if(cut(idat) < zval)        idat = idat + 1
                        if(cut(ilow) < ccut(nccut)) ilow = ilow + 1
                    
                        ! Straight linear interpolation if no data; otherwise, local linear
                        ! interpolation:
                    
                        if(idat <= 0 .OR. idat >= ncut .OR. &
                        ilow <= 0 .OR. ilow >= ncut) then
                            cdfval = powint(ccut(nccut),zmax, &
                            ccdf(nccut),1.0,zval,1.)
                        else
                            temp   = powint(cut(idat),cut(idat+1), &
                            cdf(idat),cdf(idat+1),zval,1.)
                            cdfval = powint(cdf(ilow),1.0, &
                            ccdf(nccut),1.0,temp,1.)
                        endif
                    else
                    
                        ! Computing Z value: Are there any data out in the tail?
                    
                        ilow = locate(cut,ncut,1,ncut,ccut(nccut))
                        if(cut(ilow) < ccut(nccut)) ilow = ilow + 1
                    
                        ! Straight linear interpolation if no data; otherwise, local linear
                        ! interpolation:
                    
                        if(ilow <= 0 .OR. ilow >= ncut) then
                            zval   = powint(ccdf(nccut),1.0, &
                            ccut(nccut),zmax,cdfval,1.)
                        else
                            temp = powint(ccdf(nccut),1.0, &
                            cdf(ilow),1.0,cdfval,1.)
                            idat = locate(cdf,ncut,1,ncut,temp)
                            if(cut(idat) < ccut(nccut)) idat=idat+1
                            if(idat >= ncut) then
                                zval   = powint(ccdf(nccut),1.0, &
                                ccut(nccut),zmax,cdfval,1.)
                            else
                                zval = powint(cdf(idat),cdf(idat+1), &
                                cut(idat),cut(idat+1),temp,1.)
                            endif
                        endif
                    endif
                
                    ! Fit a Hyperbolic Distribution?
                
                else if(utail == 4) then
                
                    ! Figure out "lambda" and required info:
                
                    lambda = (ccut(nccut)**utpar)*(1.0-ccdf(nccut))
                    if(zval > UNEST) then
                        cdfval = 1.0 - (lambda/(zval**utpar))
                    else
                        zval = (lambda/(1.0-cdfval))**(1.0/utpar)
                    endif
                else
                
                    ! Error situation - unacceptable option:
                
                    ierr = 2
                    return
                endif
            endif
            if(zval < zmin) zval = zmin
            if(zval > zmax) zval = zmax
        
            ! All finished - return:
        
            return
        
        end subroutine beyond


        pure real*8 function sqdist(x1,y1,z1,x2,y2,z2,ind,nrotmat,rotmat)
            !-----------------------------------------------------------------------
            !    Squared Anisotropic Distance Calculation Given Matrix Indicator
            !    ***************************************************************
            !
            !  This routine calculates the anisotropic distance between two points
            !  given the coordinates of each point and a definition of the
            !  anisotropy.
            !
            ! INPUT VARIABLES:
            !
            !   x1,y1,z1         Coordinates of first point
            !   x2,y2,z2         Coordinates of second point
            !   ind              The rotation matrix to use
            !   nrotmat          The maximum number of rotation matrices dimensioned
            !   rotmat           The rotation matrices
            !
            ! OUTPUT VARIABLES:
            !   sqdist           The squared distance accounting for the anisotropy
            !                      and the rotation of coordinates (if any).
            !
            ! NO EXTERNAL REFERENCES
            !-----------------------------------------------------------------------

            ! input variables variables
            integer, intent(in) :: nrotmat, ind 
            real*8, intent(in), dimension(nrotmat,3,3) :: rotmat
            real, intent(in) :: x1,y1,z1, x2,y2, z2

            ! internal variables
            real*8 :: cont,dx,dy,dz
            integer :: i
            ! Compute component distance vectors and the squared distance:

            dx = dble(x1 - x2)
            dy = dble(y1 - y2)
            dz = dble(z1 - z2)
            sqdist = 0.0
            do i=1,3
                cont   = rotmat(ind,i,1) * dx &
                + rotmat(ind,i,2) * dy &
                + rotmat(ind,i,3) * dz
                sqdist = sqdist + cont * cont
            end do
            return
        end function

        pure real function cova3(x1,y1,z1,x2,y2,z2,ivarg,nst,c0,it,cc,aa, irot,nrotmat,rotmat, d)
            !-----------------------------------------------------------------------
            !                    Covariance Between Two Points
            !                    *****************************
            ! This subroutine calculated the covariance associated with a variogram
            ! model specified by a nugget effect and nested varigoram structures.
            ! The anisotropy definition can be different for each nested structure.
            !
            ! INPUT VARIABLES:
            !   x1,y1,z1         coordinates of first point
            !   x2,y2,z2         coordinates of second point
            !   nst              number of nested structures (maximum of 4)
            !   ivarg            variogram number (set to 1 unless doing cokriging
            !                       or indicator kriging)
            !   c0(ivarg)        isotropic nugget constant
            !   it(i)            type of each nested structure:
            !                      1. spherical model of range a;
            !                      2. exponential model of parameter a;
            !                           i.e. practical range is 3a
            !                      3. gaussian model of parameter a;
            !                           i.e. practical range is a*sqrt(3)
            !                      4. power model of power a (a must be gt. 0  and
            !                           lt. 2).  if linear model, a=1,c=slope.
            !                      5. hole effect model
            !   cc(i)            multiplicative factor of each nested structure.
            !                      (sill-c0) for spherical, exponential,and gaussian
            !                      slope for linear model.
            !   aa(i)            parameter "a" of each nested structure.
            !   irot             index of the rotation matrix for the first nested
            !                    structure (the second nested structure will use
            !                    irot+1, the third irot+2, and so on)
            !   nrotmat          size of rotation matrix arrays
            !   rotmat           rotation matrices
            !
            ! OUTPUT VARIABLES:
            !   cmax             maximum covariance
            !   cova             covariance between (x1,y1,z1) and (x2,y2,z2)
            !
            ! EXTERNAL REFERENCES: sqdist    computes anisotropic squared distance
            !                      rotmat    computes rotation matrix for distance
            !
            ! Change from old code:
            !    -nst is now a sclar instead of an array, and all the variograms will 
            !     have the same number of structures (nst). The sise of the arrays 
            !     is nst*ivarg
            !    -size of the rotation matrix renemed to nrotmat
            !    -removed output cmax, which is not used in any other gslib program
            !    -turned cova3 into function
            !    -parameter d for damped hole effect added
            !-----------------------------------------------------------------------
            
            ! input variables
            real, intent(in) :: x1,y1,z1,x2,y2,z2          ! coordinates of the two points
            integer, intent(in) :: ivarg, irot, nst, nrotmat                  
            real, intent(in), dimension(ivarg) ::c0
            real, intent(in), dimension(nst*ivarg) :: it,cc,aa
            real*8, intent(in), dimension(nrotmat,3,3) :: rotmat
            real, intent(in), optional :: d
            
            ! internal variables
            real*8 :: hsqd
            real :: cmax, h, hr
            integer :: istart, is, ist, ir
            
        
            ! Calculate the maximum covariance value (used for zero distances and
            ! for power model covariance):
        
            istart = 1 + (ivarg-1)*nrotmat
            cmax   = c0(ivarg)
            do is=1,nst
                ist = istart + is - 1
                if(it(ist) == 4) then
                    cmax = cmax + PMX
                else
                    cmax = cmax + cc(ist)
                endif
            end do
        
            ! Check for "zero" distance, return with cmax if so:
        
            hsqd = sqdist(x1,y1,z1,x2,y2,z2,irot,nrotmat,rotmat)
            if(real(hsqd) < EPSLON) then
                cova3 = cmax
                return
            endif
        
            ! Loop over all the structures:
        
            cova3 = 0.0
            do is=1,nst
                ist = istart + is - 1
            
                ! Compute the appropriate distance:
                if(ist /= 1) then
                    ir = min((irot+is-1),nrotmat)
                    hsqd=sqdist(x1,y1,z1,x2,y2,z2,ir,nrotmat,rotmat)
                end if
                h = real(dsqrt(hsqd))
            
                ! Spherical Variogram Model?
                if(it(ist) == 1) then
                    hr = h/aa(ist)
                    if(hr < 1.) cova3=cova3+cc(ist)*(1.-hr*(1.5-.5*hr*hr))
                
                ! Exponential Variogram Model?
                else if(it(ist) == 2) then
                    cova3 = cova3 + cc(ist)*exp(-3.0*h/aa(ist))
                
                ! Gaussian Variogram Model?
                else if(it(ist) == 3) then
                    cova3 = cova3 + cc(ist)*exp(-3.*(h/aa(ist))*(h/aa(ist)))
                
                ! Power Variogram Model?
                else if(it(ist) == 4) then
                    cova3 = cova3 + cmax - cc(ist)*(h**aa(ist))
                
                ! Hole Effect Model?
                else if(it(ist) == 5) then
                    cova3 = cova3 + cc(ist)*cos(h/aa(ist)*PI)

                ! Damped Hole Effect Model?
                else if(it(ist) == 6) then

                    if(present(d)) then
                        cova3 = cova3 + cc(ist)*exp(-3.0*h/d)*cos(h/aa(ist)*PI)
                    else
                        ! d is the distance where 95% of the hole effect is dampened out. 
                        ! Here we assume d is 10 times the range.
                        cova3 = cova3 + cc(ist)*exp(-3.0*h/(10.0 * aa(ist)))*cos(h/aa(ist)*PI) 
                    end if

                ! TODO: add here Cauchy model, Matern model, Logistic model (rational quadratic model), Generalised Cauchy, Spheroidal
                
                end if

            end do
        
            ! Finished:
        
        end function cova3


        pure elemental integer function getindx(n,min,siz,loc) result (index)
            !-----------------------------------------------------------------------
            !     Gets the coordinate index location of a point within a grid
            !     ***********************************************************
            ! n       number of "nodes" or "cells" in this coordinate direction
            ! min     origin at the center of the first cell
            ! siz     size of the cells
            ! loc     location of the point being considered
            ! index   output index within [1,n]
            !-----------------------------------------------------------------------

            ! inputs 
            integer, intent(in) :: n
            real, intent (in) :: min,siz,loc
            
            ! Compute the index of "loc":
            index = int( (loc-min)/siz + 1.5 )
            
        end function getindx
            

        subroutine setrot(ang1,ang2,ang3,anis1,anis2,ind,nrotmat,rotmat)
            !-----------------------------------------------------------------------
            !              Sets up an Anisotropic Rotation Matrix
            !              **************************************
            !
            ! Sets up the matrix to transform cartesian coordinates to coordinates
            ! accounting for angles and anisotropy (see manual for a detailed
            ! definition):
            !
            ! INPUT PARAMETERS:
            !   ang1             Azimuth angle for principal direction
            !   ang2             Dip angle for principal direction
            !   ang3             Third rotation angle
            !   anis1            First anisotropy ratio
            !   anis2            Second anisotropy ratio
            !   ind              matrix indicator to initialize
            !   nrotmat          maximum number of rotation matrices dimensioned
            !   rotmat           rotation matrices
            !-----------------------------------------------------------------------

            ! inputs
            integer, intent (in) :: nrotmat, ind
            real, intent (in) :: ang1,ang2,ang3,anis1,anis2

            !output
            real*8, intent (inout), dimension(nrotmat,3,3) :: rotmat

            !internal
            real*8 :: afac1,afac2,sina,sinb,sint, cosa,cosb,cost
            real :: alpha, beta, theta
            
            ! Converts the input angles to three angles which make more
            !  mathematical sense:
            
            !         alpha   angle between the major axis of anisotropy and the
            !                 E-W axis. Note: Counter clockwise is positive.
            !         beta    angle between major axis and the horizontal plane.
            !                 (The dip of the ellipsoid measured positive down)
            !         theta   Angle of rotation of minor axis about the major axis
            !                 of the ellipsoid.
            
            if(ang1 >= 0.0 .AND. ang1 < 270.0) then
                alpha = (90.0   - ang1) * DEG2RAD
            else
                alpha = (450.0  - ang1) * DEG2RAD
            endif
            beta  = -1.0 * ang2 * DEG2RAD
            theta =        ang3 * DEG2RAD
            
            ! Get the required sines and cosines:
            
            sina  = dble(sin(alpha))
            sinb  = dble(sin(beta))
            sint  = dble(sin(theta))
            cosa  = dble(cos(alpha))
            cosb  = dble(cos(beta))
            cost  = dble(cos(theta))
            
            ! Construct the rotation matrix in the required memory:
            
            afac1 = 1.0 / dble(max(anis1,EPSLON))
            afac2 = 1.0 / dble(max(anis2,EPSLON))
            rotmat(ind,1,1) =       (cosb * cosa)
            rotmat(ind,1,2) =       (cosb * sina)
            rotmat(ind,1,3) =       (-sinb)
            rotmat(ind,2,1) = afac1*(-cost*sina + sint*sinb*cosa)
            rotmat(ind,2,2) = afac1*(cost*cosa + sint*sinb*sina)
            rotmat(ind,2,3) = afac1*( sint * cosb)
            rotmat(ind,3,1) = afac2*(sint*sina + cost*sinb*cosa)
            rotmat(ind,3,2) = afac2*(-sint*cosa + cost*sinb*sina)
            rotmat(ind,3,3) = afac2*(cost * cosb)
        
        end subroutine setrot

        subroutine ordrel(ivtype,ncut,ccdf,ccdfo,nviol,aviol,xviol)
            !-----------------------------------------------------------------------
            !                 Correct Order Relation Problems
            !                 *******************************
            ! This subroutine identifies and corrects order relation problems in a
            ! conditional distribution known at a specified number of cutoffs.
            !
            ! INPUT VARIABLES:
            !   ivtype           variable type (0=categorical, 1=continuous)
            !   ncut             number of cutoffs
            !   ccdf(i)          input ccdf values
            !
            ! OUTPUT VARIABLES:
            !   ccdfo            corrected ccdf values
            !   nviol()          number of order relation violations
            !   aviol()          average magnitude of the order relation violations
            !   xviol()          maximum magnitude of the order relation violations
            !
            ! PROGRAMMING NOTES:
            !
            !   1. the arrays ccdf1 and ccdf2 are used for temporary storage of the
            !      ccdf corrected sequentially upwards and downwards.  The program
            !      execution will be stopped if the memory allocation of these two
            !      arrays is not sufficient.
            !-----------------------------------------------------------------------

            ! inputs
            integer, intent(in) :: ncut, ivtype 
            real, intent(in), dimension(ncut) :: ccdf
            
            ! outputs
            real, intent(out), dimension(ncut) :: ccdfo, aviol, xviol
            integer, intent(out), dimension(ncut) :: nviol
            
            ! internals 
            real, dimension(ncut) :: ccdf1(ncut), ccdf2(ncut)
            integer :: i
            real :: sumcdf, viol
        
            ! Make sure conditional cdf is within [0,1]:
            
            do i=1,ncut
                if(ccdf(i) < 0.0) then
                    ccdf1(i) = 0.0
                    ccdf2(i) = 0.0
                else if(ccdf(i) > 1.0) then
                    ccdf1(i) = 1.0
                    ccdf2(i) = 1.0
                else
                    ccdf1(i) = ccdf(i)
                    ccdf2(i) = ccdf(i)
                endif
            end do
            
            ! Correct sequentially up, then down, and then average:
            
            if(ivtype == 0) then
                sumcdf = 0.0
                do i=1,ncut
                    sumcdf = sumcdf + ccdf1(i)
                end do
                if(sumcdf <= 0.0) sumcdf = 1.0
                do i=1,ncut
                    ccdfo(i) = ccdf1(i) / sumcdf
                end do
            else
                do i=2,ncut
                    if(ccdf1(i) < ccdf1(i-1)) ccdf1(i) = ccdf1(i-1)
                end do
                do i=ncut-1,1,-1
                    if(ccdf2(i) > ccdf2(i+1)) ccdf2(i) = ccdf2(i+1)
                end do
                do i=1,ncut
                    ccdfo(i) = 0.5*(ccdf1(i)+ccdf2(i))
                end do
            end if
            
            ! Accumulate error statistics:
            
            do i=1,ncut
                if(ccdf(i) /= ccdfo(i)) then
                    viol = abs(ccdf(i)-ccdfo(i))
                    nviol(i) = nviol(i) + 1
                    aviol(i) = aviol(i) + viol
                    xviol(i) = max(xviol(i),viol)
                endif
            end do
            
        end subroutine ordrel

end module gslib