module gslib

    real, parameter, private :: EPSLON=1.0e-20, UNEST=-1.0

    contains 
        
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
            implicit none

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

        end function backtr

        pure elemental real function gcum(x)
            !-----------------------------------------------------------------------
            ! Evaluate the standard normal cdf given a normal deviate x.  gcum is
            ! the area under a unit normal curve to the left of x.  The results are
            ! accurate only to about 5 decimal places.
            !-----------------------------------------------------------------------

            implicit none
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
            
            implicit none
            
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
            
            implicit none
            
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
                
                implicit none

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
            
end module gslib