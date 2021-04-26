module gslib

    real, parameter, private :: EPSLON=1.0e-20

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
            
end module gslib