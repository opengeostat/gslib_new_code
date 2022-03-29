module neighborhood

    implicit none


    contains

        ! todo: add kdtree 
        !       add efficient structure to search in grid (may be required in simulations)
        !           see for example https://link.springer.com/article/10.1007/s41060-020-00208-2
        !           consider this one: https://github.com/dongli/fortran-octree
        !           see also http://www.open3d.org/docs/release/index.html
        !           it looks like octree is the way to go 
        !         

        pure elemental logical function ball(x, y, z, x0, y0, z0, r)
            !-----------------------------------------------------------------------
            !                Naive search in sphere of radius r 
            !              ******************************************
            ! Finds points within a distance r from a target poin 
            !
            ! INPUT VARIABLES:
            !   x,y,z           arrays of reals or reals, coordinates of data points in 
            !                   isotropic space
            !   x0,y0,z0        reals (scalar), coordinates of the target point
            !   r               real,  isotropic search radius
            !
            ! OUTPUT:
            !   ball            array of booleans, True if the data point is 
            !                   within the isotropic distance r from the target 
            !                   point
            !  NOTES:
            !                   This is an elemental function and all inputs 
            !                   can be real or scalar. Users may use this function carefully
            !                   to avoid unexpected results. For example, 
            !                      
            !                      real :: x(5) = [1,2,3,4,5], y(5) = [1,2,3,4,5], z(5) = [0,0,0,0,0]
            !                      real :: x0 = 1, y0 = 1, z0 =1
            !                      real :: r =2
            !                      
            !                      ! legal, test each one of the five data points if is within r distance from a target point 
            !                      print *,  ball (  x,  y,  z, x0, y0, z0, r)   ! result  T T F F F
            !                      ! legal, test if one data point is within a distance r from five target points
            !                      print *,  ball ( x0, y0, z0,  x,  y,  z, r)   ! result  T T F F F  
            !                      ! legal, but arbitrary calls
            !                      print *,  ball (x, y, z, x, y, z, r)          ! result  T T T T T 
            !                      print *,  ball (x, y, z, x, y, z, x)          ! result  T T T T T 
            !                                                
            !-----------------------------------------------------------------------

            ! inputs
            real, intent(in) :: x, y, z, x0, y0, z0, r
            
            if ((x-x0)**2 + (y-y0)**2 + (z-z0)**2<=r**2) then
                ball = .True.
            else
                ball = .False.
            end if

        end function ball



end module neighborhood
