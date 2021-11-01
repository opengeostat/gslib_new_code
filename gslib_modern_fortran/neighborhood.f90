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
            !                Naive search in sphere or radius r 
            !              ******************************************
            ! Transform sequence of integers generated by the function lcg into uniform distribution
            !
            ! INPUT VARIABLES:
            !   x                64 bit integers generated by function lcg
            !
            ! OUTPUT
            !   real             with uniform distribution in interval ]0,1[
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
