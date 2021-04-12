    subroutine blue(value,hexrep,bfrac)
!-----------------------------------------------------------------------

! Provided with a real value ``value'' this subroutine returns the blue
! portion of the color specification.

! Note common block "color" and call to "hexa"

!-----------------------------------------------------------------------
    real ::            value
    character       hexrep*2,hexa*2
    common /color/  cmin,cmax,cint(4),cscl
    hexrep = '00'
    if(value < cint(2))then
    
    ! Scale it between (255,255):
    
        integ  = 255
    else if((value >= cint(2)) .AND. (value < cint(3)))then
    
    ! Scale it between (255,0):
    
        integ = int((cint(3)-value)/(cint(3)-cint(2))*255.)
        if(integ > 255) integ = 255
        if(integ < 0)   integ = 0
    else if(value >= cint(3))then
    
    ! Scale it between (0,0):
    
        integ  = 0
    end if

! Establish coding and return:

    bfrac  = real(integ) / 255.
    hexrep = hexa(integ)
    return
    end subroutine blue
