    subroutine red(value,hexrep,rfrac)
!-----------------------------------------------------------------------

! Provided with a real value ``value'' this subroutine returns the red
! portion of the color specification.

! Note common block "color" and call to "hexa"

!-----------------------------------------------------------------------
    real ::            value
    character       hexrep*2,hexa*2
    common /color/  cmin,cmax,cint(4),cscl
    hexrep='00'
    if(value < cint(1))then
    
    ! Scale it between (y0,0):
    
        integ=int((cint(1)-value)/(cint(1)-cmin)*cscl)
        if(integ > 255) integ = 255
        if(integ < 0)   integ = 0
    else if((value >= cint(1)) .AND. (value < cint(3)))then
    
    ! Scale it between (0,0):
    
        integ  = 0
    else if((value >= cint(3)) .AND. (value < cint(4)))then
    
    ! Scale it between (0,255):
    
        integ = int((value-cint(3))/(cint(4)-cint(3))*255.)
        if(integ > 255) integ = 255
        if(integ < 0)   integ = 0
    else if(value >= cint(4))then
    
    ! Scale it between (255,255):
    
        integ  = 255
    end if

! Establish coding and return:

    rfrac  = real(integ) / 255.
    hexrep = hexa(integ)
    return
    end subroutine red
