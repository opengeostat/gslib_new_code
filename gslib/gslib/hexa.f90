    
!-----------------------------------------------------------------------

!        Return the Hexadecimal Representation of a Number
!        *************************************************



!-----------------------------------------------------------------------
character(2) function hexa(number)
    character(1) :: hex(16)
    integer ::     number,digit1,digit2
    data hex    /'0','1','2','3','4','5','6','7','8', &
    '9','A','B','C','D','E','F'/

    if(number > 255) number = 255
    if(number < 1)   number =   1
          
    digit1 = number/16
    digit2 = number-16*digit1

    hexa=hex(digit1+1)//hex(digit2+1)

    return
end function
