! compile using the command: gfortran -static unravel.f90 -o unravel


program main
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!                                                                      %
! OpenGeostat Consulting (C) 2018                                      %
! Adrian Martinez Vargas                                               %
! License: MIT, see LICENSE for more information                       %
! This program was created modifying original GSLIB code               %
!                                                                      %
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!                                                                      %
! Copyright (C) 1996, The Board of Trustees of the Leland Stanford     %
! Junior University.  All rights reserved.                             %
!                                                                      %
! The programs in GSLIB are distributed in the hope that they will be  %
! useful, but WITHOUT ANY WARRANTY.  No author or distributor accepts  %
! responsibility to anyone for the consequences of using them or for   %
! whether they serve any particular purpose or work at all, unless he  %
! says so in writing.  Everyone is granted permission to copy, modify  %
! and redistribute the programs in GSLIB, but only under the condition %
! that this notice and the above copyright notice remain intact.       %
!                                                                      %
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-----------------------------------------------------------------------
!
!   Unravel multiple realizations in a single file
!
!-----------------------------------------------------------------------
    parameter (MAXLEN=132,MAXFIL=40,VERSION=1.000)
    
    parameter (maxnvar = 300)
    
    character str*132,datafl*500,outfl*500,fmt*36,label*20
    logical ::   testfl
    data      lin/1/,lout/2/
    real      var(500)
    integer ::   varnum, nrealiz, nrows, nvari
    integer :: test, error
    real, allocatable :: table(:,:), unraveltb(:,:)
    
    ! header parameters
       
    character(len=500) :: comment_line, tttx
    character(len=80), dimension (maxnvar) ::  varnames   
    integer  ::  nvar

    ! data parameters
    
    integer :: maxdat

    
    ! to format some integer printing 
    character(len=8) :: i_char
    
    
    ! Note VERSION number:
    write(*,*) '  '
    write(*,*) ' ***************************************** '
    write(*,*) '    OpenGeostat Consulting 2018            '
    write(*,*) ' ***************************************** '
    write(*,9999) VERSION
    9999 format(/' COLUMN Version: ',f5.3/)

    ! Get the name of the parameter file - try the default name if no input:

    write(*,*) 'Which parameter file do you want to use?'
    read (*,'(a20)') str(1:20)
    if(str(1:1) == ' ') str(1:20) = 'unravel.par        '
    inquire(file=str(1:20),exist=testfl)
    if( .NOT. testfl) then
        write(*,*) 'ERROR - the parameter file does not exist,'
        write(*,*) '        check for the file and try again  '
        write(*,*)
        stop
    endif
    open(lin,file=str(1:20),status='OLD')

    ! Find Start of Parameters:

    1 read(lin,'(a4)',end=98) str(1:4)
    if(str(1:4) /= 'STAR') go to 1

    ! Read Input Parameters:

    read(lin,'(a40)',err=98) datafl
    call chknam(datafl,MAXFIL)
    write(*,*) 'Data File:                      ',datafl

    read(lin,'(a40)',err=98) outfl
    call chknam(outfl,MAXFIL)
    write(*,*) 'Output File:                    ',outfl

    read(lin,*,err=98) varnum
    write(*,*) 'Variable to unravel:            ', varnum    
    
    read(lin,*,err=98) nrealiz
    write(*,*) 'Number or realization:          ', nrealiz
    
    read(lin,*,err=98) nrows
    write(*,*) 'Number of data rows in the file:', nrows 

    write(*,*)
    close(lin)

    ! get header parameters : comment_line, varnames and nvar 
    call read_header(datafl, comment_line, varnames, nvar, error, maxnvar)
    if (error/=0) then 
        Write (*,*) 'Error reading file header', error
        stop
    end if
    write (*,*) 'Comment_line:          ', trim(comment_line)
    write (*,*) 'Nvar:                  ', nvar
    write (*,*) 'Varnames:              '
    do i=1,nvar
        write (*,*) i, varnames(i)
    end do 
    
    ! get number of rows in the file
    call read_ndata(datafl, nvar, maxdat, error)
    if (error/=0) then 
        Write (*,*) 'Error obtaining number of raws', error
        stop
    end if
    write (*,*) 'Number of raws in file:', maxdat
    
    ! load the data in memory 
    ! first allocate memory
    allocate (table(maxdat,nvar),stat = test)
    if (test /= 0) then
        write(*,*) 'Error 1: Allocation failed due to ','insufficient memory!', test
        stop
    end if    
    ! then read the data
    call read_data(datafl, nvar, maxdat, table, error)
    if (error/=0) then 
        Write (*,*) 'Error reading data table', error
        stop
    end if
    ! print 10 first lines       
    do i=1,min(maxdat, 10)
        write (*,*) ( table(i,j), j=1,nvar )
    end do 
    write (*,*) 'End of table or first 10 rows'
  
    ! Finally unravel the data
    allocate (unraveltb(nrows,nrealiz),stat = test)
    if (test /= 0) then
        write(*,*) 'Error 2: Allocation failed due to ','insufficient memory!', test
        stop
    end if
    call unravelsim(table, maxdat, nvar, unraveltb, varnum, nrows, nrealiz, error)
    if (error/=0) then 
        if (error /= 1) then 
            Write (*,*) 'Error unstaking simulation', error
            stop
        else 
            write(*,*) 'Error: the  nrows*nrealiz<>maxdat', error
            stop
        end if
    end if
    
    ! and write results in the output file 
    open(unit=lout,file=outfl, status='UNKNOWN')
    write(lout,'(A)') trim(comment_line)
    write(tttx,'(I8)') nrealiz
    write(lout,'(A)') trim(tttx) 
    !write variable names
    do i=1,nrealiz
        write (i_char, '(i8)') i
        write(lout,'(A,A)') 'Sim_' , trim(adjustl(i_char))        
    end do 
    do i=1,nrows
        write (lout,*) ( unraveltb(i,j), j=1,nrealiz )
    end do
    
    ! close files and stop the program

    close(lin)
    close(lout)
    write(*,9998) VERSION
    9998 format(/' COLUMN Version: ',f5.3, ' Finished'/)
    stop
    98 stop ' ERROR in parameter file'
    99 stop ' ERROR in data file'


END PROGRAM


subroutine unravelsim(table, maxdat, nvar, unraveltb, varnum, nrows, nrealiz, error)
    ! takes a column from a table and unravel realizations
    integer,   intent(out) ::   error
    integer , intent(in):: maxdat, nvar, nrows, nrealiz, varnum
    real, intent(in) :: table(maxdat,nvar)

    real, intent(out) :: unraveltb (nrows, nrealiz)
    
    integer :: i,j,k
    
    error = 0
    
    ! check we have enough rows
    if (nrows*nrealiz/=maxdat) then 
        error = 1
        return
    end if 
    
    do i = 1, nrows
        do j=1, nrealiz
            k = nrows*(j-1)+i
            unraveltb (i,j) = table(k,varnum)
        end do
    end do
    
    
end subroutine


! This subroutines are from PyGSLIB code: read_gslib.f90

!*********************************************************************************
!     Subroutines for GSLIB datafile 
!*********************************************************************************
subroutine read_header(datafl, comment_line, varnames, nvar, error, maxnvar)
    ! read gslib data file header
    ! returns comment line, number of variables and variable names
    ! use this with the functions read_ndata (to get the number of lines in the file)
    ! and read_data to get an array with data
    ! warning: all values are converted to real and character variables are not admissible
    

    ! parameters:
    !--------------------------------
    ! input parameters:  
    !    datafl :  file path 
    !    maxnvar:  maximum number of variables (use a large number, ex. 500 or 500)
    !
    ! output: 
    !    comment_line : string with comments in the datafile
    !    varnames     : string array with varnames    
    !    nvar         : integer with number of variables in the file
    !    error        : integer with error code
    

    implicit none

    integer, intent(in) :: maxnvar    
    character(len=500), intent(in)  :: datafl
    character(len=500), intent(out) :: comment_line
    character(len=80), intent(out) , dimension (maxnvar) ::  varnames   
            
    integer,   intent(out)  ::  nvar
    
    character(len=80) str

    !error variable
    integer,   intent(out) ::   error
    
    integer i, j, lin

    error=0
    lin=1
    
    ! open the file and read in the header 
    open(unit=lin,file=datafl(1:500), status='old')

    read(lin,'(a500)',err=99) comment_line(1:500)
    read(lin,*,err=99)       nvar
    
    if (nvar > maxnvar) then
        goto 97
    end if
    
    do i=1,nvar
        read(lin,'(a80)',err=99) str
        varnames(i)= str
    end do

    close(lin)
    
    return
    
    97 error=97 !'nvar > maxnvar'
    close(lin)
    return
    
    99 error=99 !'error in data file!'
    close(lin)
    return
    
end subroutine read_header


subroutine read_ndata(datafl, nvar, maxdat, error)
    ! find out number of rows in a gslib file
    ! use this with the functions read_header (to get nvar and varnames)
    ! and read_data to get an array with data
    ! warning: all values are converted to real and character variables are not admissible

    ! parameters:
    !--------------------------------
    ! input parameters:  
    !    datafl :  file path 
    !    nvar   :  integer with number of variables in the file
    !
    ! output: 
    !    maxdat : integer with number of rows in the datafile. 
    !    error  : integer with error code


    implicit none
  
    character(len=250), intent(in)  :: datafl
    integer,  intent(in)  :: nvar
    integer, intent(out) ::  maxdat
    real*8, dimension (nvar) ::  var       
    
    !error variable
    integer,   intent(out) ::   error
    
    integer i, j, lin

    error=0
    lin=1
    
    ! Open the file and read number of data points 
    open(lin,file=datafl,status='OLD')
    read(lin,*,err=99) !comment line 
    read(lin,*,err=99) !nvar 
    do i=1,nvar        !varnames 
        read(lin,*)
    end do
    
    !initialize number of data to zero
    maxdat = 0
    2 read(lin,*,end=4,err=99) (var(j),j=1,nvar)
        maxdat = maxdat + 1
        go to 2
    4 continue
    
    close(lin)
    
    return
       
    99 error=99 !'ERROR in data file!'
    close(lin)
    return
    
end subroutine read_ndata


subroutine read_data(datafl, nvar, maxdat, table, error)
    ! read data in gslib file
    ! returns array with data (REAL format)
    ! use this with the functions read_header (to get nvar and varnames)
    ! and read_ndata to get the number of rows
    ! warning: all values are converted to real and character variables are not admissible


    ! parameters:
    !--------------------------------
    ! input parameters:  
    !    datafl :  file path 
    !    nvar   :  integer with number of variables in the file
    !    maxdat :  integer with number of rows in the datafile. 
    !
    ! output: 
    !    table    :  real array (maxdat,nvar) with data
    !    varnames :  string array with varnames    
    !    error        : integer with error code


    IMPLICIT NONE
  
    character(len=250), intent(in)  :: datafl
    integer,  intent(in)  :: nvar, maxdat
    real, intent(out), dimension (maxdat,nvar) ::  table
    real, dimension (nvar) ::  var       
    
    !error variable
    integer,   intent(out) ::   error
    
    integer i, j, lin

    error=0
    lin=1
    
    ! Open the file and read number of data points 
    open(lin,file=datafl,status='OLD')
    read(lin,*,err=99) !comment line 
    read(lin,*,err=99) !nvar 
    do i=1,nvar        !varnames 
        read(lin,*)
    end do


    ! Now, read the data
    do i=1, maxdat
        read(lin,*,end=94,err=99) (var(j),j=1,nvar)
        table(i,:)=var
    end do
    

    close(lin)
    
    return

    94  error=94 !'unexpected end of line'
    close(lin)
    return
    
    99 error=99 !'ERROR in data file!'
    close(lin)
    return
    
end subroutine read_data


subroutine chknam(str,len)
!-----------------------------------------------------------------------

!                   Check for a Valid File Name
!                   ***************************

! This subroutine takes the character string "str" of length "len" and
! removes all leading blanks and blanks out all characters after the
! first blank found in the string (leading blanks are removed first).



!-----------------------------------------------------------------------
    parameter (MAXLEN=132)
    character str(MAXLEN)*1

! Remove leading blanks:

    do i=1,len-1
        if(str(i) /= ' ') then
            if(i == 1) go to 1
            do j=1,len-i+1
                k = j + i - 1
                str(j) = str(k)
            end do
            do j=len,len-i+2,-1
                str(j) = ' '
            end do
            go to 1
        end if
    end do
    1 continue

! Find first blank and blank out the remaining characters:

    do i=1,len-1
        if(str(i) == ' ') then
            do j=i+1,len
                str(j) = ' '
            end do
            go to 2
        end if
    end do
    2 continue

! Return with modified file name:

    return
end subroutine chknam
