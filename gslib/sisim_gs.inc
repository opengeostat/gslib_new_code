!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!                                                                      %
! Copyright (C) 2003, Statios Software and Services Incorporated.  All %
! rights reserved.                                                     %
!                                                                      %
! This program has been modified from the one distributed in 1996 (see %
! below).  This version is also distributed in the hope that it will   %
! be useful, but WITHOUT ANY WARRANTY. Compiled programs based on this %
! code may be redistributed without restriction; however, this code is %
! for one developer only. Each developer or user of this source code   %
! must purchase a separate copy from Statios.                          %
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
! Parameters:
!
    parameter(MAXNST=4,EPSLON=1.0e-20,UNEST=-99.0,VERSION=3.000)
!
! ACORN parameters:
!
    parameter(KORDEI=12,MAXOP1=KORDEI+1,MAXINT=2**30)
!
! Static Array Dimensioning:
!
    integer :: inoct(8),ltail,middle,utail,sstrat,mults,nmult,ivtype, &
    test
    real ::    ltpar,mpar,utpar
!
! Common Blocks:
!
    common /iaco/   ixv(MAXOP1)
    
    common /grid3d/ nx,ny,nz,nxy,nxyz,xsiz,ysiz,zsiz,xmn,ymn,zmn
    
    common /data3d/ ivtype,nd,zmin,zmax,ltail,middle,test, &
    utail,ltpar,mpar,utpar,tmin,tmax
    
    common /simula/ nsim,ktype,lin,lout,idbg,ldbg,ncutv,tmpp
    
    common /srch3d/ radius,radsqd,sang1,sang2,sang3,sanis1,sanis2, &
    noct,nclose,ndmax,mik,sstrat,infoct,inoct, &
    mults,nmult,maxsec,nodmax
    
    common /cova3d/ ncutz,ng,isrot,ivrot,cbb
    
    common /clooku/ nctx,ncty,nctz,nlooku,ncnode
