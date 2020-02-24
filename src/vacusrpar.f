!##############################################################################
! include vacusrpar - nul

! This file should contain the number of PROBLEM depedent equation parameters,
! the index names for them with values neqpar+1..neqpar+nspecialpar, 
! and a string giving the names for the file header. For example:
!
! INTEGER,PARAMETER:: mass_=neqpar+1, nspecialpar=1
! CHARACTER*4,PARAMETER:: specialparname='mass'
!
! By default there are no special parameters

INTEGER,PARAMETER:: grav0_=neqpar+1
INTEGER,PARAMETER:: grav1_=neqpar+1+1,grav2_=neqpar+2+1
INTEGER,PARAMETER:: nu_=neqpar+4
INTEGER,PARAMETER:: nspecialpar=2+2

!CHARACTER*2 ,PARAMETER:: specialparname='nu'


   CHARACTER*11,PARAMETER:: specialparname='grav1 grav2'



! end include vacusrpar - nul
!##############################################################################
