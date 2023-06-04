!==============================================================================
! Module random.  
! This module defines the random number generator.
!==============================================================================
module random
 double precision lastr
contains
 subroutine init_ranLM(rseed)
 implicit none
 double precision rseed
 lastr = rseed
 return
 end subroutine init_ranLM

 double precision function ranLM()
 implicit none
 double precision, parameter :: tiny = 1e-7
 double precision pi, r
 pi = 2.0d0*asin(1.0d0)
 ranLM = 4.0d0*lastr*(1.0d0-lastr)
 if(ranLM.eq.0.0d0) ranLM = tiny
 if(ranLM.ge.1.0d0) ranLM = 1.0d0 - tiny
 lastr = ranLM
 ranLM = 2.0d0*asin(dsqrt(ranLM))/pi
 return
 end function ranLM
end module random
