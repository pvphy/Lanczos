!==============================================================================
! Module_quick_sort
! This module defines the quick sort.
!==============================================================================
module module_quick_sort
 public :: quicksort
 private :: partition
contains
 recursive subroutine quicksort(A)
 implicit none
 integer(KIND=8), intent(in out), dimension(:) :: A
 integer(KIND=8) :: iq
 if(size(A)>1)then
  call partition(A,iq)
  call quicksort(A(:iq-1))
  call quicksort(A(iq:))
 endif
 return
 end subroutine quicksort

 subroutine partition(A,marker)
 implicit none
 integer(KIND=8), intent(in out), dimension(:) :: A
 integer(KIND=8), intent(out) :: marker
 integer(KIND=8) :: i, j, temp, x
 x = A(1)
 i = 0
 j = size(A)+1
 do
  j = j -1
  do
   if(A(j) <=x) exit
   j = j-1
  enddo
  i = i + 1
  do
   if(A(i) >= x) exit
   i = i + 1
  enddo
  if(i<j)then
   temp = A(i)
   A(i) = A(j)
   A(j) = temp
  elseif(i==j)then
   marker = i + 1
   return
  else
   marker = i
   return
  endif
 enddo
 return
 end subroutine partition
end module module_quick_sort
