subroutine return_string_for_state(state,psiup,psidn,nsites)
implicit none
integer(KIND=8) state
integer(KIND=8) nsites, i
character psiup*200, psidn*200
do i = 1,200
 psiup(i:i) = ' '
 psidn(i:i) = ' '
enddo
psiup(1:1) = '|'
psidn(1:1) = '|'
psiup(nsites+2:nsites+2) = '>'
psidn(nsites+2:nsites+2) = '>'
do i = 1,nsites
 if(btest(state,i-1))then
  psidn(nsites+2-i:nsites+2-i) = '1'
 else
  psidn(nsites+2-i:nsites+2-i) = '_'
 endif
 if(btest(state,nsites+i-1))then
  psiup(nsites-i+2:nsites-i+2) = '1'
 else
  psiup(nsites-i+2:nsites-i+2) = '_'
 endif
enddo
return
end subroutine return_string_for_state

!==============================================================================
subroutine binary_search(state,idx,psi,hs_size)
implicit none
integer(KIND=8) state,hs_size,ii
integer(KIND=8) idx, mid, low, hi, maxiter, iter
integer(KIND=8), dimension(1:hs_size) :: psi
logical loop
maxiter = dlog(dfloat(hs_size))/dlog(dfloat(2))+10
low = 1
hi = hs_size
mid = (hi+low)/2
iter = 0
!did we get it right away?
if(state.eq.psi(low))then
 idx = low
 return
elseif(state.eq.psi(hi))then
 idx = hi
 return
endif

loop = .true.
do while(loop)
 iter = iter + 1
 if(state.eq.psi(mid))then
  loop = .false.
 elseif(state.lt.psi(mid))then
  hi=mid
  mid=(hi+low)/2
 else
  low = mid
  mid=(hi+low)/2
 endif
 if(iter.ge.maxiter)then
  print*, 'State not found in hilbert space after', maxiter, ' lookups.'
  print*, 'Looking for state ', state
  print*, low, mid, hi
  print*, psi(low), psi(mid), psi(hi)
  stop
 endif
enddo
idx = mid
return
end subroutine binary_search

