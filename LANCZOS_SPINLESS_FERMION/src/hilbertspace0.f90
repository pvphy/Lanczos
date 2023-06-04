module hilbertspace

contains
 include 'src/multiply.f90'

 !=================================================================================================
 subroutine remove_2particle( vecin, Nin, basisin,hs_size_in,&
                            vecout,Nout,basisout,hs_size_out, &
                            site1,site2,spin,nsites)
 use cluster, only: apply_op, annhil
 implicit none
 logical iszero
 integer(KIND=8) state, newstate, op
 integer(KIND=8) nn,i,j,k,Nin,hs_size_in,hsupin,hsdnin,spin,nsites, m
 integer(KIND=8) Nout,hs_size_out,hsupout,hsdnout,site1,site2
 double complex, dimension(1:Nin) :: vecin
 double complex, dimension(1:Nout) :: vecout
 double precision fsgn
 integer(KIND=8), dimension(1:hs_size_in) :: basisin
 integer(KIND=8), dimension(1:hs_size_out) :: basisout
 print*,site1,site2,nsites,hs_size_in,hs_size_out
 vecout(:) = 0.0d0
 do nn = 1,hs_size_in
  state = basisin(nn)
  fsgn = 1.0d0
  iszero = .false.
  newstate = state
  call apply_op(newstate,iszero,site2,spin,nsites,fsgn,annhil)
  call apply_op(newstate,iszero,site1,spin,nsites,fsgn,annhil)
  if(.not.iszero)then
   call binary_search(newstate,k,basisout,hs_size_out)
   vecout(k) = vecout(k) + fsgn*vecin(nn)
  endif
 enddo
 return
 stop
 end subroutine remove_2particle
 !=================================================================================================
 subroutine add_particle( vecin, Nin, basisin,hs_size_in,&
                            vecout,Nout,basisout,hs_size_out, &
                            site,spin,nsites)
 use cluster, only: apply_op, create
 implicit none
 logical iszero
 integer(KIND=8) state, newstate, op
 integer(KIND=8) nn,i,j,k,Nin,hs_size_in,hsupin,hsdnin,spin,nsites, m
 integer(KIND=8) Nout,hs_size_out,hsupout,hsdnout,site
 double complex, dimension(1:Nin) :: vecin
 double complex, dimension(1:Nout) :: vecout
 double precision fsgn
 integer(KIND=8), dimension(1:hs_size_in) :: basisin
 integer(KIND=8), dimension(1:hs_size_out) :: basisout
 character psiup*200, psidn*200
 vecout(:) = 0.0d0

 do nn = 1,hs_size_in
  state = basisin(nn)
  fsgn = 1.0d0
  iszero = .false.
  newstate = state
  call apply_op(newstate,iszero,site,spin,nsites,fsgn,create)
  if(.not.iszero.and..not.doubleocc(newstate,nsites))then
   call binary_search(newstate,k,basisout,hs_size_out)
   vecout(k) = vecout(k) + dcmplx(fsgn,0.0d0)*vecin(nn)
  endif
 enddo
 return
 end subroutine add_particle
 !=================================================================================================

subroutine construct_basis(nup,ndn,N, hs_size0, basis0)
  	use module_quick_sort
	implicit none
	integer(KIND=8) hs_size0
	integer(KIND=8), allocatable, dimension(:) :: basis0
	integer(KIND=8), parameter :: izero = 0

	integer(KIND=8) tmpsizeup, tmpsizedn
	double precision starttime, endtime
	integer(KIND=8) state, state1, idx, i, j, hs_size,nup, ndn, N, pow2N
	integer(KIND=8), allocatable, dimension(:) :: vup, vdn
	integer(KIND=8), allocatable, dimension(:) :: pup, pdn
	integer(KIND=8), allocatable, dimension(:) :: tmp
	character psiup*200, psidn*200
 	!First we do things for (nup, ndn)
 	!Determine how big the total hilbert space can be
	pow2N = 2**N
	
	!tmpsizeup = nchoosek(N,nup)!size_hilbert_space(nup,izero,N)
	tmpsizedn = nchoosek(N,nup)!size_hilbert_space(izero,ndn,N)

	!allocate(vup(1:nup));       
	allocate(vdn(1:ndn))

	!allocate(pup(1:tmpsizeup));
	allocate(pdn(1:tmpsizedn))


	allocate(tmp(1:tmpsizedn))

	call cpu_time(starttime)

	i = 1
	idx = 0
	!print*,pow2N,tmpsizedn

	call combinations(pdn,vdn,i,N,i,ndn,idx)
	hs_size0 = 0

  	do j = 1,tmpsizedn
   		state=0!pow2N*pdn(j)
		!print*,j,state
 		do i=1,N
			if(.not.btest(pdn(j), i-1)) then
	 			state = state+ (2**(i-1))
			endif
 		enddo
        
		 
    	hs_size0 = hs_size0 + 1
    	tmp(hs_size0) = state
		! call return_string_for_state(state,psiup,psidn,N)
		! print*,state,psiup,psidn
		!print*, "State",state, ": ", dec2bin(state, N)
  	enddo
	
	allocate(basis0(1:hs_size0))
	basis0(:) = tmp(1:hs_size0)
	call quicksort(basis0)

	!deallocate(vup)
	deallocate(vdn)
	! deallocate(pup)
	deallocate(pdn)
	deallocate(tmp)
	call cpu_time(endtime)
	print*, "Basis (nup,ndn) done at:", endtime-starttime
	

	return
end subroutine construct_basis

!=============================================================================
!This subroutine determines of a state is double occupied.
!=============================================================================
logical function doubleocc(state,N)
	implicit none
	integer(KIND=8) state, site, N
	doubleocc = .false.
	do site = 1,N
		if(btest(state,site-1).and.btest(state,site+N-1))then
			doubleocc = .true.
			return
		endif
	enddo
	return
 end function doubleocc

 !=============================================================================
 ! This routine builds a list of all possible combinations of numbers where
 ! a certain number of bits are set.
 !=============================================================================
 recursive subroutine combinations(array,v,start,n,k,maxk,idx)  !	call combinations(pdn,vdn,i,N,i,ndn,idx)
 implicit none
 integer(KIND=8) i, start, n, k, maxk, state, idx
 integer(KIND=8), dimension(1:maxk) :: V
 integer(KIND=8), intent(in out), dimension(:) :: array
 if(k.gt.maxk)then
  state = 0
  idx = idx + 1
  do i = 1,maxk
   state = ibset(state,v(i)-1)
  enddo
  array(idx) = state
  return
 endif
 do i = start, n
  v(k) = i;
  call combinations(array,v,i+1,n,k+1,maxk,idx);
 enddo
 return
 end subroutine combinations
 !=============================================================================
 ! Function size_hilbert_space
 ! This function takes the number of up and dn spins, along with the total
 ! number of sites, and returns the size of the hilbert space
 !=============================================================================
! integer(KIND=8) function size_hilbert_space(nup,ndn,nsites)
! implicit none
! integer(KIND=8) nup, ndn, nsites
! size_hilbert_space = nchoosek(nsites,nup)!*nchoosek(nsites,ndn)
! return
! end function size_hilbert_space
 !=============================================================================
 ! integer function nchoosek
 ! integer implemenation of a combination.
 !=============================================================================
 integer(KIND=8) function nchoosek(n,k)
 implicit none
 integer(KIND=8)  n, k
 integer(KIND=8) i, j, itmp2, itmp1, itmp3

 itmp1 = 1; j=1
do i = n,n-k+1,-1
  itmp1 = itmp1*i

	if(j.lt.(k+1)) then
	itmp1 = itmp1/j; j=j+1
	endif
	enddo
  nchoosek = itmp1

  return
  end function nchoosek

 !=============================================================================
 ! RECURSIVE SUBROUTINE FACTORIAL(n,fac)
 ! Recursive implementation of n factorial where n is an integer.
 !=============================================================================
 recursive subroutine factorial(n,fac)
 implicit none
 integer(KIND=8), intent(in) :: n
 integer(KIND=8), intent(out) :: fac
 integer(KIND=8) itmp
 if(n.ge.1)then
  call factorial(n-1,itmp)
  fac = n*itmp
 else
  fac = 1
 endif
 return
 end subroutine factorial

   !=================================================================================================
subroutine spinflip( vecin, Nin, basisin,hs_size_in,&
                            vecout,Nout,basisout,hs_size_out, &
                            site, N)
	use cluster, only: apply_op, create, annhil, spinup, spindn
	implicit none
	logical iszero
	integer(KIND=8) state, newstate, op, niup, nidn, njup, njdn
	integer(KIND=8) nn,i,j,k,Nin,hs_size_in,hsupin,hsdnin,spin,N, m
	integer(KIND=8) Nout,hs_size_out,hsupout,hsdnout,site, nsite
	double complex, dimension(1:Nin) :: vecin
	double complex, dimension(1:Nout) :: vecout
	double precision fsgn
	integer(KIND=8), dimension(1:hs_size_in) :: basisin
	integer(KIND=8), dimension(1:hs_size_out) :: basisout
	character psiup*200, psidn*200
	vecout(:) = 0.0d0

	do nn = 1,hs_size_in
	state = basisin(nn)
	do spin=0,1

	!do the Sz term
	niup = 0; nidn = 0;

	if(btest(state,site+spinup*N-1)) niup = niup + 1
	if(btest(state,site+spindn*N-1)) nidn = nidn + 1

	vecout(nn) = vecout(nn) + 0.5d0*dfloat(niup-nidn)*vecin(nn)  !*dfloat(njup-njdn)
	enddo
	enddo
	return
end subroutine spinflip
 !=================================================================================================
 !=================================================================================================
subroutine spinexch( vecin, Nin, basisin,hs_size_in,&
                            vecout,Nout,basisout,hs_size_out, &
                            site,nsite, N)
 	use cluster, only: apply_op, create, annhil, spinup, spindn
	implicit none
	logical iszero
	integer(KIND=8) state, newstate, op, niup, nidn, njup, njdn
	integer(KIND=8) nn,i,j,k,Nin,hs_size_in,hsupin,hsdnin,spin,N, m
	integer(KIND=8) Nout,hs_size_out,hsupout,hsdnout,site, nsite
	double complex, dimension(1:Nin) :: vecin
	double complex, dimension(1:Nout) :: vecout
	double precision fsgn
	integer(KIND=8), dimension(1:hs_size_in) :: basisin
	integer(KIND=8), dimension(1:hs_size_out) :: basisout
	character psiup*200, psidn*200
	vecout(:) = 0.0d0

	do nn = 1,hs_size_in
		state = basisin(nn)
		!do spin=0,1

   			newstate = state
			iszero = .false.
			fsgn = 1.0d0
			call apply_op(newstate,iszero,nsite,spindn,N,fsgn,annhil)
			call apply_op(newstate,iszero,site,spindn,N,fsgn,annhil)

			if(.not.iszero.and..not.doubleocc(newstate,N))then
			call binary_search(newstate,k,basisout,hs_size_out)
			vecout(k) = vecout(k) + dcmplx(fsgn,0.0d0)*vecin(nn);
			endif

			! newstate = state
			! iszero = .false.
			! fsgn = 1.0d0
		
			! call apply_op(newstate,iszero,nsite,spindn,N,fsgn,annhil)
			! call apply_op(newstate,iszero,site,spindn,N,fsgn,annhil)
		
			! if(.not.iszero.and..not.doubleocc(newstate,N))then
			! 	call binary_search(newstate,k,basisout,hs_size_out)
			! 	vecout(k) = vecout(k) + 0.5*dcmplx(fsgn,0.0d0)*vecin(nn)
			! endif


			vecout(nn) = vecout(nn) + vecin(nn)
 		!enddo
	enddo
	print*,hs_size_in,hs_size_out
	return
end subroutine spinexch
function dec2bin(decimal, digits) result(binary)
    implicit none
    integer(KIND=8), intent(in) :: decimal
    integer(KIND=8), intent(in) :: digits
    character(len=digits) :: binary
    integer(KIND=8) :: i, bit1

    binary = ''
    do i = digits, 1, -1
        IF(btest(decimal, i-1)) bit1=1
		IF(.not.btest(decimal, i-1)) bit1=0
        write(binary(i:i), '(I1)') bit1
    enddo

    return
end function dec2bin

end module hilbertspace
