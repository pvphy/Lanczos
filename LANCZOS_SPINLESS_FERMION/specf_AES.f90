include 'src/quicksort.f90'
include 'src/random.f90'
include 'src/auxilary.f90'
include 'cluster.f90'
include 'src/hilbertspace0.f90'
include 'src/lanczos.f90'

program main
    use cluster
    use hilbertspace
    use module_lanczos

    implicit none
    include 'mpif.h'
    integer rank, msize, ierror, tag, status(MPI_STATUS_SIZE)
    integer startval, endval,flag_H
    logical, parameter :: verbose = .false.
    integer(KIND=8) nistates, nfstates 
    integer(KIND=8), parameter :: num = 1
    integer(KIND=8), parameter :: kk = 1
    integer(KIND=8), parameter :: term = 0
    integer(KIND=8), parameter :: numw = 500
    integer(KIND=8) ch, i, j, nq, spin, site, nsite, nw, nch
    integer(KIND=8)	ii, mm, niup, nidn, njup, njdn !!
    integer(KIND=8), parameter :: izero = 0
    integer(KIND=8) h0len, h1len, h2len, nstates, state, newstate
    integer(KIND=8) maxh0len, maxh1len, maxh2len, hs_size0,hs_size1
    integer(KIND=8), allocatable, dimension(:) :: H0i, H0j, basis0,basis1
    integer(KIND=8), allocatable, dimension(:) :: H1i, H1j
    integer(KIND=8), allocatable, dimension(:) :: H2i, H2j
    double precision, dimension(1:nfmax) :: EE
    double precision, dimension(1:num) :: Egg
    double precision, parameter :: wscale = 4.0d0
    double precision, allocatable, dimension(:) :: gs, cgs
    double precision, allocatable, dimension(:,:) :: States
    double complex, allocatable, dimension(:,:) :: psi0, psi00
    double precision, allocatable, dimension(:) :: Hel0
    double precision, allocatable, dimension(:) :: Hel1
    double precision, allocatable, dimension(:) :: Hel2
    double precision q, pi, Eg, win, arg, delta_fn, w, JJ, s, fsgn !!
    double precision, dimension(1:numw) :: AES
    double complex, allocatable, dimension(:) :: psi1, psi2, cvec, phi0
    double complex, allocatable, dimension(:) :: fvecs, fvecs2
    double complex ctmp
    double complex, parameter :: eye = (0.0d0,1.0d0)
    double complex, dimension(1:nfmax) :: final
    logical iszero0
    character filename*200, cmdlinearg*200
    double precision starttime, endtime,dw,wmax,wmin

    wmin=-35.0d0
    wmax=8.0d0

    pi = 2.0d0*asin(1.0d0)

    CALL getarg(1, cmdlinearg)
    read(unit=cmdlinearg,fmt=*) flag_H
    WRITE (*,*) 'flag for H = ', flag_H,'','IF FLAG_H=123...t-u...FLAG_H=345...t_delta_u'    



    ! Determine the size of the Hilbert Space
    call construct_basis(nup,ndn,N, hs_size0, basis0)

    nistates = min(50, hs_size0-2)
    nfstates = min(50, hs_size0-2)

    print*, 'Number of sites: ', N
    print*, 'Number particle: ', nup
    print*, 'Number holes:    ', ndn
    print*, ' '
    print*, 'Hamiltonian size for  problem = ', hs_size0


    H0len = 0
    maxH0Len = 100*Hs_size0
    allocate(H0i(1:maxH0Len))
    allocate(H0j(1:maxH0Len))
    allocate(Hel0(1:maxH0Len))
    !construct the hamiltonian for the (nup,ndn) problem.
    ch = 0
    ! read(*,*) flag_H
    ! print*,flag_H,'123=t_u,345=t_delta_u'
    ! print*, ' '
    print*, 'Constructing H for  problem.'

    if(flag_H.eq.123) call construct_hamiltonian(H0i,H0j,Hel0,H0len,maxH0len,basis0,hs_size0,ch)
    if(flag_H.eq.345) call construct_hamiltonian_delta_U(H0i,H0j,Hel0,H0len,maxH0len,basis0,hs_size0,ch)

    

    !*********  n+2*******************************
    print*, ' '
    print*, 'Constructing H for N+2 problem'
    print*, ' '
    print*, 'Number particle: ', nup-2
    print*, 'Number holes:    ', N-(nup-2)
    call construct_basis(nup-2,N-(nup-2),N, hs_size1, basis1)
    print*, 'Hilbert space for  N+2 problem = ', hs_size1

    H1len= 0
    maxH1Len= 100*Hs_size1
    allocate(H1i(1:maxH1Len))
    allocate(H1j(1:maxH1Len))
    allocate(Hel1(1:maxH1Len))
    if(flag_H.eq.123) call construct_hamiltonian(H1i,H1j,Hel1,H1len,maxH1len,basis1,hs_size1,ch)
    if(flag_H.eq.345) call construct_hamiltonian_delta_U(H1i,H1j,Hel1,H1len,maxH1len,basis1,hs_size1,ch)

    nstates = min(50,hs_size0)
    allocate(states(1:hs_size0,1:num)); Egg = 0.0d0; states = 0.0d0;
    allocate(gs(1:hs_size0))
    ORTHOGONALIZE = .true.
    print*, 'Calling get ground state.'
    call lanczos_get_ground_state(states,Egg,num,hs_size0,hs_size0,H0i,H0j,Hel0,H0len,maxH0len,nstates)
    print*, 'Initial state energies: '
    ORTHOGONALIZE = .false.




    do i = 1,num
    print*, i, Egg(i)
    enddo

    AES = 0.0d0


    Eg = Egg(1)
    Gs(:) = states(:,1)
    ! do i=1,hs_size0
    !     write(*,*) gs(i)  !ground state vector
    ! enddo
    call MPI_INIT(ierror)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, msize, ierror)
    call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierror)
    
    !loop over core-hole positions
    startval = 11 !(N*rank)/msize+1;
    endval   = 11 !(N*(rank+1))/msize;
    !print*,rank,msize,startval,endval

    allocate(psi0(1:hs_size1,startval:endval));
    !call MPI_Barrier( MPI_COMM_WORLD, ierror)
    
   
!! Evaluation of \sum_i e^{iq.R_i} S_i . S_J|g> and \sum_i e^{iq.R_i} S_i^z S_i . S_J|g>
    
    do ch =startval,endval
        psi0(:,ch) = dcmplx(0.0d0,0.0d0)

        site = ch
        
        print*, 'Applying ch ', ch

        allocate(psi1(1:hs_size0))
        allocate(psi2(1:hs_size1))

        psi2 = 0.0d0

        ! Initialize |g>
        psi1 = dcmplx(gs,0.d0)


        nsite=site+1
        call remove_2particle(psi1,hs_size0,basis0,hs_size0,psi2,hs_size1,basis1,hs_size1,site,nsite,spin,N)
        psi0(:,ch) = psi2

     
        !enddo
        


        deallocate(psi1)
        deallocate(psi2)

    enddo

    !allocate(psi00(1:hs_size1,1:N));

    !call MPI_Barrier( MPI_COMM_WORLD, ierror)
    !call mpi_gather(psi0(:,startval:endval), size(psi0(:,startval:endval)), MPI_DOUBLE_COMPLEX, psi00, N*hs_size0/msize, &
     !               MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD, ierror)

    !call MPI_Barrier( MPI_COMM_WORLD, ierror)
    

    
    ! Evalution of factor e^{iq.R_i}
    if(rank.eq.0) then

     
        ch=startval

        allocate(phi0(1:hs_size1))
        phi0(:) = psi0(:,ch) 
        deallocate(psi0)
        !project onto the final state manifold and add the result to the AES spectra
        call project_on_final(phi0,hs_size1,EE,Final,nfstates,H1i,H1j,Hel1,H1len,maxH1len,j)

        deallocate(phi0)
        write(unit=filename,fmt="('2hspecf_L_',i2.2,'_N_',i2.2,'_U_',i2.2,'.dat')") N,nup,int(U)
        open(file=filename,unit=60,action='write')
        w=wmin
        dw=(wmax-wmin)/dble(numw)
        do nw = 1,numw
            w = w+dw
            do i = 2,j
                ctmp = final(i)
                arg = w + Eg - EE(i)
                delta_fn = ((sig/pi)/((arg)**2+((sig)**2))) !exp(-arg*arg/(2.0d0*sig*sig))/(sig*sqrt(2*pi)) 
               
                !if(abs(Eg-EE(i)).gt.0.0001d0)then
                    AES(nw) = AES(nw) + delta_fn*dreal(ctmp*conjg(ctmp))
                !endif
            enddo
            write(60,*)w, AES(nw)
        enddo

           


 
        
        close(unit=60)
    endif

    call MPI_Barrier( MPI_COMM_WORLD, ierror)
    call MPI_FINALIZE(ierror)
  
    stop

end program main


! if(rank.eq.0) then

!     do nq = 0,N/2
!         q = dfloat(nq)*pi/dfloat(N/2)

!         allocate(phi0(1:hs_size0))
!         phi0(:) = 0.0d0

!         do ch=1,N
!             if(mod(ch,3).ne.0)phi0 = phi0 + psi00(:,ch)*exp(eye*q*dfloat(ch))/sqrt(dfloat(N))
!         enddo

!         !project onto the final state manifold and add the result to the RIXS spectra
!         call project_on_final(phi0,hs_size0,EE,Final,nfstates,H0i,H0j,Hel0,H0len,maxH0len,j)

!         deallocate(phi0)
!         do nw = 0,numw-1
!             w = wscale*dfloat(nw)/dfloat(numw)
!             do i = 2,j
!             ctmp = final(i)
!             arg = w + Eg - EE(i)
!             delta =  exp(-arg*arg/(2.0d0*sig*sig))/(sig*sqrt(2*pi)) !(sig*0.5d0/pi)/(arg*arg+0.25d0*sig*sig) !
!             !RIXS( nq,nw) = RIXS( nq,nw) + delta*dreal(ctmp*conjg(ctmp))
!             if(abs(Eg-EE(i)).gt.0.0001d0)then
!             RIXS( nq,nw) = RIXS( nq,nw) + delta*dreal(ctmp*conjg(ctmp))
!             endif
!             enddo
!         enddo
!         RIXS(-nq,:) = RIXS(nq,:)
!         !RIXS2(-nq,:) = RIXS2(nq,:)
!     enddo

!     write(unit=filename,fmt="('Correlation_',i1.1,'CLen',i2.2,'Order',i1.1,'J2_',f4.2,'.dat')") channel, N, kk, J2
!     open(file=filename,unit=60,action='write')
!     do nq = -N/2,N/2
!         do nw = 0,numw-1
!             w = wscale*dfloat(nw)/dfloat(numw)
!             write(60,*) dfloat(nq)/dfloat(N/2), w, RIXS(nq,nw)
!         enddo
!     enddo
!     close(unit=60)
! endif