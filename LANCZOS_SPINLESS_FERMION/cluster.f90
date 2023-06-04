module cluster
    integer(KIND=8), parameter :: create = 1
    integer(KIND=8), parameter :: annhil = 0
    integer(KIND=8), parameter :: spinup = 1
    integer(KIND=8), parameter :: spindn = 0
    integer(KIND=8), parameter :: N = 24
    integer(KIND=8), parameter :: nup = 18
    integer(KIND=8), parameter :: ndn = N-nup
    double complex, parameter :: igam = (0.0d0,0.05d0)
    double precision, parameter :: t =-1.00d0
    double precision, parameter :: U = 8.00d0
    double precision, parameter :: delta=10.00d0
    double precision, parameter :: J1 = 1.00d0
    double precision, parameter :: J2 = 1.00d0*0.18
    double precision, parameter :: J3 = 1.00d0*0.18
    double precision, parameter :: sig = 0.10d0
    contains

    logical function doubleocc(state,Nsites)
        implicit none
        integer(KIND=8) state, site, Nsites
        doubleocc = .false.
        do site = 1,Nsites
            if(btest(state,site-1).and.btest(state,site+Nsites-1))then
            doubleocc = .true.
            return
            endif
        enddo
        return
    end function doubleocc
   

    
    subroutine apply_op(state,iszero,site,spin,N,sgn,op)
        implicit none
        logical iszero
        integer(KIND=8) k,newstate,state, site, spin, N, op
        double precision sgn
    
        newstate = state
        if(.not.iszero)then
            if(op.eq.create)then
                if(btest(newstate,site+spin*N-1))then
                    iszero = .true.
                else
                    do k = 2*N,site+1+spin*N,-1
        
                        if(btest(newstate,k-1)) sgn = sgn    !use  sgn = -sgn for fermions
                    
                    enddo
                    
                    newstate = ibset(newstate,site+spin*N-1)
                endif
            elseif(op.eq.annhil)then
                if(.not.btest(newstate,site+spin*N-1))then
                    iszero = .true.
                else
                    do k = 2*N,site+1+spin*N,-1
                        if(btest(newstate,k-1)) sgn = sgn   !use  sgn = -sgn for fermions
                    enddo
                    newstate = ibclr(newstate,site+spin*N-1)
                endif
            endif
        endif
    
        state = newstate
        
        return
    end subroutine apply_op

    
    subroutine insert_element(Hi,Hj,H,i,j,value,siz,point)
        implicit none
        integer(KIND=8) i, j, siz, point
        double precision, dimension(1:siz) :: H
        integer(KIND=8), dimension(1:siz) :: Hi, Hj
        double precision value
        point = point + 1
        if(point.ge.siz)then
            print*, 'Array is too small for storage of H, please increase the size.'
            print*, 'point = ', point, 'siz = ', siz
            stop
        else
        
            Hi(point) = i
            Hj(point) = j
            H(point) = value
        endif
        return
    end subroutine insert_element


    subroutine construct_hamiltonian(Hi,Hj,Hel,Hlen,maxHlen,basis,hs_size,ch)
        implicit none
        logical iszero
        integer(KIND=8) nn,mm, Hlen, maxHLen, hs_size, ch
        integer(KIND=8) state, newstate, i, j, jp, spin
        integer(KIND=8) niup, nidn, njup, njdn
        integer(KIND=8), dimension(1:hs_size) :: basis
        integer(KIND=8), dimension(1:maxHLen) :: Hi
        integer(KIND=8), dimension(1:maxHLen) :: Hj
        double precision, dimension(1:maxHLen) :: Hel
        double precision rtmp, rtmp2, fsgn, JJ
    
    
        do nn = 1,hs_size    !HAMILTONIAN DIMENSION
            state = basis(nn)
            
            !print*, "State",state, ": ", dec2bin(state, N)    
           
            rtmp = 0.0d0
            njdn= 0; nidn = 0;
            do i = 1,N
                j = i + 1
                if(j.eq.(N+1)) j = 1
                
            
                
                ! U * n_i * n_j
            
                if(btest(state,i-1).and.(btest(state,j+spindn*N-1))) nidn = nidn + 1
                
                !print*,i,j,btest(state,i+spindn*N-1),btest(state,j+spindn*N-1),nidn
                
                ! c_i^+ c_j^-  + c_i^- c_j^+
                newstate = state
                iszero = .false.
                fsgn = 1.0d0
                call apply_op(newstate,iszero,j,spindn,N,fsgn,annhil)
                call apply_op(newstate,iszero,i,spindn,N,fsgn,create)
                if(.not.iszero)then
                   
                    call binary_search(newstate,mm,basis,hs_size)
                   
                    call insert_element(Hi,Hj,Hel,nn,mm,fsgn*t,maxHlen,Hlen)
                endif
            
                newstate = state
                iszero = .false.
                fsgn = 1.0d0
               
                call apply_op(newstate,iszero,i,spindn,N,fsgn,annhil)
                call apply_op(newstate,iszero,j,spindn,N,fsgn,create)
                if(.not.iszero)then
                    call binary_search(newstate,mm,basis,hs_size)
                    call insert_element(Hi,Hj,Hel,nn,mm,fsgn*t,maxHlen,Hlen)
                endif
              
            enddo
            rtmp = rtmp + U*dfloat(nidn)
          
            if(rtmp.ne.0.0d0) call insert_element(Hi,Hj,Hel,nn,nn,rtmp,maxHlen,Hlen)
        enddo
   
            
        return
    end subroutine construct_hamiltonian

    subroutine construct_hamiltonian_delta_U(Hi,Hj,Hel,Hlen,maxHlen,basis,hs_size,ch)
        implicit none
        logical iszero
        integer(KIND=8) nn,mm, Hlen, maxHLen, hs_size, ch
        integer(KIND=8) state, newstate, i, j, jp, spin
        integer(KIND=8) niup, nidn, njup, njdn,n_delta
        integer(KIND=8), dimension(1:hs_size) :: basis
        integer(KIND=8), dimension(1:maxHLen) :: Hi
        integer(KIND=8), dimension(1:maxHLen) :: Hj
        double precision, dimension(1:maxHLen) :: Hel
        double precision rtmp, rtmp2, fsgn, JJ
    
    
        do nn = 1,hs_size    !HAMILTONIAN DIMENSION
            state = basis(nn)
            
            !print*, "State",state, ": ", dec2bin(state, N)    
           
            rtmp = 0.0d0
            njdn= 0; nidn = 0;
            n_delta=0
            do i = 1,N
                j = i + 1
                if(j.eq.(N+1)) j = 1
                
            
                
                ! U * n_i * n_j
            
                if(mod(i,4).eq.3)then
                    if(btest(state,i-1).and.(btest(state,j+spindn*N-1))) nidn = nidn + 1
                    if(btest(state,i-1)) n_delta=n_delta+1
                    if(btest(state,j+spindn*N-1)) n_delta=n_delta+1
                endif
                
                !print*,i,j,btest(state,i+spindn*N-1),btest(state,j+spindn*N-1),nidn
                
                ! c_i^+ c_j^-  + c_i^- c_j^+
                newstate = state
                iszero = .false.
                fsgn = 1.0d0
                call apply_op(newstate,iszero,j,spindn,N,fsgn,annhil)
                call apply_op(newstate,iszero,i,spindn,N,fsgn,create)
                if(.not.iszero)then
                   
                    call binary_search(newstate,mm,basis,hs_size)
                   
                    call insert_element(Hi,Hj,Hel,nn,mm,fsgn*t,maxHlen,Hlen)
                endif
            
                newstate = state
                iszero = .false.
                fsgn = 1.0d0
               
                call apply_op(newstate,iszero,i,spindn,N,fsgn,annhil)
                call apply_op(newstate,iszero,j,spindn,N,fsgn,create)
                if(.not.iszero)then
                    call binary_search(newstate,mm,basis,hs_size)
                    call insert_element(Hi,Hj,Hel,nn,mm,fsgn*t,maxHlen,Hlen)
                endif
              
            enddo
            rtmp = rtmp + U*dfloat(nidn)+ delta *dfloat(n_delta)
            
            if(rtmp.ne.0.0d0) call insert_element(Hi,Hj,Hel,nn,nn,rtmp,maxHlen,Hlen)
        enddo
   
           
        return
    end subroutine construct_hamiltonian_delta_U
    

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
end module cluster
   