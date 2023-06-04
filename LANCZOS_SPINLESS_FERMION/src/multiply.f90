 subroutine multiply_cvector_by_H(vec,newvec,Nhs,Hi,Hj,Hval,Hlen,maxHLen,hs_size)
 implicit none
 integer(KIND=8) nn, Nel, Hlen, maxHLen, Nhs, hs_size
 integer(KIND=8), dimension(1:maxHLen) :: Hi, Hj
 double complex, dimension(1:Nhs) :: vec, newvec
 double precision, dimension(1:maxHLen) :: Hval
 double precision rtmp
 !init the vector to zero
 newvec = dcmplx(0.0d0,0.0d0)
 !loop over the phonon subspace
 !$OMP PARALLEL DO SHARED(Hi,Hj,Hval,newvec,vec) PRIVATE(nn) 
 do nn = 1,Hlen
  !$OMP ATOMIC
  newvec(Hi(nn)) = newvec(Hi(nn)) + dcmplx(Hval(nn),0.0d0)*vec(Hj(nn))
 enddo 
 !$OMP END PARALLEL DO
 return
 end subroutine multiply_cvector_by_H 

 subroutine multiply_vector_by_H(vec,newvec,Nhs,Hi,Hj,Hval,Hlen,maxHLen,hs_size)
 implicit none
 integer(KIND=8) idx1, idx2, k, spin
 integer(KIND=8) nn, i1, i2, i3, i4, j1, j2, j3, j4
 integer(KIND=8) ii, jj, Nel, Hlen, maxHLen, Nhs, hs_size
 integer(KIND=8), dimension(1:maxHLen) :: Hi, Hj
 double precision, dimension(1:Nhs) :: vec, newvec
 double precision, dimension(1:maxHLen) :: Hval
 double precision rtmp

 !init the vector to zero
 newvec = 0.0d0
 !do the electronic part first
 do nn = 1,Hlen
  ii = Hi(nn) 
  jj = Hj(nn) 
  rtmp = Hval(nn) 
  newvec(ii) = newvec(ii) + rtmp*vec(jj)
 enddo 

 return
 end subroutine multiply_vector_by_H 
