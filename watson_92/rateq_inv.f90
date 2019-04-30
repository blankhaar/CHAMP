subroutine stateq(lF1,lF2,lbd_a,lbd_b,gam_a,gam_b,sa,Ip_r,Qp_r,Vp_r,delta_I,&
&delta_Q,delta_U,delta_V,rho_a_r,rho_b_r)
implicit none
integer a1,a2,b1,b2,F1,F2,lF1,lF2,I,J,N
integer a11,a12,a13,b11,b12,b13,INFO,ma1,mb1
integer, dimension (lF1+lF2) :: IPIV 
double precision lbd_a,lbd_b,gam_a,gam_b,om,Const,sa
double precision pi,c0,hbar,h,x0,x1,x2
double precision teller_r,teller_i,noem_r,noem_i,noemer
double precision rho_new_r,zeta_tot,zeta_r
double precision, dimension (lF1) :: rho_a_r
double precision, dimension (lF2) :: rho_b_r
double precision, dimension (lF1,3) :: Ip_r,Qp_r,Vp_r
double precision, dimension (lF1,3) :: delta_I,delta_Q,delta_U,delta_V
double precision, dimension (3) :: rat_1_r
double precision, dimension (lF1+lF2,lF1+lF2) :: Am 
double precision, dimension (lF1+lF2) :: WORK,rhovec,dvec
 
! Rate equations: N & W 1994: (Eq. A.8)

!input
!pa: density matrix *
!pb: density matrix *
!a1,a2: which states will be selected for new density
!zetp:
!zetm: integrated (gam_p/m \times Zeta) *
!lbd: pumping *
!gam: decay
!om: frequency-difference a1 and a2

!* input for a certain velocity, v
!output
!pn: new density matrix element p(a1,a2) *


!build in commons: F1,F2,lF1,lF2
!                  c,hbar,pi

c0 = 299792458.d0
h = 6.6261E-34
x0 = 0.d0
x1 = 1.d0
x2 = 2.d0
pi = dacos(-x1)
hbar = h/(x2*pi)

Const = x2*sa*pi/(c0*hbar**x2)

F1 = (lF1 - 1)/2
F2 = (lF2 - 1)/2

do a1 = 1,lF1+lF2
  do a2 = 1,lF1+lF2
    Am(a1,a2) = x0
  end do
end do

do a1 = 1,lF1
  dvec(a1) = lbd_a
end do

do b1 = 1,lF2
  dvec(b1+lF1) = lbd_b
end do

!eigenfunctions [pa(:),pb(:)]
do a1 = 1,lF1
  call rat1_a(lF1,lF2,a1,Ip_r,Qp_r,Vp_r,delta_I,delta_Q,delta_U,&
&delta_V,rat_1_r,zeta_r)

  ma1 = a1 - 1 - F1

  b11 = ma1 + 1 - 2 + F2 + 1
  b12 = ma1 + 2 - 2 + F2 + 1
  b13 = ma1 + 3 - 2 + F2 + 1

  Am(a1,a1) = Const*zeta_r + gam_a

  if (b11.gt.0.and.b11.le.lF2) then
    Am(a1,b11+lF1) = -Const*rat_1_r(1)
  endif

  if (b12.gt.0.and.b12.le.lF2) then
    Am(a1,b12+lF1) = -Const*rat_1_r(2)
  endif

  if (b13.gt.0.and.b13.le.lF2) then
  Am(a1,b13+lF1) = -Const*rat_1_r(3)
  endif
end do

do b1 = 1,lF2
  call rat1_b(lF1,lF2,b1,Ip_r,Qp_r,Vp_r,&
&delta_I,delta_Q,delta_U,delta_V,rat_1_r,zeta_r)  

  mb1 = b1 - 1 - F2

  a11 = F1 + mb1 - 1 + 3
  a12 = F1 + mb1 - 2 + 3
  a13 = F1 + mb1 - 3 + 3

  Am(b1+lF1,b1+lF1) = Const*zeta_r + gam_b
 
  if (a11.gt.0.and.a11.le.lF1) then 
    Am(b1+lF1,a11) = -Const*rat_1_r(1)
  endif
  if (a12.gt.0.and.a12.le.lF1) then  
    Am(b1+lF1,a12) = -Const*rat_1_r(2)
  endif
  if (a13.gt.0.and.a12.le.lF1) then
    Am(b1+lF1,a13) = -Const*rat_1_r(3)
  endif
end do

N = lF1 + lF2
call DGETRF(N,N,Am,N,IPIV,INFO)
if (INFO.ne.0) then
  write(*,*)'Density-matrix calculation: failed'
  stop
else
  call DGETRI(N,Am,N,IPIV,WORK,N,INFO) 
  if (INFO.ne.0) then
    write(*,*)'Density-matrix calculation: failed'
    stop
  endif
endif

call matvec_prod(Am,dvec,N,rhovec)

do a1 = 1,lF1
  rho_a_r(a1) = rhovec(a1)
end do

do b1 = 1,lF2
  rho_b_r(b1) = rhovec(b1+lF1)
end do

end subroutine

subroutine rat1_a(lF1,lF2,a1,Ip_r,Qp_r,Vp_r,delta_I,delta_Q,delta_U,delta_V,ratfaq_r,zeta_tot)
! checked 19-06-2017
! consistent with A.8
! checked numerically
implicit none
integer a1,a2,b1,b2,br1,br2,ma1,ma2,lF1,lF2,F1,F2
double precision zeta_r,zeta_i,zeta_tot
double precision Ifp_r,Qfp_r,Vfp_r
double precision, dimension (lF2) :: rho_b_r
double precision, dimension (lF1,3) :: delta_I,delta_Q,delta_U,delta_V
double precision, dimension (lF1,3) :: Ip_r,Ip_i,Qp_r,Qp_i,Up_r,Up_i,Vp_r,Vp_i
double precision, dimension (3) :: ratfaq_r

F1 = (lF1 - 1)/2
F2 = (lF2 - 1)/2

ma1 = a1 - 1 - F1
ma2 = a2 - 1 - F1

ratfaq_r(:) = 0.d0
zeta_tot = 0.d0

do br1=1,3

  b1 = ma1 + br1 - 2 + F2 + 1

  Ifp_r = Ip_r(a1,br1)
  Qfp_r = Qp_r(a1,br1)
  Vfp_r = Vp_r(a1,br1)

  zeta_r = Ifp_r*delta_I(a1,br1) - Qfp_r*delta_Q(a1,br1) + Vfp_r*delta_V(a1,br1) 

  ratfaq_r(br1) = zeta_r
  zeta_tot = zeta_tot + zeta_r

end do

end subroutine

subroutine rat1_b(lF1,lF2,b1,Ip_r,Qp_r,Vp_r,&
&delta_I,delta_Q,delta_U,delta_V,ratfaq_r,zeta_tot)
implicit none
integer a1,a2,b1,b2,br1,br2,ma1,ma2,lF1,lF2,F1,F2
integer ac1,ac2,mb1,mb2
double precision zeta_r,zeta_tot
double precision Ifp_r,Qfp_r,Vfp_r
double precision, dimension (lF1) :: rho_a_r
double precision, dimension (lF1,3) :: delta_I,delta_Q,delta_U,delta_V
double precision, dimension (lF1,3) :: Ip_r,Qp_r,Vp_r
double precision, dimension (3) :: ratfaq_r

F1 = (lF1 - 1)/2
F2 = (lF2 - 1)/2

mb1 = b1 - 1 - F2

ratfaq_r(:) = 0.d0
zeta_tot = 0.d0

do ac1=1,3

  a1 = F1 + mb1 - ac1 + 3

  Ifp_r = Ip_r(a1,ac1)
  Qfp_r = Qp_r(a1,ac1)
  Vfp_r = Vp_r(a1,ac1)

  zeta_r = Ifp_r*delta_I(a1,ac1) - Qfp_r*delta_Q(a1,ac1) + Vfp_r*delta_V(a1,ac1) 

  ratfaq_r(ac1) = zeta_r
  zeta_tot = zeta_tot + zeta_r
end do


end subroutine



subroutine matvec_prod(A,b,N,c)
!compute the product A*b = c
!A = NxN matrix
!b = N-vector
!c = N-vector
integer N,I,J
double precision, dimension (N) :: b,c
double precision, dimension (N,N) :: A

do I=1,N
  c(I) = 0.d0
end do

do I = 1,N
  do J = 1,N
    c(I) = c(I) + A(I,J)*b(J)
  end do
end do

end subroutine

