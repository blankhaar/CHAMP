subroutine stateq_aa(rho_a_r,rho_a_i,rho_b_r,rho_b_i,lF1,lF2,a1,a2,lbd,gam,om,sa,Iw,Qw,Uw,Vw,delta_I,&
&delta_Q,delta_U,delta_V,rho_new_r,rho_new_i)
implicit none
integer a1,a2,F1,F2,lF1,lF2,I,J
double precision lbd,gam,om,Const,sa
double precision pi,c0,hbar,h,x1,x2
double precision rat_1_r,rat_2_r,rat_3_r,rat_single_2_r,rat_single_3_r
double precision rat_1_i,rat_2_i,rat_3_i,rat_single_2_i,rat_single_3_i
double precision rat_4_r,rat_4_i
double precision teller_r,teller_i,noem_r,noem_i,noemer
double precision rho_new_r,rho_new_i
double precision Iw,Qw,Uw,Vw
double precision, dimension (lF1,lF1) :: rho_a_r,rho_a_i 
double precision, dimension (lF2,lF2) :: rho_b_r,rho_b_i 
double precision, dimension (3,3) :: delaaI,delaaQ,delaaU,delaaV
double precision, dimension (lF1,3,lF1,3) :: delta_I,delta_Q,delta_U,delta_V

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
x1 = 1.d0
x2 = 2.d0
pi = dacos(-x1)
hbar = h/(x2*pi)

Const = sa*pi*pi/(c0*hbar**x2) ! extra factor pi from line centre approximation


rho_new_r = 0.d0
rho_new_i = 0.d0

do I=1,3
  do J=1,3
    delaaI(I,J) = delta_I(a2,I,a1,J)
    delaaV(I,J) = delta_V(a2,I,a1,J)
    delaaQ(I,J) = delta_Q(a2,I,a1,J)
    delaaU(I,J) = delta_U(a2,I,a1,J)
  end do
end do

call rat1_a(rho_b_r,rho_b_i,lF1,lF2,a1,a2,Iw,Qw,Uw,Vw,delaaI,delaaQ,delaaU,&
&delaaV,rat_1_r,rat_1_i)

call rat2_a(rho_a_r,rho_a_i,lF1,lF2,a1,a2,Iw,Qw,Uw,Vw,delta_I,delta_Q,delta_U,&
&delta_V,rat_2_r,rat_2_i)

call rat3_a(rho_a_r,rho_a_i,lF1,lF2,a1,a2,Iw,Qw,Uw,Vw,delta_I,delta_Q,delta_U,&
&delta_V,rat_3_r,rat_3_i)

call rat_single_2_a(lF1,lF2,a1,a2,Iw,Qw,Uw,Vw,delta_I,delta_Q,delta_U,delta_V,&
&rat_single_2_r,rat_single_2_i)

call rat_single_3_a(lF1,lF2,a1,a2,Iw,Qw,Uw,Vw,delta_I,delta_Q,delta_U,delta_V,&
&rat_single_3_r,rat_single_3_i)


rat_4_r = rat_single_2_r + rat_single_3_r
rat_4_i = rat_single_2_i + rat_single_3_i

teller_r = lbd/Const + rat_1_r - rat_2_r - rat_3_r  
teller_i = rat_1_i - rat_2_i - rat_3_i

noem_r = gam/Const + rat_4_r
noem_i = om/Const + rat_4_i

noemer = noem_r**2 + noem_i**2

rho_new_r = (teller_r*noem_r + teller_i*noem_i) / (noemer)
rho_new_i = (noem_r*teller_i - teller_r*noem_i) / (noemer)

write(*,*)
write(*,*)a1,a2
write(*,*)Const*(rat_1_r-rat_2_r-rat_3_r),Const*rat_1_r,-Const*rat_2_r,-Const*rat_3_r
write(*,*)Const*(rat_1_i-rat_2_i-rat_3_i)
write(*,*)gam + Const*rat_4_r
write(*,*)om  + Const*rat_4_i

!stop

end subroutine

subroutine rat1_a(rho_b_r,rho_b_i,lF1,lF2,a1,a2,Iw,Qw,Uw,Vw,delaaI,delaaQ,delaaU,delaaV,ratfaq_r,&
&ratfaq_i)
! checked 19-06-2017
! consistent with A.8
! checked numerically
implicit none
integer a1,a2,b1,b2,br1,br2,ma1,ma2,lF1,lF2,F1,F2
double precision ratfaq_r,ratfaq_i,zeta_r,zeta_i
double precision Iw,Qw,Uw,Vw,x2 
double precision Ifp_r,Qfp_r,Ufp_r,Vfp_r 
double precision, dimension (lF2,lF2) :: rho_b_r,rho_b_i
double precision, dimension (3,3) :: delaaI,delaaQ,delaaU,delaaV
double precision, dimension (lF1,lF2) :: Ip_r,Ip_i,Qp_r,Qp_i,Up_r,Up_i,Vp_r,Vp_i


F1 = (lF1 - 1)/2
F2 = (lF2 - 1)/2

ma1 = a1 - 1 - F1
ma2 = a2 - 1 - F1

x2 = 2.d0

ratfaq_r = 0.d0
ratfaq_i = 0.d0

Ifp_r = x2*Iw
Ufp_r = x2*Uw
Qfp_r = x2*Qw
Vfp_r = x2*Vw

do b1=1,3

  b2 = b1

  br1 = ma2 + b1 - 2 + F2 + 1         ! going to b1  
  br2 = ma1 + b2 - 2 + F2 + 1         ! going to b2  

  zeta_r = Ifp_r*delaaI(b1,b2) - Qfp_r*delaaQ(b1,b2) + Vfp_r*delaaV(b1,b2) 
  zeta_i = -Ufp_r*delaaU(b1,b2) 

  ratfaq_r = ratfaq_r + rho_b_r(br2,br1)*zeta_r - rho_b_i(br2,br1)*zeta_i
  ratfaq_i = ratfaq_i + rho_b_r(br2,br1)*zeta_i + rho_b_i(br2,br1)*zeta_r

end do

b1 = 1
b2 = 3

br1 = ma2 + b1 - 2 + F2 + 1         ! going to b1  
br2 = ma1 + b2 - 2 + F2 + 1         ! going to b2  


zeta_r = Ifp_r*delaaI(b1,b2) - Qfp_r*delaaQ(b1,b2) + Vfp_r*delaaV(b1,b2) 
zeta_i = -Ufp_r*delaaU(b1,b2) 

ratfaq_r = ratfaq_r + rho_b_r(br2,br1)*zeta_r - rho_b_i(br2,br1)*zeta_i
ratfaq_i = ratfaq_i + rho_b_r(br2,br1)*zeta_i + rho_b_i(br2,br1)*zeta_r

b1 = 3
b2 = 1

br1 = ma2 + b1 - 2 + F2 + 1         ! going to b1  
br2 = ma1 + b2 - 2 + F2 + 1         ! going to b2  


zeta_r = Ifp_r*delaaI(b1,b2) - Qfp_r*delaaQ(b1,b2) + Vfp_r*delaaV(b1,b2) 
zeta_i = -Ufp_r*delaaU(b1,b2) 

ratfaq_r = ratfaq_r + rho_b_r(br2,br1)*zeta_r - rho_b_i(br2,br1)*zeta_i
ratfaq_i = ratfaq_i + rho_b_r(br2,br1)*zeta_i + rho_b_i(br2,br1)*zeta_r

b1 = 1
b2 = 2

br1 = ma2 + b1 - 2 + F2 + 1         ! going to b1  
br2 = ma1 + b2 - 2 + F2 + 1         ! going to b2  


zeta_i = Ifp_r*delaaI(b1,b2) - Qfp_r*delaaQ(b1,b2) + Vfp_r*delaaV(b1,b2) 
zeta_r = Ufp_r*delaaU(b1,b2) 

ratfaq_r = ratfaq_r + rho_b_r(br2,br1)*zeta_r - rho_b_i(br2,br1)*zeta_i
ratfaq_i = ratfaq_i + rho_b_r(br2,br1)*zeta_i + rho_b_i(br2,br1)*zeta_r

b1 = 2
b2 = 1

br1 = ma2 + b1 - 2 + F2 + 1         ! going to b1  
br2 = ma1 + b2 - 2 + F2 + 1         ! going to b2  

zeta_i = Ifp_r*delaaI(b1,b2) - Qfp_r*delaaQ(b1,b2) + Vfp_r*delaaV(b1,b2) 
zeta_r = Ufp_r*delaaU(b1,b2) 

ratfaq_r = ratfaq_r + rho_b_r(br2,br1)*zeta_r - rho_b_i(br2,br1)*zeta_i
ratfaq_i = ratfaq_i + rho_b_r(br2,br1)*zeta_i + rho_b_i(br2,br1)*zeta_r

b1 = 3
b2 = 2

br1 = ma2 + b1 - 2 + F2 + 1         ! going to b1  
br2 = ma1 + b2 - 2 + F2 + 1         ! going to b2  

zeta_i = Ifp_r*delaaI(b1,b2) - Qfp_r*delaaQ(b1,b2) + Vfp_r*delaaV(b1,b2) 
zeta_r = Ufp_r*delaaU(b1,b2) 

ratfaq_r = ratfaq_r + rho_b_r(br2,br1)*zeta_r - rho_b_i(br2,br1)*zeta_i
ratfaq_i = ratfaq_i + rho_b_r(br2,br1)*zeta_i + rho_b_i(br2,br1)*zeta_r

b1 = 2
b2 = 3

br1 = ma2 + b1 - 2 + F2 + 1         ! going to b1  
br2 = ma1 + b2 - 2 + F2 + 1         ! going to b2  

zeta_i = Ifp_r*delaaI(b1,b2) - Qfp_r*delaaQ(b1,b2) + Vfp_r*delaaV(b1,b2) 
zeta_r = Ufp_r*delaaU(b1,b2) 

ratfaq_r = ratfaq_r + rho_b_r(br2,br1)*zeta_r - rho_b_i(br2,br1)*zeta_i
ratfaq_i = ratfaq_i + rho_b_r(br2,br1)*zeta_i + rho_b_i(br2,br1)*zeta_r

end subroutine


subroutine rat2_a(rho_a_r,rho_a_i,lF1,lF2,a1,a2,Iw,Qw,Uw,Vw,delta_I,delta_Q,delta_U,delta_V,&
&ratfaq_r,ratfaq_i)
! we take the correct D&W 1990 formula A.39 for these
! checked 19-06-2017
! consistent with A.39
! not checked numerically
implicit none
integer a1,a2,a3,b1,b2,b,mb,bc1,bc2,ma1,ma2,ma3,lF1,lF2,F1,F2,dm
double precision ratfaq_r,ratfaq_i,zeta_r,zeta_i
double precision Iw,Qw,Uw,Vw
double precision, dimension (lF1,lF1) :: rho_a_r,rho_a_i
double precision, dimension (lF1,3,lF1,3) :: delta_I,delta_Q,delta_U,delta_V

!check --- also for the other programs, if dm is defined correctly

F1 = (lF1 - 1)/2
F2 = (lF2 - 1)/2

ma1 = a1 - F1 - 1
ma2 = a2 - F1 - 1

ratfaq_r = 0.d0 
ratfaq_i = 0.d0

do bc1 =1,3
  mb = ma1 + (bc1 - 2)
  b = mb + F2 + 1
  do a3 = 1,lF1
    if (a3.ne.a1) then
      ma3 = a3 - F1 - 1
      dm = mb - ma3
  
      if (abs(dm).le.1) then
        bc2 = dm + 2
  
  ! if bc1+bc2 is even, then delta is real, otherwise it is imaginary
        if (mod(bc1+bc2,2).eq.0) then
          zeta_r = Iw*delta_I(a1,bc1,a3,bc2) - Qw*delta_Q(a1,bc1,a3,bc2) + &
& Vw*delta_V(a1,bc1,a3,bc2) 
          zeta_i = Uw*delta_U(a1,bc1,a3,bc2)
        else
          zeta_i = -Iw*delta_I(a1,bc1,a3,bc2) + Qw*delta_Q(a1,bc1,a3,bc2) - &
& Vw*delta_V(a1,bc1,a3,bc2) 
          zeta_r = Uw*delta_U(a1,bc1,a3,bc2)
        endif
 
        ratfaq_r = ratfaq_r + rho_a_r(a3,a2)*zeta_r - rho_a_i(a3,a2)*zeta_i
        ratfaq_i = ratfaq_i + rho_a_r(a3,a2)*zeta_i + rho_a_i(a3,a2)*zeta_r
      endif 
    endif
      
  end do
end do

end subroutine



subroutine rat3_a(rho_a_r,rho_a_i,lF1,lF2,a1,a2,Iw,Qw,Uw,Vw,delta_I,delta_Q,delta_U,&
&delta_V,ratfaq_r,ratfaq_i)
! we take the correct D&W 1990 formula A.39 for these
! checked 19-06-2017
! consistent with A.39
! not checked numerically
implicit none
integer a1,a2,a3,b1,b2,b,mb,bc1,bc2,ma1,ma2,ma3,lF1,lF2,dm,F1,F2
double precision ratfaq_r,ratfaq_i,zeta_r,zeta_i
double precision Iw,Qw,Uw,Vw
double precision, dimension (lF1,lF1) :: rho_a_r,rho_a_i
double precision, dimension (lF1,3,lF1,3) :: delta_I,delta_Q,delta_U,delta_V

!check --- also for the other programs, if dm is defined correctly

F1 = (lF1 - 1)/2
F2 = (lF2 - 1)/2

ma1 = a1 - F1 - 1
ma2 = a2 - F1 - 1

ratfaq_r = 0.d0
ratfaq_i = 0.d0

do bc1 =1,3
  mb = ma2 + (bc1 - 2)
  b = mb + F2 + 1
  do a3 = 1,lF1
    if (a3.ne.a2) then
      ma3 = a3 - F1 - 1
      dm = mb - ma3

      if (abs(dm).le.1) then
        bc2 = dm + 2

  ! if bc1+bc2 is even, then delta is real, otherwise it is imaginary
        if (mod(bc1+bc2,2).eq.0) then
          zeta_r = Iw*delta_I(a2,bc1,a3,bc2) - Qw*delta_Q(a2,bc1,a3,bc2) + &
& Vw*delta_V(a2,bc1,a3,bc2) 
          zeta_i = -Uw*delta_U(a2,bc1,a3,bc2)
        else
          zeta_i = Iw*delta_I(a2,bc1,a3,bc2) - Qw*delta_Q(a2,bc1,a3,bc2) + &
& Vw*delta_V(a2,bc1,a3,bc2) 
          zeta_r = Uw*delta_U(a2,bc1,a3,bc2)
        endif
 
        ratfaq_r = ratfaq_r + rho_a_r(a1,a3)*zeta_r - rho_a_i(a1,a3)*zeta_i
        ratfaq_i = ratfaq_i + rho_a_r(a1,a3)*zeta_i + rho_a_i(a1,a3)*zeta_r
      endif
    endif

  end do
end do

end subroutine

subroutine rat_single_2_a(lF1,lF2,a1,a2,Iw,Qw,Uw,Vw,delta_I,delta_Q,&
&delta_U,delta_V,ratfaq_r,ratfaq_i)
! we take the correct D&W 1990 formula A.39 for these
! checked 19-06-2017
! consistent with A.39
! checked numerically
implicit none
integer a1,a2,a3,b1,b2,b,mb,bc1,bc2,ma1,ma2,ma3,lF1,lF2,dm,F1,F2
double precision ratfaq_r,ratfaq_i,zeta_r,zeta_i
double precision Iw,Qw,Uw,Vw
double precision, dimension (lF1,3,lF1,3) :: delta_I,delta_Q,delta_U,delta_V

F1 = (lF1 - 1)/2
F2 = (lF2 - 1)/2

ma1 = a1 - F1 - 1
ma2 = a2 - F1 - 1

ratfaq_r = 0.d0 
ratfaq_i = 0.d0


do bc1 =1,3
  mb = ma1 + (bc1 - 2)
  b = mb + F2 + 1
  a3 = a1 
  ma3 = a3 - F1 - 1
  dm = mb - ma3
 
  if (abs(dm).le.1) then
    bc2 = dm + 2
 
! if bc1+bc2 is even, then delta is real, otherwise it is imaginary
! zeta = <\gamma_-^{a'b} \Zeta^{ab,ab}>
    if (mod(bc1+bc2,2).eq.0) then
      zeta_r = Iw*delta_I(a1,bc1,a3,bc2) - Qw*delta_Q(a1,bc1,a3,bc2) +  &
& Vw*delta_V(a1,bc1,a3,bc2) 
      zeta_i = Uw*delta_U(a1,bc1,a3,bc2)
    else
      zeta_i = -Iw*delta_I(a1,bc1,a3,bc2) + Qw*delta_Q(a1,bc1,a3,bc2) - &
& Vw*delta_V(a1,bc1,a3,bc2) 
      zeta_r = Uw*delta_U(a1,bc1,a3,bc2)
    endif

    ratfaq_r = ratfaq_r + zeta_r 
    ratfaq_i = ratfaq_i + zeta_i 
  endif 
      
end do

end subroutine



subroutine rat_single_3_a(lF1,lF2,a1,a2,Iw,Qw,Uw,Vw,delta_I,&
&delta_Q,delta_U,delta_V,ratfaq_r,ratfaq_i)
! we take the correct D&W 1990 formula A.39 for these
! checked 19-06-2017
! consistent with A.39
! checked numerically
implicit none
integer a1,a2,a3,b1,b2,b,mb,bc1,bc2,ma1,ma2,ma3,lF1,lF2,dm,F1,F2
double precision ratfaq_r,ratfaq_i,zeta_r,zeta_i
double precision Iw,Qw,Uw,Vw
double precision, dimension (lF1,3,lF1,3) :: delta_I,delta_Q,delta_U,delta_V

F1 = (lF1 - 1)/2
F2 = (lF2 - 1)/2

ma1 = a1 - F1 - 1
ma2 = a2 - F1 - 1

ratfaq_r = 0.d0
ratfaq_i = 0.d0

do bc1 =1,3
  mb = ma2 + (bc1 - 2)
  b = mb + F2 + 1
  a3 = a2
  ma3 = a3 - F1 - 1
  dm = mb - ma3

  if (abs(dm).le.1) then
    bc2 = dm + 2

!if bc1+bc2 is even, then delta is real, otherwise it is imaginary
    if (mod(bc1+bc2,2).eq.0) then
      zeta_r = Iw*delta_I(a2,bc1,a3,bc2) - Qw*delta_Q(a2,bc1,a3,bc2) +  &
& Vw*delta_V(a2,bc1,a3,bc2) 
      zeta_i = -Uw*delta_U(a2,bc1,a3,bc2)
    else
      zeta_i = Iw*delta_I(a2,bc1,a3,bc2) - Qw*delta_Q(a2,bc1,a3,bc2) +  &
& Vw*delta_V(a2,bc1,a3,bc2) 
      zeta_r = Uw*delta_U(a2,bc1,a3,bc2)
    endif

    ratfaq_r = ratfaq_r + zeta_r 
    ratfaq_i = ratfaq_i + zeta_i 
  endif
end do

end subroutine


