subroutine stateq_bb(rho_a_r,rho_a_i,rho_b_r,rho_b_i,lF1,lF2,b1,b2,lbd,gam,om,sa,Iw,Qw,Uw,Vw,&
&delta_I,delta_Q,delta_U,delta_V,rho_new_r,rho_new_i)
implicit none
integer a1,a2,b1,b2,lF1,lF2,mb1,mb2,ma1,ma2,F1,F2,mb11,mb22
integer I,J,II,JJ
double precision lbd,gam,om,Const,sa
double precision h,x1,x2,pi,c0,hbar
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

Const = sa*pi*pi/(c0*hbar**x2)

rho_new_r = 0.d0
rho_new_i = 0.d0

! construct delaaI,delaaQ,delaaU,delaaV
F1 = (lF1 - 1)/2
F2 = (lF2 - 1)/2

mb11 = b1 - F2 - 1 
mb22 = b2 - F2 - 1

do a1 =1,lF1
  do a2 = 1,lF1
    ma1 = a1 - F1 - 1
    ma2 = a2 - F1 - 1
    do I=1,3
      do J=1,3
        mb1 = ma2 + I - 2 
        mb2 = ma1 + J - 2 

        if (mb1.eq.mb11.and.mb2.eq.mb22) then
          II = -(ma1 - mb2) + 2
          JJ = -(ma2 - mb1) + 2
          if (II.gt.0.and.II.lt.4.and.JJ.gt.0.and.JJ.lt.4) then       
            delaaI(II,JJ) = delta_I(a2,I,a1,J)
            delaaV(II,JJ) = delta_V(a2,I,a1,J)
            delaaQ(II,JJ) = delta_Q(a2,I,a1,J)
            delaaU(II,JJ) = delta_U(a2,I,a1,J)
          endif
        endif
      end do
    end do
  end do
end do


!not sure if the sums are right with regard to the real/imaginary parts.
call rat1_b(rho_a_r,rho_a_i,lF1,lF2,b1,b2,Iw,Qw,Uw,Vw,&
&delaaI,delaaQ,delaaU,delaaV,rat_1_r,rat_1_i)

call rat2_b(rho_b_r,rho_b_i,lF1,lF2,b1,b2,Iw,Qw,Uw,Vw,&
&delta_I,delta_Q,delta_U,delta_V,rat_2_r,rat_2_i)

call rat3_b(rho_b_r,rho_b_i,lF1,lF2,b1,b2,Iw,Qw,Uw,Vw,&
&delta_I,delta_Q,delta_U,delta_V,rat_3_r,rat_3_i)

! Hemel: no rho_a needed in these terms.
! destroy these out of there to make things faster.
call rat_single_2_b(lF1,lF2,b1,b2,Iw,Qw,Uw,Vw,&
&delta_I,delta_Q,delta_U,delta_V,rat_single_2_r,rat_single_2_i)

call rat_single_3_b(lF1,lF2,b1,b2,Iw,Qw,Uw,Vw,&
&delta_I,delta_Q,delta_U,delta_V,rat_single_3_r,rat_single_3_i)


rat_4_r = rat_single_2_r + rat_single_3_r
rat_4_i = rat_single_2_i + rat_single_3_i

teller_r = lbd/Const + rat_1_r - rat_2_r - rat_3_r   !Hemel: lbd is not multiplied by rho_aa'---pumping independent of the population
teller_i = rat_1_i - rat_2_i - rat_3_i

noem_r = gam/Const + rat_4_r
noem_i = om/Const + rat_4_i

noemer = noem_r**2 + noem_i**2

rho_new_r = (teller_r*noem_r + teller_i*noem_i) / (noemer)
rho_new_i = (noem_r*teller_i - teller_r*noem_i) / (noemer)

write(*,*)
write(*,*)b1,b2
write(*,*)Const*(rat_1_r-rat_2_r-rat_3_r)
write(*,*)Const*(rat_1_i-rat_2_i-rat_3_i)
write(*,*)gam+Const*rat_4_r
write(*,*)om+Const*rat_4_i

!write(*,*)
!write(*,*) 'b1,b2,rat_1_r,rat_2_r,rat_3_r,teller_r'
!write(*,*) b1,b2,rat_1_r,rat_2_r,rat_3_r,teller_r
!write(*,*) 'b1,b2,rat_1_i,rat_2_i,rat_3_i,teller_i'
!write(*,*) b1,b2,rat_1_i,rat_2_i,rat_3_i,teller_i
!write(*,*) 'rat_single_2_r,rat_single_2_i,rat_single_3_r,rat_single_3_i'
!write(*,*) rat_single_2_r,rat_single_2_i,rat_single_3_r,rat_single_3_i
!write(*,*)
!write(*,*) rho_new_r,rho_new_i
!write(*,*)rho_a_r(1,1),rho_a_r(1,2),rho_a_r(1,3)
!write(*,*)rho_a_r(2,1),rho_a_r(2,2),rho_a_r(2,3)
!write(*,*)rho_a_r(3,1),rho_a_r(3,2),rho_a_r(3,3)
!write(*,*)
!write(*,*)rho_a_i(1,1),rho_a_i(1,2),rho_a_i(1,3)
!write(*,*)rho_a_i(2,1),rho_a_i(2,2),rho_a_i(2,3)
!write(*,*)rho_a_i(3,1),rho_a_i(3,2),rho_a_i(3,3)
!
!
!stop

end subroutine

subroutine rat1_b(rho_a_r,rho_a_i,lF1,lF2,b1,b2,Iw,Qw,Uw,Vw,&
&delaaI,delaaQ,delaaU,delaaV,ratfaq_r,ratfaq_i)
implicit none
integer a1,a2,b1,b2,br1,br2,ma1,ma2,lF1,lF2,F1,F2
integer ac1,ac2,mb1,mb2
double precision ratfaq_r,ratfaq_i,zeta_r,zeta_i,x2
double precision Ifp_r,Qfp_r,Ufp_r,Vfp_r,Iw,Qw,Uw,Vw
double precision, dimension (lF1,lF1) :: rho_a_r,rho_a_i
double precision, dimension (3,3) :: delaaI,delaaQ,delaaU,delaaV
double precision, dimension (lF1,lF2) :: Ip_r,Ip_i,Qp_r,Qp_i,Up_r,Up_i,Vp_r,Vp_i


F1 = (lF1 - 1)/2
F2 = (lF2 - 1)/2

mb1 = b1 - 1 - F2
mb2 = b2 - 1 - F2

ratfaq_r = 0.d0
ratfaq_i = 0.d0

x2 = 2.d0

Ifp_r = x2*Iw 
Ufp_r = x2*Qw 
Qfp_r = x2*Uw 
Vfp_r = x2*Vw 

! real delaa
do ac1=1,3

  ac2 = ac1

  a1 = F1 + mb2 - ac1 + 3         ! going to b1  
  a2 = F1 + mb1 - ac2 + 3         ! going to b2  

  zeta_r = Ifp_r*delaaI(ac2,ac1) - Qfp_r*delaaQ(ac2,ac1) + Vfp_r*delaaV(ac2,ac1) 
  zeta_i = -Ufp_r*delaaU(ac2,ac1) 

  ratfaq_r = ratfaq_r + rho_a_r(a2,a1)*zeta_r - rho_a_i(a2,a1)*zeta_i
  ratfaq_i = ratfaq_i + rho_a_r(a2,a1)*zeta_i + rho_a_i(a2,a1)*zeta_r
end do

ac1 = 1
ac2 = 3

a1 = F1 + mb2 - ac1 + 3         ! going to b1  
a2 = F1 + mb1 - ac2 + 3         ! going to b2  

zeta_r = Ifp_r*delaaI(ac2,ac1) - Qfp_r*delaaQ(ac2,ac1) + Vfp_r*delaaV(ac2,ac1) 
zeta_i = -Ufp_r*delaaU(ac2,ac1) 

!write(*,*) a1,a2,zeta_r,zeta_i

ratfaq_r = ratfaq_r + rho_a_r(a2,a1)*zeta_r - rho_a_i(a2,a1)*zeta_i
ratfaq_i = ratfaq_i + rho_a_r(a2,a1)*zeta_i + rho_a_i(a2,a1)*zeta_r


ac1 = 3
ac2 = 1

a1 = F1 + mb2 - ac1 + 3         ! going to b1  
a2 = F1 + mb1 - ac2 + 3         ! going to b2  


zeta_r = Ifp_r*delaaI(ac2,ac1) - Qfp_r*delaaQ(ac2,ac1) + Vfp_r*delaaV(ac2,ac1) 
zeta_i = -Ufp_r*delaaU(ac2,ac1) 

!write(*,*) a1,a2,zeta_r,zeta_i

ratfaq_r = ratfaq_r + rho_a_r(a2,a1)*zeta_r - rho_a_i(a2,a1)*zeta_i
ratfaq_i = ratfaq_i + rho_a_r(a2,a1)*zeta_i + rho_a_i(a2,a1)*zeta_r

ac1 = 1
ac2 = 2

a1 = F1 + mb2 - ac1 + 3         ! going to b1  
a2 = F1 + mb1 - ac2 + 3         ! going to b2  

zeta_i = Ifp_r*delaaI(ac2,ac1) - Qfp_r*delaaQ(ac2,ac1) + Vfp_r*delaaV(ac2,ac1) 
zeta_r = Ufp_r*delaaU(ac2,ac1) 

!write(*,*) a1,a2,zeta_r,zeta_i

ratfaq_r = ratfaq_r + rho_a_r(a2,a1)*zeta_r - rho_a_i(a2,a1)*zeta_i
ratfaq_i = ratfaq_i + rho_a_r(a2,a1)*zeta_i + rho_a_i(a2,a1)*zeta_r



ac1 = 2
ac2 = 1

a1 = F1 + mb2 - ac1 + 3         ! going to b1  
a2 = F1 + mb1 - ac2 + 3         ! going to b2  

zeta_i = Ifp_r*delaaI(ac2,ac1) - Qfp_r*delaaQ(ac2,ac1) + Vfp_r*delaaV(ac2,ac1) 
zeta_r = Ufp_r*delaaU(ac2,ac1) 

!write(*,*) a1,a2,zeta_r,zeta_i

ratfaq_r = ratfaq_r + rho_a_r(a2,a1)*zeta_r - rho_a_i(a2,a1)*zeta_i
ratfaq_i = ratfaq_i + rho_a_r(a2,a1)*zeta_i + rho_a_i(a2,a1)*zeta_r


ac1 = 3
ac2 = 2

a1 = F1 + mb2 - ac1 + 3         ! going to b1  
a2 = F1 + mb1 - ac2 + 3         ! going to b2  

zeta_i = Ifp_r*delaaI(ac2,ac1) - Qfp_r*delaaQ(ac2,ac1) + Vfp_r*delaaV(ac2,ac1) 
zeta_r = Ufp_r*delaaU(ac2,ac1) 

!write(*,*) a1,a2,zeta_r,zeta_i

ratfaq_r = ratfaq_r + rho_a_r(a2,a1)*zeta_r - rho_a_i(a2,a1)*zeta_i
ratfaq_i = ratfaq_i + rho_a_r(a2,a1)*zeta_i + rho_a_i(a2,a1)*zeta_r


ac1 = 2
ac2 = 3

a1 = F1 + mb2 - ac1 + 3         ! going to b1  
a2 = F1 + mb1 - ac2 + 3         ! going to b2  


zeta_i = Ifp_r*delaaI(ac2,ac1) - Qfp_r*delaaQ(ac2,ac1) + Vfp_r*delaaV(ac2,ac1) 
zeta_r = Ufp_r*delaaU(ac2,ac1) 

!write(*,*) a1,a2,zeta_r,zeta_i

ratfaq_r = ratfaq_r + rho_a_r(a2,a1)*zeta_r - rho_a_i(a2,a1)*zeta_i
ratfaq_i = ratfaq_i + rho_a_r(a2,a1)*zeta_i + rho_a_i(a2,a1)*zeta_r

end subroutine


subroutine rat2_b(rho_b_r,rho_b_i,lF1,lF2,b1,b2,Iw,Qw,Uw,Vw,&
&delta_I,delta_Q,delta_U,delta_V,ratfaq_r,ratfaq_i)
implicit none
integer a1,a2,a3,b1,b2,b,mb,bc1,bc2,ma1,ma2,ma3,lF1,lF2,dm,F1,F2
integer ac1,ac2,b3,ma,mb1,mb2,mb3
double precision Iw,Qw,Uw,Vw
double precision ratfaq_r,ratfaq_i,zeta_r,zeta_i
double precision, dimension (lF2,lF2) :: rho_b_r,rho_b_i
double precision, dimension (lF1,3,lF1,3) :: delta_I,delta_Q,delta_U,delta_V

!check --- also for the other programs, if dm is defined correctly

F1 = (lF1 - 1)/2
F2 = (lF2 - 1)/2

mb1 = b1 - F2 - 1
mb2 = b2 - F2 - 1

ratfaq_r = 0.d0 
ratfaq_i = 0.d0

do ac1 =1,3
  ma = mb1 - (ac1 - 2) 
  a1 = ma + F1 + 1
  do b3 = 1,lF2
    if (b3.ne.b1) then
      mb3 = b3 - F2 - 1
      dm = -(ma - mb3)
  
      if (abs(dm).le.1) then
        bc1 = ac1
        bc2 = dm + 2
! check bc2 and bc1 

  ! if bc1+bc2 is even, then delta is real, otherwise it is imaginary
        if (mod(bc1+bc2,2).eq.0) then
          zeta_r = Iw*delta_I(a1,bc1,a1,bc2) - &
&Qw*delta_Q(a1,bc1,a1,bc2) + Vw*delta_V(a1,bc1,a1,bc2) 
          zeta_i = -Uw*delta_U(a1,bc1,a1,bc2)
        else
          zeta_i = Iw*delta_I(a1,bc1,a1,bc2) -& 
&Qw*delta_Q(a1,bc1,a1,bc2) + Vw*delta_V(a1,bc1,a1,bc2) 
          zeta_r = Uw*delta_U(a1,bc1,a1,bc2)
        endif
 
        ratfaq_r = ratfaq_r + rho_b_r(b3,b2)*zeta_r - rho_b_i(b3,b2)*zeta_i
        ratfaq_i = ratfaq_i + rho_b_r(b3,b2)*zeta_i + rho_b_i(b3,b2)*zeta_r
      endif 
    endif
      
  end do
end do

end subroutine



subroutine rat3_b(rho_b_r,rho_b_i,lF1,lF2,b1,b2,Iw,Qw,Uw,Vw,&
&delta_I,delta_Q,delta_U,delta_V,ratfaq_r,ratfaq_i)
implicit none
integer a1,a2,a3,b1,b2,b,mb,bc1,bc2,ma1,ma2,ma3,lF1,lF2,dm,F1,F2
integer ac1,ac2,b3,ma,mb1,mb2,mb3
double precision ratfaq_r,ratfaq_i,zeta_r,zeta_i
double precision Iw,Qw,Uw,Vw
double precision, dimension (lF2,lF2) :: rho_b_r,rho_b_i
double precision, dimension (lF1,3,lF1,3) :: delta_I,delta_Q,delta_U,delta_V

F1 = (lF1 - 1)/2
F2 = (lF2 - 1)/2

mb1 = b1 - F2 - 1
mb2 = b2 - F2 - 1

ratfaq_r = 0.d0
ratfaq_i = 0.d0

do ac1 =1,3
  ma = mb2 - (ac1 - 2)
  a1 = ma + F1 + 1
  do b3 = 1,lF2
    if (b3.ne.b2) then
      mb3 = b3 - F2 - 1
      dm = -(ma - mb3)

      if (abs(dm).le.1) then
        bc1 = ac1
        bc2 = dm + 2
! check bc2 and bc1 

  ! if bc1+bc2 is even, then delta is real, otherwise it is imaginary
        if (mod(bc1+bc2,2).eq.0) then
          zeta_r = Iw*delta_I(a1,bc2,a1,bc1) - Qw*delta_Q(a1,bc2,a1,bc1)&
& + Vw*delta_V(a1,bc2,a1,bc1) 
          zeta_i = -Uw*delta_U(a1,bc2,a1,bc1)
        else
          zeta_i = Iw*delta_I(a1,bc2,a1,bc1) - Qw*delta_Q(a1,bc2,a1,bc1)&
& + Vw*delta_V(a1,bc2,a1,bc1) 
          zeta_r = Uw*delta_U(a1,bc2,a1,bc1)
        endif

        ratfaq_r = ratfaq_r + rho_b_r(b1,b3)*zeta_r - rho_b_i(b1,b3)*zeta_i
        ratfaq_i = ratfaq_i + rho_b_r(b1,b3)*zeta_i + rho_b_i(b1,b3)*zeta_r
      endif
    endif

  end do
end do


end subroutine

subroutine rat_single_2_b(lF1,lF2,b1,b2,Iw,Qw,Uw,Vw,&
&delta_I,delta_Q,delta_U,delta_V,ratfaq_r,ratfaq_i)
implicit none
integer a1,a2,a3,b1,b2,b,mb,bc1,bc2,ma1,ma2,ma3,lF1,lF2,dm,F1,F2
integer ac1,ac2,b3,ma,mb1,mb2,mb3
double precision ratfaq_r,ratfaq_i,zeta_r,zeta_i
double precision Iw,Qw,Uw,Vw 
double precision, dimension (lF1,3,lF1,3) :: delta_I,delta_Q,delta_U,delta_V

F1 = (lF1 - 1)/2
F2 = (lF2 - 1)/2

mb1 = b1 - F2 - 1
mb2 = b2 - F2 - 1

ratfaq_r = 0.d0 
ratfaq_i = 0.d0

do ac1 =1,3
  ma = mb1 - (ac1 - 2) 
  a1 = ma + F1 + 1
  b3 = b1 
  mb3 = b3 - F2 - 1
  dm = -(ma - mb3)
  
  if (abs(dm).le.1) then
    bc1 = ac1
    bc2 = dm + 2

      if (mod(bc1+bc2,2).eq.0) then
        zeta_r = Iw*delta_I(a1,bc1,a1,bc2) - &
&Qw*delta_Q(a1,bc1,a1,bc2) + Vw*delta_V(a1,bc1,a1,bc2) 
        zeta_i = -Uw*delta_U(a1,bc1,a1,bc2)
      else
        zeta_i = Iw*delta_I(a1,bc1,a1,bc2) -&
&Qw*delta_Q(a1,bc1,a1,bc2) + Vw*delta_V(a1,bc1,a1,bc2) 
        zeta_r = Uw*delta_U(a1,bc1,a1,bc2)
      endif

 !   write(*,*)a1,b1,zeta_r,zeta_i
 
    ratfaq_r = ratfaq_r + zeta_r 
    ratfaq_i = ratfaq_i + zeta_i 
  endif 
      
end do

end subroutine



subroutine rat_single_3_b(lF1,lF2,b1,b2,Iw,Qw,Uw,Vw,&
&delta_I,delta_Q,delta_U,delta_V,ratfaq_r,ratfaq_i)
implicit none
integer a1,a2,a3,b1,b2,b,mb,bc1,bc2,ma1,ma2,ma3,lF1,lF2,dm,F1,F2
integer ac1,ac2,b3,ma,mb1,mb2,mb3
double precision ratfaq_r,ratfaq_i,zeta_r,zeta_i
double precision Iw,Qw,Uw,Vw
double precision, dimension (lF2,lF2) :: rho_b_r,rho_b_i
double precision, dimension (lF1,3,lF1,3) :: delta_I,delta_Q,delta_U,delta_V

F1 = (lF1 - 1)/2
F2 = (lF2 - 1)/2

mb1 = b1 - F2 - 1
mb2 = b2 - F2 - 1

ratfaq_r = 0.d0
ratfaq_i = 0.d0

do ac1 =1,3
  ma = mb2 - (ac1 - 2)
  a1 = ma + F1 + 1
  b3 = b2
  mb3 = b3 - F2 - 1
  dm = -(ma - mb3)

  if (abs(dm).le.1) then
    bc1 = ac1
    bc2 = dm + 2

       if (mod(bc1+bc2,2).eq.0) then
          zeta_r = Iw*delta_I(a1,bc2,a1,bc1) - Qw*delta_Q(a1,bc2,a1,bc1)&
& + Vw*delta_V(a1,bc2,a1,bc1) 
          zeta_i = -Uw*delta_U(a1,bc2,a1,bc1)
        else
          zeta_i = Iw*delta_I(a1,bc2,a1,bc1) - Qw*delta_Q(a1,bc2,a1,bc1)&
& + Vw*delta_V(a1,bc2,a1,bc1) 
          zeta_r = Uw*delta_U(a1,bc2,a1,bc1)
        endif

    ratfaq_r = ratfaq_r + zeta_r 
    ratfaq_i = ratfaq_i + zeta_i 
  endif
end do

end subroutine

