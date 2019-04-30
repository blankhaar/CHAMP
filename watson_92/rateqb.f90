subroutine stateq_bb(rho_a_r,lF1,lF2,b1,lbd,gam,sa,Ip_r,Qp_r,&
&Vp_r,delta_I,delta_Q,delta_U,delta_V,rho_new_r)
implicit none
integer a1,a2,b1,b2,lF1,lF2,mb1,mb2,ma1,ma2,F1,F2,mb11,mb22
integer I,J,II,JJ
double precision lbd,gam,om,Const,sa
double precision h,x1,x2,pi,c0,hbar
double precision rat_1_r,rat_2_r,rat_3_r,rat_single_2_r,rat_single_3_r
double precision rat_1_i,rat_2_i,rat_3_i,rat_single_2_i,rat_single_3_i
double precision rat_4_r,rat_4_i
double precision teller_r,teller_i,noem_r,noem_i,noemer
double precision rho_new_r,zeta_tot
double precision, dimension (lF1) :: rho_a_r
double precision, dimension (lF1,3) :: Ip_r,Ip_i,Qp_r,Qp_i,Up_i,Vp_r,Vp_i
double precision, dimension (lF1,3) :: delta_I,delta_Q,delta_U,delta_V

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

Const = x2*sa*pi/(c0*hbar**x2)

rho_new_r = 0.d0

F1 = (lF1 - 1)/2
F2 = (lF2 - 1)/2

call rat1_b(rho_a_r,lF1,lF2,b1,Ip_r,Qp_r,Vp_r,&
&delta_I,delta_Q,delta_U,delta_V,rat_1_r,zeta_tot)

teller_r = lbd/Const + rat_1_r 

noem_r = gam/Const + zeta_tot 


rho_new_r = teller_r / noem_r 

end subroutine

subroutine rat1_b(rho_a_r,lF1,lF2,b1,Ip_r,Qp_r,Vp_r,&
&delta_I,delta_Q,delta_U,delta_V,ratfaq_r,zeta_tot)
implicit none
integer a1,a2,b1,b2,br1,br2,ma1,ma2,lF1,lF2,F1,F2
integer ac1,ac2,mb1,mb2
double precision ratfaq_r,ratfaq_i,zeta_r,zeta_tot
double precision Ifp_r,Qfp_r,Vfp_r
double precision, dimension (lF1) :: rho_a_r
double precision, dimension (lF1,3) :: delta_I,delta_Q,delta_U,delta_V
double precision, dimension (lF1,3) :: Ip_r,Qp_r,Vp_r


F1 = (lF1 - 1)/2
F2 = (lF2 - 1)/2

mb1 = b1 - 1 - F2

ratfaq_r = 0.d0
zeta_tot = 0.d0

do ac1=1,3

  a1 = F1 + mb1 - ac1 + 3         

  Ifp_r = Ip_r(a1,ac1) 
  Qfp_r = Qp_r(a1,ac1) 
  Vfp_r = Vp_r(a1,ac1) 

  zeta_r = Ifp_r*delta_I(a1,ac1) - Qfp_r*delta_Q(a1,ac1) + Vfp_r*delta_V(a1,ac1)  

  ratfaq_r = ratfaq_r + rho_a_r(a1)*zeta_r
  zeta_tot = zeta_tot + zeta_r 
end do


end subroutine
