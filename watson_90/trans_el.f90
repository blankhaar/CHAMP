subroutine prop_stokes(A,B,C,D,E,F,G,I,Q,U,V,ds)
implicit none
integer n,nv
double precision Inew,Qnew,Unew,Vnew,ds                
double precision A,B,C,D,E,F,G
double precision I,Q,U,V 
double precision, dimension (4,4) :: eK,eeK 

! propagating the stokes parameters
! all variables are real

!input
!A,B,C,D,E,F,G: the propagation coefficients
!I,Q,U,V: the old stokes parameters

!output
!I,Q,U,V: the new stokes parameters

  eK(1,1) = A*ds 
  eK(2,2) = A*ds 
  eK(3,3) = A*ds 
  eK(4,4) = A*ds 

  eK(1,2) = B*ds 
  eK(2,1) = B*ds 

  eK(1,3) = F*ds 
  eK(3,1) = F*ds 

  eK(1,4) = C*ds 
  eK(4,1) = C*ds 

  eK(2,3) = E*ds 
  eK(3,2) = -E*ds 

  eK(2,4) = G*ds
  eK(4,2) = -G*ds

  eK(3,4) = D*ds
  eK(4,3) = -D*ds

  call exp_A(eK,8,eeK)

  I = eeK(1,1)*I + eeK(1,2)*Q + eeK(1,3)*U + eeK(1,4)*V
  Q = eeK(2,1)*I + eeK(2,2)*Q + eeK(2,3)*U + eeK(2,4)*V
  U = eeK(3,1)*I + eeK(3,2)*Q + eeK(3,3)*U + eeK(3,4)*V
  V = eeK(4,1)*I + eeK(4,2)*Q + eeK(4,3)*U + eeK(4,4)*V

end subroutine


subroutine trans_el(rho_aa_r,rho_aa_i,rho_bb_r,rho_bb_i,lF1,lF2,&
&delta_I,delta_Q,delta_U,delta_V,Aw,Bw,Cw,Dw,Ew,Fw,Gw,Aws,Bws,Cws,Fws,w0)
implicit none
integer n,m,nw,a1,a2,b1,b2,ac1,ac2,bc1,bc2
integer F1,F2,lF1,lF2,ma1,ma2,mb1,mb2
double precision transa_A_r,transa_B_r,transa_C_r,transa_D_r,transa_E_r,transa_F_r,transa_G_r 
double precision transb_A_r,transb_B_r,transb_C_r,transb_D_r,transb_E_r,transb_F_r,transb_G_r 
double precision transa_A_i,transa_B_i,transa_C_i,transa_D_i,transa_E_i,transa_F_i,transa_G_i 
double precision transb_A_i,transb_B_i,transb_C_i,transb_D_i,transb_E_i,transb_F_i,transb_G_i 
double precision Const,Constp,pi,c0,hbar,w0 
double precision Ar,Br,Cr,Dr,Er,Fr,Gr,Ai,Bi,Ci,Di,Ei,Fi,Gi
double precision Aw,Bw,Cw,Dw,Ew,Fw,Gw,Aws,Bws,Cws,Fws
double precision, dimension (lF1,lF1) :: rho_aa_r,rho_aa_i 
double precision, dimension (lF2,lF2) :: rho_bb_r,rho_bb_i 
double precision, dimension (lF1,3,lF1,3) :: delta_I,delta_Q,delta_U,delta_V

pi = dacos(-1.d0)
hbar = 1.0545718E-34
Const = 2.d0*pi*pi/hbar 
c0 = 299792458.d0 

Constp = ( 4.d0 * pi**2.d0 * hbar * w0**3.d0 ) / ( c0**2.d0 ) 

F1 = (lF1 - 1)/2
F2 = (lF2 - 1)/2


! first the sum over ab, a'

transa_A_r = 0.d0 
transa_B_r = 0.d0 
transa_C_r = 0.d0 
transa_F_r = 0.d0 

transa_A_i = 0.d0 
transa_B_i = 0.d0 
transa_C_i = 0.d0 
transa_F_i = 0.d0 

do b1 = 1,lF2
  do ac1 = 1,3
    do ac2 = 1,3
      mb1 = b1 - F2 - 1

      ma1 = mb1 - (ac1 - 2)
      ma2 = mb1 - (ac2 - 2)

      a1 = ma1 + F1 + 1
      a2 = ma2 + F1 + 1

      bc1 = (mb1 - ma1)+2  !dm = mb-ma
      bc2 = (mb1 - ma2)+2  


      if (a1.gt.0.and.a2.gt.0.and.a1.le.lF1.and.a2.le.lF1) then
  !delta even or uneven
        if (mod((bc1+bc2),2).eq.0) then 
          transa_A_r = transa_A_r +  rho_aa_r(a1,a2)*delta_I(a1,bc1,a2,bc2)
          transa_B_r = transa_B_r +  rho_aa_r(a1,a2)*delta_Q(a1,bc1,a2,bc2)
          transa_C_r = transa_C_r +  rho_aa_r(a1,a2)*delta_V(a1,bc1,a2,bc2)
          transa_F_r = transa_F_r +  rho_aa_r(a1,a2)*delta_U(a1,bc1,a2,bc2)
  
          transa_A_i = transa_A_i +  rho_aa_i(a1,a2)*delta_I(a1,bc1,a2,bc2)
          transa_B_i = transa_B_i +  rho_aa_i(a1,a2)*delta_Q(a1,bc1,a2,bc2)
          transa_C_i = transa_C_i +  rho_aa_i(a1,a2)*delta_V(a1,bc1,a2,bc2)
          transa_F_i = transa_F_i +  rho_aa_i(a1,a2)*delta_U(a1,bc1,a2,bc2)
        else
          transa_A_i = transa_A_i +  rho_aa_r(a1,a2)*delta_I(a1,bc1,a2,bc2)
          transa_B_i = transa_B_i +  rho_aa_r(a1,a2)*delta_Q(a1,bc1,a2,bc2)
          transa_C_i = transa_C_i +  rho_aa_r(a1,a2)*delta_V(a1,bc1,a2,bc2)
          transa_F_i = transa_F_i +  rho_aa_r(a1,a2)*delta_U(a1,bc1,a2,bc2)
  
          transa_A_r = transa_A_r -  rho_aa_i(a1,a2)*delta_I(a1,bc1,a2,bc2)
          transa_B_r = transa_B_r -  rho_aa_i(a1,a2)*delta_Q(a1,bc1,a2,bc2)
          transa_C_r = transa_C_r -  rho_aa_i(a1,a2)*delta_V(a1,bc1,a2,bc2)
          transa_F_r = transa_F_r -  rho_aa_i(a1,a2)*delta_U(a1,bc1,a2,bc2)
        endif 

      endif

    end do
  end do
end do

transb_A_r = 0.d0
transb_B_r = 0.d0
transb_C_r = 0.d0
transb_F_r = 0.d0

transb_A_i = 0.d0
transb_B_i = 0.d0
transb_C_i = 0.d0
transb_F_i = 0.d0

do a1 = 1,lF1
  do bc1 = 1,3
    do bc2 = 1,3
      ma1 = a1 - F1 - 1
      mb1 = ma1 + (bc1 - 2)  !dm = mb - ma
      mb2 = ma1 + (bc2 - 2)

      b1 = mb1 + F2 + 1
      b2 = mb2 + F2 + 1

      if (b1.gt.0.and.b2.gt.0.and.b1.le.lF2.and.b2.le.lF2) then

  ! delta even or uneven
        if (mod((bc1+bc2),2).eq.0) then
          transb_A_r = transb_A_r +  rho_bb_r(b2,b1)*delta_I(a1,bc1,a1,bc2)
          transb_B_r = transb_B_r +  rho_bb_r(b2,b1)*delta_Q(a1,bc1,a1,bc2)
          transb_C_r = transb_C_r +  rho_bb_r(b2,b1)*delta_V(a1,bc1,a1,bc2)
          transb_F_r = transb_F_r +  rho_bb_r(b2,b1)*delta_U(a1,bc1,a1,bc2)
  
          transb_A_i = transb_A_i +  rho_bb_i(b2,b1)*delta_I(a1,bc1,a1,bc2)
          transb_B_i = transb_B_i +  rho_bb_i(b2,b1)*delta_Q(a1,bc1,a1,bc2)
          transb_C_i = transb_C_i +  rho_bb_i(b2,b1)*delta_V(a1,bc1,a1,bc2)
          transb_F_i = transb_F_i +  rho_bb_i(b2,b1)*delta_U(a1,bc1,a1,bc2)
        else
          transb_A_i = transb_A_i +  rho_bb_r(b2,b1)*delta_I(a1,bc1,a1,bc2)
          transb_B_i = transb_B_i +  rho_bb_r(b2,b1)*delta_Q(a1,bc1,a1,bc2)
          transb_C_i = transb_C_i +  rho_bb_r(b2,b1)*delta_V(a1,bc1,a1,bc2)
          transb_F_i = transb_F_i +  rho_bb_r(b2,b1)*delta_U(a1,bc1,a1,bc2)
  
          transb_A_r = transb_A_r -  rho_bb_i(b2,b1)*delta_I(a1,bc1,a1,bc2)
          transb_B_r = transb_B_r -  rho_bb_i(b2,b1)*delta_Q(a1,bc1,a1,bc2)
          transb_C_r = transb_C_r -  rho_bb_i(b2,b1)*delta_V(a1,bc1,a1,bc2)
          transb_F_r = transb_F_r -  rho_bb_i(b2,b1)*delta_U(a1,bc1,a1,bc2)
        endif

      endif

    end do
  end do
end do


Aw = -Const*(transb_A_r - transa_A_r)
Bw = Const*(transb_B_r - transa_B_r)
Cw = -Const*(transb_C_r - transa_C_r)
Dw = 0.d0 
Ew = 0.d0 
Fw = -Const*(transb_F_i - transa_F_i)
Gw = 0.d0 

Aws = Constp*Const*(transa_A_r)
Bws = -Constp*Const*(transa_B_r)
Cws = Constp*Const*(transa_C_r)
Fws = Constp*Const*(transa_F_i)


end subroutine
