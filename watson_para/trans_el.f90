subroutine prop_stokes(A,B,C,D,E,F,G,I,Q,U,V,ds,nv)
implicit none
integer n,nv
double precision Inew,Qnew,Unew,Vnew,ds                
double precision, dimension (4,4) :: eK,eeK 
double precision, dimension (-nv:nv) :: A,B,C,D,E,F,G
double precision, dimension (-nv:nv) :: I,Q,U,V 

! propagating the stokes parameters
! all variables are real

!input
!A,B,C,D,E,F,G: the propagation coefficients
!I,Q,U,V: the old stokes parameters

!output
!I,Q,U,V: the new stokes parameters

do n=-nv,nv
!  Inew = A(n)*I(n) + B(n)*Q(n) + F(n)*U(n) + C(n)*V(n)  
!  Qnew = B(n)*I(n) + A(n)*Q(n) + E(n)*U(n) + G(n)*V(n)
!  Unew = F(n)*I(n) - E(n)*Q(n) + A(n)*U(n) + D(n)*V(n)
!  Vnew = C(n)*I(n) - G(n)*Q(n) - D(n)*U(n) + A(n)*V(n)
!  
!  I(n) = I(n) + Inew*ds 
!  Q(n) = Q(n) + Qnew*ds
!  U(n) = U(n) + Unew*ds
!  V(n) = V(n) + Vnew*ds

  eK(1,1) = A(n)*ds 
  eK(2,2) = A(n)*ds 
  eK(3,3) = A(n)*ds 
  eK(4,4) = A(n)*ds 

  eK(1,2) = B(n)*ds 
  eK(2,1) = B(n)*ds 

  eK(1,3) = F(n)*ds 
  eK(3,1) = F(n)*ds 

  eK(1,4) = C(n)*ds 
  eK(4,1) = C(n)*ds 

  eK(2,3) = E(n)*ds 
  eK(3,2) = -E(n)*ds 

  eK(2,4) = G(n)*ds
  eK(4,2) = -G(n)*ds

  eK(3,4) = D(n)*ds
  eK(4,3) = -D(n)*ds

  call exp_A(eK,8,eeK)

  I(n) = eeK(1,1)*I(n) + eeK(1,2)*Q(n) + eeK(1,3)*U(n) + eeK(1,4)*V(n)
  Q(n) = eeK(2,1)*I(n) + eeK(2,2)*Q(n) + eeK(2,3)*U(n) + eeK(2,4)*V(n)
  U(n) = eeK(3,1)*I(n) + eeK(3,2)*Q(n) + eeK(3,3)*U(n) + eeK(3,4)*V(n)
  V(n) = eeK(4,1)*I(n) + eeK(4,2)*Q(n) + eeK(4,3)*U(n) + eeK(4,4)*V(n)

end do

end subroutine


subroutine trans_el(gaa_p_r,gaa_p_i,gbb_p_r,gbb_p_i,gaa_m_r,gaa_m_i,gbb_m_r,gbb_m_i,lF1,lF2,nw,&
&delta_I,delta_Q,delta_U,delta_V,Aw,Bw,Cw,Dw,Ew,Fw,Gw,w0)
implicit none
integer n,m,nw,a1,a2,b1,b2,ac1,ac2,bc1,bc2
integer F1,F2,lF1,lF2,ma1,ma2,mb1,mb2
double precision transa_A_r,transa_B_r,transa_C_r,transa_D_r,transa_E_r,transa_F_r,transa_G_r 
double precision transb_A_r,transb_B_r,transb_C_r,transb_D_r,transb_E_r,transb_F_r,transb_G_r 
double precision transa_A_i,transa_B_i,transa_C_i,transa_D_i,transa_E_i,transa_F_i,transa_G_i 
double precision transb_A_i,transb_B_i,transb_C_i,transb_D_i,transb_E_i,transb_F_i,transb_G_i 
double precision Const,pi,c0,hbar,w0 
double precision Ar,Br,Cr,Dr,Er,Fr,Gr,Ai,Bi,Ci,Di,Ei,Fi,Gi

double precision, dimension (lF2,3,3,-nw:nw) :: gaa_p_r,gaa_p_i,gaa_m_r,gaa_m_i
double precision, dimension (lF1,3,3,-nw:nw) :: gbb_p_r,gbb_p_i,gbb_m_r,gbb_m_i
double precision, dimension (lF1,3,lF1,3) :: delta_I,delta_Q,delta_U,delta_V
double precision, dimension (-nw:nw) :: Aw,Bw,Cw,Dw,Ew,Fw,Gw

pi = dacos(-1.d0)
c0 = 299792458.d0
hbar = 1.0545718E-34
Const = pi*w0/(c0*hbar) 

F1 = (lF1 - 1)/2
F2 = (lF2 - 1)/2

!$OMP DO
do n=-nw,nw

! first the sum over ab, a'

  transa_A_r = 0.d0 
  transa_B_r = 0.d0 
  transa_C_r = 0.d0 
  transa_D_r = 0.d0 
  transa_E_r = 0.d0 
  transa_F_r = 0.d0 
  transa_G_r = 0.d0 

  transa_A_i = 0.d0 
  transa_B_i = 0.d0 
  transa_C_i = 0.d0 
  transa_D_i = 0.d0 
  transa_E_i = 0.d0 
  transa_F_i = 0.d0 
  transa_G_i = 0.d0 

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

! delta even or uneven
        if (mod((bc1+bc2),2).eq.0) then 
          transa_A_r = transa_A_r +  gaa_p_r(b1,ac1,ac2,n)*delta_I(a1,bc1,a2,bc2)
          transa_B_r = transa_B_r +  gaa_p_r(b1,ac1,ac2,n)*delta_Q(a1,bc1,a2,bc2)
          transa_C_r = transa_C_r +  gaa_p_r(b1,ac1,ac2,n)*delta_V(a1,bc1,a2,bc2)
          transa_D_r = transa_D_r +  gaa_m_r(b1,ac1,ac2,n)*delta_Q(a1,bc1,a2,bc2)
          transa_E_r = transa_E_r +  gaa_m_r(b1,ac1,ac2,n)*delta_V(a1,bc1,a2,bc2)
          transa_F_r = transa_F_r +  gaa_p_r(b1,ac1,ac2,n)*delta_U(a1,bc1,a2,bc2)
          transa_G_r = transa_G_r +  gaa_m_r(b1,ac1,ac2,n)*delta_U(a1,bc1,a2,bc2)

          transa_A_i = transa_A_i +  gaa_p_i(b1,ac1,ac2,n)*delta_I(a1,bc1,a2,bc2)
          transa_B_i = transa_B_i +  gaa_p_i(b1,ac1,ac2,n)*delta_Q(a1,bc1,a2,bc2)
          transa_C_i = transa_C_i +  gaa_p_i(b1,ac1,ac2,n)*delta_V(a1,bc1,a2,bc2)
          transa_D_i = transa_D_i +  gaa_m_i(b1,ac1,ac2,n)*delta_Q(a1,bc1,a2,bc2)
          transa_E_i = transa_E_i +  gaa_m_i(b1,ac1,ac2,n)*delta_V(a1,bc1,a2,bc2)
          transa_F_i = transa_F_i +  gaa_p_i(b1,ac1,ac2,n)*delta_U(a1,bc1,a2,bc2)
          transa_G_i = transa_G_i +  gaa_m_i(b1,ac1,ac2,n)*delta_U(a1,bc1,a2,bc2)
        else
          transa_A_i = transa_A_i +  gaa_p_r(b1,ac1,ac2,n)*delta_I(a1,bc1,a2,bc2)
          transa_B_i = transa_B_i +  gaa_p_r(b1,ac1,ac2,n)*delta_Q(a1,bc1,a2,bc2)
          transa_C_i = transa_C_i +  gaa_p_r(b1,ac1,ac2,n)*delta_V(a1,bc1,a2,bc2)
          transa_D_i = transa_D_i +  gaa_m_r(b1,ac1,ac2,n)*delta_Q(a1,bc1,a2,bc2)
          transa_E_i = transa_E_i +  gaa_m_r(b1,ac1,ac2,n)*delta_V(a1,bc1,a2,bc2)
          transa_F_i = transa_F_i +  gaa_p_r(b1,ac1,ac2,n)*delta_U(a1,bc1,a2,bc2)
          transa_G_i = transa_G_i +  gaa_m_r(b1,ac1,ac2,n)*delta_U(a1,bc1,a2,bc2)

          transa_A_r = transa_A_r -  gaa_p_i(b1,ac1,ac2,n)*delta_I(a1,bc1,a2,bc2)
          transa_B_r = transa_B_r -  gaa_p_i(b1,ac1,ac2,n)*delta_Q(a1,bc1,a2,bc2)
          transa_C_r = transa_C_r -  gaa_p_i(b1,ac1,ac2,n)*delta_V(a1,bc1,a2,bc2)
          transa_D_r = transa_D_r -  gaa_m_i(b1,ac1,ac2,n)*delta_Q(a1,bc1,a2,bc2)
          transa_E_r = transa_E_r -  gaa_m_i(b1,ac1,ac2,n)*delta_V(a1,bc1,a2,bc2)
          transa_F_r = transa_F_r -  gaa_p_i(b1,ac1,ac2,n)*delta_U(a1,bc1,a2,bc2)
          transa_G_r = transa_G_r -  gaa_m_i(b1,ac1,ac2,n)*delta_U(a1,bc1,a2,bc2)
        endif 

      end do
    end do
  end do

  transb_A_r = 0.d0
  transb_B_r = 0.d0
  transb_C_r = 0.d0
  transb_D_r = 0.d0
  transb_E_r = 0.d0
  transb_F_r = 0.d0
  transb_G_r = 0.d0

  transb_A_i = 0.d0
  transb_B_i = 0.d0
  transb_C_i = 0.d0
  transb_D_i = 0.d0
  transb_E_i = 0.d0
  transb_F_i = 0.d0
  transb_G_i = 0.d0

  do a1 = 1,lF1
    do bc1 = 1,3
      do bc2 = 1,3
        ma1 = a1 - F1 - 1
        mb1 = ma1 + (bc1 - 2)  !dm = mb - ma
        mb2 = ma1 + (bc2 - 2)
  
        b1 = mb1 + F2 + 1
        b2 = mb2 + F2 + 1

! delta even or uneven
        if (mod((bc1+bc2),2).eq.0) then
          transb_A_r = transb_A_r +  gbb_p_r(a1,bc1,bc2,n)*delta_I(a1,bc1,a1,bc2)
          transb_B_r = transb_B_r +  gbb_p_r(a1,bc1,bc2,n)*delta_Q(a1,bc1,a1,bc2)
          transb_C_r = transb_C_r +  gbb_p_r(a1,bc1,bc2,n)*delta_V(a1,bc1,a1,bc2)
          transb_D_r = transb_D_r +  gbb_m_r(a1,bc1,bc2,n)*delta_Q(a1,bc1,a1,bc2)
          transb_E_r = transb_E_r +  gbb_m_r(a1,bc1,bc2,n)*delta_V(a1,bc1,a1,bc2)
          transb_F_r = transb_F_r +  gbb_p_r(a1,bc1,bc2,n)*delta_U(a1,bc1,a1,bc2)
          transb_G_r = transb_G_r +  gbb_m_r(a1,bc1,bc2,n)*delta_U(a1,bc1,a1,bc2)

          transb_A_i = transb_A_i +  gbb_p_i(a1,bc1,bc2,n)*delta_I(a1,bc1,a1,bc2)
          transb_B_i = transb_B_i +  gbb_p_i(a1,bc1,bc2,n)*delta_Q(a1,bc1,a1,bc2)
          transb_C_i = transb_C_i +  gbb_p_i(a1,bc1,bc2,n)*delta_V(a1,bc1,a1,bc2)
          transb_D_i = transb_D_i +  gbb_m_i(a1,bc1,bc2,n)*delta_Q(a1,bc1,a1,bc2)
          transb_E_i = transb_E_i +  gbb_m_i(a1,bc1,bc2,n)*delta_V(a1,bc1,a1,bc2)
          transb_F_i = transb_F_i +  gbb_p_i(a1,bc1,bc2,n)*delta_U(a1,bc1,a1,bc2)
          transb_G_i = transb_G_i +  gbb_m_i(a1,bc1,bc2,n)*delta_U(a1,bc1,a1,bc2)
        else
          transb_A_i = transb_A_i +  gbb_p_r(a1,bc1,bc2,n)*delta_I(a1,bc1,a1,bc2)
          transb_B_i = transb_B_i +  gbb_p_r(a1,bc1,bc2,n)*delta_Q(a1,bc1,a1,bc2)
          transb_C_i = transb_C_i +  gbb_p_r(a1,bc1,bc2,n)*delta_V(a1,bc1,a1,bc2)
          transb_D_i = transb_D_i +  gbb_m_r(a1,bc1,bc2,n)*delta_Q(a1,bc1,a1,bc2)
          transb_E_i = transb_E_i +  gbb_m_r(a1,bc1,bc2,n)*delta_V(a1,bc1,a1,bc2)
          transb_F_i = transb_F_i +  gbb_p_r(a1,bc1,bc2,n)*delta_U(a1,bc1,a1,bc2)
          transb_G_i = transb_G_i +  gbb_m_r(a1,bc1,bc2,n)*delta_U(a1,bc1,a1,bc2)

          transb_A_r = transb_A_r -  gbb_p_i(a1,bc1,bc2,n)*delta_I(a1,bc1,a1,bc2)
          transb_B_r = transb_B_r -  gbb_p_i(a1,bc1,bc2,n)*delta_Q(a1,bc1,a1,bc2)
          transb_C_r = transb_C_r -  gbb_p_i(a1,bc1,bc2,n)*delta_V(a1,bc1,a1,bc2)
          transb_D_r = transb_D_r -  gbb_m_i(a1,bc1,bc2,n)*delta_Q(a1,bc1,a1,bc2)
          transb_E_r = transb_E_r -  gbb_m_i(a1,bc1,bc2,n)*delta_V(a1,bc1,a1,bc2)
          transb_F_r = transb_F_r -  gbb_p_i(a1,bc1,bc2,n)*delta_U(a1,bc1,a1,bc2)
          transb_G_r = transb_G_r -  gbb_m_i(a1,bc1,bc2,n)*delta_U(a1,bc1,a1,bc2)
        endif

      end do
    end do
  end do


  Aw(n) = -Const*(transb_A_r - transa_A_r)
  Bw(n) = Const*(transb_B_r - transa_B_r)
  Cw(n) = -Const*(transb_C_r - transa_C_r)
  Dw(n) = Const*(transb_D_i - transa_D_i)
  Ew(n) = -Const*(transb_E_i - transa_E_i)
  Fw(n) = -Const*(transb_F_i - transa_F_i)
  Gw(n) = -Const*(transb_G_r - transa_G_r)

end do      
!$OMP END DO

end subroutine


subroutine int_densities(rho_aa_r,rho_aa_i,rho_bb_r,rho_bb_i,lF1,lF2,nv,w0,dv,gam_a,gam_b,wa,wb,&
&gaa_p_r,gaa_p_i,gaa_m_r,gaa_m_i,gbb_p_r,gbb_p_i,gbb_m_r,gbb_m_i)
implicit none
integer n,m,nv,a1,a2,b1,b2,ac1,ac2,bc1,bc2
integer F1,F2,lF1,lF2,ma1,ma2,mb1,mb2,mu,md
double precision gam_a,gam_b,dv,dvt,w0,c0,dw
double precision dr_aa_r,dr_bb_r,dr_aa_i,dr_bb_i 
double precision rs_aa_r,rs_bb_r,rs_aa_i,rs_bb_i 
double precision g_r,g_i,g_w_r,g_w_i
double precision gamab_r,gamab_i,gamab_w_r,gamapb_r,gamapb_i,gamapb_w_r,gamapb_w_i,gamab_w_i

double precision, dimension (lF2,3,3,-nv:nv) :: gaa_p_r,gaa_p_i,gaa_m_r,gaa_m_i 
double precision, dimension (lF1,3,3,-nv:nv) :: gbb_p_r,gbb_p_i,gbb_m_r,gbb_m_i
double precision, dimension (lF1,lF1,-nv:nv) :: rho_aa_r,rho_aa_i 
double precision, dimension (lF2,lF2,-nv:nv) :: rho_bb_r,rho_bb_i 
double precision, dimension (lF1) :: wa 
double precision, dimension (lF2) :: wb

! integrate the densities
! product is the gaa = <\rho_ \gam>-matrices for all the omega values
! 

dvt = 2.d0*dv 

F1 = (lF1 - 1)/2
F2 = (lF2 - 1)/2

do b1=1,lF2
  do ac1 = 1,3
    do ac2 = 1,3
      do n = -nv,nv
        gaa_p_r(b1,ac1,ac2,n) = 0.d0
        gaa_p_i(b1,ac1,ac2,n) = 0.d0
        gaa_m_r(b1,ac1,ac2,n) = 0.d0
        gaa_m_i(b1,ac1,ac2,n) = 0.d0
      end do
    end do
  end do
end do

do a1=1,lF1
  do bc1 = 1,3
    do bc2 = 1,3
      do n = -nv,nv
        gbb_p_r(a1,bc1,bc2,n) = 0.d0
        gbb_p_i(a1,bc1,bc2,n) = 0.d0
        gbb_m_r(a1,bc1,bc2,n) = 0.d0
        gbb_m_i(a1,bc1,bc2,n) = 0.d0
      end do
    end do
  end do
end do


do b1 = 1,lF2
  do ac1 = 1,3
    do ac2 = 1,3
      mb1 = b1 - F2 - 1 
  
      ma1 = mb1 - (ac1 - 2)
      ma2 = mb1 - (ac2 - 2)

      a1 = ma1 + F1 + 1
      a2 = ma2 + F1 + 1
 
      if (a1.le.lF1.and.a2.le.lF1.and.a1.gt.0.and.a2.gt.0) then
  
        !$OMP DO
        do n=-nv,nv
    
    !   first element
          m = -nv
          dr_aa_r = (rho_aa_r(a1,a2,-nv+1) - rho_aa_r(a1,a2,-nv))/dv
          dr_aa_i = (rho_aa_i(a1,a2,-nv+1) - rho_aa_i(a1,a2,-nv))/dv
    
          rs_aa_r = rho_aa_r(a1,a2,m)
          rs_aa_i = rho_aa_i(a1,a2,m)

          call int_gam_v(gam_a,wa,wb,w0,lF1,lF2,dv,a1,b1,n,m,gamab_r,gamab_i,gamab_w_r,gamab_w_i) !function for eqs (A.13-16)
          call int_gam_v(gam_a,wa,wb,w0,lF1,lF2,dv,a2,b1,n,m,gamapb_r,gamapb_i,gamapb_w_r,gamapb_w_i) !function for eqs (A.13-16)

          g_r = gamab_r + gamapb_r
          g_i = gamab_i - gamapb_i

          g_w_r = gamab_w_r + gamapb_w_r
          g_w_i = gamab_w_i - gamapb_w_i
    
          gaa_p_r(b1,ac1,ac2,n) = gaa_p_r(b1,ac1,ac2,n) + g_r*rs_aa_r - g_i*rs_aa_i + g_w_r*dr_aa_r - g_w_i*dr_aa_i
          gaa_p_i(b1,ac1,ac2,n) = gaa_p_i(b1,ac1,ac2,n) + g_i*rs_aa_r + g_r*rs_aa_i + g_w_r*dr_aa_i + g_w_i*dr_aa_r

          g_r = gamab_r - gamapb_r
          g_i = gamab_i + gamapb_i

          g_w_r = gamab_w_r - gamapb_w_r
          g_w_i = gamab_w_i + gamapb_w_i

          gaa_m_r(b1,ac1,ac2,n) = gaa_m_r(b1,ac1,ac2,n) + g_r*rs_aa_r - g_i*rs_aa_i + g_w_r*dr_aa_r - g_w_i*dr_aa_i 
          gaa_m_i(b1,ac1,ac2,n) = gaa_m_i(b1,ac1,ac2,n) + g_i*rs_aa_r + g_r*rs_aa_i + g_w_r*dr_aa_i + g_w_i*dr_aa_r
    
  !   last element
          m = nv
          dr_aa_r = (rho_aa_r(a1,a2,m) - rho_aa_r(a1,a2,m-1))/dv
          dr_aa_i = (rho_aa_i(a1,a2,m) - rho_aa_i(a1,a2,m-1))/dv
    
          rs_aa_r = rho_aa_r(a1,a2,m)
          rs_aa_i = rho_aa_i(a1,a2,m)
  
          call int_gam_v(gam_a,wa,wb,w0,lF1,lF2,dv,a1,b1,n,m,gamab_r,gamab_i,gamab_w_r,gamab_w_i) !function for eqs (A.13-16)
          call int_gam_v(gam_a,wa,wb,w0,lF1,lF2,dv,a2,b1,n,m,gamapb_r,gamapb_i,gamapb_w_r,gamapb_w_i) !function for eqs (A.13-16)

          g_r = gamab_r + gamapb_r
          g_i = gamab_i - gamapb_i

          g_w_r = gamab_w_r + gamapb_w_r
          g_w_i = gamab_w_i - gamapb_w_i

          gaa_p_r(b1,ac1,ac2,n) = gaa_p_r(b1,ac1,ac2,n) + g_r*rs_aa_r - g_i*rs_aa_i + g_w_r*dr_aa_r - g_w_i*dr_aa_i
          gaa_p_i(b1,ac1,ac2,n) = gaa_p_i(b1,ac1,ac2,n) + g_i*rs_aa_r + g_r*rs_aa_i + g_w_r*dr_aa_i + g_w_i*dr_aa_r

          g_r = gamab_r - gamapb_r
          g_i = gamab_i + gamapb_i

          g_w_r = gamab_w_r - gamapb_w_r
          g_w_i = gamab_w_i + gamapb_w_i

          gaa_m_r(b1,ac1,ac2,n) = gaa_m_r(b1,ac1,ac2,n) + g_r*rs_aa_r - g_i*rs_aa_i + g_w_r*dr_aa_r - g_w_i*dr_aa_i
          gaa_m_i(b1,ac1,ac2,n) = gaa_m_i(b1,ac1,ac2,n) + g_i*rs_aa_r + g_r*rs_aa_i + g_w_r*dr_aa_i + g_w_i*dr_aa_r

  
    !   rest of elements
          do m=-nv+1,nv-1
  
            mu = m + 1
            md = m - 1
            dr_aa_r = (rho_aa_r(a1,a2,mu) - rho_aa_r(a1,a2,md))/dvt
            dr_aa_i = (rho_aa_i(a1,a2,mu) - rho_aa_i(a1,a2,md))/dvt
  
            rs_aa_r = rho_aa_r(a1,a2,m)
            rs_aa_i = rho_aa_i(a1,a2,m)
  
            call int_gam_v(gam_a,wa,wb,w0,lF1,lF2,dv,a1,b1,n,m,gamab_r,gamab_i,gamab_w_r,gamab_w_i) !function for eqs (A.13-16)
            call int_gam_v(gam_a,wa,wb,w0,lF1,lF2,dv,a2,b1,n,m,gamapb_r,gamapb_i,gamapb_w_r,gamapb_w_i) !function for eqs (A.13-16)
   
            g_r = gamab_r + gamapb_r
            g_i = gamab_i - gamapb_i
  
            g_w_r = gamab_w_r + gamapb_w_r
            g_w_i = gamab_w_i - gamapb_w_i
  
            gaa_p_r(b1,ac1,ac2,n) = gaa_p_r(b1,ac1,ac2,n) + g_r*rs_aa_r - g_i*rs_aa_i + g_w_r*dr_aa_r - g_w_i*dr_aa_i
            gaa_p_i(b1,ac1,ac2,n) = gaa_p_i(b1,ac1,ac2,n) + g_i*rs_aa_r + g_r*rs_aa_i + g_w_r*dr_aa_i + g_w_i*dr_aa_r
 
            g_r = gamab_r - gamapb_r
            g_i = gamab_i + gamapb_i
  
            g_w_r = gamab_w_r - gamapb_w_r
            g_w_i = gamab_w_i + gamapb_w_i
  
            gaa_m_r(b1,ac1,ac2,n) = gaa_m_r(b1,ac1,ac2,n) + g_r*rs_aa_r - g_i*rs_aa_i + g_w_r*dr_aa_r - g_w_i*dr_aa_i
            gaa_m_i(b1,ac1,ac2,n) = gaa_m_i(b1,ac1,ac2,n) + g_i*rs_aa_r + g_r*rs_aa_i + g_w_r*dr_aa_i + g_w_i*dr_aa_r
  
          end do
          
        end do
        !$OMP END DO

      endif
    end do
  end do
end do

!stop

do a1 = 1,lF1
  do bc1 = 1,3
    do bc2 = 1,3
      ma1 = a1 - F1 - 1
      mb1 = ma1 + (bc1 - 2)
      mb2 = ma1 + (bc2 - 2)
  
      b1 = mb1 + F2 + 1 
      b2 = mb2 + F2 + 1

      if (b1.le.lF2.and.b2.le.lF2.and.b1.gt.0.and.b2.gt.0) then

        !$OMP DO  
        do n=-nv,nv
    
    !   first element
          m = -nv
          dr_bb_r = (rho_bb_r(b2,b1,-nv+1) - rho_bb_r(b2,b1,-nv))/dv
          dr_bb_i = (rho_bb_i(b2,b1,-nv+1) - rho_bb_i(b2,b1,-nv))/dv
    
          rs_bb_r = rho_bb_r(b2,b1,m)
          rs_bb_i = rho_bb_i(b2,b1,m)
    
          call int_gam_v(gam_b,wa,wb,w0,lF1,lF2,dv,a1,b1,n,m,gamab_r,gamab_i,gamab_w_r,gamab_w_i) !function for eqs (A.13-16)
          call int_gam_v(gam_b,wa,wb,w0,lF1,lF2,dv,a1,b2,n,m,gamapb_r,gamapb_i,gamapb_w_r,gamapb_w_i) !function for eqs (A.13-16)

          g_r = gamab_r + gamapb_r
          g_i = gamab_i - gamapb_i

          g_w_r = gamab_w_r + gamapb_w_r
          g_w_i = gamab_w_i - gamapb_w_i
  
          gbb_p_r(a1,bc1,bc2,n) = gbb_p_r(a1,bc1,bc2,n) + g_r*rs_bb_r - g_i*rs_bb_i + g_w_r*dr_bb_r - g_w_i*dr_bb_i
          gbb_p_i(a1,bc1,bc2,n) = gbb_p_i(a1,bc1,bc2,n) + g_i*rs_bb_r + g_r*rs_bb_i + g_w_r*dr_bb_i + g_w_i*dr_bb_r
 
          g_r = gamab_r - gamapb_r
          g_i = gamab_i + gamapb_i

          g_w_r = gamab_w_r - gamapb_w_r
          g_w_i = gamab_w_i + gamapb_w_i

          gbb_m_r(a1,bc1,bc2,n) = gbb_m_r(a1,bc1,bc2,n) + g_r*rs_bb_r - g_i*rs_bb_i + g_w_r*dr_bb_r - g_w_i*dr_bb_i
          gbb_m_i(a1,bc1,bc2,n) = gbb_m_i(a1,bc1,bc2,n) + g_i*rs_bb_r + g_r*rs_bb_i + g_w_r*dr_bb_i + g_w_i*dr_bb_r
 
    
    !   last element
          m = nv
          dr_bb_r = (rho_bb_r(b2,b1,m) - rho_bb_r(b2,b1,m-1))/dv
          dr_bb_i = (rho_bb_i(b2,b1,m) - rho_bb_i(b2,b1,m-1))/dv
    
          rs_bb_r = rho_bb_r(b2,b1,m)
          rs_bb_i = rho_bb_i(b2,b1,m)
  
          call int_gam_v(gam_b,wa,wb,w0,lF1,lF2,dv,a1,b1,n,m,gamab_r,gamab_i,gamab_w_r,gamab_w_i) !function for eqs (A.13-16)
          call int_gam_v(gam_b,wa,wb,w0,lF1,lF2,dv,a1,b2,n,m,gamapb_r,gamapb_i,gamapb_w_r,gamapb_w_i) !function for eqs (A.13-16)

          g_r = gamab_r + gamapb_r
          g_i = gamab_i - gamapb_i

          g_w_r = gamab_w_r + gamapb_w_r
          g_w_i = gamab_w_i - gamapb_w_i

          gbb_p_r(a1,bc1,bc2,n) = gbb_p_r(a1,bc1,bc2,n) + g_r*rs_bb_r - g_i*rs_bb_i + g_w_r*dr_bb_r - g_w_i*dr_bb_i
          gbb_p_i(a1,bc1,bc2,n) = gbb_p_i(a1,bc1,bc2,n) + g_i*rs_bb_r + g_r*rs_bb_i + g_w_r*dr_bb_i + g_w_i*dr_bb_r

          g_r = gamab_r - gamapb_r
          g_i = gamab_i + gamapb_i
    
          g_w_r = gamab_w_r - gamapb_w_r
          g_w_i = gamab_w_i + gamapb_w_i 
          
          gbb_m_r(a1,bc1,bc2,n) = gbb_m_r(a1,bc1,bc2,n) + g_r*rs_bb_r - g_i*rs_bb_i + g_w_r*dr_bb_r - g_w_i*dr_bb_i
          gbb_m_i(a1,bc1,bc2,n) = gbb_m_i(a1,bc1,bc2,n) + g_i*rs_bb_r + g_r*rs_bb_i + g_w_r*dr_bb_i + g_w_i*dr_bb_r

    !   rest of elements
          do m=-nv+1,nv-1
            mu = m + 1
            md = m - 1
            dr_bb_r = (rho_bb_r(b2,b1,mu) - rho_bb_r(b2,b1,md))/dvt
            dr_bb_i = (rho_bb_i(b2,b1,mu) - rho_bb_i(b2,b1,md))/dvt
  
            rs_bb_r = rho_bb_r(b2,b1,m)
            rs_bb_i = rho_bb_i(b2,b1,m)
    
            call int_gam_v(gam_b,wa,wb,w0,lF1,lF2,dv,a1,b1,n,m,gamab_r,gamab_i,gamab_w_r,gamab_w_i) !function for eqs (A.13-16)
            call int_gam_v(gam_b,wa,wb,w0,lF1,lF2,dv,a1,b2,n,m,gamapb_r,gamapb_i,gamapb_w_r,gamapb_w_i) !function for eqs (A.13-16)
  
            g_r = gamab_r + gamapb_r
            g_i = gamab_i - gamapb_i
  
            g_w_r = gamab_w_r + gamapb_w_r
            g_w_i = gamab_w_i - gamapb_w_i
  
            gbb_p_r(a1,bc1,bc2,n) = gbb_p_r(a1,bc1,bc2,n) + g_r*rs_bb_r - g_i*rs_bb_i + g_w_r*dr_bb_r - g_w_i*dr_bb_i
            gbb_p_i(a1,bc1,bc2,n) = gbb_p_i(a1,bc1,bc2,n) + g_i*rs_bb_r + g_r*rs_bb_i + g_w_r*dr_bb_i + g_w_i*dr_bb_r
  
            g_r = gamab_r - gamapb_r
            g_i = gamab_i + gamapb_i
      
            g_w_r = gamab_w_r - gamapb_w_r
            g_w_i = gamab_w_i + gamapb_w_i 
            
            gbb_m_r(a1,bc1,bc2,n) = gbb_m_r(a1,bc1,bc2,n) + g_r*rs_bb_r - g_i*rs_bb_i + g_w_r*dr_bb_r - g_w_i*dr_bb_i
            gbb_m_i(a1,bc1,bc2,n) = gbb_m_i(a1,bc1,bc2,n) + g_i*rs_bb_r + g_r*rs_bb_i + g_w_r*dr_bb_i + g_w_i*dr_bb_r
 
          end do
        end do
        !$OMP END DO

      endif
    end do
  end do
end do

end subroutine

