subroutine int_stokes(Iw,Qw,Uw,Vw,lF1,lF2,a1,b1,nw,w0,dw,gam,wa,wb,Ipv_r,Ipv_i,Qpv_r,Qpv_i,Upv_r,Upv_i,Vpv_r,Vpv_i)
implicit none
integer n,m,nw,a1,b1,lF1,lF2,md,mu
double precision gam,gZeeman,dw,dwt,w0,x0
double precision dI,dQ,dU,dV,Is,Qs,Us,Vs
double precision gamab_r,gamab_i,gamab_w_r,gamab_w_i

double precision, dimension (-nw:nw) :: Iw,Qw,Uw,Vw
double precision, dimension (-nw:nw) :: Ipv_r,Qpv_r,Upv_r,Vpv_r
double precision, dimension (-nw:nw) :: Ipv_i,Qpv_i,Upv_i,Vpv_i
double precision, dimension (lF1) :: wa 
double precision, dimension (lF2) :: wb 

! computes:
! <\gamma_+^{a1,b1}(v_j,w) I(w)>_w = I_{a1,b1}(v_j)
! <\gamma_+^{a1,b1}(v_j,w) Q(w)>_w = Q_{a1,b1}(v_j)
! <\gamma_+^{a1,b1}(v_j,w) U(w)>_w = U_{a1,b1}(v_j)
! <\gamma_+^{a1,b1}(v_j,w) V(w)>_w = V_{a1,b1}(v_j)
! to later use in the density equations 

!input
!Stokes, frequency dependent
!a1,b1: element for gamma 
!wa,wb: frequency shifts for the a1 and b1 sublevels 
!gam: Gamma decay (= (Gam_a + Gam_b)/2
!nw: size frequency/velocity grid

!output:
!RE ( I_{a1,b1}(v_j) ) = Ipv_r
!IM ( I_{a1,b1}(v_j) ) = Ipv_i
!RE ( Q_{a1,b1}(v_j) ) = Qpv_r
!IM ( Q_{a1,b1}(v_j) ) = Qpv_i
!etc...

x0 = 0.d0

do n=-nw,nw
  Ipv_r(n) = x0 
  Ipv_i(n) = x0 
  Qpv_r(n) = x0 
  Qpv_i(n) = x0 
  Upv_r(n) = x0 
  Upv_i(n) = x0 
  Vpv_r(n) = x0 
  Vpv_i(n) = x0 
end do


dwt = 2.d0*dw ! delta w

!see Eq. (...) in watson.pdf
!<\gamma_+^{a1,b1}(v_n,w) I(w)>_w = 
!\sum_{m=-nw}^nw I(w_m) * <\gamma_+^{ab}(v_n,w_m)>_w 
!+ \sum_{m=-nw}^nw (dI/dw)_{w_m} * <w_m*\gamma_+^{ab}(v_n,w_m)>_w

do n=-nw,nw

! first element
  m = -nw
  dI = (Iw(-nw+1) - Iw(-nw))/dw
  dQ = (Qw(-nw+1) - Qw(-nw))/dw
  dU = (Uw(-nw+1) - Uw(-nw))/dw
  dV = (Vw(-nw+1) - Vw(-nw))/dw

  Is = Iw(m)
  Qs = Qw(m)
  Us = Uw(m)
  Vs = Vw(m)

  call int_gam_w(gam,wa,wb,w0,lF1,lF2,dw,a1,b1,n,m,gamab_r,gamab_i,gamab_w_r,gamab_w_i) !function for eqs (A.13-16)

  Ipv_r(n) = Ipv_r(n) + gamab_r*Is + gamab_w_r*dI
  Qpv_r(n) = Qpv_r(n) + gamab_r*Qs + gamab_w_r*dQ
  Upv_r(n) = Upv_r(n) + gamab_r*Us + gamab_w_r*dU
  Vpv_r(n) = Vpv_r(n) + gamab_r*Vs + gamab_w_r*dV
                     
  Ipv_i(n) = Ipv_i(n) + gamab_i*Is + gamab_w_i*dI
  Qpv_i(n) = Qpv_i(n) + gamab_i*Qs + gamab_w_i*dQ
  Upv_i(n) = Upv_i(n) + gamab_i*Us + gamab_w_i*dU
  Vpv_i(n) = Vpv_i(n) + gamab_i*Vs + gamab_w_i*dV

! last element
  m = nw
  dI = (Iw(nw) - Iw(nw-1))/dw
  dQ = (Qw(nw) - Qw(nw-1))/dw
  dU = (Uw(nw) - Uw(nw-1))/dw
  dV = (Vw(nw) - Vw(nw-1))/dw

  Is = Iw(m)
  Qs = Qw(m)
  Us = Uw(m)
  Vs = Vw(m)

  call int_gam_w(gam,wa,wb,w0,lF1,lF2,dw,a1,b1,n,m,gamab_r,gamab_i,gamab_w_r,gamab_w_i) 


  Ipv_r(n) = Ipv_r(n) + gamab_r*Is + gamab_w_r*dI
  Qpv_r(n) = Qpv_r(n) + gamab_r*Qs + gamab_w_r*dQ
  Upv_r(n) = Upv_r(n) + gamab_r*Us + gamab_w_r*dU
  Vpv_r(n) = Vpv_r(n) + gamab_r*Vs + gamab_w_r*dV

  Ipv_i(n) = Ipv_i(n) + gamab_i*Is + gamab_w_i*dI
  Qpv_i(n) = Qpv_i(n) + gamab_i*Qs + gamab_w_i*dQ
  Upv_i(n) = Upv_i(n) + gamab_i*Us + gamab_w_i*dU
  Vpv_i(n) = Vpv_i(n) + gamab_i*Vs + gamab_w_i*dV

! rest of elements
  do m=-nw+1,nw-1
    mu = m + 1
    md = m - 1
    dI = (Iw(mu) - Iw(md))/dwt 
    dQ = (Qw(mu) - Qw(md))/dwt 
    dU = (Uw(mu) - Uw(md))/dwt 
    dV = (Vw(mu) - Vw(md))/dwt 
    
    Is = Iw(m) 
    Qs = Qw(m) 
    Us = Uw(m) 
    Vs = Vw(m) 

    call int_gam_w(gam,wa,wb,w0,lF1,lF2,dw,a1,b1,n,m,gamab_r,gamab_i,gamab_w_r,gamab_w_i) !function for eqs (A.13-16)
  
    Ipv_r(n) = Ipv_r(n) + gamab_r*Is + gamab_w_r*dI
    Qpv_r(n) = Qpv_r(n) + gamab_r*Qs + gamab_w_r*dQ
    Upv_r(n) = Upv_r(n) + gamab_r*Us + gamab_w_r*dU
    Vpv_r(n) = Vpv_r(n) + gamab_r*Vs + gamab_w_r*dV
  
    Ipv_i(n) = Ipv_i(n) + gamab_i*Is + gamab_w_i*dI
    Qpv_i(n) = Qpv_i(n) + gamab_i*Qs + gamab_w_i*dQ
    Upv_i(n) = Upv_i(n) + gamab_i*Us + gamab_w_i*dU
    Vpv_i(n) = Vpv_i(n) + gamab_i*Vs + gamab_w_i*dV

  end do
end do

end subroutine

