subroutine int_stokes(Iw,Qw,Vw,lF1,lF2,a1,b1,nw,w0,dw,gam,wa,wb,Ipv_r,Qpv_r,Vpv_r)
implicit none
integer n,m,nw,a1,b1,lF1,lF2,md,mu
double precision pi,gam,gZeeman,dw,dwt,w0,x0
double precision dI,dQ,dU,dV,Is,Qs,Us,Vs
double precision ddE 

double precision, dimension (-nw:nw) :: Iw,Qw,Vw
double precision, dimension (-nw:nw) :: Ipv_r,Qpv_r,Vpv_r
double precision, dimension (lF1) :: wa 
double precision, dimension (lF2) :: wb 

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
pi = dacos(-1.d0)

do n=-nw,nw
  Ipv_r(n) = x0 
  Qpv_r(n) = x0 
  Vpv_r(n) = x0 
end do


dwt = 2.d0*dw ! delta w
ddE = wa(a1) - wb(b1) 

!see Eq. (15) and (18) in N&W 1992 

do n=-nw,nw

  if (n.eq.-nw) then
    dI = (Iw(-nw+1) - Iw(-nw))/dw
    dQ = (Qw(-nw+1) - Qw(-nw))/dw
    dV = (Vw(-nw+1) - Vw(-nw))/dw
  elseif (n.eq.nw) then
    dI = (Iw(nw) - Iw(nw-1))/dw
    dQ = (Qw(nw) - Qw(nw-1))/dw
    dV = (Vw(nw) - Vw(nw-1))/dw
  else 
    dI = (Iw(n+1) - Iw(n-1))/dwt
    dQ = (Qw(n+1) - Qw(n-1))/dwt
    dV = (Vw(n+1) - Vw(n-1))/dwt
  endif

  Is = Iw(n)
  Qs = Qw(n)
  Vs = Vw(n)

  Ipv_r(n) = pi * (Is + dI * ddE) 
  Qpv_r(n) = pi * (Qs + dQ * ddE)
  Vpv_r(n) = pi * (Vs + dV * ddE)
 
end do

end subroutine

