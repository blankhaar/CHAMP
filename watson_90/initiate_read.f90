subroutine initiate(lF1,lF2,F1,F2,w0,alp_Z1,alp_Z2,B,wa,wb,gam_a,gam_b,&
&Iw,Uw,Qw,Vw,Tb,lam1,del,lbd_a_r,lbd_a_i,lbd_b_r,lbd_b_i,anis,apump,p_an_u,p_an_d,inci,qin)
implicit none
integer a1,a2,b1,b2,ma1,mb1,F1,F2,lF1,lF2,I,n,nw
integer anis,inci
double precision wab,w0,mass,dw,x0,x1,x2,x3,x4,x8,pi,w,c0,Tb,vthermal,B,h
double precision conv,lam2,k_sb,Tmp,vel
double precision Iw,Uw,Qw,Vw,qin,alp,bet,gam
double precision alp_Z1,alp_Z2,lam1,del,gam_a,gam_b
double precision p_an_u,p_an_d,lbd_t_r,lbd_t_i
double precision, dimension (lF1) :: wa
double precision, dimension (lF2) :: wb
double precision, dimension (lF1,lF1) :: lbd_a_r,lbd_a_i
double precision, dimension (lF2,lF2) :: lbd_b_r,lbd_b_i
double precision, dimension (3) :: apump

x0 = 0.d0
x1 = 1.d0
x2 = 2.d0
x3 = 3.d0
x4 = 4.d0
x8 = 8.d0

c0 = 299792458.d0
k_sb = 1.38064852E-23 
pi = dacos(-x1)
h = 6.626070040E-34

if (anis.eq.1) then
  call twovec_euler(apump,alp,bet,gam)
endif



do a1=1,lF1
  ma1 = a1 - F1 - 1
  wa(a1) = B*alp_Z1*ma1
end do

do b1=1,lF2
  mb1 = b1 - F2 - 1
  wb(b1) = B*alp_Z2*mb1
end do


!initiate density matrix

lam2 = lam1 * (x1 - del/x2) / (x1 + del / x2)

do a1 = 1,lF1
  do a2 = 1,lF1

    if (anis.eq.1) then
      call anis_pump_nw(lam1,p_an_u,alp,bet,F1,lF1,a1,a2,lbd_t_r,lbd_t_i)
      lbd_a_r(a1,a2) = lbd_t_r
      lbd_a_i(a1,a2) = lbd_t_i
    else
      if (a1.eq.a2) then
        lbd_a_r(a1,a2) = lam1
        lbd_a_i(a1,a2) = x0 
      else
        lbd_a_r(a1,a2) = x0 
        lbd_a_i(a1,a2) = x0 
      endif
    endif

  end do
end do

do b1 = 1,lF2
  do b2 = 1,lF2

    if (anis.eq.1) then
      call anis_pump_nw(lam2,p_an_d,alp,bet,F2,lF2,b1,b2,lbd_t_r,lbd_t_i)

      lbd_b_r(b1,b2) = lbd_t_r
      lbd_b_i(b1,b2) = lbd_t_i
    else
      if (b1.eq.b2) then
        lbd_b_r(b1,b2) = lam2
        lbd_b_i(b1,b2) = x0 
      else
        lbd_b_r(b1,b2) = x0 
        lbd_b_i(b1,b2) = x0 
      endif
    endif

  end do
end do

!initiate incoming radiation

w = w0 
Iw = h*w**x3/(x4*pi**x3*c0*c0*(dexp(h*w/(x2*pi*k_sb*Tb)) - x1))
Qw = x0 
Vw = x0
if (inci.eq.1) then
  Uw = qin * Iw 
else 
  Uw = x0
endif

end subroutine


