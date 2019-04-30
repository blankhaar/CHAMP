subroutine stimulate(I0,Q0,U0,V0,lF1,lF2,Aij,w0,delta_I,delta_Q,delta_U,delta_V,Rtot,Tb)
implicit none
integer a1,a2,b1,b2,bc1,br2,ema1,ma2,lF1,lF2,F1,F2
double precision x1,x2,pi,h,hbar,c0,I0,Q0,U0,V0
double precision Const,Rtot,Rtot_r,Rtot_i
double precision kb,Tb,Aij,Const2,w0

double precision, dimension (lF1,3,lF1,3) :: delta_I,delta_Q,delta_U,delta_V

! N & W 1990, 354

c0 = 299792458.d0
h  = 6.626070040E-34

x1 = 1.d0
x2 = 2.d0
pi = dacos(-x1)
hbar = h/(x2*pi)

F1 = (lF1 - 1)/2
F2 = (lF2 - 1)/2

Const = x2*pi*pi/((x1*lF1)*c0*hbar*hbar)

Rtot_r = 0.d0
Rtot_i = 0.d0

do a1 = 1,lF1
  do bc1 = 1,3

    Rtot_r = Rtot_r + I0*delta_I(a1,bc1,a1,bc1) + V0*delta_V(a1,bc1,a1,bc1) - Q0*delta_Q(a1,bc1,a1,bc1) 
 
    Rtot_i = Rtot_i + U0*delta_U(a1,bc1,a1,bc1) 

  end do
end do

Rtot = Const*Rtot_r

! N & W 1992, 384 Eq. (31)
kb = 1.38064852E-23

Const2 = Aij*kb/(x2*h*w0) 

Tb = Rtot/Const2 


end subroutine


! R = Aij*k_B*Tb*sa/(h*x2*w)

! R = Aij * pi**x2 * c0**x2 * I0 / ( w**x3 * h )

subroutine stimulate_alt(I0,Aij,sa,w0,Rtot)
implicit none
double precision x1,x2,x3,pi,h,c0
double precision I0,Aij,sa,w0,Rtot

x1 = 1.d0
x2 = 2.d0
x3 = 3.d0

pi = dacos(-x1)

c0 = 299792458.d0
h  = 6.626070040E-34


! N & W 1992, 384 Eq. (31)
Rtot = Aij * x2 * pi**x3 * c0**x2 * sa * I0 / ( w0**x3 * h )

end subroutine



