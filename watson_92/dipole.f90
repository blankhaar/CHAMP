subroutine gen_delta(theta,Aij,w0,F1,lF1,F2,delta_I,delta_Q,delta_U,delta_V)
implicit none

integer F1,F2,lF1,a1,a2,I,J
double precision x0,x1,x2,xs2,fac1,fac2,fac3
double precision x3,x8,pi,c0,h,w0,e0
double precision ma,dmp,dmn,dmm
double precision theta,cth,sth
double precision Aij,ls2
double precision tel,noem
double precision, dimension (lF1,3) :: delta_I,delta_Q,delta_U,delta_V
double precision, dimension (lF1,3) :: dp,dm
!input
!theta: angle magnetic field-radiation field
!ls: line strength of the transition
!F1,F2: angular momentum upper and lower state

!output
!delta_I
!delta_Q
!delta_U
!delta_V  are either real or imaginary elements

x0 = 0.d0
x1 = 1.d0
x2 = 2.d0
xs2 = dsqrt(x2)
x3 = 3.d0
x8 = 8.d0

pi = dacos(-x1)

c0 = 299792458.d0 
h  = 6.62607004E-34
e0 = 8.854187817E-12

cth = dcos(theta)
sth = dsin(theta)
fac1 = (x1+cth)/x2
fac2 = sth/xs2
fac3 = (x1-cth)/x2

!dp dipole matrix 
!dp

do a1 = 1,lF1
  do I=1,3
    delta_I(a1,I) = 0.d0  
    delta_Q(a1,I) = 0.d0
    delta_U(a1,I) = 0.d0
    delta_V(a1,I) = 0.d0 
  end do
end do


do a1=1,lF1

  ma = x1*(a1-F1-1)

  call NED(x1*F2,x1,x1*F1,ma-x1,x1,ma,dmp) 
  call NED(x1*F2,x1,x1*F1,ma,x0,ma,dmn)
  call NED(x1*F2,x1,x1*F1,ma+x1,-x1,ma,dmm) 

  dmp = dmp/dsqrt(x1*lF1)
  dmn = dmn/dsqrt(x1*lF1)
  dmm = dmm/dsqrt(x1*lF1)

  dp(a1,1) = dmp*fac1 
  dp(a1,2) = dmn*fac2 ! imaginary
  dp(a1,3) = -x1*dmm*fac3

  dm(a1,1) = -x1*dmp*fac3
  dm(a1,2) = dmn*fac2 ! imaginary
  dm(a1,3) = dmm*fac1


end do

!Relation between Aij and ls2. D&W 1990
tel = Aij*h*x3*c0**x3
noem = x8*pi*w0**x3   

ls2 = tel/noem 

do a1 = 1,lF1
  do I=1,3
    delta_I(a1,I) = ls2*(dp(a1,I)*dp(a1,I) + dm(a1,I)*dm(a1,I)) !real
    delta_Q(a1,I) = ls2*(dp(a1,I)*dm(a1,I) + dm(a1,I)*dp(a1,I))
    delta_U(a1,I) = ls2*(dp(a1,I)*dm(a1,I) - dm(a1,I)*dp(a1,I))
    delta_V(a1,I) = ls2*(dp(a1,I)*dp(a1,I) - dm(a1,I)*dm(a1,I)) 
  end do
end do

end subroutine


