subroutine fwhm(Iw,nw,dw,fwhm_i)
implicit none
integer n,nw,fl1,fl2
double precision x1,x2,dw,flr1,flr2
double precision Imax,Im2,I,fwhm_i,dI
double precision, dimension (-nw:nw) :: Iw


x1 = 1.d0
x2 = 2.d0

Imax = maxval(Iw)
Im2 = Imax/x2

fl1 = 0 
n = -nw

do while (fl1.eq.0)
  I = Iw(n)

  if (I.gt.Im2) then
    fl1 = n 
  endif

  n = n + 1

end do

! residual in first order (flag real: flr1)

dI = ( Iw(fl1) - Iw(fl1-1) )
flr1 = x1*fl1 - ( Iw(fl1) - Im2 ) / dI

fl2 = 0 
n = nw

do while (fl2.eq.0)
  I = Iw(n)

  if (I.gt.Im2) then
    fl2 = n 
  endif

  n = n - 1

end do

! residual in first order (flag real: flr2)
dI = ( Iw(fl2) - Iw(fl2+1) )
flr2 = x1*fl2 + ( Iw(fl2) - Im2 ) / dI

!write(*,*) flr1,flr2

fwhm_i = dw * (flr2 - flr1) 

end subroutine


subroutine zeeman(Iw,Vw,nw,dw,azee,bzee)
implicit none
integer n,nw
double precision x1,x2,dw,dwt,fac
double precision azee,bzee 
double precision, dimension (2,2) :: Ap,Api
double precision, dimension (2) :: Vpr
double precision, dimension (-nw:nw) :: Iw,Vw
double precision, dimension (-nw:nw) :: dIw

! simple fitting program for the V-spectrum

x1 = 1.d0
x2 = 2.d0

! making the gradient of Iw
n = -nw
dIw(n) = (Iw(n+1) - Iw(n))/dw
n = nw
dIw(n) = (Iw(n) - Iw(n-1))/dw

dwt = x2*dw
do n=-nw+1,nw-1
  dIw(n) = (Iw(n+1) - Iw(n-1))/dwt
end do

! The matrix 
! Ap =  
! Iw'*Iw  , Iw'*dIw
! dIw'*Iw , dIw'*dIw
! is gonna be inverted later

Ap(1,1) = 0.d0
Ap(1,2) = 0.d0
Ap(2,1) = 0.d0
Ap(2,2) = 0.d0

do n=-nw,nw
  Ap(1,1) = Ap(1,1) + Iw(n)*Iw(n)
  Ap(1,2) = Ap(1,2) + Iw(n)*dIw(n)
  Ap(2,2) = Ap(2,2) + dIw(n)*dIw(n)
end do

Ap(2,1) = Ap(1,2)

!inversion of Ap
fac = x1/(Ap(1,1)*Ap(2,2) - Ap(1,2)*Ap(2,1))  ! 1/det(Ap)

Api(1,1) = fac*Ap(2,2)
Api(1,2) = -fac*Ap(2,1)
Api(2,1) = Api(1,2) 
Api(2,2) = fac*Ap(1,1)

! Vp is
! Iw'*Vw
! dIw'*Vw

Vpr(1) = 0.d0
Vpr(2) = 0.d0
do n=-nw,nw 
  Vpr(1) = Vpr(1) + Iw(n)*Vw(n)
  Vpr(2) = Vpr(2) + dIw(n)*Vw(n)
end do

! a                 Iw'*Vw
!     =  Ap^(-1) *
! b                 dIw'*Vw

azee = Api(1,1)*Vpr(1) + Api(1,2)*Vpr(2)
bzee = Api(2,1)*Vpr(1) + Api(2,2)*Vpr(2)

!write(*,*)Iw
!write(*,*)Vw
!write(*,*)
!write(*,*)dIw
!write(*,*)
!write(*,*)Ap
!write(*,*)Api
!write(*,*)
!write(*,*)Vp
!write(*,*)
!write(*,*)azee,bzee

end subroutine


subroutine vitals(Iw,Qw,Vw,nw,pQ,pV,nimax) 
implicit none
integer n,nw,nvmax,nvmin,nqmax,nimax
double precision I0,Q0,U0,Vmax,Vmin,Qmax,Imax
double precision Qtot,x2,pi,Itemp,Vtemp,Qtemp
double precision pQ,pV,pa 
double precision, dimension (-nw:nw) :: Iw,Qw,Uw,Vw
double precision, dimension (-nw:nw) :: dIw

x2 = 2.d0
pi = dacos(-1.d0)

Vmax = Vw(-nw)
Imax = Iw(-nw)
Vmin = Vw(-nw)
Qmax = Qw(-nw)

nimax = -nw
nvmax = -nw
nqmax = -nw
nvmin = -nw
do n=-nw+1,nw
  Itemp = Iw(n)
  Qtemp = dabs(Qw(n))
  Vtemp = Vw(n)
 
  if (Itemp.gt.Imax) then
    Imax = Itemp
    nimax = n
  endif

  if (Qtemp.gt.Qmax) then
    Qmax = Qtemp
    Q0 = Qw(n)
    nqmax = n
  endif

  if (Vtemp.gt.Vmax) then
    Vmax = Vtemp
    nvmax = n
  endif

  if (Vtemp.lt.Vmin) then
    Vmin = Vtemp
    nvmin = n
  endif

end do


pQ = Q0/Imax

if (nvmax.gt.nimax) then 
  pV = (Vmax - Vmin)/Imax
else
  pV = -(Vmax - Vmin)/Imax
endif

end subroutine

subroutine gaussfit(Iw,nw,dw,mu,sig,Is) 
implicit none
integer n,nw,I,J
double precision mu,sig,Is,sig2
double precision dw,w
double precision x0,x1,x2,pi 
double precision c_a,c_b,c_c 
double precision, dimension (-nw:nw) :: Iw
double precision, dimension (-nw:nw,3) :: A
double precision, dimension (3) :: C
double precision, dimension (3,3) :: B,Bi

! we are trying to solve
! g(w) = a + bw + cw^2
! g(w) = ln(I(w))
! a = ln(Is/sig*sqrt(2\pi)) - 0.5*(\mu / \sig)^2,
! b = \mu/\sig^2,
! c = -1/(2*\sig^2),

!or: g = A*[a;b;c]
!A^T g = C = B * [a;b;c]
! B-1 * C = [a;b;c]

x0 = 0.d0
x1 = 1.d0
x2 = 2.d0

pi = dacos(-x1)

do n = -nw,nw
  w = dw*n

  A(n,1) = 1.d0
  A(n,2) = w 
  A(n,3) = w**x2 
end do

do I = 1,3
  do J = 1,3
    B(I,J) = x0 
  end do
end do

do n = -nw,nw
  do I = 1,3
    do J = I,3
      B(I,J) = B(I,J) + A(n,I)*A(n,J)
      B(J,I) = B(I,J)
    end do
    C(I) = C(I) + A(n,I)*dlog(Iw(n))
  end do
end do

call matinv3(B,Bi)

c_a = x0
c_b = x0
c_c = x0

do I=1,3
  c_a = c_a + Bi(1,I)*C(I)
  c_b = c_b + Bi(2,I)*C(I)
  c_c = c_c + Bi(3,I)*C(I)
end do

!now to the gaussian paramaters
sig2 = -x1/(c_c*x2)

sig = dsqrt(sig2)
mu = c_b * sig2
Is = sig*dsqrt(x2*pi)*dexp( c_a + (mu * c_b / x2) )

end subroutine

subroutine matinv3(A,B) 
!! Performs a direct calculation of the inverse of a 3×3 matrix.
implicit none 
double precision detinv
double precision, dimension(3,3) :: A,B 

! Calculate the inverse determinant of the matrix
detinv = 1.d0/(A(1,1)*A(2,2)*A(3,3) - A(1,1)*A(2,3)*A(3,2)&
          - A(1,2)*A(2,1)*A(3,3) + A(1,2)*A(2,3)*A(3,1)&
          + A(1,3)*A(2,1)*A(3,2) - A(1,3)*A(2,2)*A(3,1))

! Calculate the inverse of the matrix
B(1,1) = +detinv * (A(2,2)*A(3,3) - A(2,3)*A(3,2))
B(2,1) = -detinv * (A(2,1)*A(3,3) - A(2,3)*A(3,1))
B(3,1) = +detinv * (A(2,1)*A(3,2) - A(2,2)*A(3,1))
B(1,2) = -detinv * (A(1,2)*A(3,3) - A(1,3)*A(3,2))
B(2,2) = +detinv * (A(1,1)*A(3,3) - A(1,3)*A(3,1))
B(3,2) = -detinv * (A(1,1)*A(3,2) - A(1,2)*A(3,1))
B(1,3) = +detinv * (A(1,2)*A(2,3) - A(1,3)*A(2,2))
B(2,3) = -detinv * (A(1,1)*A(2,3) - A(1,3)*A(2,1))
B(3,3) = +detinv * (A(1,1)*A(2,2) - A(1,2)*A(2,1))
end subroutine 

