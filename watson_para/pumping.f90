subroutine anis_pump(lbd_u,p_an_u,thp,j1,lj1,a1,a2,lbd1_m)
!function to compute the anisotropic pumping matrix
!according to lankhaar
implicit none
integer m1,m2,q,j1,j2,lj1
integer b1,b2,a1,a2

double precision x1,x2,x3,x8 
double precision fac1,st,st2,ct,ct2,s2t,f5,cg10,cg20,cg2,fas,fac 
double precision lbd_u,lbd_d,p_an_u,p_an_d,thp,lbd1_m 
double precision, dimension (-2:2) :: lbd_u_2

!input variables
!lbd: pumping rate
!p_an: anisotropic fraction
!thp: pumping angle

!output
!lbd_m: pumping matrix element lbd (a1,a2)

! the pumping matrix irreducible components 

x1 = 1.d0
x2 = 2.d0
x3 = 3.d0
x8 = 8.d0

fac1 = dsqrt(x3/x8)
st = dsin(thp)
st2 = st*st
ct = dcos(thp)
ct2 = ct*ct

s2t = dsin(x2*thp)

lbd_u_2(-2) = lbd_u*p_an_u*fac1*st2 
lbd_u_2(-1) = -x1*lbd_u*p_an_u*fac1*s2t
lbd_u_2(0) = lbd_u*p_an_u*(x3*ct2-x1)/x2 
lbd_u_2(1) = -x1*lbd_u_2(-1) 
lbd_u_2(2) = lbd_u_2(-2) 

f5 = dsqrt(5.d0)
cg10 = dsqrt(x1*lj1)

m1 = a1 - j1 - 1
m2 = a2 - j1 - 1
lbd1_m = 0.d0

if (m1.eq.m2) then
  lbd1_m = lbd_u*cg10
endif

do q = -2,2
  call NED(x1*j1,x2,x1*j1,-x1*m1,x1*q,-x1*m2,cg2)
  fas = (-x1)**(m1-m2)
  fac = f5*cg2*fas/cg10

  lbd1_m = lbd1_m + lbd_u_2(q)*fac
end do

end subroutine

subroutine anis_pump_nw(lbd_u,p_an_u,alp,thp,j1,lj1,a1,a2,lbd1_m_r,lbd1_m_i)
!function to compute the anisotropic pumping matrix
!according to: N&W and W&W 
implicit none
integer k1,k2,m1,m2,q,j1,j2,lj1
integer b1,b2,a1,a2,ma1,ma2

double precision x0,x1,x2,x3,x8
double precision fac1,fac2,st,st2,ct,ct2,s2t,f5,cg10,cg20,cg2,fas,fac 
double precision lbd_u,lbd_d,p_an_u,p_an_d,thp,lbd1_m,smalld1,smalld2,lbd_m 
double precision lbd1_m_r,lbd1_m_i,alp 

!input variables
!lbd: pumping rate
!p_an: anisotropic
!alp,thp: pumping angles

!output
!lbd_m_r: real part of the pumping matrix
!lbd_m_i: imaginary part of the pumping matrix

! the pumping matrix irreducible components 

x0 = 0.d0
x1 = 1.d0
x2 = 2.d0
x3 = 3.d0
x8 = 8.d0

ma1 = a1 - j1 - 1
ma2 = a2 - j1 - 1

lbd1_m = 0.d0

!first, create the diagonal lambda matrix
do m1 = -j1,j1
  fac1 = x3*(j1**x2 + x1*j1 - x1 + m1**x2) 
  fac2 = (x2*j1 - x1)*(x2*j1 + x3) 

  fac = fac1/fac2 - x1

  lbd_m = lbd_u*(x1 + fac*p_an_u)

  !now, rotate the matrix 
  ! L' = DLD^T
  ! L'_ij = \sum_k D_ik (D_jk)^* L_kk
  ! L'_ij = e^(-%i*(i-j)*alpha) * \sum_k d_ik d_jk L_kk
  ! if L is diagonal

  call wig_small(j1,ma1,m1,thp,smalld1)
  call wig_small(j1,ma2,m1,thp,smalld2)

  lbd1_m = lbd1_m + smalld1*smalld2*lbd_m
end do

lbd1_m_r = lbd1_m * dcos(x1*(ma2 - ma1)*alp)
lbd1_m_i = lbd1_m * dsin(x1*(ma2 - ma1)*alp)

end subroutine

subroutine wig_small(j1,m1,m2,th,smalld)
!function to compute the wigner small-d matrix element
implicit none
integer j1,m1,m2,fact
integer jpm1,jmm1,jpm2,jmm2,fpm1,fmm1,fpm2,fmm2
integer jps,jms,dms,fps,fms,fdm,fs,dm
integer mins,maxs,s

double precision th,prefac,fac1,fac2
double precision pow1,pow2,cth,sth,cp,sp,sums,smalld 
double precision x1,x2

x1 = 1.d0
x2 = 2.d0

jpm1 = j1 + m1  
jmm1 = j1 - m1
jpm2 = j1 + m2
jmm2 = j1 - m2

!write(*,*)jpm1,jmm1,jpm2,jmm2

fpm1 = fact(jpm1)  
fmm1 = fact(jmm1)
fpm2 = fact(jpm2)
fmm2 = fact(jmm2)

!write(*,*)fpm1,fmm1,fpm2,fmm2

prefac = dsqrt(x1 * fpm1 * fmm1 * fpm2 * fmm2)

dm = m1 - m2

maxs = min(jpm2,jmm1)

if (dm.lt.0) then
  mins = -dm 
else
  mins = 0
endif

sums = 0.d0

!write(*,*)

do s = mins,maxs

  fac1 = (-x1)**s

  jps = jpm2 - s
  jms = jmm1 - s
  dms = dm + s
  
!  write(*,*)jps,jms,dms,s

  fps = fact(jps)
  fms = fact(jms)
  fdm = fact(dms)
  fs  = fact(s)

!  write(*,*)fps,fms,fdm,fs

  fac2 = x1 * fps * fms * fdm * fs 

  pow1 = x1*(jps + jms)
  pow2 = x1*(dms + s)

  cth = dcos(th/x2)
  sth = dsin(th/x2)

  cp = cth**pow1
  sp = sth**pow2

  sums = sums + (fac1/fac2)*cp*sp

end do

smalld = prefac * sums

end subroutine

subroutine twovec_euler(vec,alp,bet,gam)
!subroutine to compute the euler angles between [0;0;1] and [vec(1);vec(2);vec(3)]
integer I,J
double precision x0,x1,x2,pi,lvec
double precision alp,bet,gam
double precision phi,th,psi,gpa,gma
double precision, dimension (3) :: vec,cvec 

x0 = 0.d0
x1 = 1.d0
x2 = 2.d0

pi = dacos(-x1)

! make sure that vec is a unit vector
lvec = x0
do I = 1,3
  lvec = lvec + vec(I)*vec(I)
end do

do I = 1,3
  vec(I) = vec(I)/lvec
end do

! determine axis-angle rotational part
! cvec = z \cross vec
cvec(1) = -vec(2)
cvec(2) = vec(1)
cvec(3) = x0 
! psi = acos(z \dot vec)
psi = dacos(vec(3))

! from the rotation axis n, to the two angles
phi = datan2(cvec(2),cvec(1))
th  = dacos(cvec(3))

! now, from these three rotation angles to the Euler angles
! following appendix A --- Angular Momentum Techniques in Quantum Mechanics
! Eqs. A.40, A.45 and A.47

bet = x2 * dasin(dsin(th)*dsin(psi/x2))
gpa = dacos( dcos(psi/x2) / dcos(bet/x2) )
gma = phi - pi/x2

gam = gpa + gma
alp = gpa - gma 

end subroutine

recursive function fact(i) result(j)
    integer, intent(in) :: i
    integer :: j
    if (i == 1) then
        j = 1
    elseif (i == 0) then
        j = 1
    else
        j = i * fact(i - 1)
    end if
end function fact
