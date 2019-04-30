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

subroutine turnback(d2_r,d2_i,t2_r,t2_i,th)
implicit none
integer k1,k2,k1c,k2c
double precision x0,x1,x2,pi,th,smalld,fas_r,fas_i
double precision, dimension (5) :: d2_r,d2_i,t2_r,t2_i


x0 = 0.d0
x1 = 1.d0
x2 = 2.d0

pi = dacos(-x1)

do k1c = 1,5
  t2_r(k1c) = x0
  t2_i(k1c) = x0
end do


do k1 = -2,2
  k1c = k1+2+1
  do k2 = -2,2
    k2c = k2+2+1
     
    fas_r = dcos(pi*(k2-k1)/x2)
    fas_i = dsin(pi*(k2-k1)/x2)

    call wig_small(2,k1,k2,-th,smalld)
  
    t2_r(k1c) = t2_r(k1c) + fas_r * smalld * d2_r(k2c)      
    t2_r(k1c) = t2_r(k1c) - fas_i * smalld * d2_i(k2c)      
     
    t2_i(k1c) = t2_i(k1c) + fas_i * smalld * d2_r(k2c)      
    t2_i(k1c) = t2_i(k1c) + fas_r * smalld * d2_i(k2c)      
 
  end do
end do 


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
double precision alp,bet,gam,frac
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

!it can happen, that because of numerical reasons, the frac 
!is slightly out of bounds

frac = dcos(psi/x2) / dcos(bet/x2) 
if (frac.gt.x1) then
  frac = x1
elseif (frac.lt.-x1) then
  frac = -x1
endif

gpa = dacos(frac)
gma = phi - pi/x2

gam = gpa + gma
alp = gpa - gma 

!write(*,*)th,psi,phi
!write(*,*)alp,bet,gam

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

subroutine dip_coord(th,phi,chi,theta,pump,succ)
!function to compute the coordinates for the 
!propagation-magnetic field angle and the
!pumping unit vector  
!in the correct frame of reference for the maser
!calulations

!input:
!th,phi -- angles describing the dipole vector (in frame of propagation direction)
!chi -- angle describing the place on the SiO maser ring

!output:
!theta -- propagation angle
!pump -- anisotropic pumping vector in the right frame of reference

!algorithm not always succesful -- succ = 0 succes
!				-- succ = 1 fail

implicit none
integer I,J,K,CO,succ
double precision th,phi,chi,theta,bet,nb,nn,alp
double precision mvr,x0,x1,x2,x3,pi,asp,asm,acp,acm
double precision, dimension(3) :: mv,pump,b,r,n,sp,spp,z
double precision, dimension(3,3) :: R1,R2

succ = 0

x0 = 0.d0
x1 = 1.d0
x2 = 2.d0
x3 = 3.d0

pi = dacos(-x1)

z(1) = x0
z(2) = x0
z(3) = x1

!initial frame: propagation in z-direction
mv(1) = dsin(th)*dcos(phi)
mv(2) = dsin(th)*dsin(phi)
mv(3) = dcos(th)

r(1) = dcos(chi)
r(2) = dsin(chi)
r(3) = x0         !r -- place of SiO maser

mvr = x0
do I=1,3
  mvr = mvr + r(I)*mv(I) 
end do

b(1) = x3*mvr*r(1) - mv(1)
b(2) = x3*mvr*r(2) - mv(2)
b(3) = -mv(3)

nb = dsqrt(b(1)**x2 + b(2)**x2 + b(3)**x2)
b(1) = b(1) / nb
b(2) = b(2) / nb
b(3) = b(3) / nb

!write(*,*)b

!in this frame propagtation in z-direction, so prop-mag angle is
theta = dacos(b(3))

!first, let's rotate b to its [0;0;1] position
!n = b \cross z
n(1) = b(2)
n(2) = -b(1)
n(3) = x0

nn = dsqrt(n(1)**x2 + n(2)**x2 + n(3)**x2)
n(1) = n(1)/nn 
n(2) = n(2)/nn 
!n(3) = n(3)/nn 

bet = dacos(b(3)) 
!write(*,*)n
!write(*,*)bet
call rotmat(n,bet,R1)

!write(*,*)R1(1,:)
!write(*,*)R1(2,:)
!write(*,*)R1(3,:)


do I = 1,3
  sp(I) = R1(I,3)
end do

!write(*,*)sp

!the unit vector, spp, for propagation direction, has the form
!spp = [0;sin(theta);cos(theta)]
!and will be obtained by an Rz(alp) transformation (because we want b
!in the same frame)

!now, because sp'*b' = cos(theta), we have the following sp matrix
![\pm sin(alp)*sin(theta) ; \pm cos(alp)*sin(theta) ; cos(theta]

!check for the solution that has equivalent angles
asp = dasin(sp(1)/dsin(theta))
asm = dasin(-sp(1)/dsin(theta))

acp = dacos(sp(2)/dsin(theta))
acm = dacos(-sp(2)/dsin(theta))

if (dabs(asp-acp).lt.1E-8) then
  alp = -asp
  call rotmat(z,-alp,R2)
  CO = 1
elseif (dabs(asp-acm).lt.1E-8) then
  alp = -asp
  call rotmat(z,pi+alp,R2)
  CO = 2
elseif (dabs(asm-acp).lt.1E-8) then
  alp = -asm
  call rotmat(z,alp,R2)
!  write(*,*)'3'
  CO = 3
elseif (dabs(asm-acm).lt.1E-8) then
  alp = -asm
  call rotmat(z,pi-alp,R2)
  CO = 4
else
  write(*,*)'solution to spp not found'
  stop 
endif

!write(*,*)R2(1,:)
!write(*,*)R2(2,:)
!write(*,*)R2(3,:)

!finally, check the solution
!can be turned off after confirmation

do I = 1,3
  spp(I) = x0
  do J = 1,3
    spp(I) = spp(I) + R2(I,J)*sp(J)
  end do
end do

IF (dabs(spp(2)-dsin(theta)).gt.1E-8) then
  succ = 1 !then
   
  write(*,*)CO

ENDIF
!!  write(*,*)'error'
!  if (CO.eq.3) then
!    call rotmat(z,pi-alp,R2)
!
!    do I = 1,3
!      spp(I) = x0
!      do J = 1,3
!        spp(I) = spp(I) + R2(I,J)*sp(J)
!      end do
!    end do
!  elseif (CO.eq.4) then
!    call rotmat(z,(alp),R2)
!
!    do I = 1,3
!      spp(I) = x0
!      do J = 1,3
!        spp(I) = spp(I) + R2(I,J)*sp(J)
!      end do
!    end do
!  endif
!
!!  spp(2) = -spp(2)
!  IF (dabs(spp(2)-dsin(theta)).gt.1E-8) succ = 1 
!ENDIF 

!now, put the pumping function in the right frame
do I = 1,3
  pump(I) = x0
  do J = 1,3
    do K = 1,3
      pump(I) = pump(I) + R2(I,J)*R1(J,K)*r(K) 
    end do
  end do
end do

end subroutine

subroutine rotmat(n,th,R)
!function to compute the rotation-matrix of unit vector n and angle beta
implicit none
integer I,J
double precision th,cth,sth,x1
double precision, dimension (3) :: n
double precision, dimension (3,3) :: R,nn,nx  

x1 = 1.d0
cth = dcos(th)
sth = dsin(th)

do I = 1,3
  do J = 1,3
    nn(I,J) = n(I)*n(J)
  end do
end do


nx(1,2) = -n(3)
nx(2,1) = n(3)

nx(1,3) = n(2)
nx(3,1) = -n(2)

nx(2,3) = -n(1)
nx(3,2) = n(1)

do I = 1,3
  do J = 1,3
    R(I,J) = sth*nx(I,J) + (x1-cth) * nn(I,J)
    if (I.eq.J) then
      R(I,I) = R(I,J) + cth 
    endif
  end do
end do

end subroutine
