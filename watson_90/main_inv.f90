program watned_90
implicit none
integer nF1,nF2,F1,F2,lF1,lF2,Dmm
integer MAXSTEP,STEP
integer a1,a2,bc1,bc2
integer anis,inci
double precision theta,w0,Rtot,Ttot,dR,B
double precision Tb,c0
double precision sa,ds,x1,x2,x3,pi
double precision pQ,pV,pa,Tbr,qin 
double precision Aw,Bw,Cw,Dww,Ew,Fw,Gw,Aws,Bws,Cws,Fws 
double precision Iw,Qw,Uw,Vw,gam_a,gam_b,lam1,del,p_an_u,p_an_d,Aij,alp_Z1,alp_Z2 
double precision rhoa0_r,rhob0_r,rhoa22_r,rhoa21_r,rhoa20_r,d00,rhob22_r,rhob21_r,rhob20_r
double precision rhoa2m2_r,rhoa2m2_i,rhob2m2_r,rhob2m2_i,rhoa2m1_r,rhoa2m1_i,rhob2m1_r,rhob2m1_i
double precision rhoa0_i,rhob0_i,rhoa22_i,rhoa21_i,rhoa20_i,rhob22_i,rhob21_i,rhob20_i
double precision d22_r,d22_i,d21_r,d21_i,d20_r,d20_i,wig1,wig2,d2m2_r,d2m2_i,d2m1_r,d2m1_i
double precision, dimension(:), allocatable :: wa,wb 
double precision, dimension(:,:), allocatable :: rho_aa_r,rho_aa_i,rho_bb_r,rho_bb_i 
double precision, dimension(:,:,:,:), allocatable :: delta_I,delta_Q,delta_U,delta_V
double precision, dimension(:,:), allocatable ::  lbd_a_r,lbd_b_r 
double precision, dimension(:,:), allocatable ::  lbd_a_i,lbd_b_i 
double precision, dimension(5) :: d2_r,d2_i,t2_r,t2_i
double precision, dimension(3) :: apump

Dmm = 3

c0 = 299792458.d0

!open(unit=61,file="dab0.dat")
!open(unit=62,file="p2m2.dat")
!open(unit=63,file="p2m1.dat")
!open(unit=64,file="p20.dat")
!open(unit=65,file="p21.dat")
!open(unit=66,file="p22.dat")
!open(unit=67,file="pt2m2.dat")
!open(unit=68,file="pt2m1.dat")
!open(unit=69,file="pt20.dat")
!open(unit=70,file="pt21.dat")
!open(unit=71,file="pt22.dat")

open(unit=28,file="ned_props.dat")
!open(unit=29,file="fw.dat")

open(unit=4,file="in_watned.dat")

read(4,*)

read(4,*)F1,F2,alp_Z1,alp_Z2,Aij,lam1,del,gam_a,gam_b,p_an_u,p_an_d

read(4,*)

read(4,*)w0
read(4,*)B,theta
read(4,*)Tb,ds,MAXSTEP,sa
read(4,*)anis,apump(1),apump(2),apump(3)
read(4,*)inci,qin

x1 = 1.d0
x2 = 2.d0
x3 = 3.d0
pi = dacos(-1.d0)

theta = theta*pi/180.d0
w0 = x2*pi*w0


lF1 = F1*2 + 1
lF2 = F2*2 + 1

allocate(wa(lF1))
allocate(wb(lF2))

allocate(lbd_a_r(lF1,lF1))
allocate(lbd_a_i(lF1,lF1))
allocate(lbd_b_r(lF2,lF2))
allocate(lbd_b_i(lF2,lF2))


allocate(rho_aa_r(lF1,lF1))
allocate(rho_aa_i(lF1,lF1))
allocate(rho_bb_r(lF2,lF2))
allocate(rho_bb_i(lF2,lF2))

!write(*,*)alp_Z1
call initiate(lF1,lF2,F1,F2,w0,alp_Z1,alp_Z2,B,wa,wb,gam_a,gam_b,&
&Iw,Uw,Qw,Vw,Tb,lam1,del,lbd_a_r,lbd_a_i,lbd_b_r,lbd_b_i,anis,apump,p_an_u,p_an_d,inci,qin)
!write(*,*)alp_Z1
!write(*,*)wa
!write(*,*)wb
allocate(delta_I(lF1,Dmm,lF1,Dmm))
allocate(delta_Q(lF1,Dmm,lF1,Dmm))
allocate(delta_U(lF1,Dmm,lF1,Dmm))
allocate(delta_V(lF1,Dmm,lF1,Dmm))

do a1 = 1,lF1 
  do a2 = 1,lF1
    do bc1 = 1,3
      do bc2 = 1,3
        delta_I(a1,bc1,a2,bc2) = 0.d0
        delta_Q(a1,bc1,a2,bc2) = 0.d0
        delta_U(a1,bc1,a2,bc2) = 0.d0
        delta_V(a1,bc1,a2,bc2) = 0.d0
      end do
    end do
  end do
end do

call gen_delta(theta,Aij,w0,F1,lF1,F2,delta_I,delta_Q,delta_U,delta_V)

nF1 = lF1*lF1
nF2 = lF2*lF2

STEP = 0
111 STEP = STEP + 1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!! RATE EQUATIONS !!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  call stateq(lF1,lF2,nF1,nF2,lbd_a_r,lbd_a_i,lbd_b_r,lbd_b_i,gam_a,gam_b,wa,wb,sa,Iw,&
&Qw,Uw,Vw,delta_I,delta_Q,delta_U,delta_V,rho_aa_r,rho_aa_i,rho_bb_r,rho_bb_i)
  ! paa density elements  


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!! TRANSFER ELEMENTS !!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  call trans_el(rho_aa_r,rho_aa_i,rho_bb_r,rho_bb_i,lF1,lF2,&
&delta_I,delta_Q,delta_U,delta_V,Aw,Bw,Cw,Dww,Ew,Fw,Gw,Aws,Bws,Cws,Fws,w0)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!! STOKES PROPAGATION !!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! for now, without spontaneous emission

  call prop_stokes(Aw,Bw,Cw,Dww,Ew,Fw,Gw,Iw,Qw,Uw,Vw,ds)

!write(*,*) 'step:', STEP

! if any of the densities is not a number, we might as well stop.

if (isnan(rho_aa_r(1,1)))  stop   
if (isnan(rho_aa_i(1,1)))  stop

if (isnan(rho_bb_r(1,1)))  stop
if (isnan(rho_bb_i(1,1)))  stop


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! compute rate of stimulated emission !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  call stimulate(Iw,Qw,Uw,Vw,lF1,lF2,Aij,w0,delta_I,delta_Q,delta_U,delta_V,Rtot,Ttot)
!  call stimulate_alt(Iw(0),Aij,sa,w0,Rstim)

!dR = dlog10(Rtot/gam_a)
!write(*,*)Rtot,alp_Z1,B
dR = dlog10(Rtot / (0.75d0 * B))
!write(*,*)dR
!stop
pQ = dsqrt(Qw**2.d0 + Uw**2.d0)/Iw
pa = datan2(Uw,Qw)*90.d0/pi
pV = Vw/Iw

if (mod(STEP,10).eq.0) then
  write(28,*)dR,pQ,pa
endif

if (dR.gt.3.1d0) stop

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!            and... go on             !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

ds = ds*1.0001d0

if (STEP.le.MAXSTEP) GOTO 111    

close(28)
!close(29)

!close(51)
!close(52)
!close(53)
!close(54)
!close(55)
!close(56)
!close(57)
!close(58)
!
!close(61)
!close(62)
!close(63)
!close(64)
!close(65)
!close(66)

!close(101)
!close(102)
!close(103)
!close(104)
!close(105)
!close(106)
!close(107)
!close(108)
!close(109)
!close(110)
!close(111)
!close(112)
!close(113)
!close(114)
!close(115)
!close(116)
!close(117)
!close(118)
!close(119)
!close(120)

end program
 
