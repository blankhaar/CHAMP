program watned_94
implicit none
integer F1,F2,lF1,lF2,mlF1,mlF2,Dmm
integer I,J,K,L,MAXSTEP,STEP,RAT,NBREAK
integer ma1,ma2,mb1,mb2,a1,a2,b1,b2,bc1,bc2
integer vj,loop,nw,nF,MAXRHO,THRESH
integer vvj,II
double precision theta,w0,rho_new_r,rho_new_i,om,dw,dv,Rtot,Ttot,Rstim,Rmax,dR,B
double precision tmpa1,tmpa2,tmpb1,tmpb2,Tb,gam,vthermal,lam1,lam2,c0,fwhm_i
double precision Aij,gam_a,gam_b,lbd,sa,mass,ds,x2,pi,fwhm_old,bzee,azee
double precision rata_1,rata_2,rata_3,ratb_1,ratb_2,ratb_3
double precision ratsa_2,ratsa_3,ratsb_2,ratsb_3,tthresh,Rthr
double precision pQ,pV,pa,Tbr 
double precision mu,sig,Is,x0 
double precision gpmax1,gpmax2,gpmax,delta_w
integer, dimension(:,:), allocatable :: Ftrans 
double precision, dimension(:,:), allocatable :: wall_a,wall_b 
double precision, dimension(:), allocatable :: wa,wb 
double precision, dimension(:), allocatable :: Iw,Qw,Vw
double precision, dimension(:), allocatable :: Ipv_r,Upv_r,Qpv_r,Vpv_r
double precision, dimension(:,:), allocatable :: Ip_r,Up_r,Qp_r,Vp_r
double precision, dimension(:,:,:), allocatable :: Ipabv_r,Upabv_r,Qpabv_r,Vpabv_r
double precision, dimension(:), allocatable :: Aij_all,w00,alp_Z1,alp_Z2 
double precision, dimension(:,:,:), allocatable :: all_delta_I,all_delta_Q,all_delta_U,all_delta_V
double precision, dimension(:,:,:), allocatable :: all_rho_a_r,all_rho_b_r
double precision, dimension(:,:,:), allocatable :: rho_in_a_r,rho_in_b_r
double precision, dimension(:), allocatable :: rho_old_a_r,rho_old_b_r
double precision, dimension(:,:), allocatable :: rho_aa_r,rho_bb_r
double precision, dimension(:), allocatable :: Aw,Bw,Cw
double precision, dimension(:,:), allocatable :: lbd_a,lbd_b 
double precision, dimension(:), allocatable :: Awtot,Bwtot,Cwtot
double precision, dimension(:,:,:), allocatable :: gab
double precision, dimension(:,:), allocatable :: delta_I,delta_Q,delta_U,delta_V


Dmm = 3

c0 = 299792458.d0
x0 = 0.d0

open(unit=51,file="I.dat")
open(unit=52,file="U.dat")
open(unit=53,file="Q.dat")
open(unit=54,file="V.dat")
open(unit=55,file="R.dat")
open(unit=56,file="fw.dat")
open(unit=57,file="info.dat")
open(unit=58,file="zee.dat")
open(unit=61,file="pa00.dat")
open(unit=63,file="pa11.dat")
open(unit=64,file="pam1m1.dat")
open(unit=62,file="pb.dat")

open(unit=28,file="ned_props.dat")

open(unit=4,file="in_watned.dat")

read(4,*)
read(4,*)nF
read(4,*)


allocate(Ftrans(nF,2))
allocate(w00(nF))
allocate(alp_Z1(nF))
allocate(alp_Z2(nF))
allocate(Aij_all(nF))

do I=1,nF
  read(4,*)Ftrans(I,1),Ftrans(I,2),w00(I),alp_Z1(I),alp_Z2(I),Aij_all(I)
end do

read(4,*)

read(4,*)w0
read(4,*)B,theta
read(4,*)nw,dw
read(4,*)vthermal,mass,lam1,lam2,gam_a,gam_b
read(4,*)Tb,ds,MAXSTEP,MAXRHO,THRESH
read(4,*)sa

x2 = 2.d0
pi = dacos(-1.d0)

tthresh = 10.d0**(-1.d0*THRESH)
Rthr = tthresh

theta = theta*pi/180.d0
w0 = x2*pi*w0
dv = c0*dw/w0

do I=1,nF
  w00(I) = w00(I)*x2*pi
end do

mlF1 = Ftrans(1,1)
mlF2 = Ftrans(1,2)
         
do I=2,nF 
  if (Ftrans(I,1).gt.mlF1) then 
    mlF1 = Ftrans(I,1)
  endif
  if (Ftrans(I,2).gt.mlF2) then
    mlF2 = Ftrans(I,2)
  endif
end do

mlF1 = mlF1*2 + 1
mlF2 = mlF2*2 + 1

!the adaptive parameters
do I=1,nF
  gpmax1 = (mlF1 - 1)*B*alp_Z1(I)  
  gpmax2 = (mlF2 - 1)*B*alp_Z2(I)

  if (gpmax1.ge.gpmax2) then
    gpmax = gpmax1
  else
    gpmax = gpmax2
  endif

end do

delta_w = 1000.d0 * gpmax
vthermal = c0*delta_w/(w0*1000.d0) !extra 1000.d0 because input wants it in km/s
dw = (delta_w*x2*x2)/nw

dv = c0*dw/w0 

!reading seems okay.
!!check if everything is read well...
!do I=1,nF
!  write(*,*)'Fu-Fl',Ftrans(I,1),Ftrans(I,2)
!  write(*,*)'gu-gl',alp_Z1(I),alp_Z2(I)
!  write(*,*)'A,w0',Aij_all(I),w00(I)
!end do
!!some other constants
!write(*,*)'w0,B,theta,ds',w0,B,theta,ds
!stop

allocate(Awtot(-nw:nw))
allocate(Bwtot(-nw:nw))
allocate(Cwtot(-nw:nw))

allocate(Aw(-nw:nw))
allocate(Bw(-nw:nw))
allocate(Cw(-nw:nw))

allocate(Iw(-nw:nw))
allocate(Qw(-nw:nw))
allocate(Vw(-nw:nw))

allocate(Ipv_r(-nw:nw))
allocate(Upv_r(-nw:nw))
allocate(Qpv_r(-nw:nw))
allocate(Vpv_r(-nw:nw))

allocate(lbd_a(-nw:nw,nF))
allocate(lbd_b(-nw:nw,nF))

allocate(wall_a(mlF1,nF))
allocate(wall_b(mlF2,nF))


! mlF1: max(lF1)
allocate(all_rho_a_r(mlF1,nF,-nw:nw))
allocate(all_rho_b_r(mlF2,nF,-nw:nw))

allocate(rho_in_a_r(mlF1,nF,-nw:nw))
allocate(rho_in_b_r(mlF2,nF,-nw:nw))

call initiate(all_rho_a_r,all_rho_b_r,mlF1,mlF2,nF,nw,Ftrans,w0,dw,&
&mass,w00,alp_Z1,alp_Z2,B,wall_a,wall_b,gam_a,gam_b,Iw,Qw,Vw,Tb,vthermal,lam1,lam2,lbd_a,lbd_b)


!initiate results
open(unit=13,file="initiate.dat")
do vj =-nw,nw
  write(13,*)vj*dw*c0/w0,Iw(vj),lbd_a(vj,1),lbd_a(vj,2),lbd_a(vj,3)
end do
close(unit=13)
!stop

allocate(all_delta_I(nF,mlF1,Dmm))
allocate(all_delta_Q(nF,mlF1,Dmm))
allocate(all_delta_U(nF,mlF1,Dmm))
allocate(all_delta_V(nF,mlF1,Dmm))

do I = 1,nF
  F1 = Ftrans(I,1)
  F2 = Ftrans(I,2)

  lF1 = F1*2 + 1

  do a1 = 1,lF1 
    do bc1 = 1,3
      all_delta_I(I,a1,bc1) = 0.d0
      all_delta_Q(I,a1,bc1) = 0.d0
      all_delta_U(I,a1,bc1) = 0.d0
      all_delta_V(I,a1,bc1) = 0.d0
    end do
  end do
end do

do I=1,nF
  F1 = Ftrans(I,1)
  F2 = Ftrans(I,2)
  Aij = Aij_all(I)
  lF1 = F1*2 + 1

  IF (ALLOCATED (delta_I)) DEALLOCATE (delta_I)
  IF (ALLOCATED (delta_Q)) DEALLOCATE (delta_Q)
  IF (ALLOCATED (delta_U)) DEALLOCATE (delta_U)
  IF (ALLOCATED (delta_V)) DEALLOCATE (delta_V)

  allocate(delta_I(lF1,Dmm))
  allocate(delta_Q(lF1,Dmm))
  allocate(delta_U(lF1,Dmm))
  allocate(delta_V(lF1,Dmm))


  call gen_delta(theta,Aij,w0,F1,lF1,F2,delta_I,delta_Q,delta_U,delta_V)

  do a1 = 1,lF1 
    do bc1 = 1,3
      all_delta_I(I,a1,bc1) = delta_I(a1,bc1) 
      all_delta_Q(I,a1,bc1) = delta_Q(a1,bc1) 
      all_delta_U(I,a1,bc1) = delta_U(a1,bc1) 
      all_delta_V(I,a1,bc1) = delta_V(a1,bc1) 
    end do
  end do
end do


DEALLOCATE (delta_I)
DEALLOCATE (delta_Q)
DEALLOCATE (delta_U)
DEALLOCATE (delta_V)

STEP = 0
111 STEP = STEP + 1




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!! RATE EQUATIONS !!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


do I = 1,nF

  F1 = Ftrans(I,1)
  F2 = Ftrans(I,2)
  lF1 = F1*2 + 1
  lF2 = F2*2 + 1

  IF (ALLOCATED (delta_I)) DEALLOCATE (delta_I)
  IF (ALLOCATED (delta_Q)) DEALLOCATE (delta_Q)
  IF (ALLOCATED (delta_U)) DEALLOCATE (delta_U)
  IF (ALLOCATED (delta_V)) DEALLOCATE (delta_V)
  
  allocate(delta_I(lF1,3))
  allocate(delta_Q(lF1,3))
  allocate(delta_U(lF1,3))
  allocate(delta_V(lF1,3))

  IF (ALLOCATED (wa)) DEALLOCATE (wa)
  IF (ALLOCATED (wb)) DEALLOCATE (wb)

  allocate(wa(lF1))
  allocate(wb(lF2))

  do a1=1,lF1
    wa(a1) = wall_a(a1,I)
  end do

  do b1=1,lF2
    wb(b1) = wall_b(b1,I)
  end do

  do a1 = 1,lF1
    do bc1 = 1,3
      delta_I(a1,bc1) = all_delta_I(I,a1,bc1)
      delta_Q(a1,bc1) = all_delta_Q(I,a1,bc1)
      delta_V(a1,bc1) = all_delta_V(I,a1,bc1)
    end do
  end do

  IF (ALLOCATED (Ipabv_r)) DEALLOCATE (Ipabv_r)
  IF (ALLOCATED (Qpabv_r)) DEALLOCATE (Qpabv_r)
  IF (ALLOCATED (Vpabv_r)) DEALLOCATE (Vpabv_r)

  allocate(Ipabv_r(lF1,3,-nw:nw))
  allocate(Qpabv_r(lF1,3,-nw:nw))
  allocate(Vpabv_r(lF1,3,-nw:nw))


  do a1=1,lF1
    ma1 = a1 - F1 - 1 
    do bc1=1,3
      mb1 = ma1 + (bc1 - 2)
      b1 = mb1 + F2 + 1
      
      call int_stokes(Iw,Qw,Vw,lF1,lF2,a1,b1,nw,w0,dw,gam_a,wa,wb,Ipv_r,Qpv_r,Vpv_r)

      do vj=-nw,nw
        Ipabv_r(a1,bc1,vj) = Ipv_r(vj) 
        Qpabv_r(a1,bc1,vj) = Qpv_r(vj) 
        Vpabv_r(a1,bc1,vj) = Vpv_r(vj) 
      end do
    end do
  end do

  IF (ALLOCATED (Ip_r)) DEALLOCATE (Ip_r)
  IF (ALLOCATED (Qp_r)) DEALLOCATE (Qp_r)
  IF (ALLOCATED (Vp_r)) DEALLOCATE (Vp_r)

  allocate(Ip_r(lF1,3))
  allocate(Qp_r(lF1,3))
  allocate(Vp_r(lF1,3))


  IF (ALLOCATED (rho_old_a_r)) DEALLOCATE (rho_old_a_r)
  IF (ALLOCATED (rho_old_b_r)) DEALLOCATE (rho_old_b_r)

  allocate(rho_old_a_r(lF1))
  allocate(rho_old_b_r(lF2))

  vj = -nw-1
2 vj = vj+1

  if (vj.le.nw) then  
    do a1=1,lF1
      do bc1=1,3
        Ip_r(a1,bc1) = Ipabv_r(a1,bc1,vj) 
        Qp_r(a1,bc1) = Qpabv_r(a1,bc1,vj) 
        Vp_r(a1,bc1) = Vpabv_r(a1,bc1,vj) 
      end do
    end do

    Rthr = tthresh 

    NBREAK = 0
    RAT = 0
1 RAT = RAT + 1 ! 
      if (RAT.gt.MAXRHO) then

        write(*,*)'stopped at '
        write(*,*)vj
        stop

      endif

      do a1 = 1,lF1
          rho_old_a_r(a1) = all_rho_a_r(a1,I,vj) 
      end do
  
      do b1 = 1,lF2
          rho_old_b_r(b1) = all_rho_b_r(b1,I,vj)
      end do
  
  ! paa density elements  

      do a1 = 1,lF1
        lbd = lbd_a(vj,I)

        call stateq_aa(rho_old_b_r,lF1,lF2,a1,lbd,gam_a,&
&sa,Ip_r,Qp_r,Vp_r,delta_I,delta_Q,delta_U,delta_V,rho_new_r) 

        all_rho_a_r(a1,I,vj) = rho_new_r
      end do

  ! pbb density elements 
  
      do b1 = 1,lF2
        lbd = lbd_b(vj,I)

        call stateq_bb(rho_old_a_r,lF1,lF2,b1,lbd,gam_b,&
&sa,Ip_r,Qp_r,Vp_r,delta_I,delta_Q,delta_U,delta_V,rho_new_r)
 
        all_rho_b_r(b1,I,vj) = rho_new_r
      end do

      do a1 = 1,lF1   
          tmpa1 = rho_old_a_r(a1) 
          tmpa2 = all_rho_a_r(a1,I,vj) 
          if (abs(tmpa1-tmpa2).gt.Rthr) then
            loop = loop + 1
            goto 1
          endif
      end do
    
      do b1 = 1,lF2
          tmpb1 = rho_old_b_r(b1)
          tmpb2 = all_rho_b_r(b1,I,vj)
          if (abs(tmpb1-tmpb2).gt.Rthr) then
            loop = loop + 1
            goto 1
          endif
      end do
  goto 2
  endif

end do

IF (ALLOCATED (rho_old_a_r)) DEALLOCATE (rho_old_a_r)
IF (ALLOCATED (rho_old_b_r)) DEALLOCATE (rho_old_b_r)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!! TRANSFER ELEMENTS !!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! transfer elements
do vj = -nw,nw
  Awtot(vj) = 0.d0 
  Bwtot(vj) = 0.d0 
  Cwtot(vj) = 0.d0 
end do

do I = 1,nF

  F1 = Ftrans(I,1)
  F2 = Ftrans(I,2)
  lF1 = F1*2 + 1
  lF2 = F2*2 + 1
 
  IF (ALLOCATED (wa)) DEALLOCATE (wa)
  IF (ALLOCATED (wb)) DEALLOCATE (wb)

  allocate(wa(lF1))
  allocate(wb(lF2))

  do a1=1,lF1
    wa(a1) = wall_a(a1,I)
  end do

  do b1=1,lF2
    wb(b1) = wall_b(b1,I)
  end do

  IF (ALLOCATED (rho_aa_r)) DEALLOCATE (rho_aa_r)
  IF (ALLOCATED (rho_bb_r)) DEALLOCATE (rho_bb_r)


  allocate(rho_aa_r(lF1,-nw:nw))
  allocate(rho_bb_r(lF2,-nw:nw))

  do vj = -nw,nw
    do a1 = 1,lF1
      do a2 = 1,lF1
        rho_aa_r(a1,vj) = all_rho_a_r(a1,I,vj)
      end do
    end do
  
    do b1 = 1,lF2
      do b2 = 1,lF2
        rho_bb_r(b1,vj) = all_rho_b_r(b1,I,vj)
      end do
    end do
  end do

  IF (ALLOCATED (gab)) DEALLOCATE (gab)

  allocate(gab(lF1,3,-nw:nw))

  do vj=-nw,nw
    do a1=1,lF1
      do K=1,3
        gab(a1,K,vj) = 0.d0
      end do
    end do
  end do

  call int_densities(rho_aa_r,rho_bb_r,lF1,lF2,nw,w0,dv,gam_a,gam_b,wa,wb,gab)

!  write(*,*)gab(1,3,0)
!  write(*,*)gab(2,2,0)
!  write(*,*)gab(3,1,0)
!
!  write(*,*)
!
!  write(*,*)gab(1,1,0),gab(1,2,0),gab(1,3,0)
!  write(*,*)gab(2,1,0),gab(2,2,0),gab(2,3,0)
!  write(*,*)gab(3,1,0),gab(3,2,0),gab(3,3,0)

  call trans_el(gab,lF1,lF2,nw,delta_I,delta_Q,delta_U,delta_V,Aw,Bw,Cw,w0)

!  open(unit=150,file="prop.dat")
!  do vj=-nw,nw
!    write(150,*)Aw(vj),Bw(vj),Cw(vj)
!  end do
!  close(unit=150)
!  stop

  do vj = -nw,nw
    Awtot(vj) = Awtot(vj) + Aw(vj) 
    Bwtot(vj) = Bwtot(vj) + Bw(vj) 
    Cwtot(vj) = Cwtot(vj) + Cw(vj) 
  end do


end do


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!! STOKES PROPAGATION !!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
call prop_stokes(Awtot,Bwtot,Cwtot,Iw,Qw,Vw,ds,nw)

write(*,*) 'step:', STEP

! if any of the densities is not a number, we might as well stop.

if (isnan(all_rho_a_r(1,1,1)))  stop   
if (isnan(all_rho_b_r(1,1,1)))  stop

! some further analysis and writing of the Stokes parameters

write(51,*)Iw  
write(53,*)Qw
write(54,*)Vw
!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! compute rate of stimulated emission !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
Rtot = 0.d0
Ttot = 0.d0
do I = 1,nF
  F1 = Ftrans(I,1)
  F2 = Ftrans(I,2)

  Aij = Aij_all(I) 

  lF1 = F1*2 + 1
  lF2 = F2*2 + 1

  IF (ALLOCATED (delta_I)) DEALLOCATE (delta_I)
  IF (ALLOCATED (delta_Q)) DEALLOCATE (delta_Q)
  IF (ALLOCATED (delta_U)) DEALLOCATE (delta_U)
  IF (ALLOCATED (delta_V)) DEALLOCATE (delta_V)

  allocate(delta_I(lF1,3))
  allocate(delta_Q(lF1,3))
  allocate(delta_U(lF1,3))
  allocate(delta_V(lF1,3))

!  call gen_delta(theta,Aij,w0,F1,lF1,F2,delta_I,delta_Q,delta_U,delta_V)
  do a1 = 1,lF1
    do bc1 = 1,3
      delta_I(a1,bc1) = all_delta_I(I,a1,bc1)
      delta_Q(a1,bc1) = all_delta_Q(I,a1,bc1)
      delta_U(a1,bc1) = all_delta_U(I,a1,bc1)
      delta_V(a1,bc1) = all_delta_V(I,a1,bc1)
    end do
  end do

  call stimulate(Iw(0),Qw(0),Vw(0),lF1,lF2,Aij,w0,delta_I,delta_Q,delta_U,delta_V,Rstim,Tbr)
!  call stimulate_alt(Iw(0),Aij,sa,w0,Rstim)

  Rtot = Rtot + Rstim
  Ttot = Ttot + Tbr
end do

fwhm_old = fwhm_i
!call fwhm(Iw,nw,dw,fwhm_i)
call zeeman(Iw,Vw,nw,dw,azee,bzee)
call vitals(Iw,Qw,Vw,nw,pQ,pV)

mu = x0
sig = x0
Is = x0
call gaussfit(Iw,nw,dw,mu,sig,Is)

!Rthr = tthresh / Rtot

!write(25,*)pQ
!write(26,*)pV
!write(27,*)pa

!compute the temperature

dR=dlog10(Rtot/gam_a)

write(28,*)dR,pQ,pV

write(55,*)Rtot,Ttot
write(56,*)sig,sig*c0/w0,mu,Is
write(58,*)azee,bzee

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!            and... go on             !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

if (STEP.le.MAXSTEP) GOTO 111    

! write final stokes paramaters
!open(unit=12,file="stokes.dat")
!do vj =-nw,nw 
!  write(12,*)vj*dw*c0/w0,Iw(vj),Qw(vj),Uw(vj),Vw(vj)
!end do
!close(unit=12)

write(57,*)STEP,gam_a,gam_b,lam1,lam2,sa,Rtot

!close(25)
!close(26)
!close(27)
close(28)

close(51)
close(52)
close(53)
close(54)
close(55)
close(56)
close(57)
close(58)

close(61)
close(62)
close(63)
close(64)
close(65)
close(66)
close(67)
close(68)

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
 
