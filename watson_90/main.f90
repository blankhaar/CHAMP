program watned_94
implicit none
integer F1,F2,lF1,lF2,mlF1,mlF2,Dmm
integer I,J,K,L,MAXSTEP,STEP,RAT,NBREAK
integer a1,a2,b1,b2,bc1,bc2
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
double precision Aw,Bw,Cw,Dww,Ew,Fw,Gw,Aws,Bws,Cws,Fws 
double precision Awtot,Bwtot,Cwtot,Dwtot,Ewtot,Fwtot,Gwtot 
double precision Iw,Qw,Uw,Vw 
integer, dimension(:,:), allocatable :: Ftrans 
double precision, dimension(:,:), allocatable :: wall_a,wall_b 
double precision, dimension(:), allocatable :: wa,wb 
double precision, dimension(:), allocatable :: Aij_all,w00,alp_Z1,alp_Z2 
double precision, dimension(:,:,:,:,:), allocatable :: all_delta_I,all_delta_Q,all_delta_U,all_delta_V
double precision, dimension(:,:), allocatable :: rho_aa_r,rho_aa_i,rho_bb_r,rho_bb_i 
double precision, dimension(:,:), allocatable :: rho_old_a_r,rho_old_a_i,rho_old_b_r,rho_old_b_i 
double precision, dimension(:,:,:), allocatable :: all_rho_a_r,all_rho_a_i,all_rho_b_r,all_rho_b_i 
double precision, dimension(:,:,:,:), allocatable :: delta_I,delta_Q,delta_U,delta_V
double precision, dimension(:,:), allocatable ::  lbd_a,lbd_b 


Dmm = 3

c0 = 299792458.d0
x0 = 0.d0

!open(unit=61,file="pa00.dat")
!open(unit=63,file="pa11.dat")
!open(unit=64,file="pam1m1.dat")
!open(unit=65,file="pa01.dat")
!open(unit=66,file="pa1m1.dat")
!open(unit=69,file="pa0m1.dat")
!open(unit=62,file="pb.dat")

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
read(4,*)lam1,lam2,gam_a,gam_b
read(4,*)Tb,ds,MAXSTEP,MAXRHO,THRESH
read(4,*)sa

x2 = 2.d0
pi = dacos(-1.d0)

tthresh = 10.d0**(-1.d0*THRESH)
Rthr = tthresh

theta = theta*pi/180.d0
w0 = x2*pi*w0

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

allocate(wall_a(mlF1,nF))
allocate(wall_b(mlF2,nF))

allocate(lbd_a(mlF1,nF))
allocate(lbd_b(mlF2,nF))

! mlF1: max(lF1)
allocate(all_rho_a_r(mlF1,mlF1,nF))
allocate(all_rho_a_i(mlF1,mlF1,nF))
allocate(all_rho_b_r(mlF2,mlF2,nF))
allocate(all_rho_b_i(mlF2,mlF2,nF))

call initiate(all_rho_a_r,all_rho_a_i,all_rho_b_r,all_rho_b_i,mlF1,mlF2,nF,Ftrans,w0,&
&w00,alp_Z1,alp_Z2,B,wall_a,wall_b,gam_a,gam_b,Iw,Uw,Qw,Vw,Tb,lam1,lam2,lbd_a,lbd_b)


!initiate results
!open(unit=13,file="initiate.dat")
!write(13,*)Iw,lbd_a(1,1),lbd_a(2,1),lbd_a(3,1),lbd_b(1,1)
!close(unit=13)

allocate(all_delta_I(nF,mlF1,Dmm,mlF1,Dmm))
allocate(all_delta_Q(nF,mlF1,Dmm,mlF1,Dmm))
allocate(all_delta_U(nF,mlF1,Dmm,mlF1,Dmm))
allocate(all_delta_V(nF,mlF1,Dmm,mlF1,Dmm))

do I = 1,nF
  F1 = Ftrans(I,1)
  F2 = Ftrans(I,2)

  lF1 = F1*2 + 1

  do a1 = 1,lF1 
    do a2 = 1,lF1
      do bc1 = 1,3
        do bc2 = 1,3
          all_delta_I(I,a1,bc1,a2,bc2) = 0.d0
          all_delta_Q(I,a1,bc1,a2,bc2) = 0.d0
          all_delta_U(I,a1,bc1,a2,bc2) = 0.d0
          all_delta_V(I,a1,bc1,a2,bc2) = 0.d0
        end do
      end do
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

  allocate(delta_I(lF1,Dmm,lF1,Dmm))
  allocate(delta_Q(lF1,Dmm,lF1,Dmm))
  allocate(delta_U(lF1,Dmm,lF1,Dmm))
  allocate(delta_V(lF1,Dmm,lF1,Dmm))

  call gen_delta(theta,Aij,w0,F1,lF1,F2,delta_I,delta_Q,delta_U,delta_V)

  do a1 = 1,lF1 
    do a2 = 1,lF1 
      do bc1 = 1,3
        do bc2 = 1,3
          all_delta_I(I,a1,bc1,a2,bc2) = delta_I(a1,bc1,a2,bc2) 
          all_delta_Q(I,a1,bc1,a2,bc2) = delta_Q(a1,bc1,a2,bc2) 
          all_delta_U(I,a1,bc1,a2,bc2) = delta_U(a1,bc1,a2,bc2) 
          all_delta_V(I,a1,bc1,a2,bc2) = delta_V(a1,bc1,a2,bc2) 
        end do
      end do
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
  
  allocate(delta_I(lF1,3,lF1,3))
  allocate(delta_Q(lF1,3,lF1,3))
  allocate(delta_U(lF1,3,lF1,3))
  allocate(delta_V(lF1,3,lF1,3))

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
      do a2 = 1,lF1
        do bc2 = 1,3
          delta_I(a1,bc1,a2,bc2) = all_delta_I(I,a1,bc1,a2,bc2)
          delta_Q(a1,bc1,a2,bc2) = all_delta_Q(I,a1,bc1,a2,bc2)
          delta_U(a1,bc1,a2,bc2) = all_delta_U(I,a1,bc1,a2,bc2)
          delta_V(a1,bc1,a2,bc2) = all_delta_V(I,a1,bc1,a2,bc2)
        end do
      end do
    end do
  end do


  IF (ALLOCATED (rho_old_a_r)) DEALLOCATE (rho_old_a_r)
  IF (ALLOCATED (rho_old_a_i)) DEALLOCATE (rho_old_a_i)
  IF (ALLOCATED (rho_old_b_r)) DEALLOCATE (rho_old_b_r)
  IF (ALLOCATED (rho_old_b_i)) DEALLOCATE (rho_old_b_i)

  allocate(rho_old_a_r(lF1,lF1))
  allocate(rho_old_a_i(lF1,lF1))
  allocate(rho_old_b_r(lF2,lF2))
  allocate(rho_old_b_i(lF2,lF2))

  Rthr = tthresh 

  NBREAK = 0
  RAT = 0
1 RAT = RAT + 1 ! 
      if (RAT.gt.MAXRHO) then
        write(*,*)'rho did not converge'
        stop

      endif

      do a1 = 1,lF1
        do a2 = 1,lF1
          rho_old_a_r(a1,a2) = all_rho_a_r(a1,a2,I) 
          rho_old_a_i(a1,a2) = all_rho_a_i(a1,a2,I) 
        end do
      end do
  
      do b1 = 1,lF2
        do b2 = 1,lF2
          rho_old_b_r(b1,b2) = all_rho_b_r(b1,b2,I)
          rho_old_b_i(b1,b2) = all_rho_b_i(b1,b2,I)
        end do
      end do
  
  ! paa density elements  

      do a1 = 1,lF1
        do a2 = a1,lF1   
           om = wa(a1) - wa(a2) 
           
           lbd = x0
           if (a1.eq.a2) then
             lbd = lbd_a(a1,I)
           endif

           call stateq_aa(rho_old_a_r,rho_old_a_i,rho_old_b_r,rho_old_b_i,lF1,lF2,a1,a2,lbd,gam_a,&
&om,sa,Iw,Qw,Uw,Vw,delta_I,delta_Q,delta_U,delta_V,rho_new_r,rho_new_i) 


           if (a1.eq.a2) then
             rho_new_i = x0 
           endif

           all_rho_a_r(a1,a2,I) = rho_new_r
           all_rho_a_r(a2,a1,I) = rho_new_r

           all_rho_a_i(a1,a2,I) = rho_new_i
           all_rho_a_i(a2,a1,I) = -rho_new_i
  
        end do
      end do

  ! pbb density elements 
  
      do b1 = 1,lF2
        do b2 = b1,lF2
          om = wb(b1) - wb(b2)

           lbd = 0.d0
           if (b1.eq.b2) then
             lbd = lbd_b(b1,I)
           endif

          call stateq_bb(rho_old_a_r,rho_old_a_i,rho_old_b_r,rho_old_b_i,lF1,lF2,b1,b2,lbd,gam_b,&
&om,sa,Iw,Qw,Uw,Vw,delta_I,delta_Q,delta_U,delta_V,rho_new_r,rho_new_i)
 
           if (b1.eq.b2) then
             rho_new_i = 0.d0
           endif


          all_rho_b_r(b1,b2,I) = rho_new_r
          all_rho_b_i(b1,b2,I) = rho_new_i
  
          all_rho_b_r(b2,b1,I) = rho_new_r
          all_rho_b_i(b2,b1,I) = -rho_new_i
        end do
      end do

      stop

      do a1 = 1,lF1   
!        do a2 = 1,lF1
         a2 = a1 
          tmpa1 = rho_old_a_r(a1,a2) 
          tmpa2 = all_rho_a_r(a1,a2,I) 
          if (dabs(tmpa1-tmpa2).gt.Rthr) then
            loop = loop + 1
            goto 1
          endif
!        end do
      end do
    
      do b1 = 1,lF2
!        do b2 = 1,lF2
          b2 = b1 
          tmpb1 = rho_old_b_r(b1,b2)
          tmpb2 = all_rho_b_r(b1,b2,I)
          if (dabs(tmpb1-tmpb2).gt.Rthr) then
            loop = loop + 1
            goto 1
          endif
!        end do
      end do

!     write(*,*)vj,RAT
!     write(*,*)
!     write(*,*)all_rho_a_r(1,1,1,vj),all_rho_a_r(1,2,1,vj),all_rho_a_r(1,3,1,vj)
!     write(*,*)all_rho_a_r(2,1,1,vj),all_rho_a_r(2,2,1,vj),all_rho_a_r(2,3,1,vj)
!     write(*,*)all_rho_a_r(3,1,1,vj),all_rho_a_r(3,2,1,vj),all_rho_a_r(3,3,1,vj)
!     write(*,*)
!     write(*,*)all_rho_a_i(1,1,1,vj),all_rho_a_i(1,2,1,vj),all_rho_a_i(1,3,1,vj)
!     write(*,*)all_rho_a_i(2,1,1,vj),all_rho_a_i(2,2,1,vj),all_rho_a_i(2,3,1,vj)
!     write(*,*)all_rho_a_i(3,1,1,vj),all_rho_a_i(3,2,1,vj),all_rho_a_i(3,3,1,vj)
!     write(*,*)
!     write(*,*)all_rho_b_r(1,1,1,vj)
!     write(*,*)
!    stop


end do

IF (ALLOCATED (rho_old_a_r)) DEALLOCATE (rho_old_a_r)
IF (ALLOCATED (rho_old_a_i)) DEALLOCATE (rho_old_a_i)
IF (ALLOCATED (rho_old_b_r)) DEALLOCATE (rho_old_b_r)
IF (ALLOCATED (rho_old_b_i)) DEALLOCATE (rho_old_b_i)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!! TRANSFER ELEMENTS !!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! transfer elements
Awtot = 0.d0 
Bwtot = 0.d0 
Cwtot = 0.d0 
Dwtot = 0.d0 
Ewtot = 0.d0 
Fwtot = 0.d0 
Gwtot = 0.d0 

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
  IF (ALLOCATED (rho_aa_i)) DEALLOCATE (rho_aa_i)
  IF (ALLOCATED (rho_bb_r)) DEALLOCATE (rho_bb_r)
  IF (ALLOCATED (rho_bb_i)) DEALLOCATE (rho_bb_i)


  allocate(rho_aa_r(lF1,lF1))
  allocate(rho_aa_i(lF1,lF1))
  allocate(rho_bb_r(lF2,lF2))
  allocate(rho_bb_i(lF2,lF2))

  do a1 = 1,lF1
    do a2 = 1,lF1
      rho_aa_r(a1,a2) = all_rho_a_r(a1,a2,I)
      rho_aa_i(a1,a2) = all_rho_a_i(a1,a2,I)
    end do
  end do

  do b1 = 1,lF2
    do b2 = 1,lF2
      rho_bb_r(b1,b2) = all_rho_b_r(b1,b2,I)
      rho_bb_i(b1,b2) = all_rho_b_i(b1,b2,I)
    end do
  end do
 
  call trans_el(rho_aa_r,rho_aa_i,rho_bb_r,rho_bb_i,lF1,lF2,&
&delta_I,delta_Q,delta_U,delta_V,Aw,Bw,Cw,Dww,Ew,Fw,Gw,Aws,Bws,Cws,Fws,w0)

!  open(unit=150,file="prop.dat")
!    write(150,*)Aw,Bw,Cw,Dww,Ew,Fw,Gw
!  close(unit=150)
!  stop

  Awtot = Awtot + Aw 
  Bwtot = Bwtot + Bw 
  Cwtot = Cwtot + Cw 
  Dwtot = Dwtot + Dww 
  Ewtot = Ewtot + Ew 
  Fwtot = Fwtot + Fw 
  Gwtot = Gwtot + Gw 


end do


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!! STOKES PROPAGATION !!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! for now, without spontaneous emission

call prop_stokes(Awtot,Bwtot,Cwtot,Dwtot,Ewtot,Fwtot,Gwtot,Iw,Qw,Uw,Vw,ds)

write(*,*) 'step:', STEP

! if any of the densities is not a number, we might as well stop.

if (isnan(all_rho_a_r(1,1,1)))  stop   
if (isnan(all_rho_a_i(1,1,1)))  stop

if (isnan(all_rho_b_r(1,1,1)))  stop
if (isnan(all_rho_b_i(1,1,1)))  stop

!write(51,*)Iw  
!write(52,*)Uw
!write(53,*)Qw
!write(54,*)Vw
!
!write(61,*)all_rho_a_r(2,2,1)
!write(62,*)all_rho_b_r(1,1,1)
!write(63,*)all_rho_a_r(3,3,1) 
!write(64,*)all_rho_a_r(1,1,1) 
!write(65,*)all_rho_a_r(2,3,1),all_rho_a_i(2,3,1) 
!write(66,*)all_rho_a_r(1,3,1),all_rho_a_i(1,3,1) 
!write(69,*)all_rho_a_r(1,2,1),all_rho_a_i(1,2,1) 

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

  allocate(delta_I(lF1,3,lF1,3))
  allocate(delta_Q(lF1,3,lF1,3))
  allocate(delta_U(lF1,3,lF1,3))
  allocate(delta_V(lF1,3,lF1,3))

!  call gen_delta(theta,Aij,w0,F1,lF1,F2,delta_I,delta_Q,delta_U,delta_V)
  do a1 = 1,lF1
    do bc1 = 1,3
      do a2 = 1,lF1
        do bc2 = 1,3
          delta_I(a1,bc1,a2,bc2) = all_delta_I(I,a1,bc1,a2,bc2)
          delta_Q(a1,bc1,a2,bc2) = all_delta_Q(I,a1,bc1,a2,bc2)
          delta_U(a1,bc1,a2,bc2) = all_delta_U(I,a1,bc1,a2,bc2)
          delta_V(a1,bc1,a2,bc2) = all_delta_V(I,a1,bc1,a2,bc2)
        end do
      end do
    end do
  end do

  call stimulate(Iw,Qw,Uw,Vw,lF1,lF2,Aij,w0,delta_I,delta_Q,delta_U,delta_V,Rstim,Tbr)
!  call stimulate_alt(Iw(0),Aij,sa,w0,Rstim)

  Rtot = Rtot + Rstim
  Ttot = Ttot + Tbr
end do

dR = dlog10(Rtot/gam_a)
pQ = dsqrt(Qw**2.d0 + Uw**2.d0)/Iw
pa = datan2(Uw,Qw)*90.d0/pi
pV = Vw/Iw

write(28,*)dR,pQ,Qw/Iw,Uw/Iw,pa

!write(55,*)Rtot,Ttot
!write(56,*)sig,sig*c0/w0,mu,Is
!write(58,*)azee,bzee

ds = ds*1.001d0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!            and... go on             !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

if (STEP.le.MAXSTEP) GOTO 111    

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
 
