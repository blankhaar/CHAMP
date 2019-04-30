program watned_92
implicit none
integer F1,F2,lF1,lF2,mlF1,mlF2,Dmm
integer I,J,K,L,MAXSTEP,STEP,RAT,NBREAK,nimax
integer ma1,ma2,mb1,mb2,a1,a2,b1,b2,bc1,bc2
integer vj,loop,nw,nF,MAXRHO,THRESH
integer vvj,II,inp
double precision theta,w0,rho_new_r,rho_new_i,om,dw,dv,Rtot,Ttot,Rstim,Rmax,dR,B
double precision tmpa1,tmpa2,tmpb1,tmpb2,Tb,gam,vthermal,lam2,c0,fwhm_i
double precision Aij,lbd,sa,mass,ds,x3,x2,pi,fwhm_old,bzee,azee,k_sb
double precision lbda,lbdb,delv,Ap,pVnorm,vwidth
double precision rata_1,rata_2,rata_3,ratb_1,ratb_2,ratb_3
double precision ratsa_2,ratsa_3,ratsb_2,ratsb_3,tthresh,Rthr
double precision pQ,pV,pa,Tbr,Told 
double precision mu,sig,Is,x0,maxT 
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
double precision, dimension(:), allocatable :: rho_old_a_r,rho_old_b_r
double precision, dimension(:), allocatable :: lam1,del,gam_a,gam_b 
double precision, dimension(:,:), allocatable :: rho_aa_r,rho_bb_r
double precision, dimension(:), allocatable :: Aw,Bw,Cw,omeg,vall
double precision, dimension(:,:), allocatable :: lbd_a,lbd_b 
double precision, dimension(:), allocatable :: Awtot,Bwtot,Cwtot
double precision, dimension(:,:,:), allocatable :: gab
double precision, dimension(:,:), allocatable :: delta_I,delta_Q,delta_U,delta_V


Dmm = 3

c0 = 299792458.d0
x0 = 0.d0
k_sb = 1.38064852E-23 

Told = 0.d0

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
allocate(lam1(nF))
allocate(del(nF))
allocate(gam_a(nF))
allocate(gam_b(nF))

do I=1,nF
  read(4,*)Ftrans(I,1),Ftrans(I,2),w00(I),alp_Z1(I),alp_Z2(I),Aij_all(I),lam1(I),del(I),gam_a(I),gam_b(I)
end do

read(4,*)

read(4,*)w0
read(4,*)B,theta
read(4,*)nw,vthermal,vwidth,mass
read(4,*)Tb,ds,MAXSTEP
read(4,*)sa,maxT
read(4,*)inp

if (inp.eq.1) then
  open(unit=51,file="I.dat")
  open(unit=59,file="vgrid.dat")
  open(unit=53,file="Q.dat")
  open(unit=54,file="V.dat")
endif
x2 = 2.d0
x3 = 3.d0
pi = dacos(-1.d0)

maxT = 10.d0**(maxT)

theta = theta*pi/180.d0
w0 = x2*pi*w0

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

dv = 1000.d0*vwidth / (2*nw+1)
dw = w0 * dv / c0 

Ap = (13.4d0 + 8.3d0 + 1.0d0)/(3.d0*1000.d0)

pVnorm = 10.d0**(-6.0) * x2 * Ap * B *dcos(theta*pi/180.d0) 

allocate(Awtot(-nw:nw))
allocate(Bwtot(-nw:nw))
allocate(Cwtot(-nw:nw))

allocate(Aw(-nw:nw))
allocate(Bw(-nw:nw))
allocate(Cw(-nw:nw))

allocate(Iw(-nw:nw))
allocate(Qw(-nw:nw))
allocate(Vw(-nw:nw))
allocate(omeg(-nw:nw))
allocate(vall(-nw:nw))

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

call initiate(mlF1,mlF2,nF,nw,Ftrans,w0,dw,&
&mass,w00,alp_Z1,alp_Z2,B,wall_a,wall_b,Iw,Qw,Vw,Tb,vthermal,lam1,del,lbd_a,lbd_b)

do vj = -nw,nw
  omeg(vj) = vj*dw
  vall(vj) = vj*dv
enddo

if (inp.eq.1)then
  write(59,*)vall
endif

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
      
      if (b1.gt.0.and.b1.le.lF2) then 
        call int_stokes(Iw,Qw,Vw,lF1,lF2,a1,b1,nw,w0,dw,gam_a(I),wa,wb,Ipv_r,Qpv_r,Vpv_r)
      else
        Ipv_r(:) = x0
        Qpv_r(:) = x0
        Vpv_r(:) = x0
      endif

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

  do vj = -nw,nw

    do a1=1,lF1
      do bc1=1,3
        Ip_r(a1,bc1) = Ipabv_r(a1,bc1,vj) 
        Qp_r(a1,bc1) = Qpabv_r(a1,bc1,vj) 
        Vp_r(a1,bc1) = Vpabv_r(a1,bc1,vj) 
      end do
    end do

    lbda = lbd_a(vj,I)
    lbdb = lbd_b(vj,I)
    call stateq(lF1,lF2,lbda,lbdb,gam_a(I),gam_b(I),sa,Ip_r,Qp_r,Vp_r,delta_I,&
&delta_Q,delta_U,delta_V,rho_old_a_r,rho_old_b_r)

    do a1 = 1,lF1
      all_rho_a_r(a1,I,vj) = rho_old_a_r(a1) 
    end do

    do b1 = 1,lF2
      all_rho_b_r(b1,I,vj) = rho_old_b_r(b1) 
    end do

  end do

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
      rho_aa_r(a1,vj) = all_rho_a_r(a1,I,vj)
    end do
  
    do b1 = 1,lF2
      rho_bb_r(b1,vj) = all_rho_b_r(b1,I,vj)
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

  call int_densities(rho_aa_r,rho_bb_r,lF1,lF2,nw,w0,dv,wa,wb,gab)

  IF (ALLOCATED (delta_I)) DEALLOCATE (delta_I)
  IF (ALLOCATED (delta_Q)) DEALLOCATE (delta_Q)
  IF (ALLOCATED (delta_U)) DEALLOCATE (delta_U)
  IF (ALLOCATED (delta_V)) DEALLOCATE (delta_V)

  allocate(delta_I(lF1,3))
  allocate(delta_Q(lF1,3))
  allocate(delta_U(lF1,3))
  allocate(delta_V(lF1,3))

  do a1 = 1,lF1
    do bc1 = 1,3
      delta_I(a1,bc1) = all_delta_I(I,a1,bc1)
      delta_Q(a1,bc1) = all_delta_Q(I,a1,bc1)
      delta_V(a1,bc1) = all_delta_V(I,a1,bc1)
    end do
  end do

  call trans_el(gab,lF1,lF2,nw,delta_I,delta_Q,delta_V,Aw,Bw,Cw,w0)

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


! if any of the densities is not a number, we might as well stop.

if (isnan(all_rho_a_r(1,1,1)))  stop   
if (isnan(all_rho_b_r(1,1,1)))  stop

! some further analysis and writing of the Stokes parameters

if (mod(STEP,10).eq.0.and.inp.eq.1) then
  write(51,*)Iw  
  write(53,*)Qw
  write(54,*)Vw
endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! compute rate of stimulated emission !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
call vitals(Iw,Qw,Vw,nw,pQ,pV,nimax)

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

  do a1 = 1,lF1
    do bc1 = 1,3
      delta_I(a1,bc1) = all_delta_I(I,a1,bc1)
      delta_Q(a1,bc1) = all_delta_Q(I,a1,bc1)
      delta_U(a1,bc1) = all_delta_U(I,a1,bc1)
      delta_V(a1,bc1) = all_delta_V(I,a1,bc1)
    end do
  end do

  call stimulate(Iw(nimax),Qw(nimax),Vw(nimax),lF1,lF2,Aij,w0,delta_I,delta_Q,delta_U,delta_V,Rstim,Tbr)
!  call stimulate_alt(Iw(0),Aij,sa,w0,Rstim)

  Tbr = 4.d0*(pi**3.d0)*(c0**x2)*Iw(nimax)/(k_sb*(w0**x2))

  Rtot = Rtot + Rstim
  Ttot = Ttot + Tbr
end do

call fwhm(Iw,nw,dv,delv)

dR=dlog10(Rtot)

if (mod(STEP,10).eq.0.and.Ttot.lt.maxT) then
  write(28,*)Ttot,dR,delv,pQ,pV
endif

Told = Ttot

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!            and... go on             !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

ds = ds * 1.0005d0

if (Ttot.gt.maxT) stop


if (STEP.le.MAXSTEP) GOTO 111    

close(28)

if (inp.eq.1) then
  close(51)
  close(59)
  close(53)
  close(54)
endif


end program
 
