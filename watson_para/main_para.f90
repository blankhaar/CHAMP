program watned_94
implicit none
integer nF1,nF2,F1,F2,lF1,lF2,mlF1,mlF2,Dmm
integer I,J,K,L,MAXSTEP,STEP,RAT,NBREAK
integer a1,a2,b1,b2,bc1,bc2
integer vj,loop,nw,nF,shvel,nimax
integer vvj,II,anis,inci,inp
double precision theta,w0,rho_new_r,rho_new_i,om,dw,dv,Rtot,Ttot,Rstim,Rmax,dR,B
double precision tmpa1,tmpa2,tmpb1,tmpb2,Tb,gam,vthermal,vwidth,c0,fwhm_i
double precision Aij,lbd,sa,mass,ds,x2,pi,fwhm_old,bzee,azee
double precision pQ,pV,pa,Tbr,maxTb 
double precision mu,sig,Is,x0 
double precision gpmax1,gpmax2,gpmax,delta_w
double precision qin,qvel,qga,qi
integer, dimension(:,:), allocatable :: Ftrans 
double precision, dimension(:,:), allocatable :: wall_a,wall_b 
double precision, dimension(:), allocatable :: wa,wb 
double precision, dimension(:), allocatable :: Iw,Qw,Uw,Vw
double precision, dimension(:), allocatable :: Ipv_r,Ipv_i,Upv_r,Upv_i,Qpv_r,Qpv_i,Vpv_r,Vpv_i
double precision, dimension(:,:), allocatable :: Ip_r,Ip_i,Up_r,Up_i,Qp_r,Qp_i,Vp_r,Vp_i
double precision, dimension(:,:,:), allocatable :: Ipabv_r,Ipabv_i,Upabv_r,Upabv_i,Qpabv_r,Qpabv_i,Vpabv_r,Vpabv_i
double precision, dimension(:), allocatable :: Aij_all,w00,alp_Z1,alp_Z2 
double precision, dimension(:,:,:,:,:), allocatable :: all_delta_I,all_delta_Q,all_delta_U,all_delta_V
double precision, dimension(:,:,:,:), allocatable :: all_rho_a_r,all_rho_a_i,all_rho_b_r,all_rho_b_i 
double precision, dimension(:,:), allocatable :: rho_a_r,rho_a_i,rho_b_r,rho_b_i 
double precision, dimension(:,:,:), allocatable :: rho_aa_r,rho_aa_i,rho_bb_r,rho_bb_i 
double precision, dimension(:), allocatable :: Aw,Bw,Cw,Dww,Ew,Fw,Gw
double precision, dimension(:,:,:,:), allocatable :: lbd_a_r,lbd_b_r,lbd_a_i,lbd_b_i
double precision, dimension(:,:), allocatable :: lbda_r,lbda_i,lbdb_r,lbdb_i 
double precision, dimension(:), allocatable :: Awtot,Bwtot,Cwtot,Dwtot,Ewtot,Fwtot,Gwtot
double precision, dimension(:,:,:,:), allocatable :: gaa_p_r,gaa_p_i,gbb_p_r,gbb_p_i,gaa_m_r,gaa_m_i,gbb_m_r,gbb_m_i
double precision, dimension(:,:,:,:), allocatable :: delta_I,delta_Q,delta_U,delta_V
double precision, dimension(:), allocatable :: gam_a,gam_b,lam1,lam2,del,p_an_u,p_an_d
double precision, dimension(3) :: apump 

Dmm = 3

c0 = 299792458.d0
x0 = 0.d0

!open(unit=61,file="I.dat")
!open(unit=62,file="U.dat")
!open(unit=63,file="Q.dat")
!open(unit=64,file="V.dat")
open(unit=28,file="ned_props.dat")
!open(unit=29,file="vgrid.dat")


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
allocate(p_an_u(nF))
allocate(p_an_d(nF))

do I=1,nF
  read(4,*)Ftrans(I,1),Ftrans(I,2),w00(I),alp_Z1(I),alp_Z2(I),Aij_all(I),lam1(I),del(I),gam_a(I),gam_b(I),p_an_u(I),p_an_d(I)
end do

read(4,*)

read(4,*)w0
read(4,*)B,theta
read(4,*)nw,mass,Tb
read(4,*)ds,MAXSTEP,sa,maxTb
read(4,*)anis,apump(1),apump(2),apump(3)
read(4,*)inci,qin,qvel,qga,qi
read(4,*)shvel,vthermal,vwidth
read(4,*)inp


!write(*,*)alp_Z1
!write(*,*)alp_Z2
!
!stop
!
maxTb = 10.d0 ** (maxTb)

! if anis == 1 ----> anisotropic pumping
! if inci == 1 ----> incident radiation


!apump(:) : anistropic pumping direction vector
!p_an_u(I)

!qi: intensity incoming light = 10^(-qin)
!qvel: velocity shift of the Gaussian
!qga: width of the Gaussian = qga * dw
!qin: fraction of incoming light that is polarized I = (1-qi)* qin
!                                                  Q = qi * qin

if (inp.eq.1) then
  open(unit=61,file="I.dat",recl=8824)
  open(unit=62,file="U.dat",recl=8824)
  open(unit=63,file="Q.dat",recl=8824)
  open(unit=64,file="V.dat",recl=8824)
  open(unit=29,file="vgrid.dat")
endif


x2 = 2.d0
pi = dacos(-1.d0)

theta = theta*pi/180.d0
w0 = x2*pi*w0

do I = 1,nF
  w00(I) = w00(I) * x2 * pi
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


if (shvel.ne.1) then
  delta_w = 1000.d0 * gpmax
  vthermal = c0*delta_w/(w0*1000.d0) !extra 1000.d0 because input wants it in km/s
  dw = (delta_w*x2*x2)/nw
  dv = c0*dw/w0 
else
  dv = 1000.d0 * vwidth / (2*nw + 1)
  dw = dv * w0 / c0
endif

if (inp.eq.1) then
  do vj = -nw,nw
    write(29,*)dv*vj
  end do
  close(29)
endif

qi = 10.d0**(-qi)
qga = dw * qga

allocate(Awtot(-nw:nw))
allocate(Bwtot(-nw:nw))
allocate(Cwtot(-nw:nw))
allocate(Dwtot(-nw:nw))
allocate(Ewtot(-nw:nw))
allocate(Fwtot(-nw:nw))
allocate(Gwtot(-nw:nw))

allocate(Aw(-nw:nw))
allocate(Bw(-nw:nw))
allocate(Cw(-nw:nw))
allocate(Dww(-nw:nw))
allocate(Ew(-nw:nw))
allocate(Fw(-nw:nw))
allocate(Gw(-nw:nw))

allocate(Iw(-nw:nw))
allocate(Uw(-nw:nw))
allocate(Qw(-nw:nw))
allocate(Vw(-nw:nw))

allocate(Ipv_r(-nw:nw))
allocate(Ipv_i(-nw:nw))
allocate(Upv_r(-nw:nw))
allocate(Upv_i(-nw:nw))
allocate(Qpv_r(-nw:nw))
allocate(Qpv_i(-nw:nw))
allocate(Vpv_r(-nw:nw))
allocate(Vpv_i(-nw:nw))

allocate(lbd_a_r(-nw:nw,mlF1,mlF1,nF))
allocate(lbd_a_i(-nw:nw,mlF1,mlF1,nF))
allocate(lbd_b_r(-nw:nw,mlF2,mlF2,nF))
allocate(lbd_b_i(-nw:nw,mlF2,mlF2,nF))

allocate(wall_a(mlF1,nF))
allocate(wall_b(mlF2,nF))


! mlF1: max(lF1)
allocate(all_rho_a_r(mlF1,mlF1,nF,-nw:nw))
allocate(all_rho_a_i(mlF1,mlF1,nF,-nw:nw))
allocate(all_rho_b_r(mlF2,mlF2,nF,-nw:nw))
allocate(all_rho_b_i(mlF2,mlF2,nF,-nw:nw))

call initiate(all_rho_a_r,all_rho_a_i,all_rho_b_r,all_rho_b_i,mlF1,mlF2,nF,nw,Ftrans,w0,dw,&
&mass,w00,alp_Z1,alp_Z2,B,wall_a,wall_b,gam_a,gam_b,Iw,Uw,Qw,Vw,Tb,vthermal,lam1,del,&
&lbd_a_r,lbd_a_i,lbd_b_r,lbd_b_i,&
&anis,apump,p_an_u,p_an_d,&
&inci,qin,qvel,qga,qi)

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

  IF (ALLOCATED (Ipabv_r)) DEALLOCATE (Ipabv_r)
  IF (ALLOCATED (Ipabv_i)) DEALLOCATE (Ipabv_i)
  IF (ALLOCATED (Upabv_r)) DEALLOCATE (Upabv_r)
  IF (ALLOCATED (Upabv_i)) DEALLOCATE (Upabv_i)
  IF (ALLOCATED (Qpabv_r)) DEALLOCATE (Qpabv_r)
  IF (ALLOCATED (Qpabv_i)) DEALLOCATE (Qpabv_i)
  IF (ALLOCATED (Vpabv_r)) DEALLOCATE (Vpabv_r)
  IF (ALLOCATED (Vpabv_i)) DEALLOCATE (Vpabv_i)

  allocate(Ipabv_r(lF1,lF2,-nw:nw))
  allocate(Ipabv_i(lF1,lF2,-nw:nw))
  allocate(Upabv_r(lF1,lF2,-nw:nw))
  allocate(Upabv_i(lF1,lF2,-nw:nw))
  allocate(Qpabv_r(lF1,lF2,-nw:nw))
  allocate(Qpabv_i(lF1,lF2,-nw:nw))
  allocate(Vpabv_r(lF1,lF2,-nw:nw))
  allocate(Vpabv_i(lF1,lF2,-nw:nw))


  do a1=1,lF1
    do b1=1,lF2

      call int_stokes(Iw,Qw,Uw,Vw,lF1,lF2,a1,b1,nw,w0,dw,gam_a,wa,wb,Ipv_r,Ipv_i,Qpv_r,Qpv_i,Upv_r,Upv_i,Vpv_r,Vpv_i)

      do vj=-nw,nw
        Ipabv_r(a1,b1,vj) = Ipv_r(vj) 
        Ipabv_i(a1,b1,vj) = Ipv_i(vj)

        Qpabv_r(a1,b1,vj) = Qpv_r(vj) 
        Qpabv_i(a1,b1,vj) = Qpv_i(vj) 

        Upabv_r(a1,b1,vj) = Upv_r(vj) 
        Upabv_i(a1,b1,vj) = Upv_i(vj) 

        Vpabv_r(a1,b1,vj) = Vpv_r(vj) 
        Vpabv_i(a1,b1,vj) = Vpv_i(vj) 
      end do
    end do
  end do

  IF (ALLOCATED (Ip_r)) DEALLOCATE (Ip_r)
  IF (ALLOCATED (Ip_i)) DEALLOCATE (Ip_i)
  IF (ALLOCATED (Up_r)) DEALLOCATE (Up_r)
  IF (ALLOCATED (Up_i)) DEALLOCATE (Up_i)
  IF (ALLOCATED (Qp_r)) DEALLOCATE (Qp_r)
  IF (ALLOCATED (Qp_i)) DEALLOCATE (Qp_i)
  IF (ALLOCATED (Vp_r)) DEALLOCATE (Vp_r)
  IF (ALLOCATED (Vp_i)) DEALLOCATE (Vp_i)

  allocate(Ip_r(lF1,lF2))
  allocate(Ip_i(lF1,lF2))
  allocate(Up_r(lF1,lF2))
  allocate(Up_i(lF1,lF2))
  allocate(Qp_r(lF1,lF2))
  allocate(Qp_i(lF1,lF2))
  allocate(Vp_r(lF1,lF2))
  allocate(Vp_i(lF1,lF2))


  IF (ALLOCATED (rho_a_r)) DEALLOCATE (rho_a_r)
  IF (ALLOCATED (rho_a_i)) DEALLOCATE (rho_a_i)
  IF (ALLOCATED (rho_b_r)) DEALLOCATE (rho_b_r)
  IF (ALLOCATED (rho_b_i)) DEALLOCATE (rho_b_i)

  allocate(rho_a_r(lF1,lF1))
  allocate(rho_a_i(lF1,lF1))
  allocate(rho_b_r(lF2,lF2))
  allocate(rho_b_i(lF2,lF2))

  IF (ALLOCATED (lbda_r)) DEALLOCATE (lbda_r)
  IF (ALLOCATED (lbda_i)) DEALLOCATE (lbda_i)
  IF (ALLOCATED (lbdb_r)) DEALLOCATE (lbdb_r)
  IF (ALLOCATED (lbdb_i)) DEALLOCATE (lbdb_i)

  allocate(lbda_r(lF1,lF1))
  allocate(lbda_i(lF1,lF1))
  allocate(lbdb_r(lF2,lF2))
  allocate(lbdb_i(lF2,lF2))

  !$OMP DO
  do vj = -nw,nw

    do a1=1,lF1
      do b1=1,lF2
        Ip_r(a1,b1) = Ipabv_r(a1,b1,vj) 
        Ip_i(a1,b1) = Ipabv_i(a1,b1,vj) 
        Qp_r(a1,b1) = Qpabv_r(a1,b1,vj) 
        Qp_i(a1,b1) = Qpabv_i(a1,b1,vj) 
        Up_r(a1,b1) = Upabv_r(a1,b1,vj) 
        Up_i(a1,b1) = Upabv_i(a1,b1,vj) 
        Vp_r(a1,b1) = Vpabv_r(a1,b1,vj) 
        Vp_i(a1,b1) = Vpabv_i(a1,b1,vj) 
      end do
    end do

    do a1 = 1,lF1
      do a2 = a1,lF1
         lbda_r(a1,a2) = lbd_a_r(vj,a1,a2,I) 
         lbda_i(a1,a2) = lbd_a_i(vj,a1,a2,I)
      end do
    end do
    
    do b1 = 1,lF2
      do b2 = b1,lF2
        lbdb_r(b1,b2) = lbd_b_r(vj,b1,b2,I) 
        lbdb_i(b1,b2) = lbd_b_i(vj,b1,b2,I)
      end do
    end do
 
    nF1 = lF1*lF1
    nF2 = lF2*lF2

    call stateq(lF1,lF2,nF1,nF2,lbda_r,lbda_i,lbdb_r,lbdb_i,gam_a(I),gam_b(I),wa,wb,sa,Ip_r,Ip_i,&
&Qp_r,Qp_i,Up_r,Up_i,Vp_r,Vp_i,delta_I,delta_Q,delta_U,delta_V,rho_a_r,rho_a_i,rho_b_r,rho_b_i) 

    do a1 = 1,lF1
      do a2 = a1,lF1
        all_rho_a_r(a1,a2,I,vj) = rho_a_r(a1,a2)
        all_rho_a_i(a1,a2,I,vj) = rho_a_i(a1,a2)

        all_rho_a_r(a2,a1,I,vj) = rho_a_r(a1,a2)
        all_rho_a_i(a2,a1,I,vj) = -rho_a_i(a1,a2)

        if (a1.eq.a2) then
          all_rho_a_i(a1,a2,I,vj) = 0.d0
        endif
      end do
    end do

    do b1 = 1,lF2
      do b2 = b1,lF2
        all_rho_b_r(b1,b2,I,vj) = rho_b_r(b1,b2)
        all_rho_b_i(b1,b2,I,vj) = rho_b_i(b1,b2)

        all_rho_b_r(b2,b1,I,vj) = rho_b_r(b1,b2)
        all_rho_b_i(b2,b1,I,vj) = -rho_b_i(b1,b2)

        if (b1.eq.b2) then
          all_rho_b_i(b1,b2,I,vj) = 0.d0
        endif


      end do
    end do

  end do
  !$OMP END DO

end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!! TRANSFER ELEMENTS !!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! transfer elements

do vj = -nw,nw
  Awtot(vj) = 0.d0 
  Bwtot(vj) = 0.d0 
  Cwtot(vj) = 0.d0 
  Dwtot(vj) = 0.d0 
  Ewtot(vj) = 0.d0 
  Fwtot(vj) = 0.d0 
  Gwtot(vj) = 0.d0 
end do

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


  allocate(rho_aa_r(lF1,lF1,-nw:nw))
  allocate(rho_aa_i(lF1,lF1,-nw:nw))
  allocate(rho_bb_r(lF2,lF2,-nw:nw))
  allocate(rho_bb_i(lF2,lF2,-nw:nw))

  do vj = -nw,nw
    do a1 = 1,lF1
      do a2 = 1,lF1
        rho_aa_r(a1,a2,vj) = all_rho_a_r(a1,a2,I,vj)
        rho_aa_i(a1,a2,vj) = all_rho_a_i(a1,a2,I,vj)
      end do
    end do
  
    do b1 = 1,lF2
      do b2 = 1,lF2
        rho_bb_r(b1,b2,vj) = all_rho_b_r(b1,b2,I,vj)
        rho_bb_i(b1,b2,vj) = all_rho_b_i(b1,b2,I,vj)
      end do
    end do
  end do

  IF (ALLOCATED (gaa_p_r)) DEALLOCATE (gaa_p_r)
  IF (ALLOCATED (gaa_p_i)) DEALLOCATE (gaa_p_i)
  IF (ALLOCATED (gaa_m_r)) DEALLOCATE (gaa_m_r)
  IF (ALLOCATED (gaa_m_i)) DEALLOCATE (gaa_m_i)
  IF (ALLOCATED (gbb_p_r)) DEALLOCATE (gbb_p_r)
  IF (ALLOCATED (gbb_p_i)) DEALLOCATE (gbb_p_i)
  IF (ALLOCATED (gbb_m_r)) DEALLOCATE (gbb_m_r)
  IF (ALLOCATED (gbb_m_i)) DEALLOCATE (gbb_m_i)

  allocate(gaa_p_r(lF2,3,3,-nw:nw))
  allocate(gaa_p_i(lF2,3,3,-nw:nw))
  allocate(gaa_m_r(lF2,3,3,-nw:nw))
  allocate(gaa_m_i(lF2,3,3,-nw:nw))
  allocate(gbb_p_r(lF1,3,3,-nw:nw))
  allocate(gbb_p_i(lF1,3,3,-nw:nw))
  allocate(gbb_m_r(lF1,3,3,-nw:nw))
  allocate(gbb_m_i(lF1,3,3,-nw:nw))

  do vj=-nw,nw
    do b1=1,lF2
      do K=1,3
        do L=1,3
          gaa_p_r(b1,K,L,vj) = 0.d0
          gaa_p_i(b1,K,L,vj) = 0.d0
          gaa_m_r(b1,K,L,vj) = 0.d0
          gaa_m_i(b1,K,L,vj) = 0.d0
        end do
      end do
    end do
  
    do a1=1,lF1
      do K=1,3
        do L=1,3
          gbb_p_r(a1,K,L,vj) = 0.d0
          gbb_p_i(a1,K,L,vj) = 0.d0
          gbb_m_r(a1,K,L,vj) = 0.d0
          gbb_m_i(a1,K,L,vj) = 0.d0
        end do
      end do
    end do
  end do


  call int_densities(rho_aa_r,rho_aa_i,rho_bb_r,rho_bb_i,lF1,lF2,nw,w0,dv,gam_a,gam_b,wa,wb,&
&gaa_p_r,gaa_p_i,gaa_m_r,gaa_m_i,gbb_p_r,gbb_p_i,gbb_m_r,gbb_m_i)


  call trans_el(gaa_p_r,gaa_p_i,gbb_p_r,gbb_p_i,gaa_m_r,gaa_m_i,gbb_m_r,gbb_m_i,lF1,lF2,nw,&
&delta_I,delta_Q,delta_U,delta_V,Aw,Bw,Cw,Dww,Ew,Fw,Gw,w0)

  do vj = -nw,nw
    Awtot(vj) = Awtot(vj) + Aw(vj) 
    Bwtot(vj) = Bwtot(vj) + Bw(vj) 
    Cwtot(vj) = Cwtot(vj) + Cw(vj) 
    Dwtot(vj) = Dwtot(vj) + Dww(vj) 
    Ewtot(vj) = Ewtot(vj) + Ew(vj) 
    Fwtot(vj) = Fwtot(vj) + Fw(vj) 
    Gwtot(vj) = Gwtot(vj) + Gw(vj) 
  end do


end do


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!! STOKES PROPAGATION !!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


call prop_stokes(Awtot,Bwtot,Cwtot,Dwtot,Ewtot,Fwtot,Gwtot,Iw,Qw,Uw,Vw,ds,nw)

!write(*,*) 'step:', STEP

! if any of the densities is not a number, we might as well stop.

if (isnan(all_rho_a_r(1,1,1,1)))  stop   
if (isnan(all_rho_a_i(1,1,1,1)))  stop

if (isnan(all_rho_b_r(1,1,1,1)))  stop
if (isnan(all_rho_b_i(1,1,1,1)))  stop

! some further analysis and writing of the Stokes parameters


call vitals(Iw,Qw,Uw,Vw,nw,pQ,pV,pa,nimax)

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

  call stimulate(Iw(nimax),Qw(nimax),Uw(nimax),Vw(nimax),lF1,lF2,Aij,w0,delta_I,delta_Q,delta_U,delta_V,Rstim,Tbr)
!  call stimulate_alt(Iw(0),Aij,sa,w0,Rstim)

  Rtot = Rtot + Rstim
  Ttot = Ttot + Tbr
end do


call fwhm(Iw,nw,dv,fwhm_i)
!compute the temperature

!dR=dlog10(Rtot)
dR = dlog10(Rtot/(alp_Z1(1) * B))


!if (mod(STEP,10).eq.0.and.Ttot.lt.1E13) then 
if (mod(STEP,10).eq.0) then 
  write(28,'(E16.8, E16.8, E16.8, E16.8, E16.8, E16.8)')Ttot,dR,fwhm_i,pQ,pV,pa
!  write(28,'(E16.8, E16.8, E16.8, E16.8)')dR,pQ,pV,pa

!  write(61,*)Iw  
!  write(62,*)Uw
!  write(63,*)Qw
!  write(64,*)Vw
  
endif

if (dR.gt.0.d0) then

  if (mod(STEP,10).eq.0) then
    write(28,'(E16.8, E16.8, E16.8, E16.8, E16.8, E16.8)')Ttot,dR,fwhm_i,pQ,pV,pa

    if (inp.eq.1) then
      write(61,*)Iw(-nw:nw)    
      write(62,*)Uw(-nw:nw)
      write(63,*)Qw(-nw:nw)
      write(64,*)Vw(-nw:nw)
    endif

  endif
else
  write(28,'(E16.8, E16.8, E16.8, E16.8, E16.8, E16.8)')Ttot,dR,fwhm_i,pQ,pV,pa
  if (inp.eq.1) then
    write(61,*)Iw(-nw:nw)  
    write(62,*)Uw(-nw:nw)
    write(63,*)Qw(-nw:nw)
    write(64,*)Vw(-nw:nw)
  endif
endif

ds = ds * 1.0005d0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!            and... go on             !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

if (dR.gt.3.1d0) stop
if (Ttot.gt.maxTb) stop


if (STEP.le.MAXSTEP) GOTO 111    


close(28)


end program
 
