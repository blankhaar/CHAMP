subroutine initiate(all_rho_a_r,all_rho_a_i,all_rho_b_r,all_rho_b_i,mlF1,mlF2,nF,nw,Ftrans,w0,dw,&
&mass,w00,alp_Z1,alp_Z2,B,wall_a,wall_b,gam_a,gam_b,Iw,Uw,Qw,Vw,Tb,vthermal,lam1,del,&
&lbd_a_r,lbd_a_i,lbd_b_r,lbd_b_i,&
&anis,apump,p_an_u,p_an_d,&
&inci,qin,qvel,qga,qi)
implicit none
integer a1,a2,b1,b2,ma1,mb1,F1,F2,lF1,lF2,I,n,nw,mlF1,mlF2,nF
integer anis,inci
double precision wab,w0,mass,dw,x0,x1,x2,x3,x4,x8,pi,w,c0,Tb,vthermal,B,h
double precision temp1,temp2,temp3,temp4,conv,k_sb,Tmp,vel,alp,bet,gam
double precision lbd_t_r,lbd_t_i,qin,qvel,qga,qi
double precision fac1,fac2,fac3,lam2
integer, dimension (nF,2) :: Ftrans
double precision, dimension (mlF1,mlF1,nF,-nw:nw) :: all_rho_a_r,all_rho_a_i
double precision, dimension (mlF2,mlF2,nF,-nw:nw) :: all_rho_b_r,all_rho_b_i
double precision, dimension (-nw:nw) :: Iw,Uw,Qw,Vw
double precision, dimension (mlF1,nF) :: wall_a
double precision, dimension (mlF2,nF) :: wall_b
double precision, dimension (-nw:nw,mlF1,mlF1,nF) :: lbd_a_r,lbd_a_i
double precision, dimension (-nw:nw,mlF2,mlF2,nF) :: lbd_b_r,lbd_b_i    
double precision, dimension (nF) :: alp_Z1,alp_Z2,w00,Aij_all
double precision, dimension (nF) :: lam1,del,gam_a,gam_b,p_an_u,p_an_d
double precision, dimension (3) :: apump 

x0 = 0.d0
x1 = 1.d0
x2 = 2.d0
x3 = 3.d0
x4 = 4.d0
x8 = 8.d0

c0 = 299792458.d0
k_sb = 1.38064852E-23 
pi = dacos(-x1)
h = 6.626070040E-34

Tmp=(mass/k_sb)*(0.18d0)*1.e6*(vthermal)**2.

conv = c0/w0


do I=1,nF
  F1 = Ftrans(I,1)
  F2 = Ftrans(I,2)

  lF1 = 2*F1 + 1
  lF2 = 2*F2 + 1

  do a1=1,lF1
    ma1 = a1 - F1 - 1
    wall_a(a1,I) = w00(I) + B*alp_Z1(I)*ma1
  end do

  do b1=1,lF2
    mb1 = b1 - F2 - 1
    wall_b(b1,I) = w00(I) + B*alp_Z2(I)*mb1
  end do

end do

if (anis.eq.1) then
!  write(*,*)'anisotropic pumping'
  call twovec_euler(apump,alp,bet,gam)
!else
!  write(*,*)'no anisotropic pumping'
endif


!initiate density matrix
do I =1,nF
  F1 = Ftrans(I,1)
  F2 = Ftrans(I,2)

  lF1 = 2*F1 + 1
  lF2 = 2*F2 + 1

  do n = -nw,nw
!    vel = dw*n*conv
!
!    vel=vel+w00(I)*conv
!
!    temp1=dsqrt(x2*k_sb*Tmp/mass)
!    temp2=x1/(dsqrt(pi)*temp1)
!    temp3=-x1*(vel/temp1)**x2
!    temp4=temp2*dexp(temp3)


    do a1 = 1,lF1
      do a2 = 1,lF1
        ma1 = a1 - F1 - 1
        wab = w00(I) !+ B * alp_Z1(I) * ma1

        vel = dw*n*conv
    
        vel=vel+wab*conv
    
        temp1=dsqrt(x2*k_sb*Tmp/mass)
        temp2=x1/(dsqrt(pi)*temp1)
        temp3=-x1*(vel/temp1)**x2
        temp4=temp2*dexp(temp3)


        if (anis.eq.1) then
          call anis_pump_nw(lam1(I)*temp4,p_an_u(I),alp,bet,F1,lF1,a1,a2,lbd_t_r,lbd_t_i)
 
          lbd_a_r(n,a1,a2,I) = lbd_t_r
          lbd_a_i(n,a1,a2,I) = lbd_t_i
          all_rho_a_r(a1,a2,I,n) = lbd_a_r(n,a1,a2,I)/gam_a(I)  
          all_rho_a_i(a1,a2,I,n) = lbd_a_i(n,a1,a2,I)/gam_a(I)
        else
          if (a1.eq.a2) then
            lbd_a_r(n,a1,a2,I) = lam1(I)*temp4 
            lbd_a_i(n,a1,a2,I) = x0 

            all_rho_a_r(a1,a2,I,n) = lbd_a_r(n,a1,a2,I)/gam_a(I)
            all_rho_a_i(a1,a2,I,n) = x0
          else

            lbd_a_r(n,a1,a2,I) = x0  
            lbd_a_i(n,a1,a2,I) = x0 

            all_rho_a_r(a1,a2,I,n) = x0 
            all_rho_a_i(a1,a2,I,n) = x0
          endif
        endif
      end do
    end do
  end do

  do n=-nw,nw

!    vel = dw*n*conv
!
!    vel=vel+w00(I)*conv
!
!    temp1=dsqrt(x2*k_sb*Tmp/mass)
!    temp2=x1/(dsqrt(pi)*temp1)
!    temp3=-x1*(vel/temp1)**x2
!    temp4=temp2*dexp(temp3)

    lam2 = lam1(I) * (x1 - del(I)/x2) / (x1 + del(I) / x2) 

    do b1 = 1,lF2
      do b2 = 1,lF2
        mb1 = b1 - F2 - 1
        wab = w00(I) !+ B * alp_Z2(I) * ma1

        vel = dw*n*conv
    
        vel=vel+wab*conv
    
        temp1=dsqrt(x2*k_sb*Tmp/mass)
        temp2=x1/(dsqrt(pi)*temp1)
        temp3=-x1*(vel/temp1)**x2
        temp4=temp2*dexp(temp3)

        if (anis.eq.1) then
          call anis_pump_nw(lam2*temp4,p_an_d(I),alp,bet,F2,lF2,b1,b2,lbd_t_r,lbd_t_i)

          lbd_b_r(n,b1,b2,I) = lbd_t_r
          lbd_b_i(n,b1,b2,I) = lbd_t_i
          all_rho_b_r(b1,b2,I,n) = lbd_b_r(n,b1,b2,I)/gam_b(I)
          all_rho_b_i(b1,b2,I,n) = lbd_b_i(n,b1,b2,I)/gam_b(I)
        else
 
          if (b1.eq.b2) then
            lbd_b_r(n,b1,b2,I) = lam2*temp4
            lbd_b_i(n,b1,b2,I) = x0 

            all_rho_b_r(b1,b2,I,n) = lbd_b_r(n,b1,b2,I)/gam_b(I) 
            all_rho_b_i(b1,b2,I,n) = x0 
          else
            lbd_b_r(n,b1,b2,I) = x0 
            lbd_b_i(n,b1,b2,I) = x0 

            all_rho_b_r(b1,b2,I,n) = x0 
            all_rho_b_i(b1,b2,I,n) = x0
          endif
        endif
      end do
    end do
  end do
end do

!initiate incoming radiation

if (inci.eq.1) then
!  write(*,*) 'inci'
!  fac1 = (qga * dsqrt(x2*pi))**(-x1)
!
!  do n = -nw,nw
!    w = dw*n - (qvel / conv)
!  
!    fac2 = -((w/qga)**x2)/x2
!    fac3 = qi * fac1 * dexp(fac2)
!  
!    Iw(n) = (x1 ) * fac3
!    Uw(n) = qin * fac3
!    Vw(n) = 0.d0
!    Qw(n) = 0.d0
!  end do
  do n = -nw,nw
    w = w0 + dw*n
!    Iw(n) = h*w**x3/(x4*pi**x3*c0*c0*(dexp(h*w/(x2*pi*k_sb*Tb)) - x1))
    Iw(n) = Tb * (k_sb*w*w)/(x4*(pi**x3)*c0*c0)
    Qw(n) = 0.d0
    Vw(n) = 0.d0 
    Uw(n) = qin * Iw(n) 
  end do

else
!  write(*,*) 'no inci'
  do n = -nw,nw
    w = w0 + dw*n
!    Iw(n) = h*w**x3/(x4*pi**x3*c0*c0*(dexp(h*w/(x2*pi*k_sb*Tb)) - x1))
    Iw(n) = Tb * (k_sb*w*w)/(x4*(pi**x3)*c0*c0)
    Qw(n) = 0.d0 
    Vw(n) = 0.d0 
    Uw(n) = 0.d0 
  end do
endif

end subroutine


