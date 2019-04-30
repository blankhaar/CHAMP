subroutine initiate(mlF1,mlF2,nF,nw,Ftrans,w0,dw,&
&mass,w00,alp_Z1,alp_Z2,B,wall_a,wall_b,Iw,Qw,Vw,Tb,vthermal,lam1,del,lbd_a,lbd_b)

implicit none
integer a1,a2,b1,b2,ma1,mb1,F1,F2,lF1,lF2,I,n,nw,mlF1,mlF2,nF
double precision wab,w0,mass,dw,x0,x1,x2,x3,x4,x8,pi,w,c0,Tb,vthermal,B,h
double precision temp1,temp2,temp3,temp4,conv,lam2,k_sb,Tmp,vel
integer, dimension (nF,2) :: Ftrans
double precision, dimension (-nw:nw) :: Iw,Qw,Vw
double precision, dimension (mlF1,nF) :: wall_a
double precision, dimension (mlF2,nF) :: wall_b
double precision, dimension (-nw:nw,nF) :: lbd_a
double precision, dimension (-nw:nw,nF) :: lbd_b
double precision, dimension (nF) :: alp_Z1,alp_Z2,w00,Aij_all,del,lam1

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
!write(*,*)nF

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


!initiate density matrix
do I =1,nF
  F1 = Ftrans(I,1)
  F2 = Ftrans(I,2)

  lF1 = 2*F1 + 1
  lF2 = 2*F2 + 1

  wab = w00(I)
!  write(*,*)I,wab*conv

!  do a1 = 1,lF1
!    ma1 = a1 - F1 - 1
!    wab = w00(I) + B * alp_Z1(I) * ma1

    lam2 = lam1(I) * (x1 - del(I)/x2) / (x1 + del(I) / x2)
!    write(*,*)I,lam1(I),lam2
    do n=-nw,nw      
         
      vel = dw*n*conv
    
      vel=vel+wab*conv
   
      temp1=dsqrt(x2*k_sb*Tmp/mass)
      temp2=x1/(dsqrt(pi)*temp1)
      temp3=-x1*(vel/temp1)**x2
      temp4=temp2*dexp(temp3)

      lbd_a(n,I) = lam1(I)*temp4 
      lbd_b(n,I) = lam2*temp4      !at first didn't work, but makes no diff 

    end do
!  end do
end do

!stop

!  do b1 = 1,lF2
!    mb1 = b1 - F2 - 1
!    wab = w00(I) + B * alp_Z2(I) * mb1

!    do n=-nw,nw
!      vel = dw*n*conv
!
!      vel=vel+wab*conv
!
!      temp1=dsqrt(x2*k_sb*Tmp/mass)
!      temp2=x1/(dsqrt(pi)*temp1)
!      temp3=-x1*(vel/temp1)**x2
!      temp4=temp2*dexp(temp3)
! 
!      lbd_b(n,I) = Lam2*temp4
!    end do
!  end do
!end do

!initiate incoming radiation

do n = -nw,nw
  w = w0 + dw*n
!  Iw(n) = h*w**x3/(x4*pi**x3*c0*c0*(dexp(h*w/(x2*pi*k_sb*Tb)) - x1))
  Iw(n) = Tb * (k_sb*w*w)/(x4*(pi**x3)*c0*c0)
  Qw(n) = 0.d0 
  Vw(n) = 0.d0 

end do

end subroutine


