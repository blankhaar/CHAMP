subroutine prop_stokes(A,B,C,I,Q,V,ds,nv)
implicit none
integer n,nv
double precision Inew,Qnew,Unew,Vnew,ds,x0                
double precision, dimension (4,4) :: eK,eeK 
double precision, dimension (-nv:nv) :: A,B,C
double precision, dimension (-nv:nv) :: I,Q,V 

!input
!A,B,C,D,E,F,G: the propagation coefficients
!I,Q,U,V: the old stokes parameters

!output
!I,Q,U,V: the new stokes parameters

x0 = 0.d0

do n=-nv,nv
  Inew = A(n)*I(n) + B(n)*Q(n) + C(n)*V(n)  
  Qnew = B(n)*I(n) + A(n)*Q(n) 
  Vnew = C(n)*I(n) + A(n)*V(n)
  
  I(n) = I(n) + Inew*ds 
  Q(n) = Q(n) + Qnew*ds
  V(n) = V(n) + Vnew*ds

!  eK(1,1) = A(n)*ds 
!  eK(2,2) = A(n)*ds 
!  eK(3,3) = A(n)*ds 
!  eK(4,4) = A(n)*ds 
!
!  eK(1,2) = B(n)*ds 
!  eK(2,1) = B(n)*ds 
!
!  eK(1,3) = x0 
!  eK(3,1) = x0 
!
!  eK(1,4) = C(n)*ds 
!  eK(4,1) = C(n)*ds 
!
!  eK(2,3) = x0
!  eK(3,2) = x0
!
!  eK(2,4) = x0 
!  eK(4,2) = x0 
!
!  eK(3,4) = x0 
!  eK(4,3) = x0 
!
!  call exp_A(eK,8,eeK)
!
!  if (n.eq.0) then
!    write(*,*)eeK(1,:)
!    write(*,*)eeK(2,:)
!    write(*,*)eeK(3,:)
!    write(*,*)eeK(4,:)
!  
!    write(*,*)eeK(1,1),1.d0+A(n)*ds
!    write(*,*)eeK(1,2),B(n)*ds
!    write(*,*)eeK(1,4),C(n)*ds
!  
!    stop 
!  endif
!  I(n) = eeK(1,1)*I(n) + eeK(1,2)*Q(n) + eeK(1,4)*V(n)
!  Q(n) = eeK(2,1)*I(n) + eeK(2,2)*Q(n) + eeK(2,4)*V(n)
!  V(n) = eeK(4,1)*I(n) + eeK(4,2)*Q(n) + eeK(4,4)*V(n)

end do

end subroutine


subroutine trans_el(gab,lF1,lF2,nw,delta_I,delta_Q,delta_V,Aw,Bw,Cw,w0)
implicit none
integer n,m,nw,a1,a2,b1,b2,ac1,ac2,bc1,bc2
integer F1,F2,lF1,lF2,ma1,ma2,mb1,mb2
double precision Const,pi,c0,hbar,w0,x0 

double precision, dimension (lF1,3,-nw:nw) :: gab
double precision, dimension (lF1,3) :: delta_I,delta_Q,delta_V
double precision, dimension (-nw:nw) :: Aw,Bw,Cw

pi = dacos(-1.d0)
c0 = 299792458.d0
hbar = 1.0545718E-34
!Const = 2.d0*pi*w0/(c0*hbar) 
Const = 2.d0*pi*pi/hbar 

F1 = (lF1 - 1)/2
F2 = (lF2 - 1)/2

x0 = 0.d0

do n=-nw,nw

  Aw(n) = x0
  Bw(n) = x0
  Cw(n) = x0

  do a1 = 1,lF1
    do bc1 = 1,3
      ma1 = a1 - F1 - 1
      mb1 = ma1 + bc1 - 2
      b1 = mb1 + F2 + 1
      
      if (b1.le.lF2.and.b1.gt.0) then
        Aw(n) = Aw(n) + gab(a1,bc1,n)*delta_I(a1,bc1)
        Bw(n) = Bw(n) + gab(a1,bc1,n)*delta_Q(a1,bc1)
        Cw(n) = Cw(n) + gab(a1,bc1,n)*delta_V(a1,bc1)
!        if (n.eq.0) then
!          write(*,*)Ma1,Mb1,gab(a1,bc1,n)
!          write(*,*)Ma1,Mb1,delta_I(a1,bc1),gab(a1,bc1,n)*delta_I(a1,bc1),Aw(n)
!          write(*,*)Ma1,Mb1,delta_Q(a1,bc1),gab(a1,bc1,n)*delta_Q(a1,bc1),Bw(n)
!          write(*,*)Ma1,Mb1,delta_V(a1,bc1),gab(a1,bc1,n)*delta_V(a1,bc1),Cw(n)
!        endif
      endif

    end do
  end do

  
!  stop

  Aw(n) = Const * Aw(n)
  Bw(n) = -Const * Bw(n)
  Cw(n) = Const * Cw(n)

end do      

end subroutine


subroutine int_densities(rho_aa_r,rho_bb_r,lF1,lF2,nv,w0,dv,wa,wb,gab)
implicit none
integer n,m,nv,a1,a2,b1,b2,ac1,ac2,bc1,bc2
integer F1,F2,lF1,lF2,ma1,ma2,mb1,mb2,mu,md
double precision dv,dvt,w0,c0,dw
double precision dr_aa_r,dr_bb_r,dr_aa_i,dr_bb_i 
double precision rs_aa_r,rs_bb_r,rs_aa_i,rs_bb_i 
double precision nat,nbt,ddE,pi 

double precision, dimension (lF1,3,-nv:nv) :: gab
double precision, dimension (lF1,-nv:nv) :: rho_aa_r
double precision, dimension (lF2,-nv:nv) :: rho_bb_r
double precision, dimension (lF1) :: wa 
double precision, dimension (lF2) :: wb

! integrate the densities
! product is the gaa = <\rho_ \gam>-matrices for all the omega values
! 

dvt = 2.d0*dv 
c0 = 299792458.d0 
pi = dacos(-1.d0)

F1 = (lF1 - 1)/2
F2 = (lF2 - 1)/2

do a1 = 1,lF1
  do bc1 = 1,3
    ma1 = a1 - F1 - 1 
    mb1 = ma1 + bc1 - 2 
    b1 = mb1 + F2 + 1 
 
    if (b1.le.lF2.and.b1.gt.0) then
     
      ddE = wa(a1) - wb(b1) 
 
      do n=-nv,nv

        if (n.eq.-nv) then    
          dr_aa_r = (rho_aa_r(a1,-nv+1) - rho_aa_r(a1,-nv))/dv
          dr_bb_r = (rho_bb_r(b1,-nv+1) - rho_bb_r(b1,-nv))/dv
        elseif (n.eq.nv) then
          dr_aa_r = (rho_aa_r(a1,nv) - rho_aa_r(a1,nv-1))/dv
          dr_bb_r = (rho_bb_r(b1,nv) - rho_bb_r(b1,nv-1))/dv
        else
          dr_aa_r = (rho_aa_r(a1,n+1) - rho_aa_r(a1,n-1))/dvt
          dr_bb_r = (rho_bb_r(b1,n+1) - rho_bb_r(b1,n-1))/dvt
        endif
  
!        nat = ( pi * c0 / w0 ) * ( rho_aa_r(a1,n) - dr_aa_r * (c0 / w0) * ddE )
!        nbt = ( pi * c0 / w0 ) * ( rho_bb_r(b1,n) - dr_bb_r * (c0 / w0) * ddE ) 
 
        nat = ( rho_aa_r(a1,n) - dr_aa_r * (c0 / w0) * ddE )
        nbt = ( rho_bb_r(b1,n) - dr_bb_r * (c0 / w0) * ddE ) 

        gab(a1,bc1,n) = nat - nbt 

      end do
    endif  
  end do
end do


end subroutine

