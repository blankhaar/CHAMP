subroutine stateq(lF1,lF2,nF1,nF2,lbd_a_r,lbd_a_i,lbd_b_r,lbd_b_i,gam_a,gam_b,wa,wb,sa,Ip_r,Ip_i,&
&Qp_r,Qp_i,Up_r,Up_i,Vp_r,Vp_i,delta_I,delta_Q,delta_U,delta_V,rho_a_r,rho_a_i,rho_b_r,rho_b_i)
implicit none
integer a1,a2,a3,a4,b1,b2,b3,b4,F1,F2,lF1,lF2,I,J,N,nF1,nF2,INFO
integer, dimension (nF1+nF2) :: IPIV
double precision gam_a,gam_b,Const,sa,om
double precision pi,c0,hbar,h,x0,x1,x2
double precision ratfaq_r,ratfaq_i 
double precision ratfaq1_r,ratfaq1_i,ratfaq2_r,ratfaq2_i,ratfaq3_r,ratfaq3_i 
double precision rat_4_r,rat_4_i
double precision teller_r,teller_i,noem_r,noem_i,noemer
double precision, dimension (lF1) :: wa 
double precision, dimension (lF2) :: wb 
double precision, dimension (lF1,lF1) :: rho_a_r,rho_a_i 
double precision, dimension (lF2,lF2) :: rho_b_r,rho_b_i 
double precision, dimension (lF1,lF1) :: lbd_a_r,lbd_a_i 
double precision, dimension (lF2,lF2) :: lbd_b_r,lbd_b_i 
double precision, dimension (lF1,lF2) :: Ip_r,Ip_i,Qp_r,Qp_i,Up_r,Up_i,Vp_r,Vp_i
double precision, dimension (nF1+nF2,nF1+nF2) :: Am,Bm
double precision, dimension (nF1+nF2) :: rhovec,dvec
double precision, dimension (nF1+nF2) :: WORK 
double precision, dimension (3,3) :: delaaI,delaaQ,delaaU,delaaV
double precision, dimension (lF1,3,lF1,3) :: delta_I,delta_Q,delta_U,delta_V

! Rate equations: N & W 1994: (Eq. A.8)

!input
!a1,a2: which states will be selected for new density
!zetp:
!zetm: integrated (gam_p/m \times Zeta) *
!lbd: pumping *
!gam: decay
!om: frequency-difference a1 and a2

!* input for a certain velocity, v
!output
!rho_a_r,rho_b_r


!build in commons: F1,F2,lF1,lF2
!                  c,hbar,pi

c0 = 299792458.d0
h = 6.6261E-34
x0 = 0.d0
x1 = 1.d0
x2 = 2.d0
pi = dacos(-x1)
hbar = h/(x2*pi)

Const = sa*pi/(c0*hbar**x2)

!let's set up dvec
!times minus, because: -d = A*p
I = 0
do a1 = 1,lF1
  I = I + 1
  dvec(I) = -x1*lbd_a_r(a1,a1)
end do

do a1 = 1,lF1
  do a2 = a1+1,lF1
    !real part
    I = I + 1
    dvec(I) = -x1*lbd_a_r(a1,a2)

    !imaginary part
    I = I + 1
    dvec(I) = -x1*lbd_a_i(a1,a2)
  end do
end do


do b1 = 1,lF2
  I = I + 1
  dvec(I) = -x1*lbd_b_r(b1,b1)
end do

do b1 = 1,lF2
  do b2 = b1+1,lF2
    I = I + 1
    dvec(I) = -x1*lbd_b_r(b1,b2)

    I = I + 1
    dvec(I) = -x1*lbd_b_i(b1,b2)
  end do
end do

! now we'll set up the A-matrix 
! matrix with eigenstates
! [diag(pa),real/im(offdidag(pa)),diag(pb)real/im(offdiag(pb))
! so A will be a completely real matrix
! this way, we ensure the hermicity of the end-solutions 
! and decrease the searching space for the inverse of the matrix A
do I = 1,nF1+nF2
  do J = 1,nF1+nF2
    Am(I,J) = x0
  end do
end do


I = 0
do a1 = 1,lF1
  I = I + 1
  J = 0
  do a3 = 1,lF1
    J = J + 1

    if (a1.eq.a3) then
      call rat2_a(lF1,lF2,a1,a1,a3,Ip_r,Ip_i,Qp_r,Qp_i,Up_r,Up_i,Vp_r,Vp_i,delta_I,delta_Q,delta_U,delta_V,&
  &ratfaq2_r,ratfaq2_i)
  
      call rat3_a(lF1,lF2,a1,a1,a3,Ip_r,Ip_i,Qp_r,Qp_i,Up_r,Up_i,Vp_r,Vp_i,delta_I,delta_Q,delta_U,&
  &delta_V,ratfaq3_r,ratfaq3_i)
  !    write(*,*)ratfaq2_r,ratfaq3_r
      Am(I,J) = -Const*(ratfaq2_r + ratfaq3_r)
    endif
  end do

  do a3 = 1,lF1
    do a4 = a3+1,lF1 
      J = J + 1

! real part
      if (a1.eq.a3) then

        call rat2_a(lF1,lF2,a1,a1,a4,Ip_r,Ip_i,Qp_r,Qp_i,Up_r,Up_i,Vp_r,Vp_i,delta_I,delta_Q,delta_U,delta_V,&
&ratfaq2_r,ratfaq2_i)

        call rat3_a(lF1,lF2,a1,a1,a4,Ip_r,Ip_i,Qp_r,Qp_i,Up_r,Up_i,Vp_r,Vp_i,delta_I,delta_Q,delta_U,&
&delta_V,ratfaq3_r,ratfaq3_i)

        Am(I,J) = -Const*(ratfaq2_r + ratfaq3_r)

      endif

      if (a1.eq.a4) then

        call rat2_a(lF1,lF2,a1,a1,a3,Ip_r,Ip_i,Qp_r,Qp_i,Up_r,Up_i,Vp_r,Vp_i,delta_I,delta_Q,delta_U,delta_V,&
&ratfaq2_r,ratfaq2_i)

        call rat3_a(lF1,lF2,a1,a1,a3,Ip_r,Ip_i,Qp_r,Qp_i,Up_r,Up_i,Vp_r,Vp_i,delta_I,delta_Q,delta_U,&
&delta_V,ratfaq3_r,ratfaq3_i)

        Am(I,J) = -Const*(ratfaq2_r + ratfaq3_r)

      endif

! imaginary part
      J = J + 1  
      if (a1.eq.a3) then
        Am(I,J) = -Const*(ratfaq2_i - ratfaq3_i)
      endif  
 
      if (a1.eq.a4) then
        Am(I,J) = -Const*(ratfaq3_i - ratfaq2_i)
      endif  

    end do
  end do

! now the b-part
  do b1 = 1,lF2
    J = J + 1
   
    call rat1_a(lF1,lF2,a1,a1,b1,b1,Ip_r,Ip_i,Qp_r,Qp_i,Up_r,Up_i,Vp_r,Vp_i,delta_I,delta_Q,delta_U,delta_V,ratfaq_r,&
&ratfaq_i) 

    Am(I,J) = Const*ratfaq_r
  end do
  do b1 = 1,lF2
    do b2 = b1+1,lF2 
      J = J + 1
  
      call rat1_a(lF1,lF2,a1,a1,b1,b2,Ip_r,Ip_i,Qp_r,Qp_i,Up_r,Up_i,Vp_r,Vp_i,delta_I,delta_Q,delta_U,delta_V,ratfaq1_r,&
&ratfaq1_i)

      call rat1_a(lF1,lF2,a1,a1,b2,b1,Ip_r,Ip_i,Qp_r,Qp_i,Up_r,Up_i,Vp_r,Vp_i,delta_I,delta_Q,delta_U,delta_V,ratfaq2_r,&
&ratfaq2_i)

      Am(I,J) = Const*(ratfaq1_r + ratfaq2_r)

      J = J + 1
 
      Am(I,J) = -Const*(ratfaq1_i - ratfaq2_i)
    end do
  end do 
end do

do a1 = 1,lF1
  do a2 = a1+1,lF1
    I = I + 1
    J = 0
! bra --- real part 
    do a3 = 1,lF1
 
      J = J + 1   

      if (a1.eq.a3) then
        call rat3_a(lF1,lF2,a1,a2,a1,Ip_r,Ip_i,Qp_r,Qp_i,Up_r,Up_i,Vp_r,Vp_i,delta_I,delta_Q,delta_U,&
&delta_V,ratfaq_r,ratfaq_i)
   
        Am(I,J) = -Const*ratfaq_r 
      endif

      if (a2.eq.a3) then
        call rat2_a(lF1,lF2,a1,a2,a2,Ip_r,Ip_i,Qp_r,Qp_i,Up_r,Up_i,Vp_r,Vp_i,delta_I,delta_Q,delta_U,&
&delta_V,ratfaq_r,ratfaq_i) !switched it to rat2_a
   
        Am(I,J) = -Const*ratfaq_r 
      endif

    end do

    do a3 = 1,lF1
      do a4 = a3+1,lF1
        J = J + 1

!       real part
        if (a3.eq.a1) then
          call rat3_a(lF1,lF2,a1,a2,a4,Ip_r,Ip_i,Qp_r,Qp_i,Up_r,Up_i,Vp_r,Vp_i,delta_I,delta_Q,delta_U,&
&delta_V,ratfaq3_r,ratfaq3_i)       
 
          Am(I,J) = Am(I,J) - Const*ratfaq3_r 
        endif 

        if (a4.eq.a1) then
          call rat3_a(lF1,lF2,a1,a2,a3,Ip_r,Ip_i,Qp_r,Qp_i,Up_r,Up_i,Vp_r,Vp_i,delta_I,delta_Q,delta_U,&
&delta_V,ratfaq3_r,ratfaq3_i)

          Am(I,J) = Am(I,J) - Const*ratfaq3_r
        endif

        if (a2.eq.a4) then
          call rat2_a(lF1,lF2,a1,a2,a3,Ip_r,Ip_i,Qp_r,Qp_i,Up_r,Up_i,Vp_r,Vp_i,delta_I,delta_Q,delta_U,&
&delta_V,ratfaq2_r,ratfaq2_i)       
 
          Am(I,J) = Am(I,J) - Const*ratfaq2_r 
        endif 

        if (a2.eq.a3) then
          call rat2_a(lF1,lF2,a1,a2,a4,Ip_r,Ip_i,Qp_r,Qp_i,Up_r,Up_i,Vp_r,Vp_i,delta_I,delta_Q,delta_U,&
&delta_V,ratfaq2_r,ratfaq2_i)       
 
          Am(I,J) = Am(I,J) - Const*ratfaq2_r 
        endif 

        J = J + 1

!       imaginary part
        if (a3.eq.a1) then
           Am(I,J) = Am(I,J) + Const*ratfaq3_i
        endif 

        if (a4.eq.a1) then
           Am(I,J) = Am(I,J) - Const*ratfaq3_i
        endif 

        if (a3.eq.a2) then
           Am(I,J) = Am(I,J) - Const*ratfaq2_i
        endif 

        if (a4.eq.a2) then
           Am(I,J) = Am(I,J) + Const*ratfaq2_i
        endif 
      end do
    end do
   
    do b1 = 1,lF2
      J = J + 1
      call rat1_a(lF1,lF2,a1,a2,b1,b1,Ip_r,Ip_i,Qp_r,Qp_i,Up_r,Up_i,Vp_r,Vp_i,delta_I,delta_Q,delta_U,delta_V,ratfaq_r,&
  &ratfaq_i)
  
      Am(I,J) = Const*ratfaq_r
    end do

    do b1 = 1,lF2
      do b2 = b1+1,lF2
        J = J + 1

        call rat1_a(lF1,lF2,a1,a2,b1,b2,Ip_r,Ip_i,Qp_r,Qp_i,Up_r,Up_i,Vp_r,Vp_i,delta_I,delta_Q,delta_U,delta_V,ratfaq1_r,&
  &ratfaq1_i)
        call rat1_a(lF1,lF2,a1,a2,b2,b1,Ip_r,Ip_i,Qp_r,Qp_i,Up_r,Up_i,Vp_r,Vp_i,delta_I,delta_Q,delta_U,delta_V,ratfaq2_r,&
  &ratfaq2_i) 

        Am(I,J) = Const*(ratfaq1_r + ratfaq2_r) 
 
        J = J + 1

        Am(I,J) = -Const*(ratfaq1_i - ratfaq2_i) 
      end do
    end do

    I = I + 1
    J = 0
!   bra---imaginary part
    do a3 = 1,lF1
      J = J + 1   

      if (a1.eq.a3) then
        call rat3_a(lF1,lF2,a1,a2,a1,Ip_r,Ip_i,Qp_r,Qp_i,Up_r,Up_i,Vp_r,Vp_i,delta_I,delta_Q,delta_U,&
&delta_V,ratfaq_r,ratfaq_i)
   
        Am(I,J) = -Const*ratfaq_i 
      endif

      if (a2.eq.a3) then
        call rat2_a(lF1,lF2,a1,a2,a2,Ip_r,Ip_i,Qp_r,Qp_i,Up_r,Up_i,Vp_r,Vp_i,delta_I,delta_Q,delta_U,&
&delta_V,ratfaq_r,ratfaq_i) !changed it to rat2_a
   
        Am(I,J) = -Const*ratfaq_i 
      endif

    end do

    do a3 = 1,lF1
      do a4 = a3+1,lF1
        J = J + 1

!       ket --- real part
        if (a3.eq.a1) then
          call rat3_a(lF1,lF2,a1,a2,a4,Ip_r,Ip_i,Qp_r,Qp_i,Up_r,Up_i,Vp_r,Vp_i,delta_I,delta_Q,delta_U,&
&delta_V,ratfaq3_r,ratfaq3_i)       
 
          Am(I,J) = Am(I,J) - Const*ratfaq3_i 
        endif 

        if (a4.eq.a1) then
          call rat3_a(lF1,lF2,a1,a2,a3,Ip_r,Ip_i,Qp_r,Qp_i,Up_r,Up_i,Vp_r,Vp_i,delta_I,delta_Q,delta_U,&
&delta_V,ratfaq3_r,ratfaq3_i)

          Am(I,J) = Am(I,J) - Const*ratfaq3_i
        endif

        if (a2.eq.a4) then
          call rat2_a(lF1,lF2,a1,a2,a3,Ip_r,Ip_i,Qp_r,Qp_i,Up_r,Up_i,Vp_r,Vp_i,delta_I,delta_Q,delta_U,&
&delta_V,ratfaq2_r,ratfaq2_i)       
 
          Am(I,J) = Am(I,J) - Const*ratfaq2_i 
        endif 

        if (a2.eq.a3) then
          call rat2_a(lF1,lF2,a1,a2,a4,Ip_r,Ip_i,Qp_r,Qp_i,Up_r,Up_i,Vp_r,Vp_i,delta_I,delta_Q,delta_U,&
&delta_V,ratfaq2_r,ratfaq2_i)       
 
          Am(I,J) = Am(I,J) - Const*ratfaq2_i 
        endif 

        J = J + 1

!       imaginary part
        if (a3.eq.a1) then
           Am(I,J) = Am(I,J) - Const*ratfaq3_r
        endif 

        if (a4.eq.a1) then
           Am(I,J) = Am(I,J) + Const*ratfaq3_r
        endif 

        if (a3.eq.a2) then
           Am(I,J) = Am(I,J) + Const*ratfaq2_r
        endif 

        if (a4.eq.a2) then
           Am(I,J) = Am(I,J) - Const*ratfaq2_r
        endif 
      end do
    end do
   
    do b1 = 1,lF2
      J = J + 1
      call rat1_a(lF1,lF2,a1,a2,b1,b1,Ip_r,Ip_i,Qp_r,Qp_i,Up_r,Up_i,Vp_r,Vp_i,delta_I,delta_Q,delta_U,delta_V,ratfaq_r,&
  &ratfaq_i)
  
      Am(I,J) = Const*ratfaq_i
    end do

    do b1 = 1,lF2
      do b2 = b1+1,lF2
        J = J + 1

        call rat1_a(lF1,lF2,a1,a2,b1,b2,Ip_r,Ip_i,Qp_r,Qp_i,Up_r,Up_i,Vp_r,Vp_i,delta_I,delta_Q,delta_U,delta_V,ratfaq1_r,&
  &ratfaq1_i)
        call rat1_a(lF1,lF2,a1,a2,b2,b1,Ip_r,Ip_i,Qp_r,Qp_i,Up_r,Up_i,Vp_r,Vp_i,delta_I,delta_Q,delta_U,delta_V,ratfaq2_r,&
  &ratfaq2_i) 

        Am(I,J) = Const*(ratfaq1_i + ratfaq2_i) 
 
        J = J + 1

        Am(I,J) = Const*(ratfaq1_r - ratfaq2_r) 
      end do
    end do
  end do 
end do

! bra -- b-part
do b1 = 1,lF2
  I = I + 1
  J = 0
  do a1 = 1,lF1
    J = J + 1

    call rat1_b(lF1,lF2,a1,a1,b1,b1,Ip_r,Ip_i,Qp_r,Qp_i,Up_r,Up_i,Vp_r,Vp_i,&
&delta_I,delta_Q,delta_U,delta_V,ratfaq_r,ratfaq_i)         

    Am(I,J) = Const*ratfaq_r
  end do

  do a1 = 1,lF1
    do a2 = a1+1,lF1
      J = J + 1
      
      call rat1_b(lF1,lF2,a1,a2,b1,b1,Ip_r,Ip_i,Qp_r,Qp_i,Up_r,Up_i,Vp_r,Vp_i,&
&delta_I,delta_Q,delta_U,delta_V,ratfaq1_r,ratfaq1_i)         
 
      call rat1_b(lF1,lF2,a2,a1,b1,b1,Ip_r,Ip_i,Qp_r,Qp_i,Up_r,Up_i,Vp_r,Vp_i,&
&delta_I,delta_Q,delta_U,delta_V,ratfaq2_r,ratfaq2_i)         

      Am(I,J) = Const*(ratfaq1_r + ratfaq2_r) 

      J = J + 1

      Am(I,J) = -Const*(ratfaq1_i - ratfaq2_i)

    end do
  end do

  do b3 = 1,lF2
    J = J + 1
 
    if (b3.eq.b1) then
      call rat2_b(lF1,lF2,b1,b1,b1,Ip_r,Ip_i,Qp_r,Qp_i,Up_r,Up_i,Vp_r,Vp_i,&
&delta_I,delta_Q,delta_U,delta_V,ratfaq2_r,ratfaq2_i) 

      call rat3_b(lF1,lF2,b1,b1,b1,Ip_r,Ip_i,Qp_r,Qp_i,Up_r,Up_i,Vp_r,Vp_i,&
&delta_I,delta_Q,delta_U,delta_V,ratfaq3_r,ratfaq3_i) 

      Am(I,J) = -Const*(ratfaq2_r + ratfaq3_r)
    endif 
    
  end do

  do b3 = 1,lF2
    do b4 = b3+1,lF2
      J = J + 1 
      if (b3.eq.b1) then
        call rat2_b(lF1,lF2,b1,b1,b4,Ip_r,Ip_i,Qp_r,Qp_i,Up_r,Up_i,Vp_r,Vp_i,&
&delta_I,delta_Q,delta_U,delta_V,ratfaq2_r,ratfaq2_i)

        call rat3_b(lF1,lF2,b1,b1,b4,Ip_r,Ip_i,Qp_r,Qp_i,Up_r,Up_i,Vp_r,Vp_i,&
&delta_I,delta_Q,delta_U,delta_V,ratfaq3_r,ratfaq3_i)
 
        Am(I,J) = Am(I,J) - Const*(ratfaq2_r + ratfaq3_r)
      endif

      if (b4.eq.b1) then
        call rat2_b(lF1,lF2,b1,b1,b3,Ip_r,Ip_i,Qp_r,Qp_i,Up_r,Up_i,Vp_r,Vp_i,&
&delta_I,delta_Q,delta_U,delta_V,ratfaq2_r,ratfaq2_i)

        call rat3_b(lF1,lF2,b1,b1,b3,Ip_r,Ip_i,Qp_r,Qp_i,Up_r,Up_i,Vp_r,Vp_i,&
&delta_I,delta_Q,delta_U,delta_V,ratfaq3_r,ratfaq3_i)
 
        Am(I,J) = Am(I,J) - Const*(ratfaq2_r + ratfaq3_r)
      endif

      J = J + 1      

      if (b3.eq.b1) then 
        Am(I,J) = Am(I,J) + Const*(ratfaq3_i - ratfaq2_i)
      endif

      if (b4.eq.b1) then 
        Am(I,J) = Am(I,J) + Const*(ratfaq2_i - ratfaq3_i)
      endif

    end do
  end do
end do

do b1 = 1,lF2
  do b2 = b1+1,lF2
    I = I + 1 
    J = 0
    do a1 = 1,lF1
      J = J + 1
  
      call rat1_b(lF1,lF2,a1,a1,b1,b2,Ip_r,Ip_i,Qp_r,Qp_i,Up_r,Up_i,Vp_r,Vp_i,&
  &delta_I,delta_Q,delta_U,delta_V,ratfaq_r,ratfaq_i)         
  
      Am(I,J) = Const*ratfaq_r
    end do
  
    do a1 = 1,lF1
      do a2 = a1+1,lF1
        J = J + 1
        
        call rat1_b(lF1,lF2,a1,a2,b1,b2,Ip_r,Ip_i,Qp_r,Qp_i,Up_r,Up_i,Vp_r,Vp_i,&
  &delta_I,delta_Q,delta_U,delta_V,ratfaq1_r,ratfaq1_i)         
   
        call rat1_b(lF1,lF2,a2,a1,b1,b2,Ip_r,Ip_i,Qp_r,Qp_i,Up_r,Up_i,Vp_r,Vp_i,&
  &delta_I,delta_Q,delta_U,delta_V,ratfaq2_r,ratfaq2_i)         
  
        Am(I,J) = Const*(ratfaq1_r + ratfaq2_r) 
  
        J = J + 1
  
        Am(I,J) = -Const*(ratfaq1_i - ratfaq2_i)
  
      end do
    end do
  
    do b3 = 1,lF2
      J = J + 1
   
      if (b3.eq.b2) then
        call rat2_b(lF1,lF2,b1,b2,b2,Ip_r,Ip_i,Qp_r,Qp_i,Up_r,Up_i,Vp_r,Vp_i,&
  &delta_I,delta_Q,delta_U,delta_V,ratfaq2_r,ratfaq2_i) 
  
        Am(I,J) = -Const*ratfaq2_r
      endif

      if (b3.eq.b1) then
        call rat3_b(lF1,lF2,b1,b2,b1,Ip_r,Ip_i,Qp_r,Qp_i,Up_r,Up_i,Vp_r,Vp_i,&
  &delta_I,delta_Q,delta_U,delta_V,ratfaq3_r,ratfaq3_i) 
  
        Am(I,J) = -Const*ratfaq3_r
      endif 
      
    end do
  
    do b3 = 1,lF2
      do b4 = b3+1,lF2
        J = J + 1  
        if (b3.eq.b2) then
          call rat2_b(lF1,lF2,b1,b2,b4,Ip_r,Ip_i,Qp_r,Qp_i,Up_r,Up_i,Vp_r,Vp_i,&
  &delta_I,delta_Q,delta_U,delta_V,ratfaq2_r,ratfaq2_i)
  
          Am(I,J) = Am(I,J) - Const*ratfaq2_r
        endif

        if (b4.eq.b2) then
          call rat2_b(lF1,lF2,b1,b2,b3,Ip_r,Ip_i,Qp_r,Qp_i,Up_r,Up_i,Vp_r,Vp_i,&
  &delta_I,delta_Q,delta_U,delta_V,ratfaq2_r,ratfaq2_i)
  
          Am(I,J) = Am(I,J) - Const*ratfaq2_r
        endif
  
        if (b3.eq.b1) then
          call rat3_b(lF1,lF2,b1,b2,b4,Ip_r,Ip_i,Qp_r,Qp_i,Up_r,Up_i,Vp_r,Vp_i,&
  &delta_I,delta_Q,delta_U,delta_V,ratfaq3_r,ratfaq3_i)
   
          Am(I,J) = Am(I,J) - Const*ratfaq3_r 
        endif
  
        if (b4.eq.b1) then
          call rat3_b(lF1,lF2,b1,b2,b3,Ip_r,Ip_i,Qp_r,Qp_i,Up_r,Up_i,Vp_r,Vp_i,&
  &delta_I,delta_Q,delta_U,delta_V,ratfaq3_r,ratfaq3_i)
   
          Am(I,J) = Am(I,J) - Const*ratfaq3_r 
        endif
 
        J = J + 1
 
        if (b3.eq.b1) then 
          Am(I,J) = Am(I,J) + Const*ratfaq3_i 
        endif
  
        if (b4.eq.b1) then 
          Am(I,J) = Am(I,J) - Const*ratfaq3_i 
        endif
  
        if (b3.eq.b2) then 
          Am(I,J) = Am(I,J) - Const*ratfaq2_i
        endif
  
        if (b4.eq.b2) then 
          Am(I,J) = Am(I,J) + Const*ratfaq2_i 
        endif

      end do
    end do

!   bra --- imaginary
    I = I + 1
    J = 0

    do a1 = 1,lF1
      J = J + 1
  
      call rat1_b(lF1,lF2,a1,a1,b1,b2,Ip_r,Ip_i,Qp_r,Qp_i,Up_r,Up_i,Vp_r,Vp_i,&
  &delta_I,delta_Q,delta_U,delta_V,ratfaq_r,ratfaq_i)         
  
      Am(I,J) = Const*ratfaq_i
    end do
  
    do a1 = 1,lF1
      do a2 = a1+1,lF1
        J = J + 1
        
        call rat1_b(lF1,lF2,a1,a2,b1,b2,Ip_r,Ip_i,Qp_r,Qp_i,Up_r,Up_i,Vp_r,Vp_i,&
  &delta_I,delta_Q,delta_U,delta_V,ratfaq1_r,ratfaq1_i)         
   
        call rat1_b(lF1,lF2,a2,a1,b1,b2,Ip_r,Ip_i,Qp_r,Qp_i,Up_r,Up_i,Vp_r,Vp_i,&
  &delta_I,delta_Q,delta_U,delta_V,ratfaq2_r,ratfaq2_i)         
  
        Am(I,J) = Const*(ratfaq1_i + ratfaq2_i) 
  
        J = J + 1
  
        Am(I,J) = Const*(ratfaq1_r - ratfaq2_r)
  
      end do
    end do
  
    do b3 = 1,lF2
      J = J + 1
   
      if (b3.eq.b2) then
        call rat2_b(lF1,lF2,b1,b2,b2,Ip_r,Ip_i,Qp_r,Qp_i,Up_r,Up_i,Vp_r,Vp_i,&
  &delta_I,delta_Q,delta_U,delta_V,ratfaq2_r,ratfaq2_i) 
  
        Am(I,J) = -Const*ratfaq2_i
      endif

      if (b3.eq.b1) then
        call rat3_b(lF1,lF2,b1,b2,b1,Ip_r,Ip_i,Qp_r,Qp_i,Up_r,Up_i,Vp_r,Vp_i,&
  &delta_I,delta_Q,delta_U,delta_V,ratfaq3_r,ratfaq3_i) 
  
        Am(I,J) = -Const*ratfaq3_i
      endif 
      
    end do
  
    do b3 = 1,lF2
      do b4 = b3+1,lF2
        J = J + 1  
        if (b3.eq.b2) then
          call rat2_b(lF1,lF2,b1,b2,b4,Ip_r,Ip_i,Qp_r,Qp_i,Up_r,Up_i,Vp_r,Vp_i,&
  &delta_I,delta_Q,delta_U,delta_V,ratfaq2_r,ratfaq2_i)
  
          Am(I,J) = Am(I,J) - Const*ratfaq2_i
        endif

        if (b4.eq.b2) then
          call rat2_b(lF1,lF2,b1,b2,b3,Ip_r,Ip_i,Qp_r,Qp_i,Up_r,Up_i,Vp_r,Vp_i,&
  &delta_I,delta_Q,delta_U,delta_V,ratfaq2_r,ratfaq2_i)
  
          Am(I,J) = Am(I,J) - Const*ratfaq2_i
        endif
  
        if (b3.eq.b1) then
          call rat3_b(lF1,lF2,b1,b2,b4,Ip_r,Ip_i,Qp_r,Qp_i,Up_r,Up_i,Vp_r,Vp_i,&
  &delta_I,delta_Q,delta_U,delta_V,ratfaq3_r,ratfaq3_i)
   
          Am(I,J) = Am(I,J) - Const*ratfaq3_i 
        endif
  
        if (b4.eq.b1) then
          call rat3_b(lF1,lF2,b1,b2,b3,Ip_r,Ip_i,Qp_r,Qp_i,Up_r,Up_i,Vp_r,Vp_i,&
  &delta_I,delta_Q,delta_U,delta_V,ratfaq3_r,ratfaq3_i)
   
          Am(I,J) = Am(I,J) - Const*ratfaq3_i 
        endif
 
        J = J + 1
 
        if (b3.eq.b1) then 
          Am(I,J) = Am(I,J) - Const*ratfaq3_r 
        endif
  
        if (b4.eq.b1) then 
          Am(I,J) = Am(I,J) + Const*ratfaq3_r 
        endif
  
        if (b3.eq.b2) then 
          Am(I,J) = Am(I,J) + Const*ratfaq2_r
        endif
  
        if (b4.eq.b2) then 
          Am(I,J) = Am(I,J) - Const*ratfaq2_r 
        endif

      end do
    end do

  end do
end do

!...add the constant, diagonal parts and the omega term, that is not diagonal in this symmetry adapted basis...
I = 0
do a1 = 1,lF1
  I = I + 1
  Am(I,I) = Am(I,I) - gam_a    
end do
do a1 = 1,lF1
  do a2 = a1+1,lF1
    I = I + 1

    om = wa(a1) - wa(a2)

    Am(I,I) =  Am(I,I) - gam_a
    Am(I,I+1) = Am(I,I+1) - om
    I = I + 1
    Am(I,I) = Am(I,I) - gam_a
    Am(I,I-1) = Am(I,I-1) + om

  end do
end do

do b1 = 1,lF2
  I = I + 1
  Am(I,I) = Am(I,I) - gam_b
end do

do b1 = 1,lF2
  do b2 = b1+1,lF2
    I = I + 1
 
    om = wb(b1) - wb(b2)

    Am(I,I) = Am(I,I) - gam_b
    Am(I,I+1) = Am(I,I+1) - om
    
    I = I + 1
  
    Am(I,I) = Am(I,I) - gam_b
    Am(I,I-1) = Am(I,I-1) + om
  end do
end do
 
!write Am
!write(*,*)
!do I = 1,nF1+nF2
!  write(*,*)Am(I,:)
!end do
!write(*,*)
!stop
!now invert Am
N = nF1 + nF2 
call DGETRF(N,N,Am,N,IPIV,INFO)
if (INFO.ne.0) then
  write(*,*)'Density-matrix calculation: failed'
  stop
else
!  call ZGETRI(N,Am,N,IPIV,WORK,-1,INFO)
  call DGETRI(N,Am,N,IPIV,WORK,N,INFO) !this one works better
  if (INFO.ne.0) then
    write(*,*)'Density-matrix calculation: failed'
    stop
  endif
endif 


call matvec_prod(Am,dvec,nF1+nF2,rhovec)

I = 0
do a1 = 1,lF1
  I = I + 1
  rho_a_r(a1,a1) = rhovec(I)
end do

do a1 = 1,lF1
  do a2 = a1+1,lF1
    I = I + 1
    rho_a_r(a1,a2) = rhovec(I)
    rho_a_r(a2,a1) = rhovec(I)

    I = I + 1
    rho_a_i(a1,a2) = rhovec(I)
    rho_a_i(a2,a1) = -rhovec(I)
  end do
end do

do b1 = 1,lF2
  I = I + 1
  rho_b_r(b1,b1) = rhovec(I)
end do

do b1 = 1,lF2
  do b2 = b1+1,lF2
    I = I + 1 
    rho_b_r(b1,b2) = rhovec(I)
    rho_b_r(b2,b1) = rhovec(I)

    I = I + 1
    rho_b_i(b1,b2) = rhovec(I)
    rho_b_i(b2,b1) = -rhovec(I)
  end do
end do

end subroutine

subroutine rat1_a(lF1,lF2,a1,a2,b1,b2,Ip_r,Ip_i,Qp_r,Qp_i,Up_r,Up_i,Vp_r,Vp_i,delta_I,delta_Q,delta_U,delta_V,ratfaq_r,&
&ratfaq_i)
! checked 19-06-2017
! consistent with A.8
! checked numerically
implicit none
integer a1,a2,b1,b2,bc1,bc2,ma1,ma2,lF1,lF2,F1,F2
double precision ratfaq_r,ratfaq_i,zeta_r,zeta_i
double precision Ifp_r,Qfp_r,Ufp_r,Vfp_r
double precision Ifp_i,Qfp_i,Ufp_i,Vfp_i
double precision, dimension (lF2,lF2) :: rho_b_r,rho_b_i
double precision, dimension (lF1,3,lF1,3) :: delta_I,delta_Q,delta_U,delta_V
double precision, dimension (lF1,lF2) :: Ip_r,Ip_i,Qp_r,Qp_i,Up_r,Up_i,Vp_r,Vp_i


F1 = (lF1 - 1)/2
F2 = (lF2 - 1)/2

ma1 = a1 - 1 - F1
ma2 = a2 - 1 - F1

ratfaq_r = 0.d0
ratfaq_i = 0.d0

bc1 = b1 - ma1 + 2 - F2 - 1
bc2 = b2 - ma2 + 2 - F2 - 1

if (bc1.gt.0.and.bc1.le.3.and.bc2.gt.0.and.bc2.le.3) then

  Ifp_r = Ip_r(a1,b2) + Ip_r(a2,b1)
  Ifp_i = Ip_i(a1,b2) - Ip_i(a2,b1) !gamma plus and gamma minus contributions
  Ufp_r = Up_r(a1,b2) + Up_r(a2,b1)
  Ufp_i = Up_i(a1,b2) - Up_i(a2,b1) !gamma plus and gamma minus contributions
  Qfp_r = Qp_r(a1,b2) + Qp_r(a2,b1)
  Qfp_i = Qp_i(a1,b2) - Qp_i(a2,b1) !gamma plus and gamma minus contributions
  Vfp_r = Vp_r(a1,b2) + Vp_r(a2,b1)
  Vfp_i = Vp_i(a1,b2) - Vp_i(a2,b1) !gamma plus and gamma minus contributions

  if (mod(bc1+bc2,2).eq.0) then

    zeta_r = Ifp_r*delta_I(a2,bc2,a1,bc1) - Qfp_r*delta_Q(a2,bc2,a1,bc1) + Vfp_r*delta_V(a2,bc2,a1,bc1)  &
&+ Ufp_i*delta_U(a2,bc2,a1,bc1)
    zeta_i = -Ufp_r*delta_U(a2,bc2,a1,bc1) + Ifp_i*delta_I(a2,bc2,a1,bc1) - Qfp_i*delta_Q(a2,bc2,a1,bc1) &
&+ Vfp_i*delta_V(a2,bc2,a1,bc1)

  else

    zeta_i = Ifp_r*delta_I(a2,bc2,a1,bc1) - Qfp_r*delta_Q(a2,bc2,a1,bc1) + Vfp_r*delta_V(a2,bc2,a1,bc1) &
&+ Ufp_i*delta_U(a2,bc2,a1,bc1) ! only switched zeta_i and zeta_r for these two
    zeta_r = Ufp_r*delta_U(a2,bc2,a1,bc1) - Ifp_i*delta_I(a2,bc2,a1,bc1) + Qfp_i*delta_Q(a2,bc2,a1,bc1) &
&- Vfp_i*delta_V(a2,bc2,a1,bc1) !and this one minus
  
  endif
  
  ratfaq_r = zeta_r 
  ratfaq_i = zeta_i 
  
endif
 
end subroutine


subroutine rat2_a(lF1,lF2,a1,a2,a3,Ip_r,Ip_i,Qp_r,Qp_i,Up_r,Up_i,Vp_r,Vp_i,delta_I,delta_Q,delta_U,delta_V,&
&ratfaq_r,ratfaq_i)
! we take the correct D&W 1990 formula A.39 for these
! checked 19-06-2017
! consistent with A.39
implicit none
integer a1,a2,a3,b1,b2,b,mb,bc1,bc2,ma1,ma2,ma3,lF1,lF2,F1,F2,dm
double precision ratfaq_r,ratfaq_i,zeta_r,zeta_i
double precision, dimension (lF1,lF1) :: rho_a_r,rho_a_i
double precision, dimension (lF1,lF2) :: Ip_r,Ip_i,Qp_r,Qp_i,Up_r,Up_i,Vp_r,Vp_i
double precision, dimension (lF1,3,lF1,3) :: delta_I,delta_Q,delta_U,delta_V


F1 = (lF1 - 1)/2
F2 = (lF2 - 1)/2

ma1 = a1 - F1 - 1
ma2 = a2 - F1 - 1

ratfaq_r = 0.d0 
ratfaq_i = 0.d0

do bc1 =1,3
  mb = ma1 + (bc1 - 2)
  b = mb + F2 + 1
  zeta_r = 0.d0
  zeta_i = 0.d0

      ma3 = a3 - F1 - 1
      dm = mb - ma3
  
      if (abs(dm).le.1) then
        bc2 = dm + 2
 
        if (b.gt.0.and.b.le.lF2) then 
! if bc1+bc2 is even, then delta is real, otherwise it is imaginary
          if (mod(bc1+bc2,2).eq.0) then
            zeta_r = Ip_r(a2,b)*delta_I(a1,bc1,a3,bc2) - Qp_r(a2,b)*delta_Q(a1,bc1,a3,bc2) + &
& Vp_r(a2,b)*delta_V(a1,bc1,a3,bc2) + Up_i(a2,b)*delta_U(a1,bc1,a3,bc2)
      !\gam_p is conjugated---it is a minus
            zeta_i = -Ip_i(a2,b)*delta_I(a1,bc1,a3,bc2) + Qp_i(a2,b)*delta_Q(a1,bc1,a3,bc2) - &
& Vp_i(a2,b)*delta_V(a1,bc1,a3,bc2) + Up_r(a2,b)*delta_U(a1,bc1,a3,bc2)
          else
            zeta_i = -Ip_r(a2,b)*delta_I(a1,bc1,a3,bc2) + Qp_r(a2,b)*delta_Q(a1,bc1,a3,bc2) - &
& Vp_r(a2,b)*delta_V(a1,bc1,a3,bc2) - Up_i(a2,b)*delta_U(a1,bc1,a3,bc2)
            zeta_r = -Ip_i(a2,b)*delta_I(a1,bc1,a3,bc2) + Qp_i(a2,b)*delta_Q(a1,bc1,a3,bc2) - &
& Vp_i(a2,b)*delta_V(a1,bc1,a3,bc2) + Up_r(a2,b)*delta_U(a1,bc1,a3,bc2)
          endif

        endif 
        ratfaq_r = ratfaq_r + zeta_r 
        ratfaq_i = ratfaq_i + zeta_i 
      endif 
!    endif
      
!  end do
end do

end subroutine



subroutine rat3_a(lF1,lF2,a1,a2,a3,Ip_r,Ip_i,Qp_r,Qp_i,Up_r,Up_i,Vp_r,Vp_i,delta_I,delta_Q,delta_U,&
&delta_V,ratfaq_r,ratfaq_i)
! we take the correct D&W 1990 formula A.39 for these
! checked 19-06-2017
! consistent with A.39
implicit none
integer a1,a2,a3,b1,b2,b,mb,bc1,bc2,ma1,ma2,ma3,lF1,lF2,dm,F1,F2
double precision ratfaq_r,ratfaq_i,zeta_r,zeta_i
double precision, dimension (lF1,lF1) :: rho_a_r,rho_a_i
double precision, dimension (lF1,lF2) :: Ip_r,Ip_i,Qp_r,Qp_i,Up_r,Up_i,Vp_r,Vp_i
double precision, dimension (lF1,3,lF1,3) :: delta_I,delta_Q,delta_U,delta_V

F1 = (lF1 - 1)/2
F2 = (lF2 - 1)/2

ma1 = a1 - F1 - 1
ma2 = a2 - F1 - 1

ratfaq_r = 0.d0
ratfaq_i = 0.d0

do bc1 =1,3
  mb = ma2 + (bc1 - 2)
  b = mb + F2 + 1

  zeta_r = 0.d0
  zeta_i = 0.d0

      ma3 = a3 - F1 - 1
      dm = mb - ma3

      if (abs(dm).le.1) then
        bc2 = dm + 2

        if (b.gt.0.and.b.le.lF2) then
 
    ! if bc1+bc2 is even, then delta is real, otherwise it is imaginary
          if (mod(bc1+bc2,2).eq.0) then
            zeta_r = Ip_r(a1,b)*delta_I(a2,bc1,a3,bc2) - Qp_r(a1,b)*delta_Q(a2,bc1,a3,bc2) + &
& Vp_r(a1,b)*delta_V(a2,bc1,a3,bc2) + Up_i(a1,b)*delta_U(a2,bc1,a3,bc2)
            zeta_i = Ip_i(a1,b)*delta_I(a2,bc1,a3,bc2) - Qp_i(a1,b)*delta_Q(a2,bc1,a3,bc2) + &
& Vp_i(a1,b)*delta_V(a2,bc1,a3,bc2) - Up_r(a1,b)*delta_U(a2,bc1,a3,bc2)
          else
            zeta_i = Ip_r(a1,b)*delta_I(a2,bc1,a3,bc2) - Qp_r(a1,b)*delta_Q(a2,bc1,a3,bc2) + &
& Vp_r(a1,b)*delta_V(a2,bc1,a3,bc2) + Up_i(a1,b)*delta_U(a2,bc1,a3,bc2)
            zeta_r = -Ip_i(a1,b)*delta_I(a2,bc1,a3,bc2) + Qp_i(a1,b)*delta_Q(a2,bc1,a3,bc2) - &
& Vp_i(a1,b)*delta_V(a2,bc1,a3,bc2) + Up_r(a1,b)*delta_U(a2,bc1,a3,bc2)
          endif

        endif
 
        ratfaq_r = ratfaq_r + zeta_r 
        ratfaq_i = ratfaq_i + zeta_i 
      endif
!    endif

!  end do
end do

end subroutine

subroutine rat1_b(lF1,lF2,a1,a2,b1,b2,Ip_r,Ip_i,Qp_r,Qp_i,Up_r,Up_i,Vp_r,Vp_i,&
&delta_I,delta_Q,delta_U,delta_V,ratfaq_r,ratfaq_i)
implicit none
integer a1,a2,b1,b2,bc1,bc2,ma1,ma2,lF1,lF2,F1,F2
integer ac1,ac2,mb1,mb2
double precision ratfaq_r,ratfaq_i,zeta_r,zeta_i
double precision Ifp_r,Ifp_i,Qfp_r,Qfp_i,Ufp_r,Ufp_i,Vfp_r,Vfp_i
double precision, dimension (lF1,lF1) :: rho_a_r,rho_a_i
double precision, dimension (lF1,3,lF1,3) :: delta_I,delta_Q,delta_U,delta_V
double precision, dimension (lF1,lF2) :: Ip_r,Ip_i,Qp_r,Qp_i,Up_r,Up_i,Vp_r,Vp_i


F1 = (lF1 - 1)/2
F2 = (lF2 - 1)/2

ma1 = a1 - 1 - F1
ma2 = a2 - 1 - F1

mb1 = b1 - 1 - F2
mb2 = b2 - 1 - F2

ratfaq_r = 0.d0
ratfaq_i = 0.d0

bc1 = mb1 - ma1 + 2 
bc2 = mb2 - ma2 + 2

if (bc1.gt.0.and.bc1.le.3.and.bc2.gt.0.and.bc2.le.3) then

  Ifp_r = Ip_r(a1,b2) + Ip_r(a2,b1)
  Ifp_i = Ip_i(a1,b2) - Ip_i(a2,b1) !gamma plus and gamma minus contributions
  Ufp_r = Up_r(a1,b2) + Up_r(a2,b1)
  Ufp_i = Up_i(a1,b2) - Up_i(a2,b1) !gamma plus and gamma minus contributions
  Qfp_r = Qp_r(a1,b2) + Qp_r(a2,b1)
  Qfp_i = Qp_i(a1,b2) - Qp_i(a2,b1) !gamma plus and gamma minus contributions
  Vfp_r = Vp_r(a1,b2) + Vp_r(a2,b1)
  Vfp_i = Vp_i(a1,b2) - Vp_i(a2,b1) !gamma plus and gamma minus contributions

  if (mod(bc1+bc2,2).eq.0) then
    zeta_r = Ifp_r*delta_I(a1,bc1,a2,bc2) - Qfp_r*delta_Q(a1,bc1,a2,bc2) + Vfp_r*delta_V(a1,bc1,a2,bc2) &
&+ Ufp_i*delta_U(a1,bc1,a2,bc2)
    zeta_i = -Ufp_r*delta_U(a1,bc1,a2,bc2) + Ifp_i*delta_I(a1,bc1,a2,bc2) - Qfp_i*delta_Q(a1,bc1,a2,bc2)&
& + Vfp_i*delta_V(a1,bc1,a2,bc2)
  else
    zeta_i = Ifp_r*delta_I(a1,bc1,a2,bc2) - Qfp_r*delta_Q(a1,bc1,a2,bc2) + Vfp_r*delta_V(a1,bc1,a2,bc2) &
&+ Ufp_i*delta_U(a1,bc1,a2,bc2)
    zeta_r = Ufp_r*delta_U(a1,bc1,a2,bc2) - Ifp_i*delta_I(a1,bc1,a2,bc2) + Qfp_i*delta_Q(a1,bc1,a2,bc2) &
&- Vfp_i*delta_V(a1,bc1,a2,bc2)
  endif

  ratfaq_r = zeta_r
  ratfaq_i = zeta_i

endif

end subroutine


subroutine rat2_b(lF1,lF2,b1,b2,b3,Ip_r,Ip_i,Qp_r,Qp_i,Up_r,Up_i,Vp_r,Vp_i,&
&delta_I,delta_Q,delta_U,delta_V,ratfaq_r,ratfaq_i)
implicit none
integer a1,a2,a3,b1,b2,b,mb,bc1,bc2,ma1,ma2,ma3,lF1,lF2,dm,F1,F2
integer ac1,ac2,b3,ma,mb1,mb2,mb3
double precision ratfaq_r,ratfaq_i,zeta_r,zeta_i
double precision, dimension (lF2,lF2) :: rho_b_r,rho_b_i
double precision, dimension (lF1,lF2) :: Ip_r,Ip_i,Qp_r,Qp_i,Up_r,Up_i,Vp_r,Vp_i
double precision, dimension (lF1,3,lF1,3) :: delta_I,delta_Q,delta_U,delta_V



F1 = (lF1 - 1)/2
F2 = (lF2 - 1)/2

mb1 = b1 - F2 - 1
mb2 = b2 - F2 - 1

ratfaq_r = 0.d0 
ratfaq_i = 0.d0

do ac1 =1,3
  ma = mb1 - (ac1 - 2) 
  a1 = ma + F1 + 1
      mb3 = b3 - F2 - 1
      dm = -(ma - mb3)
  
      if (abs(dm).le.1.and.a1.gt.0.and.a1.le.lF1) then
        bc1 = ac1
        bc2 = dm + 2
! check bc2 and bc1 

  ! if bc1+bc2 is even, then delta is real, otherwise it is imaginary
        if (mod(bc1+bc2,2).eq.0) then
          zeta_r = Ip_r(a1,b2)*delta_I(a1,bc1,a1,bc2) - &
&Qp_r(a1,b2)*delta_Q(a1,bc1,a1,bc2) + Vp_r(a1,b2)*delta_V(a1,bc1,a1,bc2) + Up_i(a1,b2)*delta_U(a1,bc1,a1,bc2)
          zeta_i = Ip_i(a1,b2)*delta_I(a1,bc1,a1,bc2) - &
&Qp_i(a1,b2)*delta_Q(a1,bc1,a1,bc2) + Vp_i(a1,b2)*delta_V(a1,bc1,a1,bc2) - Up_r(a1,b2)*delta_U(a1,bc1,a1,bc2)
        else
          zeta_i = Ip_r(a1,b2)*delta_I(a1,bc1,a1,bc2) -& 
&Qp_r(a1,b2)*delta_Q(a1,bc1,a1,bc2) + Vp_r(a1,b2)*delta_V(a1,bc1,a1,bc2) + Up_i(a1,b2)*delta_U(a1,bc1,a1,bc2)
          zeta_r = -Ip_i(a1,b2)*delta_I(a1,bc1,a1,bc2) + &
&Qp_i(a1,b2)*delta_Q(a1,bc1,a1,bc2) - Vp_i(a1,b2)*delta_V(a1,bc1,a1,bc2) + Up_r(a1,b2)*delta_U(a1,bc1,a1,bc2)
        endif
 
        ratfaq_r = ratfaq_r + zeta_r 
        ratfaq_i = ratfaq_i + zeta_i 
      endif 
!    endif
!      
!  end do
end do

end subroutine



subroutine rat3_b(lF1,lF2,b1,b2,b3,Ip_r,Ip_i,Qp_r,Qp_i,Up_r,Up_i,Vp_r,Vp_i,&
&delta_I,delta_Q,delta_U,delta_V,ratfaq_r,ratfaq_i)
implicit none
integer a1,a2,a3,b1,b2,b,mb,bc1,bc2,ma1,ma2,ma3,lF1,lF2,dm,F1,F2
integer ac1,ac2,b3,ma,mb1,mb2,mb3
double precision ratfaq_r,ratfaq_i,zeta_r,zeta_i
double precision, dimension (lF2,lF2) :: rho_b_r,rho_b_i
double precision, dimension (lF1,lF2) :: Ip_r,Ip_i,Qp_r,Qp_i,Up_r,Up_i,Vp_r,Vp_i
double precision, dimension (lF1,3,lF1,3) :: delta_I,delta_Q,delta_U,delta_V

F1 = (lF1 - 1)/2
F2 = (lF2 - 1)/2

mb1 = b1 - F2 - 1
mb2 = b2 - F2 - 1

ratfaq_r = 0.d0
ratfaq_i = 0.d0

do ac1 =1,3
  ma = mb2 - (ac1 - 2)
  a1 = ma + F1 + 1
      mb3 = b3 - F2 - 1
      dm = -(ma - mb3)

      if (abs(dm).le.1.and.a1.gt.0.and.a1.le.lF1) then
        bc1 = ac1
        bc2 = dm + 2
! check bc2 and bc1 

  ! if bc1+bc2 is even, then delta is real, otherwise it is imaginary
        if (mod(bc1+bc2,2).eq.0) then
          zeta_r = Ip_r(a1,b1)*delta_I(a1,bc2,a1,bc1) - Qp_r(a1,b1)*delta_Q(a1,bc2,a1,bc1)&
& + Vp_r(a1,b1)*delta_V(a1,bc2,a1,bc1) - Up_i(a1,b1)*delta_U(a1,bc2,a1,bc1)
    !\gam_p is conjugated---it is a minus
          zeta_i = -Ip_i(a1,b1)*delta_I(a1,bc2,a1,bc1) + Qp_i(a1,b1)*delta_Q(a1,bc2,a1,bc1)&
& - Vp_i(a1,b1)*delta_V(a1,bc2,a1,bc1) - Up_r(a1,b1)*delta_U(a1,bc2,a1,bc1)
        else
          zeta_i = Ip_r(a1,b1)*delta_I(a1,bc2,a1,bc1) - Qp_r(a1,b1)*delta_Q(a1,bc2,a1,bc1)&
& + Vp_r(a1,b1)*delta_V(a1,bc2,a1,bc1) - Up_i(a1,b1)*delta_U(a1,bc2,a1,bc1)
          zeta_r = Ip_i(a1,b1)*delta_I(a1,bc2,a1,bc1) - Qp_i(a1,b1)*delta_Q(a1,bc2,a1,bc1)&
& + Vp_i(a1,b1)*delta_V(a1,bc2,a1,bc1) + Up_r(a1,b1)*delta_U(a1,bc2,a1,bc1)
        endif

        ratfaq_r = ratfaq_r + zeta_r 
        ratfaq_i = ratfaq_i + zeta_i 
      endif
!    endif
!
!  end do
end do

end subroutine

subroutine matvec_prod(A,b,N,c)
!compute the product A*b = c
!A = NxN matrix
!b = N-vector
!c = N-vector
integer N,I,J
double precision, dimension (N) :: b,c 
double precision, dimension (N,N) :: A

do I=1,N
  c(I) = 0.d0
end do

do I = 1,N
  do J = 1,N
    c(I) = c(I) + A(I,J)*b(J)
  end do
end do

end subroutine


