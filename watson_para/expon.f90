subroutine exp_A(A,ith,C)
!subroutine to compute the exponent of a 4x4 matrix
!C = exp(A);
implicit none
integer ith,i,j,N
double precision thresh,fac,diff,x1
double precision, dimension(4,4) :: A,B,C,D

x1 = 1.d0
thresh = A(1,1)*10.d0**(-x1*ith) ! threshold is 10^-ith times the A-factor 
!write(*,*)A(1,1),thresh

!write(*,*)A(1,1),A(1,2),A(1,3),A(1,4)
!write(*,*)A(2,1),A(2,2),A(2,3),A(2,4)
!write(*,*)A(3,1),A(3,2),A(3,3),A(3,4)
!write(*,*)A(4,1),A(4,2),A(4,3),A(4,4)
!
!write(*,*)
!write(*,*)

!initiate B and C
do i = 1,4
  do j = 1,4
    if (i.eq.j) then
      B(i,j) = x1
      C(i,j) = x1
    else
      B(i,j) = 0.d0
      C(i,j) = 0.d0
    endif
  end do
end do


N = 0
fac = x1

101 N = N+1
  
call prod_44(A,B,D) !product of A * B = D 
                    !A = original matrix
                    !B = A^(N-1) 
                    !D = A*B= B^N


!write(*,*)N
!
!write(*,*)D(1,1),D(1,2),D(1,3),D(1,4)
!write(*,*)D(2,1),D(2,2),D(2,3),D(2,4)
!write(*,*)D(3,1),D(3,2),D(3,3),D(3,4)
!write(*,*)D(4,1),D(4,2),D(4,3),D(4,4)
!

fac = fac*(x1*N)

!write(*,*)
!write(*,*)fac
!write(*,*)

do i = 1,4
  do j = 1,4 

    C(i,j) = C(i,j) + D(i,j)/fac
    B(i,j) = D(i,j)

  end do
end do

!write(*,*)C(1,1),C(1,2),C(1,3),C(1,4)
!write(*,*)C(2,1),C(2,2),C(2,3),C(2,4)
!write(*,*)C(3,1),C(3,2),C(3,3),C(3,4)
!write(*,*)C(4,1),C(4,2),C(4,3),C(4,4)
!
!write(*,*)
!write(*,*)


do i = 1,4
  do j = 1,4

    diff = D(i,j)/fac 
    if (dabs(diff).gt.dabs(thresh)) then
      goto 101
    endif

  end do
end do      


end subroutine


subroutine prod_44(A,B,C) 
!subroutine to compute the product of A^T * B = C for two 4x4 matrices
implicit none
integer i,j,k
double precision x0
double precision, dimension(4,4) :: A,B,C

x0 = 0.d0

do i = 1,4

  do j = 1,4
    C(i,j) = x0

    do k = 1,4
      C(i,j) = C(i,j) + A(i,k)*B(k,j) 
    end do

  end do
end do

end subroutine

