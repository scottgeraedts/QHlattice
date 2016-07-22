subroutine compactified_gaussian(norb,u,norm)
  implicit none
  integer, parameter :: dp = kind(1.0d0)
  integer, intent (in) :: norb

  real(kind=dp), intent(out) :: u(norb,norb), norm
  
  !----------------------------
  ! U(q(i,j)) =   sum_{m,n}  exp ((-0.5*|(q(i,j) + Norb*q(m,n))|**2) ) * 
  !             (-1)**(i*n) * (-1)**(j*m) + eta(m,n)**norb
  !
  ! eta(m,n) = (-1)**(m + n + m*n).
  ! 
  !   exp (i (q(i,j) + Norb*q(m,n)) . R) =  (eta(m,n)**norb) *(-1)*(m*j)* (-1)**(n*i)
  !    times exp (i q(i,j) . R)
  
  complex(kind=dp) :: q, q1, q2_plus,q2_minus, g1,g2
  complex (kind=dp) :: l1,l2  
  
  real(kind=dp) :: abs_u, test_two_rows, row, two_rows, u_plus ,u_minus, u_sum
  integer :: m,n, i,j, sign, i1, j1, norb_z
  
  
  
  
  call get_l(norb_z,l1,l2)
  if(norb_z /= norb) then
     write(6,'("COMPACIFIED GAUSSIAN: nirn mismatch",2i5)') norb, norb_z
  endif

  
  g1 = cmplx(0,-1,kind=dp)*l2
  g2 = cmplx(0,1,kind=dp)*l1
  
  ! conjg(g1)*l1 + g1*conjg(l1) =  2*pi*norb
  ! conjg(g1)*l2 + g1*conjg(l2) =  0
  ! conjg(g2)*l1 + g2*conjg(l1) =  0
  ! conjg(g2)*l2 + g2*conjg(l2) =  2*pi*norb
  
  
  
  u = cmplx(0,kind=dp)
  
  do i = 1,norb
     do j = 1,norb
        q =  (i*g1 + j*g2)/norb
        abs_u = 0_dp
        m = 0 
        do    ! loop over m = 0,1,2,.....
           two_rows= 0_dp
           test_two_rows = 0_dp
           do sign = 1,-1,-2    ! drop sign = -1 case if m = 0
              q1 = q + sign*m*g1           ! n = 0 case
              row = exp(-real(conjg(q1)*q1)/2)
              if(mod(norb*m + j*m,2) /=0) row = -row
              n = 1
              do    ! loop over n = 1 , 2,.....
                 q2_plus  = q1 + n*g2
                 q2_minus = q1 - n*g2
                 u_plus =  exp(-real(conjg(q2_plus)*q2_plus)/2)
                 u_minus =  exp(-real(conjg(q2_minus)*q2_minus)/2)
                 u_sum = u_plus + u_minus
                 if(n > 2 .and. abs_u + u_sum == abs_u) exit
                 if(mod(norb*(m + n + m*n) + i*n + j*m,2) == 0) then
                    row = row + u_sum
                 else
                    row = row - u_sum
                 endif
                 n = n + 1
              enddo
              two_rows = two_rows + row
              test_two_rows = test_two_rows + abs(row)
              if(m == 0) exit
           enddo
           if(m > 2 .and. abs_u + test_two_rows == abs_u) exit
           u(i,j) = u(i,j) + two_rows
           abs_u = abs(u(i,j))
           m = m + 1
        enddo
     enddo
  enddo
  
  norm = 1/abs(u(norb,norb))
  u = norm*u
  
  
  
  return
end subroutine compactified_gaussian



subroutine make_gaussian_projector(norb,l1,l2,projector)
  implicit none
  integer, parameter :: dp = kind(1.0d0)
  integer, intent(in) :: norb
  complex (kind=dp), intent(in) :: l1,l2
  complex (kind=dp), intent(out) :: projector(norb,norb,norb,norb)
  

  complex(kind=dp), allocatable :: omega(:)
  real(kind=dp), allocatable :: u(:,:)
  real(kind=dp) :: pin, norm
  integer :: i1,j1,i2,j2,m,n,i

  allocate(u(-norb:norb,-norb:norb),omega(0:2*norb-1))
  call compactified_gaussian(norb,l1,l2,u,norm)

  do i1 = -norb,norb
     do j1 = -norb,norb
        write(6,'(3i5,f25.16)') i1,j1,norb,u(i1,j1)
     enddo
  enddo
        



  pin = 1_dp
  pin = 2*asin(pin)/norb
  do i = 0,2*norb-1
     omega(i) = cmplx(cos(i*pin),sin(i*pin),kind=dp)
  enddo


  projector = cmplx(0,kind=dp)
  do i1 = 1,norb
     do j1 = 1,norb
        do i2 = 1,norb
           do j2 = 1,norb
              m = j1 -j2
              n = i2 -i1
              projector(i1,j1,i2,j2) = u(m,n)*omega(modulo(m*i1+n*j1,2*norb))
           enddo
        enddo
     enddo
  enddo

  deallocate(u,omega)
  return
  
end subroutine make_gaussian_projector
