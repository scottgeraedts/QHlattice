! module for the coulomb energy for lattice-based
! monte carlo on the torus
! to use:
! link with z_function_m.f90
! (calls get_l(norb,l1,l2) to get geometry)
! 
! call coulomb_setup to initialize table
! then call 
! coulomb_energy(nel,norb,x,u)
! where (x(1,i), x(2,i)), i = 1,.., NEL
! are the lattice coordinates of the particles
! to get the energy u.
!
module new_coulomb_m
integer, parameter :: dp = kind(1.0d0)
integer :: norb_coulomb = 0
complex (kind=dp) :: l1_coulomb, l2_coulomb
real (kind=dp), allocatable :: coulomb(:,:)

end module new_coulomb_m

function new_coulomb(norb,m,n) result(v)
  use new_coulomb_m
  integer, intent(in) :: norb,m,n
  real(kind=dp) :: v
  
  if(norb /= norb_coulomb) then
     write(6,'(" coulomb: norb mismatch:",2i10)') norb, norb_coulomb
     stop
  endif
  
  v = coulomb(1 + modulo(m-1,norb), 1 + modulo(n-1,norb))
  return
end function new_coulomb

subroutine coulomb_energy(nel,norb,x,u)
  use new_coulomb_m
  implicit none
  integer, intent(in) :: nel, norb
  integer, intent(in) :: x(2,nel)
  real(kind=dp), intent(out) :: u


  integer :: i,j,k,m,n
  u = 0_dp
  if(norb /= norb_coulomb) then
     write(6,'(" COULOMB_ENERGY: mismatched norb",2i5)') norb,norb_coulomb
     stop
  endif
! madelung energy
  u = nel*coulomb(norb,norb)/2
  do i = 1,nel
     do j = 1, i-1
        m = 1 + modulo(x(1,i) - x(1,j) -1,norb)
        n = 1 + modulo(x(2,i) - x(2,j) -1,norb)
        if(m == norb .and. n == norb) then
           write(6,'("infinite coulomb energy: two particles",&
                &" on same site")')
           do k = 1,nel
              write(6,'(2i5)') x(:,k)
           enddo
           stop
        endif
        u = u + coulomb(m,n)
     enddo
  enddo
  return
end subroutine coulomb_energy



function new_v_coulomb(norb,m,n, l1,l2) result(v)
  implicit none
  integer, parameter :: dp = kind(1.0d0)
  integer, intent(in) :: norb,m,n
  complex (kind=dp), intent(in)  :: l1,l2
  real(kind=dp) :: v


  real(kind=dp):: v1,v2
  real(kind=dp), parameter :: a = 1.0_dp
  !write(*,*) "a=", a
 
  call coulomb_1(m,n,norb,l1,l2,a,v1,v2)
  v = v1+v2
  return
end function new_v_coulomb


subroutine coulomb_setup
  use new_coulomb_m
  implicit none
  real(kind=dp) :: a, v1,v2
  complex(kind=dp) :: l1,l2, new_tau,l1a,l2a
  integer :: norb, m, n, sl2z(2,2), m1,n1
  
  call get_l(norb,l1,l2)
  if(allocated(coulomb)) deallocate(coulomb)
  allocate (coulomb(norb,norb))
  l1_coulomb = l1
  l2_coulomb = l2
  norb_coulomb = norb

  call optimize_tau(l2/l1,new_tau,sl2z)
  write(6,'(2i5)') sl2z(1,1), sl2z(1,2),sl2z(2,1),sl2z(2,2)


  a   =real(1,kind=dp)
  do  m = 1,norb_coulomb
     do n = 1,norb_coulomb
        m1 = sl2z(2,2)*m - sl2z(2,1)*n
        n1 = -sl2z(1,2)*m + sl2z(1,1)*n
        l1a = sl2z(1,1)*l1 + sl2z(1,2)*l2
        l2a = sl2z(2,1)*l1 + sl2z(2,2)*l2
        call coulomb_1(m1,n1,norb,l1a,l2a,a,v1,v2)
        coulomb(m,n) = v1 + v2
!        call coulomb_1(m,n,norb,l1,l2,a,v1,v2)
!        write(6,'(2f25.15,e12.3)') coulomb(m,n), v1 + v2, v1+v2 -coulomb(m,n)

     enddo
  enddo
  return
end subroutine coulomb_setup


 subroutine coulomb_1(m,n,norb,l10,l20,a,v1,v2)
    implicit none
    integer, parameter :: dp = kind(1.0d0)
    integer, intent(in) :: m,n,norb
    complex(kind=dp), intent(in) :: L10,L20
    real(kind=dp),intent(in) :: a
    real(kind=dp), intent(out) :: v1,v2
!-----------------------------------------
! v1 + v2 is the 
! periodic Coulomb interaction on the lattice
!  z(m,n)  (m*L1 + n*L2)/norb, with periodicity
! under (m,n) -> (m + i*norb, n + j*norb)
! and the metric  d(m,n) = |z(m,n)|
! length units are scaled so the unit cell area is
! 2*pi*norb
!
! For (m,n) = (0,0) modulo norb, the madelung
! constant is returned.
! in exact arithmetic  v1 + v2 is independent of a
!------------------------------
    logical, parameter :: test = .false.

! length units are scaled so area = 2*pi*norb    
! (area of pbc unit cell)

    complex (kind=dp) :: l1,l2,g1,g2,z,z1,z2
    real(kind=dp) :: area, scale, rtpi, rt2
    real(kind=dp) :: alpha ,x,x1,x2,x3,x4,q,q1,q2,vprev
    real(kind=dp) :: factor, factor1, factor2
    integer :: m1,n1,i,j
    real(kind=dp), allocatable, save :: cosine(:)
    integer, save :: norb_prev = 0
    real(kind=dp), save :: pi



    if (norb /= norb_prev) then
       pi = 2*asin(real(1,kind=dp))
       if(allocated(cosine)) deallocate(cosine)
       allocate(cosine(0:norb-1))
       cosine(0) = real(1,kind=dp)
       if(mod(norb,2) == 0) cosine(norb/2) = -cosine(0)
       do i = 1,(norb-1)/2
          cosine(i) = cos(2*i*pi/norb)
          cosine(norb -i) = cosine(i)
       enddo
       norb_prev = norb
    endif

!    do i = 0,norb -1
!       write(6,'(i5,f25.16)') i, cosine(i)
!    enddo
       
    area = abs(aimag(conjg(L10)*L20))
    scale =  sqrt(area/(2*pi*norb))
    area= 2*pi*norb
    l1 = l10/scale
    l2 = l20/scale
    g1 = cmplx(0,-1,kind=dp)*l2/norb
    g2 = cmplx(0,1,kind=dp)*l1/norb
    
!    write(6,'(2f25.16)') l1,l2,g1,g2

!    write(6,'(2f25.16)') conjg(g1)*l1
!    write(6,'(2f25.16)') conjg(g1)*l2
!    write(6,'(2f25.16)') conjg(g2)*l1
!    write(6,'(2f25.16)') conjg(g2)*l2



    rtpi = sqrt(pi)
    rt2 = sqrt(real(2,kind=dp))
    
    m1 = modulo(m,norb)
    if (2*m1 > norb) m1 = m1 - norb
    n1 = modulo(n,norb)
    if (2*n1 > norb) n1 = n1 - norb
    z = (m1*L1 + n1*L2)/norb
    
    ! |x(m,n)| =  |L1|*abs(m + n*tau)/norb
    !           = s1*abs(m + n*tau)           
    
    alpha = 1/(rt2*sqrt(norb*a**2))
    if(m1 == 0 .and. n1 == 0) then
       v1 = -(2/rtpi)*alpha
    else
       x = abs(z) 
       v1 =  erfc(alpha*x)/x
    endif
    i = 0
    do
       if(test) write(6,'("x: i=",i12)') i
       i = i + 1
       x1 = abs(z + i*L1)
       x2 = abs(z - i*L1)
       vprev = v1
       v1 = v1 + erfc(alpha*x1)/x1  + erfc(alpha*x2)/x2
       if(v1 == vprev) exit
    enddo
    j= 0
    do
       j = j + 1
       z1 = z + j*L2
       z2 = z - j*L2
       x1 = abs(z1)
       x2 = abs(z2)
       vprev = v1
       v1 = v1 + erfc(alpha*x1)/x1  + erfc(alpha*x2)/x2
       if(v1 == vprev) exit
       i = 0
       do
          if(test) write(6,'("x: j,i=",2i12)') j,i
          i =  i  + 1
          x1 = abs(z1 + i*L1)
          x2 = abs(z1 - i*L1)
          x3 = abs(z2 + i*L1)
          x4 = abs(z2 - i*L1)
          vprev = v1
          v1 = v1 + erfc(alpha*x1)/x1  + erfc(alpha*x2)/x2 &
               & + erfc(alpha*x3)/x3  + erfc(alpha*x4)/x4 
          if(v1 == vprev) exit
       enddo
    enddo
    
    alpha = sqrt(norb*a**2)/rt2
    v2 = - (2/rtpi)*alpha
    
    i = 0
    do
       if(test) write(6,'("q: i=",i12)') i
       i = i + 1
       q = abs(i*G1)
       factor = 2*erfc(alpha*q)/q
       if(v2 + factor == v2) exit
       v2 = v2 + factor*cosine(modulo(i*m1,norb))
    enddo
    j = 0
    do 
       j = j + 1
       q = abs(j*G2)
       factor = 2*erfc(alpha*q)/q     
       if(v2 + factor == v2) exit
       v2 = v2 + factor*cosine(modulo(j*n1,norb))
       i = 0
       do 
          if(test) write(6,'("q: j,i=",2i12)') j,i
          i = i + 1
          q1 = abs(j*G2 + i*G1)
          q2 = abs(j*G2 - i*G1)
          factor1 = 2*erfc(alpha*q1)/q1     
          factor2 = 2*erfc(alpha*q2)/q2     
          if(v2 + factor1 + factor2 == v2) exit
          v2 = v2 + factor1*cosine(modulo(j*n1+ i*m1,norb)) &
               &  + factor2*cosine(modulo(j*n1 -i*m1,norb)) 
       enddo
    enddo
    v2 = v2*(2*pi/area)
    
    write(*,*) v1
    write(*,*) v2
    
    return
  end subroutine coulomb_1

