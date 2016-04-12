! coulomb_m.f90
! v0.01 2016-03-20
! F. Duncan. M. Haldane, haldane@princeton.edu




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

module coulomb_m
integer, parameter :: dp = kind(1.0d0)
integer :: norb_coulomb = 0
complex (kind=dp) :: l1_coulomb, l2_coulomb
real (kind=dp), allocatable :: coulomb(:,:)

end module coulomb_m

subroutine coulomb_setup
  use coulomb_m
  implicit none
  integer :: m, n, norb

  real (kind=dp) ::  v_coulomb
  complex(kind=dp) :: l1,l2



  call get_l(norb,l1,l2)
  if(allocated(coulomb)) deallocate(coulomb)
  allocate (coulomb(norb,norb))


  do m = 1,norb
     do n = 1,norb
        coulomb(m,n) = v_coulomb(norb,m,n, l1,l2)
     enddo
  enddo
  norb_coulomb = norb
  l1_coulomb = l1
  l2_coulomb = l2
  return
end subroutine coulomb_setup

subroutine coulomb_energy(nel,norb,x,u)
  use coulomb_m
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



function v_coulomb(norb,m,n, l1,l2) result(v)
  implicit none
  integer, parameter :: dp = kind(1.0d0)
  integer, intent(in) :: norb,m,n
  complex (kind=dp), intent(in)  :: l1,l2
  real(kind=dp) :: v


  real(kind=dp):: v1,v2
  real(kind=dp), parameter :: alpha = 3_dp
 
  call coulomb_parts(norb,m,n,l2/l1,alpha,v1,v2)
  v = v1+v2
  return
end function v_coulomb

subroutine coulomb_parts(norb,m,n,tau,alpha,v1,v2)
  implicit none
  integer, parameter :: dp = kind(1.0d0)
  integer, intent(in) :: norb,m,n
  complex (kind=dp), intent(in)  :: tau
  real(kind=dp), intent(in) :: alpha
  real(kind=dp),intent(out) :: v1,v2
  
  !------------------------------------
  !  periodic version of coulomb interaction
  !  v =  v1 + v2  = ell_B/|x|,   area of cell = twopi*norb*(ell_B**2)
  !
  !  x is restricted to the lattice
  ! x(m,n) = (m*L1 + n*L2)/norb
  !  where L1 and L1 are primitive pbc translations.
  !
  !   
  !  alpha parametrizes a split into short and long range
  !  parts v1 + v2, which  are calulated separately.
  !  In exact arithmetic, v1 + v2 is independent
  !  of alpha.     The best alpha seems to be around pi.
  !
  !   the madelung energy is returned for (m,n) = (0,0) modulo norb
  !--------------------------------------
  real(kind=dp) ::  pin,a,x,q, factor,a0
  real (kind=dp), save :: pi = 1_dp,  twopi, rootpi, rt2
  logical, save  :: have_pi = .false.
  integer, save :: norb_table = 0
  real (kind=dp), allocatable, save :: cos_factor(:)
  integer :: m1,n1, i,j
  real (kind=dp) :: vprev1, vprev2, increment, vft

  if(.not. have_pi) then
     pi = 2*asin(pi)
     have_pi = .true.
     twopi = 2*pi
     rootpi = sqrt(pi)
     rt2 = sqrt(real(2,kind=dp))
  endif
  

  
  if(norb /= norb_table) then
     if(allocated(cos_factor)) deallocate(cos_factor)
     allocate (cos_factor(0:norb-1))
     norb_table = norb
     pin = twopi/norb
     forall (i = 0:norb-1) cos_factor(i) = cos(i*pin)      
  endif
  

  
  m1 = modulo(m,norb)
  n1 = modulo(n,norb)
  if( m1 > norb/2) m1 = m1 - norb
  if( n1 > norb/2) n1 = n1 - norb

  

  a0 = sqrt(pi/(norb*abs(aimag(tau))))
  factor = 1/(rt2*a0)

  a= a0*alpha
  if(m1 == 0 .and. n1 == 0) then
     v1 = -2*a/rootpi
  else
     x = abs(cmplx(m1,kind=dp) + n1*tau)
     v1 = erfc(a*x)/x
  endif

  i = 0
  do
     vprev1 = v1
     j = 0
     do
        vprev2 = v1
        if(i == 0 .and. j == 0) then
          j = j + norb
          cycle
       endif
       x = abs(cmplx(m1 + i,kind=dp) + (n1 + j)*tau)
       v1 = v1 + erfc(a*x)/x
        if(j /= 0 ) then
           x = abs(cmplx(m1 + i,kind=dp) + (n1 - j)*tau)
           v1 = v1 + erfc(a*x)/x
        endif        
        if(i /= 0 ) then
           x = abs(cmplx(m1 - i,kind=dp) + (n1 + j)*tau)
           v1 = v1 + erfc(a*x)/x
        endif
        
        if(i /= 0 .and. j /= 0) then
           x = abs(cmplx(m1 - i,kind=dp) + (n1 - j)*tau)
           v1 = v1 + erfc(a*x)/x
        endif
        if (v1 == vprev2) exit
        j = j + norb
     enddo
     if(v1 == vprev1) exit
     i = i + norb
  enddo
  v1 = factor*v1


  factor= factor/norb
  a = a0/alpha
 v2  = -2*a/rootpi
  i = 0
  do
     vprev1 = v2
     j = 0
     do 
        vprev2 = v2 
        if(i == 0 .and. j == 0) then
           j = j + 1
           cycle
        endif
        q = abs(cmplx(i,kind=dp)  + j*tau)
        vft = erfc(a*q)/q
        increment = vft
        v2 = v2 + vft*cos_factor(modulo(m1*j - n1*i,norb))
        if(i /= 0) then
           q = abs(cmplx(-i,kind=dp)  + j*tau)
           vft = erfc(a*q)/q
           increment = increment + vft
           v2 = v2 + vft*cos_factor(modulo(m1*j + n1*i,norb))
        endif
        
        if(j /= 0) then
           q = abs(cmplx(i,kind=dp)  - j*tau)
           vft = erfc(a*q)/q
           increment = increment + vft
           v2 = v2 + vft*cos_factor(modulo(-m1*j -n1*i,norb))
        endif
        
        if(i /= 0 .and. j /= 0) then
           q = abs(cmplx(-i,kind=dp)  - j*tau)
           vft = erfc(a*q)/q
           increment = increment + vft
           v2 = v2 + vft*cos_factor(modulo(-m1*j + n1*i,norb))
        endif
        if(increment + vprev2 == vprev2) exit
        j = j + 1
     enddo
     if(v2 == vprev1) exit
     i = i + 1
  enddo
  v2 = factor*v2

  return
end subroutine coulomb_parts

