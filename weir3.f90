!program main
!implicit none

!!call test_weierstrass_sigma

!!call test_jacobi
!call test_landau_basis
!stop
!end program main


subroutine test_landau_basis
implicit none
integer, parameter :: dp = kind(1.0d0)

integer :: norb, m, j, n, jk, k,kk
complex (kind=dp) :: l1,l2, x, z, wf,u,tau,w
real (kind=dp) :: l1x,l1y,l2x,l2y, norm, pi, theta
complex (kind=dp), allocatable :: w0(:), w1(:), factor(:),w2(:)

pi = 1_dp
pi = 2*asin(pi)

write(6,'("give norb, l1,l2")')
read(5,*) norb, l1x, l1y,l2x,l2y
l1 = cmplx(l1x,l1y,kind=dp)
l2 = cmplx(l2x,l2y,kind=dp)




allocate (w0(norb),w1(norb),w2(norb))
w0 = cmplx(0,kind=dp)
w1 = 0
w2 = w0
tau = -l1/(norb*l2)
do j = 1,norb
   u = j*pi/norb
   n = 1
   if(mod(norb,2)==0) n = 4
   call jacobi_theta(n,u,tau,w0(j),.true.)
   if(mod(j,2)/= 0) w0(j) = -w0(j)
enddo
norm = real(sum(conjg(w0)*w0))
w0 = w0/sqrt(norm)

do j = 1,norb
   write(6,'(2i5,4f25.18,es12.3)') j, norb,w0(j)
enddo

allocate (factor(0:norb))
do j = 0,norb
   theta = (pi*j)/norb
   factor(j) = cmplx(cos(theta),sin(theta),kind=dp)
enddo

w1 = cmplx(0,kind=dp)
do j = 1,norb
   do k = 1,norb
      jk = mod(j*k,norb)
      w1(j) = w1(j) + w0(k)*factor(jk)
   enddo
enddo
norm = real(sum(conjg(w1)*w1))
w1 = w1/sqrt(norm)



tau = l2/(l1*norb)
w2 = cmplx(0,kind=dp)
do j = 1,norb
   u = j*pi/norb
   n = 2
   if(mod(norb,2) == 0) n= 3
   call jacobi_theta(n,u,tau,w2(j),.true.)
   if(mod(j,2)/= 0) w2(j) = -w2(j)
enddo
norm = real(sum(conjg(w2)*w2))
w2 = w2/sqrt(norm)


write(6,'(20("-"))')
do j = 1,norb
   write(6,'(2i5,6f25.16,6es25.16)') j, norb, w0(j), w1(j), w2(j)
enddo




return
end subroutine test_landau_basis










subroutine test_weierstrass_sigma
implicit none
integer, parameter :: dp = kind(1.0d0)
complex (kind=dp) :: l1,l2,l,z,sigma, u, sigma0
real (kind=dp) :: a, pi, scale,area
integer :: i,j

pi = 1_dp
pi = 2*asin(pi)

a = 1_dp
l1 = cmplx(a,a/5,kind=dp)
l2 = cmplx(a/7,a,kind=dp)
area = real(cmplx(0,1,kind=dp)*(conjg(l1)*l2 - conjg(l2)*l1))

write(6,'("area:",f20.10)') area

scale = (2*pi)/area
if(scale < 0) then
   l2 = -l2
   scale = -scale
endif
z = (l1/5) + (l2/3)
do i = -3,3
   do j = -3,3
      l = i*l1 + j*l2
!      call weierstrass_sigma(z-l,l1,l2,sigma)
!      call weierstrass_sigma(z-l,l1,l2,sigma0,.true.)
      write(6,'("ratio",2f25.16)') sigma/sigma0
      if(mod(i,2)/=0 .or. mod(j,2)/= 0) sigma = -sigma
      u = scale*conjg(l)*(z-(l/2)) 
      write(6,'(2i5,2f20.12)') i,j,sigma*exp(u)
   enddo
enddo
return
end subroutine test_weierstrass_sigma
     


subroutine test_jacobi
implicit none

integer :: n,nn,k
integer, parameter :: dp = kind(1.0d0)
complex (kind=dp) :: theta1,theta2,tau,z
real (kind=dp) :: tau_x, tau_y, x, y, pi

pi = 1_dp
pi = 2*asin(pi)
1 write(6,'(" give tau, n")')
read(5,*) tau_x,tau_y,nn
tau = cmplx(tau_x,tau_y,kind=dp)
z = cmplx(x,y,kind=dp)


do n = 1,4
   do k = 0, 2*nn
      z = k*cmplx(pi,kind=dp)/nn
      call jacobi_theta(n,z, tau,theta1,.true.)
      call jacobi_theta(n,z, tau,theta2,.false.)
      write(6,'(3i5,4es20.12)') n,k,nn,theta1,theta2
   enddo
enddo
goto 1
end subroutine test_jacobi




subroutine weierstrass_sigma(z,l1,l2,sigma,reduce)  
  implicit none
  integer, parameter :: dp = kind(1.0d0)
  complex (kind=dp), intent(in):: z,l1,l2
  complex(kind=dp), intent (out) :: sigma
  logical, intent(in), optional :: reduce
! computes Weierstrass sigma function with periods L = {m*L1 + n*L2}
!
!   sigma(z) = sigma(z;{L})  = z *  prod_{L \= 0} f(z/L),    f(u) = (1-u)exp (u +  u**2/2)
!
!   where {L(m,n)} ={m*L1 + n*L2},   and  L1*conjg(L2) - L2*conjg(L1) = 2*A*cmplx(0,1) \= 0.
!
!    L1*G2 - L2*G1 = 2*pi*cmplx(0,1)
!
!    G1 = conjg(L1)*pi/A
!    G2 = conjg(L2)*pi/A
!
!  for the lattice  L = m*L1 + n*L2:
!  the mapping to the reciprocal lattce is  G(L) = m*G1 + n*G2, 
!  the parity of a lattice point is eta(L) = (-1)**(m+n+m*n)
!  parity of L  is even (eta(L) = 1)  if L/2 is a lattice point, 
!  parity is odd odd (eta(L) = -1)  otherwise.
!
!  the sigma function is odd:
!  sigma(z) = -sigma(-z)
!
! The sigma function is quasiperiodic  and entire:   
!   sigma(z + L)  = eta(L)* exp (G(L)*(z + L/2)) * sigma(z)
!
! the zeroes of the sigma function are at the lattice points
!  sigma(L) = 0
!
!  the optional logical argument REDUCE (if present, and .true.)
!  calculates with z -> z0 = z-L as close as possible to z=0, thet uses quasperiodicity
!  to move back to z.  
!-----------------------------------------------
  real(kind=dp) :: area, pi
  complex (kind=dp) scale, dtheta, u, tau, q_4,q2,lshift,l10,l20,z0,sigma0
  logical :: shift
  integer :: mm(2), nn(2), n(2),sign
  shift = .false.
  if(present(reduce)) then
     shift = reduce
  endif


  pi = 1_dp
  pi = 2*asin(pi)
  tau = l2/l1
  area = -aimag(l1*conjg(l2))
  if (area == 0_dp) then
     write(6,'("weierstrass_sigma invalid L1,L2: L1*L2 is real")')
     write(6,'("L1    ",2d25.16)') l1
     write(6,'("L2    ",2d25.16)') l2
     write(6,'("L1*L2 ",d25.16)') l1*l2
     stop
  endif

  call reduce_periods(l1,l2,mm,nn)
  l10 = mm(1)*l1 + nn(1)*l2
  l20 = mm(2)*l1 + nn(2)*l2



  area = aimag(conjg(l10)*l20)
  if (area < 0_dp) then
      area = -area
      l20 = -l20
   endif
   tau = l20/l10

   q_4 = exp(cmplx(0,pi/4,kind=dp)*tau)
   q2 = (q_4)**2
   q2 = q2**2
   q2 = q2**2

   sign = 1
   z0 = z
   if(shift) then

      n(1) = -nint(aimag(conjg(z)*l20)/area)
      n(2) = nint(aimag(conjg(l10)*z)/area)

      lshift = n(1)*l10 + n(2)*l20
      if(mod(n(1),2) /= 0 .or. mod(n(2),2) /= 0) sign = -1
      z0 = z - lshift
   endif

   call theta_1(z0*pi/l10,sigma)
   call dtheta_1_0(dtheta)
   dtheta =  dtheta*pi/l10
   sigma = sigma/dtheta

   u = (z0**2)*conjg(l10)/(2*l10)
   if (shift) then
      u = u  +   conjg(lshift)*(z + z0)/2
   endif
   u  = u*pi/area
   sigma = sigma*exp(u)
   if(sign == -1) sigma = -sigma

   return
contains
  subroutine theta_1(z, theta)
    implicit none
    complex (kind=dp), intent(in) :: z
    complex (kind=dp), intent(out) :: theta
    ! evaluate elliptic theta function theta_1(z|tau) 
    ! theta_1 =  2*sum_{k>=0} (-1)**k (exp(i*pi*tau*(k + 1/2)**2) * sin( (2*k+1)*z)
    integer :: c,m
    complex (kind=dp) :: theta_prev, twoz,zz,qq,q2n
    twoz = 2*z
    zz = z
    qq = q_4
    q2n = cmplx(1,kind=dp)
    c = 2
    theta = c*q_4*sin(zz)
    do
       zz = zz + twoz
       c = -c
       q2n = q2 * q2n
       qq = qq *q2n
       theta_prev = theta
       theta  =  theta  + (c * qq * sin(zz)) 
       if(theta == theta_prev) exit
    enddo
    return
  end subroutine theta_1  
  subroutine dtheta_1_0( dtheta)
    complex (kind=dp), intent(out) :: dtheta
    ! evaluate derivative of elliptic theta function dtheta_1(z|tau)/dz at z = 0 
    integer :: c,m
    complex (kind=dp) :: dtheta_prev,qq,q2n
    m = 1
    c = 2
    qq = q_4
    q2n = cmplx(1,kind=dp)
    dtheta = 2 * qq
    do
       m = m + 2
       c = -c
       q2n = q2n * q2
       qq = qq *q2n
       dtheta_prev = dtheta
       dtheta  =  dtheta  + (m * c  * qq) 
       if(dtheta == dtheta_prev) exit
    enddo
    return
  end subroutine dtheta_1_0
end subroutine weierstrass_sigma



subroutine jacobi_theta(n,z, tau,theta,sum)
  implicit none
  integer, parameter :: dp = kind(1.0d0)
  integer, intent(in) :: n
  complex (kind=dp), intent(in) :: z, tau
  logical, intent(in) :: sum
  complex (kind=dp), intent(out) :: theta
  !----------------------------------------------------------
  ! evaluate elliptic theta function theta_n(z|tau), n = 1,2,3,or 4 
  ! theta_1 =  2*sum_{k>=0} (-1)**k (exp(i*pi*tau*(k + 1/2)**2)) * sin( (2*k+1)*z)
  ! theta_2 =  2*sum_{k>=0}  (exp(i*pi*tau*(k + 1/2)**2)) * cos( (2*k+1)*z)
  ! theta_3 =  1+  2*sum_{k>=0} (exp(i*pi*tau*(k**2)) * cos (2*k*z)    
  ! theta_1 =  2*sum_{k>=0} (-1)**k (exp(i*pi*tau*(k**2) )* cos(2*k*z)
  !----------------------------------------------------------
  
  integer :: c,m
  complex (kind=dp) ::  theta_prev,q_4,q_2,q,q2,qq,twoz,zz,q2n,b,factor
  complex (kind=dp),  parameter :: one = (1_dp,0_dp)
  real (kind=dp) :: pi, a

  a = 1_dp
  a = asin(a)/2
  q_4 = exp( cmplx(0,a,kind=dp)*tau)   ! q**(1/4)
  q = (q_4)**4
  q2 = q**2
  twoz = 2*z
  c = 2

  select case (n)
  case (1)
     theta = 2*q_4*sin(z)
     if (sum) then
        q2n = one
        qq = q_4
        zz = z
        do
           q2n = q2n * q2
           qq = qq * q2n    
           zz = zz + twoz
           c = -c
           theta_prev = theta
           theta  =  theta  + (c * qq  * sin(zz)) 
           if(theta == theta_prev) exit
        enddo
     else
        b = -2*cos(twoz)
        q2n = one
        do
           q2n = q2n * q2
           factor = (one - q2n) * (one  + q2n*(b + q2n))
           if(factor == one) exit
           theta = theta*factor
        enddo
     endif
  case (2)
     theta = 2*q_4*cos(z)
     if (sum) then
        q2n = one
        qq = q_4
        zz = z
        do
           q2n = q2n * q2
           qq = qq * q2n    
           zz = zz + twoz
           theta_prev = theta
           theta  =  theta  + (c * qq  * cos(zz)) 
           if(theta == theta_prev) exit
        enddo
     else
        b = 2*cos(twoz)
        q2n = one
        do
           q2n = q2n * q2
           factor = (one - q2n) * (one  + q2n*(b + q2n))
           if(factor == one) exit
           theta = theta*factor
        enddo
     endif
  case (3)
     theta =  one
     if (sum) then
        zz = cmplx(0,kind=dp)
        qq = one
        q2n = one
        do 
           zz = zz + twoz
           qq = qq * q * q2n
           q2n = q2n * q2
           theta_prev = theta
           theta  =  theta  + (c * qq * cos(zz)) 
           if(theta == theta_prev) exit
        enddo
     else
        b = 2*cos(twoz)
        qq = q
        q2n = one
        do
           factor = (one + qq*(b + qq))*(one - q*qq)
           if(factor == one) exit
           theta = theta*factor
           qq = qq * q2
        enddo
     endif
  case (4)
     theta = one
     if (sum) then
        zz = cmplx(0,kind=dp)
        qq = one
        q2n = one
        do 
           zz = zz + twoz
           c = -c
           qq = qq * q * q2n
           q2n = q2n * q2
           theta_prev = theta
           theta  =  theta  + (c * qq * cos(zz)) 
           if(theta == theta_prev) exit
        enddo
 




    else
        b = -2*cos(twoz)
        qq = q
        do
           factor = (one + (qq * (b + qq)) ) * (one - q*qq)
           if(factor == one) exit
           theta = theta * factor
           qq = qq * q2
        enddo
     endif
  case default
     write(6,'("error: Jacobi theta must have n = 1,2,3, or 4; n  was ",i10)') n
     stop
  end select
  return
end subroutine jacobi_theta

subroutine reduce_periods(l1,l2,m,n)
  implicit none
  integer, parameter :: dp = kind(1.0d0)
  complex (kind=dp), intent(in) :: l1,l2
  integer, intent(out) :: m(2), n(2)

  !  l1_new = m(1)*l1 + n(1)*l2
  !  l2_new = m(2)*l1 + n(2)*l2
  !  m(1)*n(2) - n(1)*m(2) = 1
  !  l1_new is shortest period
  !  l2_new  is next shortest period

  real (kind=dp) :: a, b, c, aa, bb, cc, dd
  integer::  k, mm, nn
  logical :: changed

  a = real(l1*conjg(l1))
  b = real(l1*conjg(l2))
  c = real(l2*conjg(l2))
  
  m(1) = 1
  n(1) = 0

  m(2) = 0
  n(2) = 1
  
  aa = a
  bb = b
  cc = c

  changed = .true.
  do while (changed)
     changed = .false. 
     k = nint(bb/cc)
     if (k /= 0) then    ! reduce l1 if possible
        mm = m(1) - k*m(2)
        nn = n(1) - k*n(2)
        dd = mm*mm*a + 2*mm*nn*b + nn*nn*c
        if(dd < aa) then
           changed = .true.
           m(1) = mm
           n(1) = nn
           aa = dd
           bb = m(1)*m(2)*a + (m(1)*n(2) + n(1)*m(2))*b + n(1)*n(2)*c
        endif
     endif
     k = nint(bb/aa)
     if (k /= 0) then  ! reduce l2 if possible
        mm = m(2) - k*m(1)
        nn = n(2) - k*n(1)
        dd = mm*mm*a + 2*mm*nn*b + nn*nn*c
        if (dd < cc) then
           changed = .true.
           m(2) = mm
           n(2) = nn
           bb = m(1)*m(2)*a + (m(1)*n(2) + n(1)*m(2))*b + n(1)*n(2)*c
           cc = dd
        endif
     endif
     if (cc < aa) then   ! swap l1, l2
        k = m(1)
        m(1) = m(2)
        m(2) = -k
        k = n(1)
        n(1) = n(2)
        n(2) = -k
        dd  = aa
        aa = cc
        cc = dd
        bb = -bb
     endif
  enddo
  return
end subroutine reduce_periods


subroutine landau_basis_state(norb,m,z,l1,l2,wf)
implicit none
integer, parameter :: dp= kind(1.0d0)
integer, intent(in) :: norb, m
complex (kind=dp), intent(in) :: z,l1,l2
complex (kind=dp), intent(out) :: wf


real (kind=dp) :: area, pi,scale2,scale
integer :: j
complex (kind=dp):: u,z0,l10,l20,w,beta, sigma, alpha

write(6,'(2f12.5)') l1,l2

pi = 1_dp
pi = 2*asin(pi)


area  = aimag(conjg(l1)*l2)

scale2 = (pi*norb)/area    ! unite with 2 * (ell_B)**2 = 1
scale = dsqrt(scale2)

z0 = scale*z
l10 = scale*l1
l20 = scale*l2


beta = m*l20/norb
alpha = l10/(2*norb)
u = conjg(beta)*(z0 - beta/2)
u = u - (z0*conjg(z0)/2)
wf = exp(u)

write(6,'(2f12.5)') l10,l20

do j = -(norb-1),(norb-1),2
   w = j*alpha  + beta
!   call weierstrass_sigma(z0 - w,l10,l20,sigma)
   wf = wf*sigma
enddo
if(mod(m,2)/=0) wf = -wf
return

end subroutine landau_basis_state

