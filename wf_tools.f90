! wf_tools.f90
! (c) F. Duncan M. Haldane, March 2016
! haldane@princeton.edu
!   v0.1   2016-03-08


! tools for lattice representation of guiding-center
! states on the torus for fractional qunatum all calculations


! contents:
!   module laughlin_m
!   subroutine get_m_laughlin(m)
!   subroutine setup_laughlin_state(nel,m,sl2z,k)
!   subroutine get_laughlin_cm (x,wf)
!   subroutine get_flux_factor(norb,x_center,exclude,nel,x,power,flux_factor)
!   subroutine get_vandermonde_to_power(norb,nel,x,power,vandermonde_power)
!   subroutine create_one_particle_landau_basis(norb,sl2z,basis)
!   function unnormalized_laughlin_wf(m,nel,x) result(wf)
!   subroutine make_lattice_landau_state(norb,k,sl2z,landau,nzeroes,k1,k2)
!   subroutine get_next_occupation_state(nel,norb,fermion,full,k,config,restart,done)
!   subroutine set_geometry(norb,metric)

!==================================
module laughlin_m
  integer, parameter :: dp = kind(1.0d0)
  complex (kind=dp), allocatable :: laughlin_cm(:,:),omega(:)
  integer :: nbar,cm_zeroes,list, k_cm, sl2z_cm(2,2),norb = 0
end module laughlin_m



function lattice_phase_factor(norb,m,n) result(factor)
  implicit none
  integer, parameter :: dp = kind(1.0d0)
  integer, intent(in) :: norb,m,n
  complex(kind=dp) :: factor
  
  real (kind=dp), save :: pi = 1.0_dp
  real (kind=dp) :: pin
  integer :: i
  integer, save :: norb0 = 0, list
  complex (kind=dP), allocatable, save :: omega(:)
  
  
  if(norb /= norb0) then
     if (allocated(omega)) deallocate(omega)
     list = 2*norb
     allocate (omega(0:list-1))     
     if(pi == 1_dp) pi = 2*asin(pi)
     pin = pi /norb
     forall ( i = 0:list-1) omega(i) = cmplx(cos(i*pin),sin(i*pin),kind=dp)
     norb0 = norb
  endif
  
  factor = omega(modulo(m*n,list))
  if(mod(m + n,2) /= 0) factor = -factor

  return
end function lattice_phase_factor
  

!=================================
subroutine get_m_laughlin(m)
  use laughlin_m
  integer, intent(out) :: m
  m = cm_zeroes
  return
end subroutine get_m_laughlin



!=================================
subroutine setup_laughlin_state(nel,m,sl2z,k)
  use laughlin_m
  implicit none
  integer, intent(in) :: nel,m,sl2z(2,2),k

  real (kind=dp) :: pin
  integer :: i,j

  write(6,'("set up Laughlin state with  m = ",i5)') m
  norb  = m*nel
  nbar = nel
  cm_zeroes = m
  list =  2*nbar
  k_cm = k
  sl2z_cm = sl2z

  if(allocated(laughlin_cm)) deallocate(laughlin_cm)
  if(allocated(omega)) deallocate(omega)

  allocate (laughlin_cm(norb,norb),omega(0:list-1))
  

  call make_lattice_landau_state(norb,k,sl2z,laughlin_cm,m,0,0)


  pin = 1_dp
  pin = 2*asin(pin)/nbar
  do i = 0,list-1
     omega(i) = cmplx(cos(i*pin),sin(i*pin),kind=dp)
  enddo
  
  return
end subroutine setup_laughlin_state

!=================================
subroutine get_laughlin_cm (x,wf)
  use laughlin_m
  implicit none
  integer, intent (in) :: x(2)
  complex (kind=dp), intent(out) :: wf

  integer :: m,n,j,k

  
  m = 1 + modulo(x(1)-1,norb)
  n = 1 + modulo(x(2)-1,norb)

  wf = laughlin_cm(m,n)

  j = (x(1)-m)/norb
  k = (x(2)-n)/norb
  wf = wf*omega(modulo(j*x(2)-k*x(1),list))
  
  if(mod(cm_zeroes,2) /= 0) then
     if(.not. (mod(j,2) == 0 .and.  mod(k,2) == 0)) wf = -wf
  endif


  return
end subroutine get_laughlin_cm



!=================================
subroutine get_flux_factor(norb,x_center,exclude,nel,x,power,flux_factor)
  implicit none
  integer, parameter :: dp = kind(1.0d0)
  integer, intent(in) :: norb, nel, power
  integer, intent(in) :: x_center(2),x(2,nel)
  logical, intent(in) :: exclude(nel)
  complex(kind=dp),intent(out) :: flux_factor
  !-----------------------------------------
  ! returns the "flux attachment factor"
  !
  !   prod_{j(not excluded} Z(x_center  - x(j))**power  
  !
  !   the flux is centered at x_center, and particles excluded from the product
  !    can be positioned  relative to x_center
  !  
  !   if exclude(i) = .true. the product omits particle i
  !------------------------------------
  integer :: i,m,n, k
  complex (kind=dp) :: logfactor, logz
  real(kind=dp) :: normalization
  
  do i = 1,nel
     if(exclude(i)) cycle
     if(modulo(x(1,i) - x_center(1),norb) /= 0) cycle
     if(modulo(x(2,i)-x_center(2),norb) == 0) then
        flux_factor = cmplx(0,kind=dp)
        return
     endif
  enddo

  logfactor = cmplx(0,kind=dp)
  do i = 1,nel
     if(exclude(i)) cycle
     m =  x_center(1) - x(1,i)
     n =  x_center(2) - x(2,i)
     call get_log_lattice_z(norb,m,n,logz,normalization)
     if(normalization == 0_dp) then
        write(6,'("ERROR in GET_FLUX_FACTOR")')
        write(6,'("X_center = ",2i5)') x_center(:)
        do k = 1,nel
           if(exclude(k)) cycle
           write(6,'("X(",i5,") = ",2i5)') k, x(:,k)
        enddo
        write(6,'("anomaly for particle  i  =",i5)') i
        stop
     endif
     logfactor = logfactor + logz
  enddo
  flux_factor = exp (power*logfactor)      
  return
end subroutine get_flux_factor

!=================================
subroutine get_vandermonde_to_power(norb,nel,x,power,vandermonde_power)
implicit none
integer, parameter :: dp = kind(1.0d0)
integer, intent(in) :: norb,nel,power
integer, intent(in) :: x(2,nel)
complex (kind=dp), intent (out) :: vandermonde_power
!-----------------------------
!
!for a set of nel particles on distict  sites  x(i) = (x(1,i),x(2,i))
! so that for i /= j
!  (modulo(x(1,i)-x(1,j),norb) == 0 .and. x(2,i)-x(2,j) == 0) is
!   .false. 
!
!  this returns  vandermonde_power =  prod_{i <j}  Z(x(j) - x(j))**power
! m is any integer, positive, negative or zero.
!
!  cmplx(0) is returned if any particles are on the same site
!-----------------------------------
integer ::  i,j,m,n,k
complex (kind=dp) :: logfactor, logz
real(kind=dp) :: normalization


vandermonde_power = cmplx(0,kind=dp)
do j = 1,nel
   do i =  1,j-1
      if(modulo(x(1,i)-x(1,j),norb) /= 0) cycle
      if(modulo(x(2,i)-x(2,j),norb) == 0) return
   enddo
enddo
logfactor = cmplx(0,kind=dp)
do j = 1,nel
   do  i = 1, j-1
      m = x(1,i) - x(1,j)
      n = x(2,i) - x(2,j)
      call get_log_lattice_z(norb,m,n,logz,normalization)
      if(normalization == 0_dp) then
         write(6,'("ERROR in GET_VANDERMONDE_TO_POWER")')
         do k = 1,nel
            write(6,'("X(",i5,") = ",2i5)') k, x(:,k)
         enddo
         write(6,'("anomaly for  pair  i,j =",2i5)') i,j
         stop
      endif
      logfactor = logfactor + logz
   enddo
enddo
vandermonde_power = exp (power*logfactor)

return
end subroutine get_vandermonde_to_power

!=================================

subroutine make_quantum_distance_table(norb,s2)
implicit none
integer,  parameter :: dp = kind(1.0d0)
integer, intent(in) :: norb
real(kind=dp), intent(out) :: s2(norb,norb)
!----------------------
! computes important one-particle lattice properties
! that derive from the choice of metric
!
!  s2(m,n) = |<psi(m,n)|psi0,0)||**2
!  where |psi(m,n)> is the coherent state on lattice site
! x =(m,n)
!
!--------------------------------
integer :: sl2z(2,2), i,j

complex (kind=dp), allocatable :: basis(:,:,:)

!-----------------------------------------
   
! the result is modular invariant, choice of sl2z doesnt matter
sl2z(1,1) = 1
sl2z(1,2) = 0
sl2z(2,1) = 0
sl2z(2,2) = 1
allocate (basis(norb,norb,norb))
call create_one_particle_landau_basis(norb,sl2z,basis)
forall ( i = 1:norb, j=1:norb) s2(i,j) = &  
     & abs(sum( conjg(basis(norb,norb,:)) * basis(i,j,:)  ))**2
s2(norb,norb) = 1

return
end subroutine make_quantum_distance_table


subroutine make_gaussian_form_factor(norb,u,norm)
  implicit none
  integer, parameter :: dp = kind(1.0d0)
  integer, intent (in) :: norb

  real(kind=dp), intent(out) :: u(norb,norb), norm
  
  !--------------------------------------------------
  ! U(q(i,j)) =   sum_{m,n}  exp ((-0.5*|(q(i,j) + Norb*q(m,n))|**2) ) * 
  !             (-1)**(i*n) * (-1)**(j*m) + eta(m,n)**norb
  !
  ! eta(m,n) = (-1)**(m + n + m*n).
  ! 
  !   exp (i (q(i,j) + Norb*q(m,n)) . R) =  (eta(m,n)**norb) *(-1)*(m*j)* (-1)**(n*i)
  !    times exp (i q(i,j) . R)
  !-----------------------------------------------------
 
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
end subroutine make_gaussian_form_factor


subroutine get_omega(norb,omega)
  implicit none
  integer, parameter :: dp = kind(1.0d0)
  integer, intent(in) :: norb
  complex(kind=dp), intent(out) :: omega(0:norb-1)
!-----------------------------------------
!
!  produce and save the table omega(i = 0:norb-1)) = exp (2*pi*i/norb)
!     
!-----------------------------------------

  real (kind=dp) :: theta, two_pi
  integer :: i

  two_pi  = 4*asin(real(1,kind=dp))
  theta = two_pi/norb
  omega(0) = real(1,kind=dp)
  if (mod(norb,2) == 0) omega(norb/2) = real(-1,kind=dp)
  forall (i = 1:(norb-1)/2)
     omega(i) = cmplx(cos(i*theta),sin(i*theta), kind=dp)
     omega(norb-i) = conjg(omega(i))
  end forall


  return
end subroutine get_omega




subroutine create_one_particle_landau_basis(norb,sl2z,basis)
  implicit none
  integer, parameter :: dp = kind(1.0d0)
  integer, intent(in) :: norb, sl2z(2,2)
  complex(kind=dp), intent(out) :: basis(norb,norb,norb)

  
  integer, parameter :: bigint = 100000000
  integer:: k1,k2
  complex (kind=dp) :: s12
  logical :: test
  
  if(sl2z(1,1)*sl2z(2,2) - sl2z(1,2)*sl2z(2,1) /= 1) then
     write(6,'("error: CREATE_ONE_PARTICLE_LANDAU_BASIS:")')
     write(6,'(" called with invalid SL(2,Z) matrix:")')
     write(6,'(2i5)') sl2z(1,:), sl2z(2,:)
     stop
  endif
  
  
  do k1 = 1,norb
     call make_lattice_landau_state(norb,k1,sl2z,basis(:,:,k1),0,0,0)
  enddo
  
  test = .true.
  do k1 = 1,norb
     do k2 = 1,norb
        s12 = sum(conjg(basis(:,:,k1))*basis(:,:,k2))/norb
        if(k1 /= k2 .and. nint(bigint*abs(s12)) /= 0) then
           test = .false.
           write(6,'(2i5,2f25.16)') k1,k2,s12
        else if (k1 == k2 .and. nint(bigint*abs(s12)) /= bigint) then
           test = .false.
           write(6,'(2i5,2f25.16)') k1,k2,s12
        endif
     enddo
  enddo
  if (test) return
  
  write(6,'("orthonormality test of one-particle basis failed")')
  stop
end subroutine create_one_particle_landau_basis

!=================================

function unnormalized_laughlin_wf(m,nel,x) result(wf)
  implicit none
  integer, parameter :: dp = kind(1.0d0)
  integer, intent(in) :: m, nel
  integer, intent(in) :: x(2,nel)
  complex (kind=dp) :: wf

  complex (kind=dp) :: vandermonde
  integer :: x_cm(2), i, norb

  norb = m*nel

  forall (i = 1:2) x_cm(i) = sum(x(i,:))
  call get_laughlin_cm (x_cm,wf)
  if(abs(wf) == 0_dp) return

  call get_vandermonde_to_power(norb,nel,x,m,vandermonde)
  wf = wf * vandermonde

  return
end function unnormalized_laughlin_wf

!=================================
subroutine make_lattice_landau_state(norb,k,sl2z,landau,nzeroes,k1,k2)
  implicit none
  integer, parameter :: dp = kind(1.0d0)
  integer, intent(in) :: norb,  k, sl2z(2,2)
  integer, intent(in):: nzeroes, k1,k2
  complex(kind=dp) , intent(out) :: landau(norb,norb)
  
  ! when called with NZEROES <= 0, or NZEROES = NORB,  K1 and K2
  ! are ignored, and this provides the Landau basis state K of the 
  ! basis with 
  !  L1 = sl2z(1,1)*l(1)  + sl2z(1,2)*l(2)
  !  L2 = sl2z(2,1)*l(1)  + sl2z(2,2)*l(2)
  ! evaluated on the lattice
  !   z(m,n) = (m*l(1) + n*l(2))/norb,    m,n = 1:norb
  !
  !  otherwise, NZEROES must be a positive factor of NORB, 
  !  this can be used for center-of-mass factors of Laughlin or CFL
  !  states at fliing 1/NZEROES.   For Laughlin states, K1 = K2 = 0,
  !  while for CFL states, they fix the many-body translational quantum number
  !
  !   The sum of the NZEROES cm_wavefunction zeroes is   (k1*l(1) + k2*l(2))/nbar, modulo L,
  !   where nbar * NZEROES = NORB.
  !
  !   here the basis  L(1:2) for the lattice is taken from the 
  !   values set set by the  last call to  SET_GEOMETRY
  !
  !  for NZEROES = NORB:
  !   <i,j| k> = landau_basis(i,j,k)  is state k of the Landau basis
  !   defined by t1 = -t(L1), t2 = -t(L(2) )
  !   where t1|k} =  exp 2 pi i k/norb)|k>
  !       t2|k> = |k+1>
  !       |k+Norb> = |k>
  !  and L1 x L2 = 1
  !  
  !
  
  integer :: m,n,i,j,p,list, det, count, nflux, nbar
  complex (kind=dp) :: factor, l(2), l1_new, l2_new,zij, z, w, tau,new_tau,zz
  integer, parameter :: bigint = 10000
  real (kind=dp) :: area, pin, xij, yij, x, y 
  integer :: sl2z0(2,2),norb_z
  logical :: rationalize = .true.
  real (kind=dp), allocatable ::  zeroes(:,:)
  complex (kind=dp), allocatable ::  omega(:)
  real (kind=dp), save :: pi = 1_dp

  if(pi == 1_dp) pi = 2*asin(pi)

! get l(1:2) that is stored in module z_function_m
  call get_l(norb_z,l(1),l(2))

  area = aimag(conjg(l(1))*l(2))/pi
  if(nint(bigint*area)/= bigint*norb) then
     write(6,'("MAKE_LATTICE_LANDAU_STATE: area mismatch:",i5,f25.16)') norb,area
     stop
  endif
  
  det = sl2z(1,1)*sl2z(2,2) - sl2z(1,2)*sl2z(2,1)
  if(det /= 1) then
     write(6,'("invalid SL(2,z) matrix, det /= 1")')
     write(6,'(2i5)') sl2z(1,:)
     write(6,'(2i5)') sl2z(2,:)
     stop
  endif
  
  if (nzeroes > 0)  then
     if(mod(norb,nzeroes)/= 0) then
        write(6,'("MAKE_LANDAU__BASIS_STATE: 0 <  NZEROES  = ",i5)') nzeroes
        write(6,'("must be a factor of NORB = ",i5)') norb
        stop
     endif
     nflux = nzeroes
     nbar = norb/nflux
  else   
     nflux = norb
     nbar = 1
  endif
  
  
  list = 2*norb*nbar
  allocate (omega(0:list-1))
  pin = pi/(norb*nbar)
  do i = 0, list - 1
     omega(i) = cmplx(cos(i*pin),sin(i*pin),kind=dp)
  enddo
  


  
  allocate (zeroes(2,nflux))
  count = 0
  m = modulo(k1,nbar) +  nbar*k*sl2z(2,1) 
  n = modulo(k2,nbar) +  nbar*k*sl2z(2,2)
  do i = -(nflux-1),(nflux-1),2
     count = count + 1
     zeroes(1,count) =  real(2*m  + i*nbar*sl2z(1,1),kind=dp)/(2*norb) 
     zeroes(2,count) =  real(2*n  + i*nbar*sl2z(1,2),kind=dp)/(2*norb) 
  enddo
  
  
  tau = l(2)/l(1)
  call optimize_tau(tau,new_tau,sl2z0)
  do i = 1,norb
     do j= 1,norb
        landau(i,j) = omega(modulo(m*j-n*i,list))      
        if(mod(k,2)/=0) landau(i,j)  = -landau(i,j)
        xij = real(i,kind=dp)/norb
        yij = real(j,kind=dp)/norb
        do p = 1,nflux
           x  = xij - zeroes(1,p)
           y  = yij - zeroes(2,p)
 !          call z_function_with_modular_transform(x,y,l(1),l(2),&
 !               &rationalize,2*norb,zz,sl2z0)
           call z_function(x,y,l(1),l(2),rationalize,2*norb,zz)
           landau(i,j) = landau(i,j) *  zz
        enddo
     enddo
  enddo
!  landau = landau/sqrt(real(sum(conjg(landau)*landau))/norb )
  
  deallocate (zeroes, omega)
  return
end subroutine make_lattice_landau_state



!==================================
subroutine set_geometry(norb,metric)
  implicit none
  integer, parameter :: dp = kind(1.0d0)
  real(kind=dp), intent(inout) :: metric(2,2)
  integer, intent(in) :: norb

!---------------------------------------------------------
!
!   D= sqrt(det metric)
! 0.5 * ( METRIC(1,1)              METRIC(1,2) + CMPLX(0,D) )
!       ( METRIC(2,1) -CMPLX(0,D)  METRIC(2,2)              )
!
!     D* (conjg(omega(1))*omega(1)  conjg(omega(1))*omega(2) )
!        (conjg(omega(2))*omega(1)  conjg(omega(2))*omega(2) )
!
!---------------------------------------------------------
!
! metric is  proportional to      L1.L1    L1.L2  
!                                 L2.L1    L2.L2
!
!   Lattice distances are  propootional to
!   
!      metric(1,1) m*m + 2*metric(1,2)*m*n + metric(2,2)*n*n 
!   
!   The aspect ratio/angle model sets
!   L1 =   C (1,0)    L2 =  C(aspct*cosa, aspct*sina) with sina = sqrt(1-cosa**2)
!   and L1 x L2 = pi*norb  (units where  magnetic area = pi)
!
!   The unimodular metric is then      (1/sina)  x     1/aspct       cosa
!                                                        cosa         aspct   
!
!      
!   this routine rescales any valid input metric to unimodular form
!
!-------------------------------------------------------

  complex (kind=dp) :: omega(2),l(2)
  real(kind=dp) :: scaled_metric(2,2), epsilon(2,2), det, area,pi,l1,l2 , &
      & scale,  cosa,aspct
  logical :: test

  integer :: i,j, m(2),n(2)

  epsilon(1,1) = 0_dp
  epsilon(1,2) = 1_dp
  epsilon(2,1) = -1_dp
  epsilon(2,2) = 0_dp
  
  test = .true.
  if(metric(1,1) <= 0_dp) test= .false.
  if(metric(2,2) <= 0_dp) test= .false.
  if(metric(1,2) /= metric(2,1)) test= .false.
  det = metric(1,1)*metric(2,2) -  metric(1,2)*metric(2,1) 
  if(det <= 0_dp) test = .false.

  if(.not. test) then
     write(6,'("creating new metric:")')
     write (6,'("give COSINE = L1.L2/|L1||L1} and ASPECT ratio |L2|/|L1|")')
     write(6,'(" The absolute value of COSINE  MUST be less than 1")')
     read(5,*) cosa,  aspct
     if (abs(cosa) >= 1_dp .or. aspct <= 0_dp ) then
        write(6,'("ïnvalid cosine or aspect ratio:")'), cosa, aspct
        stop
     endif
     metric(1,1) = 1/aspct
     metric(1,2) = cosa
     metric(2,1) = cosa
     metric(2,2) = aspct
     det = 1_dp - cosa**2 
  endif

  metric = metric/sqrt(det)

!  now produce the factorization
!
!  (1/2)  *    cmplx( metric(1,1))  cmplx(metric(1,2), 1)  
!              cmplx(metric(2,1) -1)   cmplx(metric(2,2))
!
!          =  conjg(omega(1)) * omega(1)        conjg(omega(1)) *  omega(2) 
!             conjg(omega(2)) * omega(1)       conjg(omega(2)) *  omega(2) 
! 



  omega(1) = cmplx(metric(1,1),kind=dp)

  omega(2) = cmplx(metric(1,2), 1,kind=dp)

  omega = omega/sqrt(2*metric(1,1))

 



  write(6,'("resolution of metric")')
  do i = 1,2
     do j = 1,2
        write(6,'(2i5,4f25.12)') i,j, conjg(omega(i))*omega(j), &
           &  cmplx(metric(i,j),epsilon(i,j),kind=dp)/2
     enddo
  enddo
  
  pi = 1_dp
  pi = 2*asin(pi)
  area = pi*norb
  scale = sqrt(area)
  
  l = sqrt(2*area)*omega

! place geometry  in z_function_m
  call set_l(norb,l(1),l(2))

 
  write(6,'("L1: (",f25.16,",",f25.16,")")') l(1)
  write(6,'("L2: (",f25.16,",",f25.16,")")') l(2)

  area = aimag(conjg(l(1))*l(2))/pi

!  write(6,'(" test: norb = ",i5,f25.16)') norb, area
  return
end subroutine set_geometry

