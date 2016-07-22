program single_sample
  implicit none
  integer, parameter :: dp = kind(1.0d0)

  integer :: sl2z(2,2), nel, m_laughlin, norb
  logical :: fermionlike, changed
  
  real(kind=dp) :: metric (2,2)

  complex (kind=dp) :: wf
  
  real (kind=dp) :: weight, test
  
  integer (kind=8) :: n_reject, n_accept, steps, history, step
  integer :: units, nswap, strip_width, unit


  integer, allocatable :: x(:,:)
  logical, allocatable :: swap_mask(:,:)

  integer :: seed
  real::dd, r4_uni


  integer ::  i, j, k, m, n
  
  integer, allocatable :: strip_count(:,:)
  real(kind=dp), allocatable ::strip_data(:,:,:,:), strip(:,:,:)
  real(kind=dp), allocatable ::energy_data(:,:), running_energy(:,:)
  real(kind=dp) :: energy(2), uc, variance, step_cutoff

  real (kind=dp), allocatable :: q_distance(:,:), indexd(:,:), sq(:,:), sq_data(:,:,:)


  logical, allocatable :: k_space_mask(:,:)



!  write(6,'(" nu = 1/m Laughlin : give nel, m ")')
!  read (5,*) nel, m_laughlin
  nel = 20
!  nel = 4
  m_laughlin = 3
  norb = nel*m_laughlin
  fermionlike = .true.
  
  
  allocate(swap_mask(1,1))

  allocate (k_space_mask(norb,norb))
  k_space_mask = .true.
!  do i = -norb/2,norb/2
!!     do j = -norb/2,norb/2
!        if(i**2 > 4*norb .or. j**2 > 4*norb) cycle
!        k_space_mask(1 + modulo(i-1,norb),1 + modulo(j-1,norb)) = .true.
!     enddo
!  enddo
 


  ! choose the holomorphic geometry
  ! if an invalid metric is used as input,
  ! SET_GEOMETRY will request cosa and aspct
  !  (|L2|/L1 = aspct, L1 .L2 = |L1|L2| cos a
  ! in either case, the metric will be returned
  ! in unimodular (det = 1) form
  metric = 0_dp
  call set_geometry(norb,metric)     !from wf_tools.f90
  
  ! set up Weierstrass for the vandermonde
  ! z_function is the Weierstrass sigma function
  ! times a gaussian.  Its absolute value is periodic
  call setup_z_function_table        !from z_function_m.f90
  
  ! select center-of-mass quantum numbers for Laughlin state  
  sl2z(1,1) = 1
  sl2z(1,2) = 0
  sl2z(2,1) = 0
  sl2z(2,2) = 1
  k = 0
  call setup_laughlin_state(nel,m_laughlin,sl2z,k)  !from wf_tools.f90
  write(6,'("Laughlin state on lattice with nel = ",i5)') nel
  call coulomb_setup


! chose the step length
  step_cutoff =  0.9d0
  call step_table(step_cutoff)
  


  write(6,'("give number of data units, and mc steps per unit")')
  read(5,*) units,steps
  history = steps
  allocate(running_energy(2,0:history-1))
  allocate (energy_data(2,units),sq(norb,norb),sq_data(norb,norb,units))
  allocate (strip_data(0:nel,norb,norb,units),strip(0:nel,norb,norb),strip_count(norb,norb))

  running_energy = 0.0_dp
  energy_data = 0_dp
  strip_data = 0_dp
  sq_data = 0_dp
  
! lattice coordinates for the  system
  allocate (x(2,nel))






  write(6,'("give seed for random number generator")')
  read(5,*) seed

  

  !  initialize the state
  call initialize(nel,norb,fermionlike,x,seed,.false.,nswap,swap_mask)

  call coulomb_energy(nel,norb,x,uc)

  running_energy(1,:) = uc/nel
  running_energy(2,:) = (uc/nel)**2
  call get_strip_count(nel,norb,x,strip_count)
  call get_sq(nel,norb,x,k_space_mask,sq)
  call get_wf(nel,x,wf)
  weight = abs(conjg(wf)*wf)

! now start the walkers
  do unit = 1,units
     n_reject = 0
     n_accept = 0
     do step = 1,steps
        !$OMP PARALLEL !SHARED(NUMTHREADS), PRIVATE(OMP_THREAD_NUM)
!$OMP SECTIONS
        !$OMP SECTION
        call metropolis_step(nel,norb,fermionlike,x,seed,weight,wf,changed)
        if(.not. changed) then
           n_reject = n_reject + 1
        else
           n_accept = n_accept + 1
           call coulomb_energy(nel,norb,x,uc)
           call get_strip_count(nel,norb,x,strip_count)
           call get_sq(nel,norb,x,k_space_mask,sq)
        endif
     
        forall (i = 1:norb, j = 1:norb, k_space_mask(i,j))  &
             &  sq_data(i,j,unit) = sq_data(i,j,unit) + sq(i,j)

        energy_data(1,unit) = energy_data(1,unit) + (uc/nel)
        energy_data(2,unit) = energy_data(2,unit) + (uc/nel)**2
        forall (i = 1:norb, j= 1:norb, i /= j) & 
             & strip_data(strip_count(i,j),i,j,unit) = &
             &  strip_data(strip_count(i,j),i,j,unit) + 1
     enddo
!$OMP END SECTIONS
!$OMP END PARALLEL
     energy_data(:,unit)  = energy_data(:,unit)/steps
     variance = energy_data(2,unit) - energy_data(1,unit)**2
     strip_data(:,:,:,unit)  = strip_data(:,:,:,unit)/steps
     sq_data(:,:,unit) = sq_data(:,:,unit)/stepS
     write(12,'(3i6,3i16)') unit,nel,norb,steps, n_accept, n_reject
      do i = 1,nel
        write(12,'(2i5)') x(1,i),x(2,i)
     enddo
     call print_sq(norb,sq_data(:,:,unit), k_space_mask)
     write(12,'(2f25.16)') energy_data(1,unit),energy_data(2,unit),variance
     do strip_width = 1,norb-1
        do n = 0,nel
           test = 0_dp
           do i = 1,norb
              test = test + strip_data(n,i,1 + modulo(i + strip_width -1,norb),unit)
           enddo
           test = test/norb
           if (test < 1.0e-7_dp) cycle
 
           write(16,'("strip width:",i5," n = ",i5)') strip_width, n
           write(16,'(10f10.7)') (strip_data(n,i,1 + modulo(i + strip_width -1,norb),unit), i=1,norb)
        enddo
     enddo
    write(6,'(3i6,i16," accepted =",i16," rejected =",i16)') unit,nel,norb,steps, n_accept, n_reject
    write(6,'("energy per particle:",f25.16, &
          &" variance:",f25.16)') energy_data(1,unit),variance
 enddo
  
  ! final summary
  energy(1) = sum(energy_data(1,:))/units
  energy(2) = sum(energy_data(2,:))/units
  forall (i = 1:norb, j = 1:norb, k_space_mask(i,j))  &
       &  sq(i,j) = sum(sq_data(i,j,:))/units
   variance = energy(2) - energy(1)**2
  forall (n = 0:nel,i = 1:norb, j=1:norb) strip(n,i,j) = sum(strip_data(n,i,j,:))/units
  write(13,'(3i6,3i16)') nel,norb,steps*units
  do i = 1,nel
     write(12,'(2i5)') x(1,i),x(2,i)
  enddo
  write(13,'(2f25.16)') energy(1),energy(2),variance
  do strip_width = 1,norb-1
     do n = 0,nel
           test = 0_dp
           do i = 1,norb
              test = test + strip(n,i,1 + modulo(i + strip_width -1,norb))
           enddo
           test = test/norb
           if (test < 1.0e-7_dp) cycle
           write(6,'("strip width:",i5," n = ",i5)') strip_width, n
        write(6,'(10f10.7)') (strip(n,i,1 + modulo(i + strip_width -1,norb)), i=1,norb)
     enddo
  enddo
  call print_sq(norb,sq, k_space_mask)
  write(6,'("energy per particle:",f25.16, &
       &" variance:",f25.16)') energy(1),variance
contains
     subroutine  print_sq(norb,sq, k_space_mask)
       integer :: norb
       real(kind=dp) :: sq(norb,norb)
       logical :: k_space_mask(norb,norb)
       integer :: i,j,i1,j1     
       do i = -(norb+1)/2, (norb+1)/2
          do j = -(norb+1)/2, (norb+1)/2
             i1 = 1 + modulo(i-1,norb)
             j1 = 1 + modulo(j-1,norb)
             if (.not. k_space_mask(i1,j1)) cycle
             write(6,'(2i5,f25.16)')   i,j, sq(i1,j1)
          enddo
       enddo
       return
     end subroutine print_sq
end program single_sample



subroutine get_wf(nel,x,wf)
  implicit none
  integer, parameter :: dp = kind(1.0d0)
  integer, intent(in) :: nel
  integer, intent(in) :: x(2,nel)
  complex(kind=dp), intent(out) :: wf

  complex (kind=dp) :: unnormalized_laughlin_wf
  integer :: m
  call get_m_laughlin(m)
  wf = unnormalized_laughlin_wf(m,nel,x) 
  return
end subroutine get_wf
  
subroutine get_random(r4_random,seed)
  implicit none
  integer, intent(inout) :: seed
  real, intent(out) :: r4_random
!------------------------------
! random default real number with a 
! uniform distribution in the range
! 0 <= uni < 1
! seed is a default integer
! r4_uni is a random number generator
! based on the ziggurat algorithm
!----------------------------
  real :: r4_uni
  r4_random = r4_uni(seed)
  return
end subroutine get_random
