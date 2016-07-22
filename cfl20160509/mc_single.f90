module sampling_m
  logical ::  do_strip = .false.
  logical ::  do_sq = .true. 
  logical ::  do_energy = .true.

contains
  subroutine sampling_setup(nel,norb)
    implicit none
    integer, intent(in) :: nel, norb
    
    write(6,'(1x)')

    if (do_strip) then
       write(6,'("setting up strips for entanglement calculation")')
       call strip_setup(nel,norb)
    end if
    
    if (do_sq) then
       write(6,'("setting up structure factor calculation")')
       call sq_setup(nel,norb)
    endif
    
    
    if (do_energy) then
       write(6,'("setting up coulomb energy calculation")')
       call energy_setup(nel,norb)
    endif
    return
  end subroutine sampling_setup

  subroutine  data_initialize(nel,norb,n_unit)
    implicit none
    integer, intent(in) :: nel, norb, n_unit

    if (do_strip) call strip_data_init(nel,norb,n_unit)
    if (do_sq) call sq_data_init(nel,norb,n_unit)
    if (do_energy) call energy_data_init(nel,norb,n_unit)
    return
  end subroutine data_initialize
  
  subroutine sample(nel,x,changed_x,unit)
    implicit none
    integer, intent(in) :: nel, x(2,nel), unit
    logical, intent(in) :: changed_x(nel)

    if (do_strip) call strip_sample(nel,x,changed_x,unit)
    if (do_sq)  call  sq_sample(nel,x, changed_x,unit)
    if (do_energy) call energy_sample(nel,x,changed_x,unit)
    return
  end subroutine sample

   subroutine data_combine(unit,n_steps)
     implicit none
     integer, intent(in) :: unit, n_steps

     if (do_strip) call strip_data_combine(unit,n_steps)     
     if (do_sq) call sq_data_combine(unit,n_steps)
     if (do_energy) call energy_data_combine(unit,n_steps)
     return
   end subroutine data_combine

   subroutine data_finalize(n_unit)
     implicit none
     integer, intent(in) :: n_unit

     if (do_strip) call strip_data_finalize(n_unit)
     if (do_sq) call sq_data_finalize(n_unit)
     if (do_energy) call energy_data_finalize(n_unit)
     return
   end subroutine data_finalize

end module sampling_m


program single_sample
  use sampling_m
  implicit none
  integer, parameter :: dp = kind(1.0d0)

  integer :: sl2z(2,2), nel, m_laughlin, norb
  logical :: fermionlike,  updated
  
  real(kind=dp) :: metric (2,2)
  logical, allocatable :: swap_mask(:,:), changed_x(:)
  

  
  integer :: n_reject, n_accept, n_steps, step, equilibration, ring
  integer :: n_unit, nswap,  unit, norb2,id, nsq_max, count
  integer :: seed
  integer, allocatable :: x(:,:), previous_x(:,:)
  complex (kind=dp) :: wf
  real (kind=dp) :: weight,prev, norm, mean_weight, max_weight
  integer :: current, history, sampling_interval, unchanged_x
  integer :: max_changes, training_period, target
  integer, allocatable :: update_list(:)

  integer ::  unchanged, changed
! quantum distance
  real (kind=dp), allocatable :: s2(:,:)

  real (kind=dp), allocatable :: weight_history(:)


  complex (kind=dp) :: l1,l2


  real (kind=dp) :: ratio
  integer :: nbreaks,breakpoint,stepcount
  logical :: success


  logical, allocatable :: bzb_list(:)







  real(kind=dp) ::  step_cutoff



!  real(kind=dp) :: small_ff

  integer ::  i, j, k, m, n, w, k1,k2
  integer :: i0,j0

  complex (kind=dp) :: q


  write(6,'(" nu = 1/m Laughlin : give nel, m ")')
  read (5,*) nel, m_laughlin
  m_laughlin = 3
  norb = nel*m_laughlin
  fermionlike = .true.
  
  
  ! choose the holomorphic geometry
  ! if an invalid metric is used as input,
  ! SET_GEOMETRY will request cosa and aspct
  !  (|L2|/L1 = aspct, L1 .L2 = |L1|L2| cos a
  ! in either case, the metric will be returned
  ! in unimodular (det = 1) form
  metric = 0_dp
  call set_geometry(norb,metric)     !from wf_tools.f90

  ! set up "modified Weierstrass sigma function"for the vandermonde
  ! z_function is the modified Weierstrass sigma function
  ! times a gaussian.  Its absolute value is periodic

  write(6,'("setting up z_function table")')
  call setup_z_function_table        !from z_function_m.f90

  ! make quantum distance and form factor tables



  write(6,'("# setting up form-factor table")')
  call setup_formfactor_table


  


  !  setup the wavefunction, int this case Laughlin
  ! select center-of-mass quantum numbers for Laughlin state  
  sl2z(1,1) = 1
  sl2z(1,2) = 0
  sl2z(2,1) = 0
  sl2z(2,2) = 1
  k = 0

  write(6,'("setting up model state wavefunction")')
  call setup_laughlin_state(nel,m_laughlin,sl2z,k)  !from wf_tools.f90
  write(6,'("Laughlin state on lattice with nel = ",i5)') nel



! chose the step length
  write(6,'("setting up quantum distance table for step lengths")')
  allocate (s2(norb,norb))
  call make_quantum_distance_table(norb, s2) 

  write(6,'("give seed for random number generator")')
  read(5,*) seed




! max_changes is the maximum number of particles
! that may move in a single MC step.   This is currently
! restructed to max_changes = 1
  max_changes = 1
  allocate (update_list(max_changes))
  allocate (x(2,nel), changed_x(nel))
  


!----------------------------------------------
! equilibration:

  
  write (6,'("give initial step_cutoff d < 1,  number of initial equilibration steps and  history_length")')
  read(5,*)   step_cutoff,equilibration, history
  call step_table(norb,step_cutoff,s2)


  !  initialize the state
  allocate(swap_mask(1,1))   ! not used here
  call initialize(nel,norb,fermionlike,x,seed,.false.,nswap,swap_mask)

  ring = history
  allocate (weight_history(0:ring-1))
! compute  wavefunction amplitude
  call get_wf(nel,x,wf)
  weight = abs(conjg(wf)*wf)
  weight_history  = weight
  n_reject = 0
  n_accept = 0
 
  max_weight = 0
  do step = 1,equilibration
     call metropolis_step(nel,norb,fermionlike,x,seed,weight,wf,max_changes,update_list)
     changed_x(:)  = .false.
     if(update_list(1) == 0) then
        n_reject = n_reject + 1
     else
        n_accept = n_accept + 1
     endif
     current = modulo(step,ring)
     weight_history(current) =  weight
     mean_weight = sum(weight_history(:))/ring 
     if(weight > max_weight) max_weight = weight
     write(6,'("n=",i10,"  a=",i10,"  r=",i10," weight =",e20.10," history mean=",e20.10," max=",e20.10)') &
             &  step, n_accept, n_reject, weight, mean_weight, max_weight
  enddo

  write(6,'("equilibration period over")')
  deallocate (weight_history)


!-----------------------------------------------------
! tune acceptance ration
  write(6,'(" give desired acceptance ratio")')
  read(5,*) ratio
  call set_step_break_point(1,success)
  if(.not.success) then
     write(6,'(" set_step_break_point failed")')
     stop
  endif
  training_period = 2000
  do
     call get_step_break_point(breakpoint,nbreaks,stepcount)
     write(6,'("step_count =",i12)') stepcount
     n_accept = 0
     n_reject = 0
     do  step = 1, training_period
        call metropolis_step(nel,norb,fermionlike,x,seed,weight,wf,max_changes,update_list)
        if(update_list(1) == 0) then
           n_reject = n_reject + 1
        else
           n_accept = n_accept + 1
        endif
     enddo
     write(6,'("steps",i12," accept",i12," reject:",i12)') training_period,n_accept,n_reject
     if (real(n_accept,kind=dp)/n_reject < ratio) exit
     call set_step_break_point(breakpoint+1,success)
     if(.not.success) exit
  enddo
!------------------------------------------------------




! determine sampling interval
  write(6,'(" gve target number (up to ",i3,")  of particles &
       &that move per sampling interval")') nel
  read(5,*) target
  target = min(target,nel)
  sampling_interval = 0
  changed = 0
  training_period = 500
  do while (changed < training_period*target)
     sampling_interval = sampling_interval + 1
     write(6,'("proposed sampling_interval = ",i10)') sampling_interval
     changed = 0
     do j = 1,training_period
        changed_x(:)  = .false.
        do k = 1,sampling_interval
           call metropolis_step(nel,norb,fermionlike,x,seed,weight,wf,max_changes,update_list)
           if(update_list(1) > 0)  then
              do i = 1,max_changes
                 if(update_list(i) == 0) exit
                 changed_x(update_list(i)) = .true.
              enddo
           endif
        enddo
        changed = changed+count(changed_x)
     enddo
  enddo

!----------------------------------------------------------
!  prepare for data sampling
  call sampling_setup(nel,norb)
  


  write(6,'("give number of data units, and mc steps per unit")')
  read(5,*) n_unit,n_steps
 


  


  
! initial data values
  call data_initialize(nel,norb,n_unit)
  changed_x = .true.
  call sample(nel,x,changed_x,0)



! now start the metropolis walkers   
! data is accumulated in n_unit units of n_steps metropolis steps 
! in each unit of the calculation, the data from n_steps of MC steps is
! accumulated, and the average values ove the n_step steps is save,
! n_unit units of such data is  produced and stored.


  do unit = 1,n_unit

     ! initalise data stores for this unit
     n_reject = 0
     n_accept = 0
!     energy_data(:,unit) = 0_dp
!     strip_data(:,:,:,unit) = 0_dp
!     sq_data(:,:,unit) = cmplx(0,kind=dp)
!     sq2_data(:,:,unit) = real(0,kind=dp)

     ! start the next unit of n_steps steps
     do step = 1,n_steps
!$OMP PARALLEL !SHARED(NUMTHREADS), PRIVATE(OMP_THREAD_NUM)
!$OMP SECTIONS
!$OMP SECTION

        changed_x(:)  = .false.
!        write(6,'("step",i12)') step
        do k = 1,sampling_interval
           call metropolis_step(nel,norb,fermionlike,x,seed,weight,wf,max_changes,update_list)
           if(update_list(1) == 0)  then
              n_reject = n_reject + 1
           else
              n_accept = n_accept + 1
              do i = 1,max_changes
                 if(update_list(i) == 0) exit
                 changed_x(update_list(i)) = .true.
              enddo
           endif
        enddo
        changed = count(changed_x)
!        write(6,'(i3,1x,100L2)') changed, changed_x

!   sample the data
!        write(6,'("sample",i5)') nel
        call sample(nel,x, changed_x,unit)

     enddo
!$OMP END SECTIONS
 !$OMP END PARALLEL
     
     ! unit is completed, process data

      do i = 1,nel
        write(12,'(2i5)') x(1,i),x(2,i)
     enddo

! combine all data in this unit
     call data_combine(unit,n_steps)

      write(12,'(3i6,3i16)') unit,nel,norb,n_steps, n_accept, n_reject
      write(6,'(3i6,i16," accepted =",i16," rejected =",i16)')  &
          & unit,nel,norb,n_steps, n_accept, n_reject

  enddo
  
  ! final summary
  write(13,'(3i6,3i16)') nel,norb,n_steps*n_unit
  do i = 1,nel
     write(12,'(2i5)') x(1,i),x(2,i)
  enddo
  call data_finalize(n_unit)
  stop
end program single_sample



subroutine get_wf(nel,x,wf)
  implicit none
  integer, parameter :: dp = kind(1.0d0)
  integer, intent(in) :: nel
  integer, intent(in) :: x(2,nel)
  complex(kind=dp), intent(out) :: wf

  complex (kind=dp) :: unnormalized_laughlin_wf
  integer, save :: m = 0
  if( m == 0) call get_m_laughlin(m)
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


