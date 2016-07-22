program single_sample
  implicit none
  integer, parameter :: dp = kind(1.0d0)

  integer :: sl2z(2,2), nel, m_laughlin, norb, npair
  logical :: fermionlike, changed, updated
  
  real(kind=dp) :: metric (2,2)
  logical, allocatable :: swap_mask(:,:), changed_x(:)
  
  complex (kind=dp), allocatable :: omega(:)
  integer,  allocatable :: num(:)
  
  integer (kind=8) :: n_reject, n_accept, n_steps, step, equilibration, ring
  integer :: n_unit, nswap, strip_width, unit, norb2,id, nsq, n_sample, nsq_max, count
  integer :: seed
  integer, allocatable :: x(:,:), x_prev(:,:)
  complex (kind=dp) :: wf
  real (kind=dp) :: weight,prev, norm
  integer :: current, history, sampling_interval, unchanged_x
  integer :: max_changes, training_period, target
  integer, allocatable :: update_list(:)


! quantum distance
  real (kind=dp), allocatable :: s2(:,:),s2_temp(:,:,:,:)



  complex (kind=dp), allocatable ::  rho(:,:), q_list(:,:,:)
  integer, allocatable :: q_set(:,:),q2_index(:,:)
  real (kind=dp), allocatable :: sq2_data(:,:,:), sq2(:,:), q2(:),  &
       & form_factor(:,:)
  complex (kind=dp),allocatable  :: sq_data(:,:,:),sq(:,:)
  complex (kind=dp), allocatable:: sq_combined(:)
  real(kind=dp), allocatable :: sq2_combined(:)



  real (kind=dp) :: ratio
  integer :: nbreaks,breakpoint,stepcount
  logical :: success


  logical, allocatable :: bzb_list(:)




  integer, allocatable :: strip_count(:,:), strip_start(:)
  integer :: min_width,max_width, unchanged
  logical, allocatable ::  strip_mask(:,:,:)
  real(kind=dp), allocatable ::strip_data(:,:,:,:), strip(:,:,:)

  logical :: do_strip = .true. , do_sq = .true. , do_energy = .true.
  real (kind=dp), allocatable  :: sq2_current(:,:)



  real(kind=dp), allocatable ::energy_data(:,:), running_energy(:,:)
  real(kind=dp) :: energy(2), uc, variance, step_cutoff, mean

  logical, parameter :: exclude_inverse_q = .true., exclude_small_form_factors = .true.
  real(kind=dp) :: small_ff

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
  allocate (form_factor(norb,norb))


  write(6,'("setting up form-factor table")')
  call make_gaussian_form_factor(norb,form_factor,norm)
  form_factor  = form_factor**2
 





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



!---------structures needed for s(q)------------------------

 
! set up table of Fourier factors

  allocate (omega(0:norb-1))
  call get_omega(norb,omega)
 







  ! make ordered list of q in the Brillouin zone
  ! q_list(1,i,j) is the reduced-to-the-Brillouin-zone
  ! complex q with indices i,j
  ! if its on the BZ boundary (up tp 4 equivalent points)
  ! are listed as non-zero q_list(2:4,i,j)
  ! q2_index(i,j) are the indices of the k'th q in list of
  ! increasing abs(q)

  allocate (q_list(4,norb,norb),q2_index(2,norb**2))
  call  get_ordered_q_list(norb,q_list,q2_index)
  
!  rescale q to units with magnetic length = 1
  q_list = q_list/sqrt(real(2,kind=dp))

  ! select sq_set (q points to compute sq at)
  ! for the moment do all NON-INVERSION RELATED for testing purposes


  ! up to norb**2 q values are possible
  ! only nsq of them will be used
  allocate (q_set(2,norb**2))
  

  write(6,'("give smallest values of form-factor f(q)",&
& /" for which s(q) should be  obtained (0 for all q)")')
  read(5,*)  small_ff
  if ( small_ff > 0_dp .and. exclude_small_form_factors) & 
       & write(6,'("s(q): excluding q with form_factor < ",es12.4)') small_ff
     
  nsq = 0
  do k = 1,norb**2
     i = q2_index(1,k)
     j = q2_index(2,k)
     
     if (exclude_inverse_q) then
        i0 = modulo(i,norb)
        if(i0 > norb/2) i0 = i0 - norb
        if (i0 < 0) cycle
        if ( modulo(2*i0,norb) == 0) then
           j0  = modulo(j,norb)
           if (j0 > norb/2) j0 = j0 - norb
           if (j0 < 0) cycle
        endif
     endif
     
     if (exclude_small_form_factors) then
        if (form_factor(i,j)  < small_ff) cycle
     endif

     nsq = nsq + 1
     q_set(:,nsq) = q2_index(:,k)
  end do






!--------structures for renyi entropy  of strips

!  define strip masks, for computing renyi entropy of vertical
!  strip decompositions
  min_width = 1
  max_width = norb - 1

!  the selected strip has width w from 1: max_width < norb, starting at x(1) = i
  max_width = min(max_width,norb-1)

! we will maximise the strip data harvested (n_sample = norb for each desired width)
  allocate(strip_start(norb))
  n_sample = norb
  forall (i=1:norb) strip_start(i) = i


! set up the strip_mask array
  allocate (strip(0:nel,n_sample,min_width:max_width))
  allocate (strip_mask(norb,n_sample,min_width:max_width),strip_count(n_sample,min_width:max_width))
  do w = min_width,max_width
     do k = 1,n_sample
!   the selected region is from  i to j-1, (modulo norb) inclusive
        i = strip_start(k)
        j =  i + w 
        j = 1 +  modulo(j - 1, norb)
        if (j >  i) then
           strip_mask(:,k,w) = .false.
           strip_mask(i:j-1,k,w) = .true.
        else
           strip_mask(:,k,w) = .true.
           strip_mask(j:i-1,k,w) = .false.
        endif
     enddo
  enddo

!-----------------------------------------


  call coulomb_setup

! chose the step length

  write(6,'("setting up quantum distance table for step lengths")')

! could be parallelized over norb processors,  or even norb**2
! if s2_temp was  extended to s2_temp(norb,norb,norb,norb)
! and do loop over n was made a forall loop
  allocate (s2(norb,norb),s2_temp(norb,norb,norb,1))
  do n = 1,norb
     forall  (m = 1:norb)
        forall(i = 1:norb,j = 1:norb) s2_temp(i,j,m,1) = &
             & real(omega(modulo(m*i + n*j, norb)))
        s2(m,n) = sum(form_factor(:,:)*s2_temp(:,:,m,1))/norb
     end forall
  end do
  deallocate(s2_temp)












  max_changes = 1
  allocate (update_list(max_changes))
  allocate (x(2,nel), changed_x(nel),x_prev(2,nel))
  



!--------------MC calculations begin----------------------
  call step_table(norb,0.8_dp,s2)

  write(6,'("give seed for random number generator")')
  read(5,*) seed

  
  write (6,'("give number of initial equilibration steps and  history_length")')
  read(5,*)   equilibration, history

  !  initialize the state
  allocate(swap_mask(1,1))   ! not used here
  call initialize(nel,norb,fermionlike,x,seed,.false.,nswap,swap_mask)

! compute  wavefunction ampliitude
  call get_wf(nel,x,wf)
  weight = abs(conjg(wf)*wf)
  call coulomb_energy(nel,norb,x,uc)

  ring = max(1,history)
  allocate (running_energy(2,0:ring-1))
  running_energy(1,:) = uc/nel
  running_energy(2,:) = (uc/nel)**2

  n_reject = 0
  n_accept = 0
  do step = 1,equilibration
     call metropolis_step(nel,norb,fermionlike,x,seed,weight,wf,max_changes,update_list)
     if(update_list(1) == 0) then
        n_reject = n_reject + 1
     else
        n_accept = n_accept + 1
        ! update state properties
        call coulomb_energy(nel,norb,x,uc)
     endif
     current = modulo(step,ring)
     running_energy(1,current) = uc/nel
     running_energy(2,current) = (uc/nel)**2

     if(modulo(step,ring) == 0) then
        mean = sum(running_energy(1,:))/ring
        variance = (sum(running_energy(2,:))/ring) - mean**2
        write(6,'("n=",i10,"  a=",i10,"  r=",i10,"  e=",f20.10,"  v=",f20.10)') &
             &  step, n_accept, n_reject, mean,variance
     endif
  enddo

  write(6,'("equilibration period over")')
  deallocate (running_energy)


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



  write(6,'("give number of data units, and mc steps per unit")')
  read(5,*) n_unit,n_steps
  allocate (energy_data(2,n_unit))






! initial data values
  if(do_sq) then
     ! structure factor
     npair = (nel*(nel-1))/2
     allocate (sq2_current(npair,nsq), rho(nel,nsq),sq(nel,nsq), sq2(npair,nsq))  
     allocate (sq_data(nel,nsq,n_unit),sq2_data(npair,nsq,n_unit))
     allocate (sq_combined(nsq), sq2_combined(nsq))


     forall (i = 1:nel, k = 1:nsq)  rho(i,k) = &
          &  omega(modulo(sum(x(:,i)*q_set(:,k)),norb))
     forall (i = 2:nel, j = 1:nel-1, k = 1:nsq, &
          & j < i) &
          & sq2_current(((i-1)*(i-2))/2+j,k) = &
          & real(conjg(rho(i,k))*rho(j,k))
  endif

  if (do_strip) then
     ! strip occupation count for entanglement calculation
     allocate (num(norb))
     allocate (strip_data(0:nel,n_sample,min_width:max_width,n_unit))
     num = 0
     do i = 1,nel
        num(x(1,i)) = num(x(1,i)) + 1
     enddo
  endif
  


! now start the metropolis walkers   
! data is accumulated in n_unit units of n_steps metropolis steps 
! in each unit of the calculation, the data from n_steps of MC steps is
! accumulated, and the average values ove the n_step steps is save,
! n_unit units of such data is  produced and stored.


  do unit = 1,n_unit

     ! initalise data stores for this unit
     n_reject = 0
     n_accept = 0
     energy_data(:,unit) = 0_dp
     strip_data(:,:,:,unit) = 0_dp
     sq_data(:,:,unit) = cmplx(0,kind=dp)
     sq2_data(:,:,unit) = real(0,kind=dp)

     ! start the next unit of n_steps steps
     x_prev(:,:) = x(:,:)
     do step = 1,n_steps
!$OMP PARALLEL !SHARED(NUMTHREADS), PRIVATE(OMP_THREAD_NUM)
!$OMP SECTIONS
!$OMP SECTION

        changed_x(:)  = .false.
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

        ! update state properties
        ! this could be modified to just compute change in coulomb energy
        if(do_energy) then
           call coulomb_energy(nel,norb,x,uc)
           energy_data(1,unit) = energy_data(1,unit) + (uc/nel)
           energy_data(2,unit) = energy_data(2,unit) + (uc/nel)**2
        endif

       ! update strip counts
        if (do_strip) then
           do i = 1,nel
              if(.not.changed_x(i)) cycle
                 num(x_prev(1,i)) = num(x_prev(1,i)) - 1
                 num(x(1,i)) = num(x(1,i)) + 1
           enddo
           ! these can be paralelized
           forall (i=1:n_sample, w = min_width:max_width) & 
                & strip_count(i,w)  = sum(pack(num(:),strip_mask(:,i,w)))
           forall (i = 1:n_sample, w = min_width:max_width) & 
             &  strip_data(strip_count(i,w),i,w,unit) = &
             &  strip_data(strip_count(i,w),i,w,unit) + 1        
        endif

        forall (i = 1:nel, changed_x(i)) x_prev(:,i) = x(:,i)

        if (do_sq) then
           ! update rho and sq2_current for particles that moved 
           forall (i = 1:nel,changed_x(i))
              forall (k = 1:nsq) rho(i,k) = &
                &  omega(modulo(sum(x(:,i)*q_set(:,k)),norb))
           end forall
           forall (i = 2:nel, j = 1:nel-1,  &
                & ( j < i .and. (changed_x(i) .or. changed_x(j)))) &
                & sq2_current(((i-1)*(i-2))/2+j,:) = &
                & real(conjg(rho(i,:))*rho(j,:))
           
           sq_data(:,:,unit) =  sq_data(:,:,unit) + rho(:,:)
           sq2_data(:,:,unit) =  sq2_data(:,:,unit) + sq2_current(:,:)
        endif

     enddo
!$OMP END SECTIONS
 !$OMP END PARALLEL
     
     ! unit is completed, process data

      do i = 1,nel
        write(12,'(2i5)') x(1,i),x(2,i)
     enddo

     if (do_strip) then
        strip_data(:,:,:,unit)=strip_data(:,:,:,unit)/n_steps
        do w = min_width,max_width
           call print_strip(12,nel,n_sample,w,strip_data(:,:,w,unit))
        enddo

        ! print  average of strip  data so far accumulated:
        forall (n = 0:nel, i = 1:n_sample, w=min_width:max_width) &
             & strip(n,i,w) = sum(strip_data(n,i,w,1:unit))/unit     
        do w = min_width,max_width
           call print_strip( 6,nel,n_sample,w,strip(:,:,w))
        enddo
     endif


     if (do_sq) then
        sq_data(:,:,unit) = sq_data(:,:,unit)/n_steps
        sq2_data(:,:,unit) = sq2_data(:,:,unit)/n_steps

        ! compute s(q) using combined data computed so far

        forall (i = 1:nel,k = 1:nsq) sq(i,k) = sum(sq_data(i,k,1:unit))/unit
        forall (i = 1:npair,k = 1:nsq) sq2(i,k)= sum(sq2_data(i,k,1:unit))/unit

        forall (i = 2:nel, j = 1:nel-1, k = 1:nsq,  j < i)  sq2(((i-1)*(i-2))/2 + j,k) &
             & = sq2(((i-1)*(i-2))/2 + j,k)  - real(conjg(sq(i,k))*sq(j,k))
        
        forall (k = 1:nsq)  sq_combined(k) = sum(sq(:,k))/nel
        forall (k = 1:nsq) sq2_combined(k) =  2*sum(sq2(:,k))/norb

        call print_sq(6,norb,nsq,sq_combined,sq2_combined,q_set)
     endif

     write(12,'(3i6,3i16)') unit,nel,norb,n_steps, n_accept, n_reject

     if (do_energy) then
        energy_data(:,unit)  = energy_data(:,unit)/n_steps
        variance = energy_data(2,unit) - energy_data(1,unit)**2
        energy(1) = sum(energy_data(1,1:unit))/unit
        energy(2) = sum(energy_data(2,1:unit))/unit
        variance = energy(2) - energy(1)**2
        write(12,'(2f25.16)') &  energy_data(1,unit),energy_data(2,unit),variance
     end if
     

     write(6,'(3i6,i16," accepted =",i16," rejected =",i16)')  &
          & unit,nel,norb,n_steps, n_accept, n_reject
     if (do_energy) then
        write(6,'("energy per particle:",f25.16, &
             &" variance:",f25.16)') energy(1),variance
     endif
  enddo
  
  ! final summary

  if (do_energy) then
     energy(1) = sum(energy_data(1,:))/n_unit
     energy(2) = sum(energy_data(2,:))/n_unit
     variance = energy(2) - energy(1)**2
  endif


  forall (n = 0:nel, i = 1:n_sample, w=min_width:max_width) &
       & strip(n,i,w) = sum(strip_data(n,i,w,:))/n_unit


  write(13,'(3i6,3i16)') nel,norb,n_steps*n_unit
  do i = 1,nel
     write(12,'(2i5)') x(1,i),x(2,i)
  enddo
  write(13,'(2f25.16)') energy(1),energy(2),variance

  do w = min_width,max_width
     call print_strip(13,nel,n_sample,w,strip(:,:,w))
     call print_strip( 6,nel,n_sample,w,strip(:,:,w))
  enddo

  call print_sq(13,norb,nsq,sq_combined,sq2_combined,q_set)
  call print_sq( 6,norb,nsq,sq_combined,sq2_combined,q_set)

  
  if (do_energy) then
     write(6,'("energy per particle:",f25.16, &
          &" variance:",f25.16)') energy(1),variance
  endif
  stop
contains
  subroutine  print_sq(print_unit,norb,nsq,sq,sq2,q_set)
    implicit none
    integer, intent(in) :: norb, print_unit, nsq
    integer :: q_set(2,nsq)
    real(kind=dp), intent(in) :: sq2(nsq)
    complex(kind=dp),intent(in) :: sq(nsq)
    
    complex(kind=dp) :: sbar
    real(kind=dp) ::  abs_q, ff, nu, sbar2 
    logical :: bzb
    
    integer :: i,j,i1,j1     
    nu =  real(nel,kind=dp)/norb 
    write(print_unit,'(25x,"|q|",12x,"form_factor",15x,"s(q)",15x,"s(q)-nu",10x, &
         &"(complex) <rho(q)>/Nel",6x,"|<rho(q)>|/nel")')
    do i = 1,nsq
       q =  q_list(1,q_set(1,i),q_set(2,i))
       abs_q = abs(q)
       ff = form_factor(q_set(1,i),q_set(2,i))
       sbar = sq(i)
       sbar2 = sq2(i)
       write(print_unit,'(3i5,f20.12,es20.12,2f20.12," (",es15.8,",",es15.8,")",es15.8)') &
            &  i,q_set(:,i),abs_q,ff,sbar2+nu,sbar2,sbar, abs(sbar)
    enddo
    return
  end subroutine print_sq

  subroutine print_strip(print_unit,nel,n_sample,strip_width,strip_data)
    implicit none
    integer, intent(in) :: nel,n_sample, print_unit, strip_width
    real(kind=dp), intent(in) ::strip_data(0:nel,n_sample)
    integer ::  n, i
    real(kind=dp) :: test
    real(kind=dp),parameter :: small = 1.0e-7_dp
    
    do n = 0,nel
       test = real(sum(strip_data(n,:)))/n_sample
      if (test < 1.0e-7_dp) cycle
       write(print_unit,'(20("-"))')
       write(print_unit,'("strip width:",i5," n = ",i5,e12.3)') strip_width, n, test
       write(print_unit,'(10f13.10)') (strip_data(n,i), i=1,n_sample)
    enddo
    return
  end subroutine print_strip

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
