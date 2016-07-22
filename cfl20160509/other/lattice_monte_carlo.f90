module mc_data_m

  integer, parameter, private:: dp = kind(1.0d0)
  type  mcdata_t
     integer (kind=8) :: nf1,nf2,nswaps, count
     real (kind=dp) :: dswaps
     complex (kind=dp) :: zswaps 
  end type mcdata_t

contains
  subroutine initialize_mc_data(nel,data)
    type (mcdata_t) :: data(0:nel)
    integer, intent(in) :: nel
    
    integer :: i
    
    forall (i = 0:nel) data(i)%count = 0
    forall (i = 0:nel) data(i)%nf1 = 0
    forall (i = 0:nel) data(i)%nf2 = 0 
    forall (i = 0:nel) data(i)%nswaps = 0
    forall (i = 0:nel) data(i)%dswaps = 0_dp
    forall (i = 0:nel) data(i)%zswaps = cmplx(0,kind=dp)
    return
  end subroutine initialize_mc_data
end module mc_data_m
  
program full_mag
  use mc_data_m
  implicit none
  integer, parameter :: dp = kind(1.0d0)

  integer :: sl2z(2,2), nel, m_laughlin, norb
  logical :: fermionlike, changed
  
  real(kind=dp) :: metric (2,2)

  complex (kind=dp) :: wf1, wf2, wf1_swap, wf2_swap, phase
  
  real (kind=dp) :: weight1,weight2, amplitude
  
  integer (kind=8) :: nswap1, nswap2, nswap,  n_reject1, n_reject2
  
  integer, allocatable :: x1(:,:),  x2(:,:),  x1_swap(:,:), x2_swap(:,:)
  logical, allocatable :: swap_mask(:,:), swap1(:), swap2(:)
  
  integer :: seed1, seed2
  real::dd, r4_uni
  real(kind=dp) :: Renyi_2, Reni
  integer(kind=8) :: mcsteps, loop, counter
 
  integer(kind=8), allocatable :: Nf1(:),Nf2(:),Nswaps(:)
  complex(kind=dp), allocatable :: Zswaps(:)
  real(kind=dp), allocatable:: Dswaps(:)


  complex (kind=dp) , allocatable :: zrenyi_2(:)
  real (kind=dp), allocatable :: drenyi_2(:), factor1(:),factor2(:)




  type (mcdata_t), allocatable:: full_data(:)
  type (mcdata_t), allocatable:: running_data(:)
  real (kind=dp) , allocatable :: udata(:,:), energy_data(:,:)

  real (kind=dp) :: u,u_runing, u2_running, uc1,uc2, uav,var,uav1,var1, uav2,var2
  integer :: segment, range
  integer ::  i, j, k, m, n
  
!  write(6,'(" nu = 1/m Laughlin : give nel, m ")')
!  read (5,*) nel, m_laughlin
  nel = 20
!  nel = 4
  m_laughlin = 3
  norb = nel*m_laughlin
  fermionlike = .true.
  
  u = 0_dp
  u_runing = 0
  
  range  = 1000000
  allocate (udata(0:range-1,2))
  allocate(energy_data(2,2))
  energy_data = 0.0_dp


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




! setup a logical array  swap_mask to define the swap region
  allocate (swap_mask(norb,norb), swap1(nel),swap2(nel))
  swap_mask = .false.
  forall (i = 1 : norb/2) swap_mask(i,:) = .true.
  
! lattice coordinates for the two  systems
  allocate (x1(2,nel),x2(2,nel))
  ! coordinates for the swaps
  allocate (x1_swap(2,nel),x2_swap(2,nel))
  
  
  ! two independent monte_carlo evaluations are run simultaneeously
  ! for computation of 2nd renyi entangement entropy
  
  
  seed1 = 1763528
  seed2 = 9165426
  
 
  ! dd is the MC step length empirically chosen to
  ! equalize the accept/reject rates
  ! this changes with particle number
  !dd=0.94300 ! N=4 fermi LJ
  dd=0.77050    ! N=20 Lauglin NJ
  call set_step_range(dd)
  
  
  
  
  ! loop is kind = 8 integer that sets the
  ! number requested MC stes=ps for this run
  loop=1000000000
  segment = 1000000
  
  ! full_data(i)%nf1, full_data(i)%nf2(i) are total numbers of
  ! configurations sampled   where i particle were present  in
  !  copies 1 and 2 of the swap region.   full_data(i)%nswaps is 
  ! the number of samples in which i electrons were found
  ! in both copies of the swap region, and a swap was computed

  ! full_data(i)%zswaps isthe total i-particle swap amplitude.
  !  = 1 for i = 0 or nel

  ! full_data(i)%sswaps isthe total i-particle swap absolute amplitude.
 
  allocate (full_data(0:nel),running_data(0:nel))
  call initialize_mc_data(nel,full_data)
  call initialize_mc_data(nel,running_data)


  allocate (zrenyi_2(0:nel),drenyi_2(0:nel), factor1(nel), factor2(nel))
  

  
  n_reject1 = 0
  n_reject2 = 0

  !  initialize the two copies of the state
  call initialize(nel,norb,fermionlike,x1,seed1,.false.,nswap,swap_mask)
  nswap = 0
  do i = 1,nel
     if(swap_mask(x1(1,i),x1(2,i))) nswap = nswap + 1
  enddo
  call initialize(nel,norb,fermionlike,x2,seed2,.true.,nswap,swap_mask)

  call coulomb_energy(nel,norb,x1,uc1)
  call coulomb_energy(nel,norb,x2,uc2)
  uc1 = uc1/nel
  uc2 = uc2/nel
  udata(:,1) = uc1
  udata(:,2) = uc2
  


! random number seeds for the two copies
  seed1 = 1763592
  seed2 = 9165428

  call get_wf(nel,x1,wf1)
  weight1 = abs(conjg(wf1)*wf1)
  call get_wf(nel,x2,wf2)
  weight2 = abs(conjg(wf2)*wf2)


  forall (i = 1:nel) swap1(i) = swap_mask(x1(1,i),x1(2,i))
  nswap1 = count(swap1)
  forall (i = 1:nel) swap2(i) = swap_mask(x2(1,i),x2(2,i))



  nswap2 = count(swap2)


  full_data(nswap1)%nf1 = 1
  full_data(nswap2)%nf2 = 1
  running_data(nswap1)%nf1 = 1
  running_data(nswap2)%nf2 = 1
  


  if(nswap1/=nswap2) then
     write(6,'("intialization failure",2i5)') nswap1,nswap2
     stop
  endif
  !numbers in swap region match, do swap 
  nswap = nswap1 



  call swap(nel,nswap,swap1,swap2,x1,x2,x1_swap,x2_swap)
  call get_wf(nel,x1_swap,wf1_swap)
  call get_wf(nel,x2_swap,wf2_swap)
  call swap_amplitude(wf1,wf2,wf1_swap,wf2_swap,phase,amplitude)

  full_data(nswap)%zswaps = amplitude*phase
  full_data(nswap)%dswaps = amplitude
  full_data(nswap)%nswaps = 1

  running_data(nswap)%zswaps = amplitude*phase
  running_data(nswap)%dswaps = amplitude
  running_data(nswap)%nswaps = 1


  mcsteps = 0
  counter = 0
! now start the walkers
  

     do while (mcsteps <= loop)

        mcsteps = mcsteps + 1
        counter = counter + 1



!$OMP PARALLEL !SHARED(NUMTHREADS), PRIVATE(OMP_THREAD_NUM)
!$OMP SECTIONS
     
!$OMP SECTION
     call metropolis_step(nel,norb,fermionlike,x1,seed1,weight1,wf1,changed)
     if(.not. changed) n_reject1 = n_reject1 + 1
     if(changed)call coulomb_energy(nel,norb,x1,uc1)
     udata(modulo(mcsteps,range),1) = uc1/nel
     energy_data(1,1) = energy_data(1,1) + uc1/nel
     energy_data(2,1) = energy_data(2,1) + (uc1/nel)**2

     !   now repeat this for the second copy
!$OMP SECTION
     call metropolis_step(nel,norb,fermionlike,x2,seed2,weight2,wf2,changed)
     if(.not. changed) n_reject2 = n_reject2 + 1
     if(changed)call coulomb_energy(nel,norb,x2,uc2)
     udata(modulo(mcsteps,range),2) = uc2/nel
     energy_data(1,2) = energy_data(1,2) + uc2/nel
     energy_data(2,2) = energy_data(2,2) + (uc2/nel)**2


!$OMP END SECTIONS
!$OMP END PARALLEL
! 




     if(mod(mcsteps,segment) == 0)then
         write(6,'(20("-"))')

         uav = sum(udata(:,1))/range
         var = (sum(udata(:,1)**2)/range) - uav**2
         write(6,'("U1, varU1",2f16.10)') uav,var
         uav = sum(udata(:,2))/range
         var = (sum(udata(:,2)**2)/range) - uav**2
         write(6,'("U2, varU2",2f16.10)') uav,var



         write(6,'(4e25.16)') wf1,wf2
         print*,mcsteps,n_reject1,mcsteps-n_reject1
         print*,mcsteps,n_reject2,mcsteps-n_reject2
         call sample_renyi(counter,nel,running_data,factor1,factor2,drenyi_2,zrenyi_2)
         counter = 0
         call initialize_mc_data(nel,running_data)
         
         write(6,'("  nm     nf1/nswp    nf2/nswp        drenyi(nm)            zrenyi(nm)")')
         do i = 0,nel
            write(6,'(i3,2f12.5,f25.16," (",f25.16,",",f25.16,")")') i,&
                 &  factor1(i), factor2(i),drenyi_2(i), zrenyi_2(i)
         enddo
         write(6,'("exp -R_2:",2f25.16)') sum(zrenyi_2)

         write(6,'(1x)')
         write(6,'(" full data so so far")')
         call sample_renyi(mcsteps,nel,full_data,factor1,factor2,drenyi_2,zrenyi_2)

         write(6,'("  nm     nf1/nswp    nf2/nswp        drenyi(nm)            zrenyi(nm)")')
         do i = 0,nel
            write(6,'(i3,2f12.5,f25.16," (",f25.16,",",f25.16,")")') i,&
                 &  factor1(i), factor2(i),drenyi_2(i), zrenyi_2(i)
         enddo
         write(6,'("exp -R_2:",2f25.16)') sum(zrenyi_2)

         uav1 = energy_data(1,1)/mcsteps
         var1 = energy_data(2,1)/mcsteps - uav1**2
         uav2 = energy_data(1,2)/mcsteps
         var2 = energy_data(2,2)/mcsteps - uav2**2
         write(6,'("energy: u_av1, variance1 =",2f25.16)') uav1,var1
         write(6,'("energy: u_av2, variance2 =",2f25.16)') uav2,var2

      endif


!----------------------------------
!     now do the swap calculation
! swap_mask(1:norb,1:norb) = .true.  defines the swap region

      forall (i=1:nel) swap1(i) = swap_mask(x1(1,i),x1(2,i))
      forall (i=1:nel) swap2(i) = swap_mask(x2(1,i),x2(2,i))
      nswap1 = count(swap1)
      nswap2 = count(swap2)
 
      full_data(nswap1)%nf1 = full_data(nswap1)%nf1 + 1
      full_data(nswap2)%nf2 = full_data(nswap2)%nf2 + 1

      running_data(nswap1)%nf1 = running_data(nswap1)%nf1 + 1
      running_data(nswap2)%nf2 = running_data(nswap2)%nf2 + 1

      if(nswap1==nswap2) then     !numbers in swap region match, do swap 
         nswap = nswap1
         call swap(nel,nswap,swap1,swap2,x1,x2,x1_swap,x2_swap)
         call get_wf(nel,x1_swap,wf1_swap)
         call get_wf(nel,x2_swap,wf2_swap)

         call swap_amplitude(wf1,wf2,wf1_swap,wf2_swap,phase,amplitude)

         full_data(nswap)%zswaps = full_data(nswap)%zswaps + amplitude*phase
         full_data(nswap)%dswaps = full_data(nswap)%dswaps + amplitude
         full_data(nswap)%nswaps = full_data(nswap)%nswaps + 1

         running_data(nswap)%zswaps = running_data(nswap)%zswaps + amplitude*phase
         running_data(nswap)%dswaps = running_data(nswap)%dswaps + amplitude
         running_data(nswap)%nswaps = running_data(nswap)%nswaps + 1


      endif
   enddo

   write(6,'("monte-carlo steps:",i20)') mcsteps
   write(6,'("system 1: accepts,rejects =",2i20)') mcsteps - n_reject1, n_reject1
   write(6,'("system 2: accepts,rejects =",2i20)') mcsteps - n_reject2, n_reject2
   uav1 = energy_data(1,1)/mcsteps
   var1 = energy_data(2,1)/mcsteps - uav1**2
   uav2 = energy_data(1,2)/mcsteps
   var2 = energy_data(2,2)/mcsteps - uav2**2
   write(6,'("full run: uav1, var1 =",2f25.16)') uav1,var1
   write(6,'("full run: uav2, varw =",2f25.16)') uav2,var2
   
   
   call sample_renyi(mcsteps,nel,full_data,factor1,factor2,drenyi_2,zrenyi_2)



!   do i = 0,nel
!      write(22,*) i,Zswaps(i)/Nswaps(i),Dswaps(i)/Nswaps(i)
!      write(32,*) i,real(Nf1(i),kind=dp)/mcsteps
!      write(42,*) i,real(Nf2(i),kind=dp)/mcsteps
!      if(Nswaps(i).ne.0)then
!         Renyi_2=(real(Zswaps(i))/Nswaps(i))*(real(Nf1(i),kind=dp)/mcsteps)**2 + Renyi_2
!         Reni=(aimag(Zswaps(i))/Nswaps(i))*(real(Nf1(i),kind=dp)/mcsteps)**2 + Reni
!         Renyi_2d(i) = (dswaps(i)/Nswaps(i))*(real(Nf1(i),kind=dp)/mcsteps)**2 + Renyi_2d
!      endif
!   enddo      

   do i = 0,nel
      write(6,'(i3,2f12.5,f25.16," (",f25.16,",",f25.16,")")') i,&
           & factor1(i), factor2(i),drenyi_2(i), zrenyi_2(i)
   enddo
   write(6,'("exp -R_2:",2f25.16)') sum(zrenyi_2)
   


   print*,'S2  = ', -log(abs(sum(zrenyi_2)))
   print*,'imaginary part of exp(-S_2) = ',aimag(sum(zrenyi_2))
  stop
end program full_mag


subroutine metropolis_step(nel,norb,fermionlike,x,seed,weight,wf,changed)
  implicit none
  integer, parameter :: dp = kind(1.0d0)
  integer, intent(in) :: nel, norb
  logical, intent(in) :: fermionlike
  integer, intent(inout) :: seed,x(2,nel)
  complex (kind=dp), intent(inout) :: wf
  real(kind=dp), intent(inout) :: weight
  logical, intent(out):: changed
!----------------------------------------  
! FERMIONLIKE= .true. means no site can be 
! multiply-occupied
!
!  metropolis algorithm
!  with one-particle moves by a single "walker"
!---------------------------------------


  integer:: walker, d(2), x_prev(2), i, nel1
  complex (kind = dp) :: new_wf
  real(kind=dp) :: new_weight
  logical :: acceptance, accept
  real :: ratio

! select a step displacement
  call get_step(d,seed)

! null moves that dont change anything return changed = .false.
  if (modulo(d(1),norb) == 0 .and. (modulo(d(2),norb) == 0)) then
     changed = .false.
     return
  endif

  accept = .true.
! select a particle to be the "walker"  
  call select_particle(nel,walker,seed)
  
! store the old position of the walker
  x_prev(:) = x(:,walker)

! displace the walker
  x(:,walker) = 1 + modulo(x(:,walker) + d(:) - 1,norb)
  
  if(fermionlike) then
     ! check for double occupancy, if found reject move
     do i = 1,nel
        if(i == walker) cycle
        if(x(1,i) /= x(1,walker)) cycle
        if(x(2,i) /= x(2,walker)) cycle
        accept = .false.
        exit
     enddo
  endif
  if(accept) then
     nel1 = nel
     call get_wf(nel1,x,new_wf)
     new_weight = abs(conjg(new_wf)*new_wf)
     if(new_weight== 0_dp) then
        accept = .false.
     else if (new_weight < weight) then
        ratio =  real(new_weight/weight)
        accept = acceptance(ratio, seed)
     endif
  endif
  if(accept) then
     wf = new_wf
     weight = new_weight
     changed = .true.
  else
     x(:,walker) = x_prev(:)  ! put walker back in old position
     changed = .false.
  end if
  return
end subroutine metropolis_step

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
  
subroutine select_particle(nel,n,seed)
  implicit none
  integer, intent(in) :: nel
  integer, intent(inout) :: seed
  integer, intent(out) :: n
  
  real :: r4_uni

  n = 1 + floor(nel*r4_uni(seed))
  return
end subroutine select_particle

function acceptance(ratio, seed) result(accept)
  implicit none
  logical :: accept
  real, intent(in)  :: ratio
  integer, intent(inout) :: seed
  real :: r4_uni


  if (ratio > r4_uni(seed)) then
     accept = .true.
  else
     accept = .false.
  endif
  return
end function acceptance
    
subroutine initialize(nel,norb,fermionlike,x,seed,partition,nswap,swap_mask)
implicit none
integer, intent(in) :: nel, norb, nswap
integer, intent(inout) :: seed
integer, intent(out) :: x(2,nel)
logical, intent(in) :: partition,swap_mask(norb,norb)
logical, intent(in) :: fermionlike

!--------------------------------------
! initialize x1(2,nel) randomly,
! with no multiple site occupation if
! FERMIONLIKE = .true.
!
! if partition = ,true. initialize with nswap particles
! in the swap regiion defined by swap_mask(i,j) = .true.
!-----------------------------------------------
integer num, i,  j, swap_size, n, n0
logical,allocatable :: occupied(:,:)
real :: r4_uni

If (nel < 1 .or. norb < 1) then
   write(6,'("invalid input to INITIALIZE: nel, norb=",2i5)') nel, norb
   stop
endif

if (fermionlike .and. nel > norb**2) then
   write(6,'("cannot place NEL = ",i5,&
        & " fermionlike particles on ",i5, &
        &" x ",i5," sites ")') nel, norb, norb
   stop
endif
if (partition) then
   swap_size = count(swap_mask)
   if (nswap < 0 .or. nswap > nel) then
      write(6,'("invalid partition input to INITIALIZE: nel, nswap",2i5)') &
           & nel, nswap
      stop
   endif
   if (fermionlike .and. (nswap > swap_size .or.  &
        & nel - nswap > norb**2 - swap_size)) then
      write(6,'("cannot place nswap = ",i5,&
           & " fermionlike particles on ",i5," swap sites ")') nswap,swap_size
      write(6,'("or cannot place nel - nswap = ",i5,&
           & " fermionlike particles on ",i5," non-swap sites ")') &
           & nel - nswap, norb**2 - swap_size
      stop
   endif
endif


allocate (occupied(norb,norb))
if(partition) then
   occupied  = swap_mask
   n = nel - nswap
else
   occupied = .false.
   n = nel
endif

! if partition = .true., put n  particles in the non-swap region,
! and the rest in the swap region.
num = 0
do while (num < nel)
   if(partition  .and. num == n)  occupied = (.not. swap_mask)
   n0 = num
   do while (num == n0)
      i = 1 + floor(r4_uni(seed)*norb)
      j = 1 + floor(r4_uni(seed)*norb)
      if(occupied(i,j)) cycle
      num = num + 1
      x(1,num) = i
      x(2,num) = j
      if(fermionlike) occupied(i,j) = .true.
   enddo
enddo
return
end subroutine initialize

subroutine swap(nel,nswap,swap1,swap2,x1,x2,x1_swap,x2_swap)
  implicit none
  integer, intent(in):: nel,nswap
  logical, intent(in) :: swap1(nel),swap2(nel)
  integer, intent(in) :: x1(2,nel),x2(2,nel)
  integer, intent(out) :: x1_swap(2,nel),x2_swap(2,nel) 
  
  integer :: nswap1, nswap2, i, j
  
  integer, allocatable, save :: swapper_1(:), swapper_2(:)
  integer, save :: nel_swap = 0

  if (nel /= nel_swap) then
     if(allocated(swapper_1)) deallocate(swapper_1)
     if(allocated(swapper_2)) deallocate(swapper_2)
     allocate (swapper_1(nel),swapper_2(nel))
     nel_swap = nel
  endif

  nswap1 = 0
  nswap2 = 0
  do i = 1,nel
     if(swap1(i)) then
        nswap1 = nswap1 + 1
        swapper_1(nswap1) = i
     endif
     if(swap2(i)) then
        nswap2 = nswap2 + 1
        swapper_2(nswap2) = i
     endif
  enddo

  if(nswap1 /= nswap2) then
     write(6,'("swap count error",3i5)') nel,nswap1,nswap2
     stop
  endif
  
  forall (i = 1:nel , .not. swap1(i)) x1_swap(:,i) = x1(:,i)
  forall (i = 1:nel , .not. swap2(i)) x2_swap(:,i) = x2(:,i)
  forall (j = 1:nswap)  x1_swap(:,swapper_1(j)) = x2(:,swapper_2(j))
  forall (j = 1:nswap)  x2_swap(:,swapper_2(j)) = x1(:,swapper_1(j))
  
  return
end subroutine swap



subroutine swap_amplitude(wf1,wf2,wf1_swap,wf2_swap,phase,amplitude)
  implicit none
  integer, parameter :: dp = kind(1.0d0)
  complex(kind=dp), intent(in) :: wf1,wf2,wf1_swap,wf2_swap
  complex(kind=dp), intent(out) :: phase
  real(kind=dp), intent(out) :: amplitude
  
  real(kind=dp) :: amp1,amp2,amp1_swap,amp2_swap
  
  amp1 = abs(wf1)
  amp2 = abs(wf2)
  amp1_swap = abs(wf1_swap)
  amp2_swap = abs(wf2_swap)
  if(amp1_swap == 0_dp .or. amp2_swap == 0_dp) then
     amplitude = 0_dp
     phase = 1_dp
     return
  endif

  phase = (wf1_swap/amp1_swap) * (wf2_swap/amp2_swap) 
  phase = phase *(conjg(wf1)/amp1) * (conjg(wf2)/amp2) 

  amplitude = (max(amp1_swap,amp2_swap)/max(amp1,amp2)) * &
       & (min(amp1_swap,amp2_swap)/min(amp1,amp2)) 

  ! amp = abs(phase)
  ! phase = phase/amp
  ! amplitude = amplitude * amp

  return
end subroutine swap_amplitude

module step_m
real dd
end module step_m

subroutine set_step_range(dd1)
  use step_m
  implicit none
  real, intent(in) :: dd1
  dd = dd1
  return
end subroutine set_step_range

subroutine get_step(d,seed)
use step_m
implicit none
integer, intent(inout) :: seed
integer, intent(out) :: d(2)
!-------------------------------------------
! step chooser:
! this follow ed's code : should it be  improved ?
! DD (from module step_m)  must be tuned to equalize 
! acceptance and rejection rates
!------------------------------------------
real::  r4_uni

d(1) = nint((2*r4_uni(seed) - 1.0)*dd)
d(2) = nint((2*r4_uni(seed) - 1.0)*dd)

return
end subroutine get_step



subroutine sample_renyi(mcsteps,nel,data,factor1,factor2,drenyi_2,zrenyi_2)
  use mc_data_m
  implicit none
  integer, parameter :: dp= kind(1.0d0)
  integer, intent(in) :: nel
  integer (kind=8) :: mcsteps
  type (mcdata_t), intent(in):: data(0:nel)
  
  real(kind=dp), intent(out) :: factor1(0:nel), factor2(0:nel),drenyi_2(0:nel)
  complex (kind=dp), intent(out) :: zrenyi_2(0:nel)

  integer :: i

  zrenyi_2 = cmplx(0,kind=dp)
  drenyi_2 = real(0,kind=dp)
  do i = 0,nel
     factor1(i) = real(data(i)%nf1,kind=dp)/mcsteps
     factor2(i) = real(data(i)%nf2,kind=dp)/mcsteps
     if(data(i)%nswaps == 0) cycle 
     zrenyi_2(i) = (data(i)%zswaps/data(i)%nswaps)*factor1(i)*factor2(i)
     drenyi_2(i) = (data(i)%dswaps/data(i)%nswaps)*factor1(i)*factor2(i)
  enddo
  return
end subroutine sample_renyi
