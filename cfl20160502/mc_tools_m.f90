subroutine metropolis_step(nel,norb,fermionlike,x,seed,weight,wf,max_changes,update_list)
  implicit none
  integer, parameter :: dp = kind(1.0d0)
  integer, intent(in) :: nel, norb, max_changes
  logical, intent(in) :: fermionlike
  integer, intent(inout) :: seed,x(2,nel)
  complex (kind=dp), intent(inout) :: wf
  real(kind=dp), intent(inout) :: weight
  integer, intent(out) :: update_list(max_changes)
!----------------------------------------  
! FERMIONLIKE= .true. means no site can be 
! multiply-occupied!
!  metropolis algorithm
!  with one-particle moves by a single "walker"
! max_changes is the limit to the number of particles that coan
! be moved in a step.   This is currently set as 1
!---------------------------------------
  integer:: walker, d(2), x_prev(2), i, nel1
  complex (kind = dp) :: new_wf
  real(kind=dp) :: new_weight
  logical :: acceptance, accept
  real :: ratio

! currently coded for just one walker
  update_list = 0


! select a step displacement
  call get_step(d,seed)

! null moves that dont change anything return changed = .false.
  if (modulo(d(1),norb) == 0 .and. (modulo(d(2),norb) == 0)) then
     update_list(:) = 0
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
     update_list(1) = walker
  else
     x(:,walker) = x_prev(:)  ! put walker back in old position
  end if
  return
end subroutine metropolis_step



subroutine select_particle(nel,n,seed)
  implicit none
  integer, intent(in) :: nel
  integer, intent(inout) :: seed
  integer, intent(out) :: n
  
  real :: r4_random
  
  call get_random(r4_random,seed)
  n = 1 + floor(nel*r4_random)
  return
end subroutine select_particle

function acceptance(ratio, seed) result(accept)
  implicit none
  logical :: accept
  real, intent(in)  :: ratio
  integer, intent(inout) :: seed
  real :: r4_random

  call get_random(r4_random,seed)
  if (ratio > r4_random) then
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
! initialize x(2,nel) randomly,
! with no multiple site occupation if
! FERMIONLIKE = .true.
!
! if partition = ,true. initialize with nswap particles
! in the swap regiion defined by swap_mask(i,j) = .true.
!-----------------------------------------------
integer num, i,  j, swap_size, n, n0
logical,allocatable :: occupied(:,:)
real :: r4_random

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
      call get_random(r4_random,seed)
      i = 1 + floor(r4_random*norb)
      call get_random(r4_random,seed)
      j = 1 + floor(r4_random*norb)
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
! stores the allowed steps   step(1:2,1:step_count)
!
! also the ordering by quantum distance of the table
! s2(norb,norb) used to set up the step list

  integer, allocatable :: step(:,:)

  integer, allocatable, target :: list_order(:,:) 

  integer :: n_steps, n_breaks
  integer, allocatable :: break_point(:)
  integer :: step_count, break



  ! identification stamps for documenting what geometry
  !  was used to setup the above tables
  integer :: norb_s = 0
end module step_m
subroutine get_step_break_point(breakpoint,nbreaks,stepcount)
  use step_m
  implicit none
  integer, intent(out) :: breakpoint,nbreaks,stepcount
  nbreaks = n_breaks
  breakpoint = break
  stepcount = step_count
  return
end subroutine get_step_break_point
subroutine set_step_break_point(breakpoint,success)
  use step_m
  implicit none
  integer, intent(in) :: breakpoint
  logical, intent(out) :: success
  if (breakpoint < 1 .or. breakpoint > n_breaks) then
     success = .false.
     return
  endif
  success = .true.
  break = breakpoint
  step_count = break_point(break)
  return
end subroutine set_step_break_point
  

subroutine get_distance_list_order(norb,distance_list_order_p)
use step_m
implicit none
integer, intent(in) :: norb
integer, pointer :: distance_list_order_p(:,:)
!---------------------------------------
! returns a pointer to the stored list of the
! displacements ordered in increasing quantum distance.
!  distance_list_order(1:2,1:norb**2)
!----------------------------------------


 if(norb /= norb_s) then
    write(6,'("GET_DISTANCE_LIST_ORDER: norb mismatch",2i5)') norb, norb_s
    stop
 endif
 if (.not. allocated(list_order)) then
    write(6,'("GET_DISTANCE_LIST_ORDER: list_order not allocated")')
    stop
 endif
 distance_list_order_p => list_order

return
end subroutine get_distance_list_order


subroutine get_step(d,seed)
use step_m
implicit none
integer, intent(inout) :: seed
integer, intent(out) :: d(2)
!-------------------------------------------
! step chooser:
!------------------------------------------
real::  r4_random
call get_random(r4_random,seed)
d = step(:,1 + floor(r4_random*step_count))
return
end subroutine get_step


subroutine step_table(norb,distance_cutoff,s2)
  use step_m
  implicit none
  integer, parameter :: dp = kind(1.d0)
  integer, intent(in) :: norb
  real(kind=dp), intent(in) :: distance_cutoff
  real(kind=dp) , intent(in) :: s2(norb,norb)
!------------------------------------------
! distance_cutoff is a cutoff on the squared
! quantum distance  (hilbert-schmidt distance)
! of lattice steps.   It must be smaller than 1.
! all steps shorter than this are allowed with
! equal probabiility
!
! it should be tuned to approximately equalize 
! step acceptance and rejection rates.
!----------------------------------------
 real(kind=dp), parameter :: small =  0.001_dp
 real(kind=dp), allocatable ::  d(:), distance(:)
 integer :: i, j, ij,p, q, k, norb2, pos, count
 complex(kind=dp) :: l1,l2
 integer, allocatable :: index(:)


 if(distance_cutoff > 1_dp - small) then
    write(6,'("STEP_TABLE: called with DISTANCE_CUTOFF =",f25.16)') distance_cutoff
    write(6,'(" this must be smaller than 1 - small =")') 1_dp - small
    stop
 endif

 if(allocated(step)) deallocate(step)
 if(allocated(list_order)) deallocate(list_order)
 if(allocated(break_point)) deallocate(break_point)
 norb_s = norb

! order the possible displacements by the quantum distcnce
! between the coherent states on the initial and final sites

 norb2 = norb**2
 allocate (index(norb2),d(norb2))

 d = pack(-s2,.true.)
 call indexd(norb2,d,index)

 count =  0
 do k = 1,norb2
    if (d(index(k)) + small > 0) exit
    count = count + 1
 enddo

 d = d + 1
 n_steps = count
 allocate (step(2,n_steps))
 do k = 1,n_steps
    i = 1 + mod(index(k)-1,norb)
    j = 1 + (index(k)-i)/norb 
    if (2*i > norb) i = i - norb
    if (2*j > norb) j = j - norb
    step(1,k) = i
    step(2,k) = j 
 enddo

 pos = 1
 count = 0
 do 
    if(pos == n_steps) exit
    pos = pos + 1
    if(d(index(pos+1)) - d(index(pos)) < small) cycle
    count = count + 1
 enddo

 n_breaks = count + 1
 allocate (break_point(n_breaks))
 pos = 1
 count = 0
 do 
    if(pos == n_steps) exit
    pos = pos + 1
    if (d(index(pos+1)) - d(index(pos)) < small) cycle
    count = count + 1
    break_point(count) = pos
 enddo
 break_point(n_breaks) = n_steps

do i = 1,n_breaks
   if (i == 1)  then
      k = 0
   else
      k = break_point(i-1)
   endif
   do j = k + 1,break_point(i)
      write(6,'(3i5,f25.16)') i, step(1,j),step(2,j), d(index(j))      
   enddo
enddo


! initial choice of break point in step list set by distance_cutoff
 break = 1
 do i = 2,n_breaks
    if(d(index(break_point(i))) > distance_cutoff) exit 
    break = i
 enddo
write(6,'("break_point = ",i5)') break
step_count = break_point(break) 
 
 deallocate (index,d)
 return
end subroutine step_table

subroutine indexd(n,list,order)
  implicit none
  integer, parameter :: dp = kind(1.0d0)
  integer n
  real (kind=dp), intent(in) :: list(n)
  real (kind=dp) :: item
  integer, intent(out) ::order(n)
  !-------------------------------------
  ! sort a  list, list(1),...,list(n) so
  ! list(order(i)) <=list(order(i+1)) for i = 1,...,n-1
  ! based on  indexx  from "Numerical Recipes" by Press et al.
  !------------------------------------------
  integer low,high,id
  integer i,j
  
  
  ! simple returns, do nothing
  if(n <= 0) return
  
  
  ! initial ordering  
  forall (j=1:n) order(j) = j
  if(n==1) return
    
    
  low = n/2 + 1
  high = n
    
10 continue
  if(low > 1) then
     low = low - 1
     id = order(low)
     item = list(id)
  else
     id = order(high)
     item = list(id)
     order(high) = order(1)
     high = high - 1
     if(high==1) then
        order(1) = id
        return
     endif
  endif
    
  i = low
  j = low + low
  do while (j <= high)
     if(j /= high)then
        if(list(order(j)) < list(order(j+1))) j=j+1
     endif
     if(item < list(order(j))) then
        order(i) = order(j)
        i = j
        j = j + j
     else
        j = high + 1
     endif
  enddo
  order(i) = id
  
  goto 10  
end subroutine indexd



pure function rho_q(nel,norb,x,q,omega) result (rho)
  implicit none
  integer, parameter :: dp = kind(1.0d0)
  integer, intent(in) :: norb, nel
  integer, intent(in) :: x(2,nel), q(2)
  complex (kind=dp) , intent(in) :: omega(0:norb-1)
  complex (kind=dp) :: rho

  rho =  sum(omega(modulo(q(1)*x(1,:) + q(2)*x(2,:),norb)))
  return
end function rho_q


subroutine get_ordered_q_list(norb,q_list,q2_index)
  implicit none
  integer, parameter :: dp = kind(1.0d0)
  integer, intent(in) :: norb
  complex(kind=dp), intent(out) :: q_list(4,norb,norb)
  integer, intent(out) :: q2_index(2,norb**2)
!------------------------------------------------
! produces a list of the norb**2 reciprocal vectors in the Brillouin zone
! defined by the metric (wigner-seitz or voronoi cell), ordered by magnitude.
!
!   q = q_list(k), k = 1:norb**2) is the complex reciprocal vector
!
!   q_index(1:2,k) are its integer indices in the range 1:norb, 
!   so q.x = (2*pi/norb)*(q(1)*x(1) + q(2)*x(2))b
!
!    bzb_list(k) is .true. if q is on the boundary of the Brillouin zone
!-----------------------------------
!

  integer, allocatable :: index(:)
  real(kind=dp), allocatable :: q2(:)
  integer :: i,j,k,id
  complex (kind=dp) :: q4(4)
  logical :: bzb



  allocate (q2(norb**2),index(norb**2))
  do i = 1,norb
     do j = 1,norb
        call reduced_q(norb,i,j,q_list(:,i,j))
        id = i + (j-1)*norb
        q2(id) = abs(q_list(1,i,j))
     enddo
  enddo

! make an ordered list by magnitude of allowed reciprocal vectors in the Brillouin zone

  call indexd(norb**2,q2,index)
  do k = 1,norb**2
     id = index(k)
     i =  1 + mod(id-1,norb)
     j  = 1 + ((id - i)/norb)
     q2_index(1,k) = i
     q2_index(2,k) = j
  enddo
  deallocate (index,q2)

  return
end subroutine get_ordered_q_list


subroutine  reduced_q(norb,i,j,q)
  implicit none
  integer, parameter :: dp = kind(1.0d0) 
  integer, intent(in) :: norb,i,j
  complex(kind=dp), intent(out)  :: q(4) 
  !----------------------------------------------                                                                                       
  !  for q = i*g(1) + j*g(2) 
  !  return  the value of q reduced to the wigner-Seitz
  ! (voronoi cell ) Brillouin zone
  !  also return up to 4 equivalent points if
  ! q is on the  Brillouin zone boundary
  !--------------------------------
  complex(kind=dp), save :: g1_prev = (0_dp,0_dp), g2_prev = (0_dp,0_dp)
  integer:: sl2z(2,2), i0, j0, norb_w,s,sa,sb,i0a,i0b,j0a,j0b,ma,na,mb,nb,k, nq
  integer, save :: m(3),n(3)
  logical :: change(3)
  real (kind=dp), save  :: a,b,c
  complex (kind=dp) :: new_tau, tau, g1,g2,qa,qb,q1,l1,l2
  real (kind=dp) :: x,y, shift,test, pi 
  real(kind=dp) :: tiny = 100*epsilon(1.0_dp)
  ! get geometry parameters
                                                                                                                                               
  call get_g(norb_w,g1,g2)
  if(norb /= norb_w) then
     write(6,'("REDUCED_Q: call with mismatched norb:",2i5)') norb,norb_w
     stop
  endif
  call get_l(norb_w,l1,l2)
  

  q = cmplx(0,kind=dp)
  nq = 0

  i0 = modulo(i,norb)
  j0 = modulo(j,norb)
  if(2*i0 > norb) i0 = i0 - norb
  if(2*j0 > norb) j0 = j0 - norb
  if(i0 == 0 .and. j0 == 0) then
     q = cmplx(0,kind=dp)  ; return
  endif


  if(g1 /= g1_prev .or. g2 /= g2_prev) then
     tau = g2/g1
     call optimize_tau(tau,new_tau,sl2z)
     m(1:2) = norb*sl2z(:,1)
     n(1:2) = norb*sl2z(:,2)
     ma = m(1) + m(2) ;   na = n(1) + n(2)
     mb = m(1) - m(2) ;   nb = n(1) - n(2)
     if (abs(ma*g1 + na*g2) < abs(mb*g1+nb*g2)) then
        m(3) = ma
        n(3) = na
     else
        m(3) = mb
        n(3) = nb
     endif
     a = real(conjg(g1)*g1)
     b = real(conjg(g1)*g2)
     c = real(conjg(g2)*g2)
     g1_prev = g1
     g2_prev = g2
  endif


  ! the three optimum step directions are
  !  (i,j) ->  ( i + k*m(:), j + k*n(:)) 

  change(:) = .true.
  do while (any(change))
     change(:) = .true.
     do k = 1,3
        y = a*m(k)*i0 + b*(n(k)*i0 + m(k)*j0) + c*n(k)*j0
        x = a*m(k)**2   + b*2*m(k)*n(k) + c*n(k)**2
        shift = y/x
    ! is shift close to a half integer ?
        test = shift  - 0.5_dp*(2*floor(shift) + 1)
        if(abs(test) > tiny) then
           s = nint(shift)
        else  ! shift is close to a half-integer, need to break a tie
           sa = floor(shift)
           sb = sa + 1
           i0a = i0 - sa*m(k) ; j0a = j0 - sa*n(k) ; qa = i0a*g1 + j0a*g2
           i0b = i0 - sb*m(k) ; j0b = j0 - sb*n(k) ; qb = i0b*g1 + j0b*g2
           if(real(qa) <  real(qb)) then
              s = sb
           else if(real(qa) > real(qb)) then
              s = sa
           else if (aimag(qa) >  aimag(qb)) then
              s = sa
           else
              s = sb
           endif
        endif
        if(s == 0) then
           change(k) = .false.
        else
           i0 = i0 - s*m(k) ; j0 = j0 - s*n(k)
        endif
     enddo
  enddo

  nq = 1
  q(1) = i0*g1 + j0*g2
  ! BZB TEST
  do k = 1,3
     do s = -1,1,2
        q1 = (i0 + s*m(k))*g1 + (j0 + s*n(k))*g2
        if (abs(abs(q(1))-abs(q1)) > tiny ) cycle
        nq = nq + 1
        if(nq > 4) then
           write(6,'("ERROR in REDUCED Q, NQ > 4!")')
           stop
        endif
        q(nq) = q1
     enddo
  enddo

  return
end subroutine reduced_q

