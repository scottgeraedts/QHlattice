module sq_m
integer, parameter, private :: dp1 = kind(1.0d0)

real (kind=dp1), allocatable :: form_factor(:,:)
integer, allocatable :: q2_index(:,:), q_set(:,:)
integer :: nsq, norb_sq, npair, nel_sq

  complex (kind=dp1), allocatable ::  rho(:,:),q_list(:,:,:)
  real (kind=dp1), allocatable :: sq2_data(:,:,:), sq2(:,:), q2(:)
   complex (kind=dp1),allocatable  :: sq_data(:,:,:),sq(:,:)
  complex (kind=dp1), allocatable:: sq_combined(:)
  real(kind=dp1), allocatable :: sq2_combined(:)
  real (kind=dp1), allocatable  :: sq2_current(:,:)
  complex (kind=dp1), allocatable :: omega(:)

  logical, parameter :: exclude_inverse_q = .true., exclude_small_form_factors = .true.


end module sq_m





subroutine sq_setup(nel,norb)
  use sq_m
  implicit none
  integer, parameter :: dp = kind(1.0d0)
  integer, intent(in) :: nel, norb
  
  real (kind=dp) :: small_ff
  integer :: i,j,k,i0,j0

  complex (kind=dp) :: l1,l2
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

  allocate (form_factor(norb,norb))
  call get_form_factor (norb,l1,l2,form_factor)  



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

  
  allocate (omega(0:norb-1))
  call get_omega(norb,omega)
  norb_sq = norb
  nel_sq = nel

  return
end subroutine sq_setup

subroutine sq_data_init(nel, norb, n_unit)
  use sq_m
  implicit none
  integer, parameter :: dp = kind(1.0d0)
  integer, intent(in) :: nel, norb, n_unit


  integer :: i,j,k

  ! structure factor
  npair = (nel*(nel-1))/2
  allocate (sq2_current(npair,nsq), rho(nel,nsq),sq(nel,nsq), sq2(npair,nsq))  
  allocate (sq_data(nel,nsq,n_unit),sq2_data(npair,nsq,n_unit))
  allocate (sq_combined(nsq), sq2_combined(nsq))
  sq_data = cmplx(0,kind=dp)
  sq2_data = real(0,kind=dp)
  sq2_current = real(0,kind=dp)
  rho = cmplx(0,kind=dp)
  sq = cmplx(0,kind=dp)
  sq2 = real(0,kind=dp)
  sq_combined = cmplx(0,kind=dp)
  sq2_combined = real(0,kind=dp)
end subroutine sq_data_init

subroutine sq_data_combine(unit,n_steps)
  use sq_m
  implicit none
  integer, intent(in) :: unit,n_steps

  integer :: i,j,k, nel, norb
  
  nel = nel_sq
  norb = norb_sq

  sq_data(:,:,unit) = sq_data(:,:,unit)/n_steps
  sq2_data(:,:,unit) = sq2_data(:,:,unit)/n_steps

  ! compute s(q) using combined data computed so far
  
  forall (i = 1:nel,k = 1:nsq) sq(i,k) = sum(sq_data(i,k,1:unit))/unit
  forall (i = 1:npair,k = 1:nsq) sq2(i,k)= sum(sq2_data(i,k,1:unit))/unit
  
  forall (i = 2:nel, j = 1:nel-1, k = 1:nsq,  j < i)  sq2(((i-1)*(i-2))/2 + j,k) &
       & = sq2(((i-1)*(i-2))/2 + j,k)  - real(conjg(sq(i,k))*sq(j,k))
  
  forall (k = 1:nsq)  sq_combined(k) = sum(sq(:,k))/nel
  forall (k = 1:nsq) sq2_combined(k) =  2*sum(sq2(:,k))/norb
  
  call print_sq(6,nsq,sq_combined,sq2_combined)

  return
end subroutine sq_data_combine




subroutine sq_sample(nel,x,changed_x,unit)
  use sq_m
  implicit none
  integer, intent(in) :: nel,unit
  integer, intent(in) :: x(2,nel)
  logical, intent(in) :: changed_x(nel)
  integer :: i,j,k

  ! update rho and sq2_current for particles that moved 
  forall (i = 1:nel,changed_x(i))
     forall (k = 1:nsq) rho(i,k) = &
          &  omega(modulo(sum(x(:,i)*q_set(:,k)),norb_sq))
  end forall

  forall (i = 2:nel, j = 1:nel-1,  &
       & ( j < i  .and. (changed_x(i) .or. changed_x(j)))) &
       & sq2_current(((i-1)*(i-2))/2+j,:) = &
       & real(conjg(rho(i,:))*rho(j,:))
 
  if(unit /= 0) then
     sq_data(:,:,unit) =  sq_data(:,:,unit) + rho(:,:)
     sq2_data(:,:,unit) =  sq2_data(:,:,unit) + sq2_current(:,:)     
  endif

 return
end subroutine sq_sample



subroutine  print_sq(print_unit,nsq1,sq_print,sq2_print)
  use sq_m
  implicit none
  integer, parameter :: dp = kind(1.0d0)
  integer, intent(in) :: print_unit, nsq1
  real(kind=dp), intent(in) :: sq2_print(nsq1)
  complex(kind=dp),intent(in) :: sq_print(nsq1)
  
  complex(kind=dp) :: sbar, q
  real(kind=dp) ::  abs_q, ff, nu, sbar2, s0, ss 
  logical :: bzb
  
  integer :: i,j,i1,j1,nel, norb

  nel = nel_sq
  norb = norb_sq

  if(nsq1 /= nsq) then
     write(6,'(" print_sq: incompatible nsq:",2i12)') nsq, nsq1
     stop
  endif


  nu =  real(nel,kind=dp)/norb 
  write(print_unit,'(25x,"|q|",12x,"form_factor",15x,"s(q)",15x,"sbar(q)",15x, "s0(q)",10x, &
       &"(complex) <rho(q)>/Nel",6x,"|<rho(q)>|/nel")')
  do i = 1,nsq
     q =  q_list(1,q_set(1,i),q_set(2,i))
     abs_q = abs(q)
     ff = form_factor(q_set(1,i),q_set(2,i))
     sbar = sq_print(i)
     sbar2 = sq2_print(i)
     ss = sbar2 + nu
     if (abs(q) /= 0_dp) sbar2 = sbar2  + ff*nu
     s0 = sbar2/ff
     write(print_unit,'(3i5,f20.12,es20.12,3f20.12," (",es15.8,",",es15.8,")",es15.8)') &
          &  i,q_set(:,i),abs_q,ff,ss,sbar2,s0,sbar, abs(sbar)
  enddo
  return
end subroutine print_sq



subroutine sq_data_finalize(n_unit)
  use sq_m
  implicit none
  integer, intent(in) :: n_unit
  call print_sq(13,nsq,sq_combined,sq2_combined)
  call print_sq( 6,nsq,sq_combined,sq2_combined)
end subroutine sq_data_finalize
