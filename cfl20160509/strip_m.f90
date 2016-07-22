
module strip_m
  integer, parameter, private :: dp2=kind(1.0d0)
  integer, allocatable :: strip_count(:,:), strip_start(:)
  integer :: min_width,max_width
  logical, allocatable ::  strip_mask(:,:,:)
  real(kind=dp2), allocatable ::strip_data(:,:,:,:), strip(:,:,:)
  integer :: strip_width
  integer :: n_sample, nel_s
  integer, allocatable :: x_prev(:,:), num(:)
end module strip_m






  subroutine strip_setup(nel,norb)
    use strip_m
    implicit none
    integer, parameter :: dp = kind(1.0d0)
    integer, intent(in) ::nel, norb
    
!--------structures for renyi entropy  of strips

    integer :: i,j, k, w




!  define strip masks, for computing renyi entropy of vertical
!  strip decompositions
  min_width = 1
  max_width = norb - 1

  nel_s = nel



!  the selected strip has width w from 1: max_width < norb, starting at x(1) = i
  max_width = min(max_width,norb-1)

! we will maximise the strip data harvested (n_sample = norb for each desired width)
  allocate(strip_start(norb))
  n_sample = norb
  forall (i=1:norb) strip_start(i) = i


! set up the strip_mask array
  allocate (strip(0:nel,n_sample,min_width:max_width))
  strip = real(0,kind=dp)
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

  allocate (x_prev(2,nel), num(norb))
  num = 0

  return
end subroutine strip_setup


subroutine strip_data_init(nel, norb,n_unit)
  use strip_m
  implicit none
  integer, parameter :: dp = kind(1.0d0)
  integer, intent(in) :: nel, norb, n_unit
  allocate (strip_data(0:nel,n_sample,min_width:max_width,n_unit))
  strip_data = 0_dp
  
  return
end subroutine strip_data_init

subroutine strip_data_combine(unit,n_steps)
  use strip_m
  implicit none
  integer, intent(in) :: unit,n_steps
  integer :: w, i, n, nel

  nel = nel_s
  strip_data(:,:,:,unit) = strip_data(:,:,:,unit)/n_steps


  do w = min_width,max_width
     call print_strip(12,nel,n_sample,w,strip_data(:,:,w,unit))
  enddo
  
  ! print  average of strip  data so far accumulated:
  forall (n = 0:nel, i = 1:n_sample, w=min_width:max_width) &
       & strip(n,i,w) = sum(strip_data(n,i,w,1:unit))/unit     
  do w = min_width,max_width
     call print_strip( 6,nel,n_sample,w,strip(:,:,w))
  enddo  
  return
end subroutine strip_data_combine




subroutine  strip_sample(nel,x,changed_x,unit)
  use strip_m
  implicit none
  integer, intent(in) :: nel,unit
  integer, intent(in) :: x(2,nel)
  logical, intent(in) :: changed_x(nel)
 
  integer :: i, w


  if(all(changed_x)) then
     num = 0
     do i = 1,nel
        num(x(1,i)) = num(x(1,i)) + 1
     enddo
  else
     do i = 1,nel
        if(.not.changed_x(i)) cycle
        num(x_prev(1,i)) = num(x_prev(1,i)) - 1
        num(x(1,i)) = num(x(1,i)) + 1
     enddo
  endif
  forall (i = 1:nel, changed_x(i)) x_prev(:,i) = x(:,i) 

  ! these can be paralelized
  forall (i=1:n_sample, w = min_width:max_width) & 
       & strip_count(i,w)  = sum(pack(num(:),strip_mask(:,i,w)))

  
  if(unit /= 0) then
     forall (i = 1:n_sample, w = min_width:max_width) & 
          &  strip_data(strip_count(i,w),i,w,unit) = &
          &  strip_data(strip_count(i,w),i,w,unit) + 1        
  endif
  

  return
end subroutine strip_sample


subroutine print_strip(print_unit,nel,n_sample,strip_width,strip_data)
  implicit none
  integer, parameter :: dp = kind(1.0d0)
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


subroutine strip_data_finalize(n_unit)
  use strip_m
  implicit none
  integer, intent(in) :: n_unit
  integer :: w, n, i

  forall (n = 0:nel_s, i = 1:n_sample, w=min_width:max_width) &
       & strip(n,i,w) = sum(strip_data(n,i,w,:))/n_unit

  do w = min_width,max_width
     call print_strip(13,nel_s,n_sample,w,strip(:,:,w))
     call print_strip( 6,nel_s,n_sample,w,strip(:,:,w))
  enddo
end subroutine strip_data_finalize



