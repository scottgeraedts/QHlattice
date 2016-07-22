program twobody_test
  implicit none
  integer, parameter :: dp = kind(1.0d0)

  integer :: norb, sl2z(2,2),i
  complex (kind=dp) :: l1,l2
  complex(kind=dp), allocatable :: ham(:,:)
  real(kind=dp), allocatable ::eval(:)
  real (kind=dp) :: metric(2,2)


  write(6,'("give norb")')
  read(5,*) norb
  metric = 0_dp
  call set_geometry(norb,metric)     !from wf_tools.f90

  write(6,'("setting up z_function table")')
  call setup_z_function_table        !from z_function_m.f90


! the Coulomb metric  corresponds to
! |L| = abs(L), where L1 and L1 are the complex representaion of
! the pbc translations and The Landa basis is an eigenstate of
! translation by L1/norb 

  call get_l(norb,l1,l2)  
  write(6,'("L1=",2f25.16)') l1
  write(6,'("L2=",2f25.16)') l2


  call coulomb_setup    ! from new_coulomb_m.f90
  write(6,'("make_landau_coulomb")')
  call make_landau_coulomb

  
  allocate (ham(2*norb,2*norb),eval(2*norb))
  call two_particle_hamiltonian(norb,ham)
  write(6,'("diagonalize")')
  call diagonalize_hermitian(2*norb,ham,eval)
  write(6,'("two-body coulomb energies")')
  write(6,'(4f25.16)')  (eval(i), i  = 2*norb,1,-1)
  
  
  
  stop
end program twobody_test


subroutine two_particle_hamiltonian(norb,ham)
  implicit none
  integer, parameter :: dp = kind(1.0d0)
  integer, intent(in) :: norb
  complex(kind=dp) , intent(out) :: ham(2*norb, 2*norb)
  !-------------------------------
  ! returns the complex norb**2 x norb**2 Coulomb
  ! hamiltonian for two particles without neutralising
  ! background
  !----------------------------
  real(kind=dp) ::  shift
  complex(kind=dp) :: landau_coulomb
  real(kind=dp) :: madelung
  integer, allocatable :: basis2(:,:)
  
  ! if norb is odd, keep mod(k1 + k2, norb) = 0 or 1
  !  this has two copies of norb distinct eigenvalues
  ! if norb there 2*norb distinct eigenvalues
  
  integer :: k1,k2,k3,k4,i,j, count
  
  allocate (basis2(2,2*norb))
  count = 0
  do k1 = 1,norb
     do k2 = 1,norb
        if(modulo(k1 + k2, norb) > 1) cycle
        count = count + 1
        if(count > 2*norb) then
           write(6,'("count error")')
           stop
        endif
        basis2(1,count) = k1
        basis2(2,count) = k2
     enddo
  enddo
  if(count /= 2*norb) then
     write(6,'("error",2i5)') count ,2*norb
     stop
  endif
  
  
  ham = cmplx(0,kind=dp)
  do i = 1,count
     do j = 1,count
        k1 = basis2(1,i)
        k2 = basis2(2,i)
        k4 = basis2(1,j)
        k3 = basis2(2,j)
        ham(i,j) = landau_coulomb(k1,k2,k3,k4)
     enddo
  enddo
  deallocate (basis2)
  ! remove background neutralization  (to epose pseudopotentials)
  ! this is not the usual use of madelung:
  ! in the neutralized system, the background  energy is
  !  nel * madeling()/2

  shift = -madelung()
  write(6,'("background shift",f25.10)') shift
  forall (i = 1:2*norb) ham(i,i) = ham(i,i) +  cmplx(shift,kind=dp)
  return
end subroutine two_particle_hamiltonian



subroutine diagonalize_hermitian(n,matrix,eval)
  implicit none
  integer, parameter :: dp = kind(1.0d0)
  integer, intent(in) :: n
  complex(kind=dp), intent(in) :: matrix(n,n)
  real(kind=dp), intent(out) :: eval(n)

  complex (kind=dp), allocatable :: evec(:,:), work(:)
  real(kind=dp), allocatable  :: rwork(:)
  integer:: lwork, info, i
  
  if (n < 1) stop
  if (n == 1) then
     eval(1) = matrix(1,1)
     return
  endif
  
  lwork = 2*n
  allocate(evec(n,n),work(lwork), rwork(3*n))
  
  evec = matrix
  call zheev('n','u',n,evec,n,eval,work,lwork,rwork,info)
  if(info > 0) then
     write(6,'("ZHEEV: INFO =",i12)') info
     stop
  endif
  deallocate(evec,work,rwork)
  return
end subroutine diagonalize_hermitian
