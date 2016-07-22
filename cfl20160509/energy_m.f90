
module energy_m
  integer, parameter :: dp3 = kind(1.0d0)
  real (kind=dp3), allocatable :: energy_data(:,:),coulomb(:,:)
  real(kind=dp3) :: emadelung
  real(kind=dp3) :: energy(2), variance, mean
  integer :: norb_e
  integer, allocatable :: x_prev1(:,:)
end module energy_m


  subroutine  energy_setup(nel,norb)
    use energy_m
    implicit none
    integer, intent(in) :: nel, norb
    
    allocate (x_prev1(2,nel))

    norb_e = norb
    call coulomb_setup
    return
  end subroutine energy_setup


  subroutine  energy_data_init(nel,norb,n_unit)
    use energy_m
    implicit none
    integer, parameter :: dp = kind(1.0d0)
    integer, intent(in) :: nel, norb, n_unit
    integer :: i,j
    real (kind=dp3) :: coulomb_table

    allocate (energy_data(2,n_unit))
    energy_data = real(0,kind=dp)


    allocate (coulomb(norb,norb))
    do i = 1,norb
       do j = 1,norb
          coulomb(i,j) = coulomb_table(norb,i,j)
       enddo
    enddo
    emadelung = -nel*coulomb(norb,norb)/(2*norb)
    return
  end subroutine energy_data_init

  subroutine energy_data_combine(unit,n_steps)
    use energy_m
    implicit none
    integer, intent(in) :: n_steps, unit
    energy_data(:,unit) = energy_data(:,unit)/n_steps

    variance = energy_data(2,unit) - energy_data(1,unit)**2
    write(12,'(3f25.16)')  energy_data(1,unit),energy_data(2,unit),variance

    energy(1) = sum(energy_data(1,1:unit)/unit)
    energy(2) = sum(energy_data(2,1:unit)/unit)
    variance = energy(2) - energy(1)**2
    write(6,'("energy per particle:",f25.16, &
          &" variance:",f25.16)') energy(1),variance
    return
  end subroutine energy_data_combine

  subroutine  energy_sample(nel,x,changed_x,unit)
    use energy_m
    implicit none
    integer,parameter :: dp = kind(1.0d0)
    integer, intent(in) :: nel,unit
    integer, intent(in) :: x(2,nel)
    logical, intent(in) :: changed_x(nel)

    real(kind=dp), save :: uc
    integer :: i,j, mn(2)
    if (all(changed_x)) then
       uc = emadelung
       do i = 1,nel
          do j = 1,i-1
             mn(:)   = 1 + modulo(x(:,i) - x(:,j) - 1,norb_e)
             uc = uc + coulomb(mn(1),mn(2))
          enddo
       enddo
       x_prev1 = x
    else
       do i = 1,nel
          do j = 1,i-1
             if (changed_x(i) .or. changed_x(j)) then
                mn = 1 + modulo (x(:,i) - x(:,j) -1, norb_e)
                uc = uc + coulomb(mn(1), mn(2))
                mn = 1 + modulo (x_prev1(:,i) - x_prev1(:,j) -1, norb_e)
                uc = uc - coulomb(mn(1), mn(2))
             endif
          enddo
       enddo
       do i = 1,nel
          if( changed_x(i)) x_prev1(:,i) = x(:,i)
       enddo
    end if
    if (unit /= 0) then
       energy_data(1,unit) = energy_data(1,unit) + (uc/nel)
       energy_data(2,unit) = energy_data(2,unit) + (uc/nel)**2
    endif
    return
  end subroutine energy_sample

  subroutine energy_data_finalize(n_unit)
    use energy_m
    implicit none
    integer, intent(in) :: n_unit
    energy(1) = sum(energy_data(1,:))/n_unit
    energy(2) = sum(energy_data(2,:))/n_unit
    variance = energy(2) - energy(1)**2
    write(13,'(2f25.16)') energy(1),energy(2),variance
    write(6,'("energy per particle:",f25.16, &
          &" variance:",f25.16)') energy(1),variance
    return
  end subroutine energy_data_finalize

