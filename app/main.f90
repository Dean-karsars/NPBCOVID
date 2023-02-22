program main
  use iso_fortran_env, only:real32, int32
  use simple_rnd, only: random_mvn => mnorm_smp
  use stdlib_stats, only:var
  use mod_matrix
  implicit none
  real(real32) :: cov22(2, 2) = reshape([2.0, 1.0, 0.0, 1.0], [2, 2]), &
                  mean(2) = [-113.012342, 200.022], res(100000, 2), &
                  AAA(2, 2) = reshape([2, 0, 0, 2], [2, 2]), &
                  xxx(2, 2), ah
  integer :: i, info
  call mat_mul(2, cov22, AAA, xxx)
  do i = 1, 2
    write(*,'(*(f12.6,3x))') cov22(i,:)
  end do
  
  !call inv(2, cov22, info)
  

  call det(cov22, ah, info)
  print * , "det = ", ah
  
end program
