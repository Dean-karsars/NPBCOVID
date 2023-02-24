module mod_data
    use iso_fortran_env, only:real32, int32, real64
    use mod_prior, only:prior
    implicit none
    real(real32), parameter :: beta(1, 2) = reshape([1.0, 2.0, 3.0, 4.0], [1, 2]), &
                               p_gamma(1, 2) = reshape([1.0, 2.0, 3.0, 4.0], [1, 2]), &
                               gamma(1, 2) = reshape([1.0, 2.0, 3.0, 4.0], [1, 2]), &
                               alpha(1, 2) = reshape([1.0, 2.0, 3.0, 4.0], [1, 2]), & 
                               tau(1, 2) = reshape([1.0, 2.0, 3.0, 4.0], [1, 2]), &  
                               omega(1, 2) = reshape([1.0, 2.0, 3.0, 4.0], [1, 2]), &
                               eta(1, 2) = reshape([1.0, 2.0, 3.0, 4.0], [1, 2])

    real(real32) :: p_theta(2) = [1, 1], theta(2) = [1, 1]
    type(prior) :: p_prior
    integer(int32), parameter :: N(1) = [100]
    real(real32) :: shape_scale(1) = [0.0006], variance_scale = 20
    
end module mod_data