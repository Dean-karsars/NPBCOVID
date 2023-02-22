program main
    use mod_events, only:event, rand_int, &
                         Delete_Event, Update_Event, Insert_Event, sample, gibbs_sampling, &
                         Update_Gaus
    use iso_fortran_env, only:real32, int32
    use mod_prior, only:prior
    use functional, only: subscript, put => insert, arange, last, sort
    use stdlib_stats_distribution_uniform, only: random_uniform => rvs_uniform
    use stdlib_stats, only: mean

    implicit none
    
    real(real32) :: time(10), new
    integer(int32) :: outcome(10)
    integer(int32) :: group(10)
    integer(int32) :: N(2) = [50, 50], i
    type(event) :: test
    real(real32) :: shape_scale(2) = [0.2, 0.2], variance_scale = 2
    !priors
    real(real32) :: beta(2, 2) = reshape([1.0, 2.0, 3.0, 4.0], [2, 2]), &
                    p_gamma(2, 2) = reshape([1.0, 2.0, 3.0, 4.0], [2, 2]), &
                    gamma(2, 2) = reshape([1.0, 2.0, 3.0, 4.0], [2, 2]), &
                    alpha(2, 2) = reshape([1.0, 2.0, 3.0, 4.0], [2, 2]), & 
                    tau(2, 2) = reshape([1.0, 2.0, 3.0, 4.0], [2, 2]), &  
                    omega(2, 2) = reshape([1.0, 2.0, 3.0, 4.0], [2, 2]), &
                    eta(2, 2) = reshape([1.0, 2.0, 3.0, 4.0], [2, 2])

    real(real32) :: p_theta(2) = [1, 1], theta(2) = [1, 1]
    type(prior) :: p_prior
    logical, allocatable :: mask(:)
    real(real32), allocatable :: result(:)
    p_prior = prior(beta = beta, p_gamma = p_gamma, gamma = gamma, eta = eta, &
                    alpha = alpha, tau = tau, omega = omega, p_theta = p_theta, theta = theta)

    time = sort(random_uniform(-10.0, 10.0, 10))
    outcome = [1, 1, 1, 1, 2, 3, -1, -1, 1, 11]                
    group = [1, 2, 1, 2, 1, 1, 1, 1, 1, 0]
    test = event(time, outcome, group, N, 1, p_prior, variance_scale, shape_scale)
    call test%output_time()
    print *, "********************************************************************************************************************"
    call test%output_parms()
    print * , "********************************************************************************************************************"
    allocate(result(1000000))
    do i = 1, 1000000
      new = random_uniform(test%time(1), last(test%time) - test%time(1))
      test = Insert_Event(new, 1, 1, test)
      mask = test%outcome == 1 .and. test%group == 1 .and. test%time /= test%time(1)
      new = sample(pack(test%time, mask))
      test = Delete_Event(new, 1, 1, test)
      new = random_uniform(test%time(1), last(test%time) - test%time(1))
      mask = test%outcome == 1 .and. test%group == 1 .and. test%time /= test%time(1)
      test = Update_Event(sample(pack(test%time, mask)), new, 1, 1, test)
      test = Update_Gaus(test, 1)
      call gibbs_sampling(test)
      result(i) = test%beta(1)
    end do
    
    print * , "FINALL:", mean(result(5000:1000000))
    call test%output_time
  end program main