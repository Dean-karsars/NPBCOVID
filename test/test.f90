program check
    use mod_events, only:event, rand_int, &
                         Delete_Event, Update_Event, Insert_Event, sample, gibbs_sampling, &
                         Update_Gaus, Update_shape_scale
    use iso_fortran_env, only:real32, int32, real64
    use mod_prior, only:prior
    use functional, only: subscript, put => insert, arange, last, sort
    use stdlib_stats_distribution_uniform, only: random_uniform => rvs_uniform
    use stdlib_stats, only: mean
    use mod_matrix
    use forbear, only : bar_object

    implicit none
    
    real(real32) :: time(10), new
    integer(int32) :: outcome(10)
    integer(int32) :: group(10), loop_group
    integer(int32) :: N(2) = [50, 50], i, j
    integer(int32), parameter :: step(7) = [-1, 3, 5, 6, 7, 8, 9]
    type(event) :: test
    real(real32) :: shape_scale(2) = [1.0, 0.2], variance_scale = 2
    real(real32), allocatable :: sample_vector(:)
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
    type(bar_object) :: bar
    call bar%initialize(filled_char_string='+', prefix_string='progress |', suffix_string='| ', add_progress_percent=.true., &
                        max_value = real(1000000, real64))
    call bar%start

    p_prior = prior(beta = beta, p_gamma = p_gamma, gamma = gamma, eta = eta, &
                    alpha = alpha, tau = tau, omega = omega, p_theta = p_theta, theta = theta)

    time = sort(random_uniform(-10.0, 10.0, 10))
    outcome = [1, 1, 1, 1, 2, 3, 1, 1, 1, 11]                
    group = [1, 2, 1, 2, 1, 1, 1, 1, 1, 0]
    test = event(time, outcome, group, N, 1, p_prior, variance_scale, shape_scale)
    call test%output_time()
    print *, "********************************************************************************************************************"
    call test%output_parms()
    print * , "********************************************************************************************************************"
    allocate(result(1000000))
    do i = 1, 1000000
      call bar%update(real(i, real64))
      !group
      loop_group = random_uniform(1, test%num_groups - 1)
      !Infection Insert or Update
      if ( random_uniform(1) == 0 ) then
        ! Insert
        new = random_uniform(test%time(1), last(test%time) - test%time(1))
        test = Insert_Event(new, 1, loop_group, test)
      else 
        ! Delete
        mask = test%outcome == 1 .and. test%group == loop_group .and. test%time /= test%time(1)
        sample_vector = pack(test%time, mask)
        if(size(sample_vector) >= 1) then
          new = sample(sample_vector)
          test = Delete_Event(new, 1, loop_group, test)
        end if
      end if
   
      !Infection Update
      new = random_uniform(test%time(1), last(test%time) - test%time(1))
      mask = test%outcome == 1 .and. test%group == loop_group .and. test%time /= test%time(1)
      sample_vector = pack(test%time, mask)
      if(size(sample_vector) >= 1) then
        test = Update_Event(sample(sample_vector), new, 1, loop_group, test)
      end if

      do j = 1, 7 !M-H for other events

        !Insert or Update
        if ( random_uniform(1) == 0 ) then
          ! Insert
          mask = test%outcome == 1 .and. test%group == loop_group !Find Ii1
          sample_vector = pack(test%time, mask)
          if(size(sample_vector) < 1) cycle
          new = random_uniform(minval(pack(test%time, mask)), last(test%time) - minval(pack(test%time, mask))) !sample from[Ii1, T]
          test = Insert_Event(new, step(j), loop_group, test)
        else 
          ! Delete
          mask = test%outcome == step(j) .and. test%group == loop_group .and. test%time /= test%time(1)
          sample_vector = pack(test%time, mask)
          if(size(sample_vector) < 1) cycle
          new = sample(sample_vector)
          test = Delete_Event(new, step(j), loop_group, test)
        end if

        !Update
        mask = test%outcome == 1 .and. test%group == loop_group !Find Ii1
        new = random_uniform(minval(pack(test%time, mask)), last(test%time) - minval(pack(test%time, mask))) !sample from[Ii1, T]
        mask = test%outcome == step(j) .and. test%group == loop_group .and. test%time /= test%time(1)
        sample_vector = pack(test%time, mask)
        if(size(sample_vector) < 1) cycle
        test = Update_Event(sample(sample_vector), new, step(j), loop_group, test)

      end do

      test = Update_Gaus(test, loop_group)
      test = Update_shape_scale(1.0, loop_group, test)

      call gibbs_sampling(test)

      result(i) = test%beta(1)
    end do
    
    print * , "FINALL:", mean(result(5000:1000000))
    call test%output_parms
    call test%output_time
  end program check