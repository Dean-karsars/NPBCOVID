program main
  use iso_fortran_env, only:real32, int32
  use mod_events, only:event, Delete_Event, Update_Event, Insert_Event, sample, gibbs_sampling, &
                       Update_Gaus, Update_shape_scale
  use mod_io, only:read_event_data, write_event_data
  use mod_prior, only:prior
  use mod_data
  use forbear, only : bar_object
  use stdlib_stats_distribution_uniform, only: random_uniform => rvs_uniform
  use functional, only: subscript, put => insert, arange, last, sort
  use stdlib_stats, only: mean

  implicit none
    !Data Readin
    type(event) :: sim
    character(*), parameter :: filename = "data/Fortran.csv"
    real(real32), allocatable :: time(:)
    integer(int32), allocatable :: outcome(:)
    integer(int32), allocatable :: group(:)
    !Loop Variables
    real(real32), allocatable :: sample_vector(:)
    logical, allocatable :: mask(:) 
    type(bar_object) :: bar
    integer(int32) :: loop_group, i, j
    integer(int32), parameter :: step(7) = [-1, 3, 5, 6, 7, 8, 9]
    real(real32) :: new
    integer(int32), parameter :: n_iter = 10

    call read_event_data(filename, time, outcome, group)
    
    p_prior = prior(beta = beta, p_gamma = p_gamma, gamma = gamma, eta = eta, &
                    alpha = alpha, tau = tau, omega = omega, p_theta = p_theta, theta = theta)           
    sim = event(time, outcome, group, N, 1, p_prior, variance_scale, shape_scale)
    
    call sim%output_time()
    print *, "********************************************************************************************************************"
    call sim%output_parms()
    print * , "********************************************************************************************************************"
    call bar%initialize(filled_char_string='+', prefix_string='progress |', suffix_string='| ', add_progress_percent=.true., &
                        max_value = real(n_iter, real64), add_date_time = .true., add_progress_speed = .true.)
    call bar%start
    allocate(beta_sample(n_iter, sim%num_groups), pgamma_sample(n_iter, sim%num_groups), gamma_sample(n_iter, sim%num_groups), &
             eta_sample(n_iter, sim%num_groups), alpha_sample(n_iter, sim%num_groups), tau_sample(n_iter, sim%num_groups), &
             omega_sample(n_iter, sim%num_groups), ptheta_sample(n_iter, 1), theta_sample(n_iter,1))
    do i = 1, n_iter
      call bar%update(real(i, real64))
      !group
      if(sim%num_groups > 1) then
        loop_group = random_uniform(1, sim%num_groups - 1)
      else 
        loop_group = 1
      end if
      
      !Infection Insert or Update
      if ( random_uniform(1) == 0 ) then
        ! Insert
        new = random_uniform(sim%time(1), last(sim%time) - sim%time(1))
        sim = Insert_Event(new, 1, loop_group, sim)
      else 
        ! Delete
        mask = sim%outcome == 1 .and. sim%group == loop_group .and. sim%time /= sim%time(1)
        sample_vector = pack(sim%time, mask)
        if(size(sample_vector) >= 1) then
          new = sample(sample_vector)
          sim = Delete_Event(new, 1, loop_group, sim)
        end if
      end if
   
      !Infection Update
      new = random_uniform(sim%time(1), last(sim%time) - sim%time(1))
      mask = sim%outcome == 1 .and. sim%group == loop_group .and. sim%time /= sim%time(1)
      sample_vector = pack(sim%time, mask)
      if(size(sample_vector) >= 1) then
        sim = Update_Event(sample(sample_vector), new, 1, loop_group, sim)
      end if

      do j = 1, 7 !M-H for other events
        !Insert or Update
        if ( random_uniform(1) == 0 ) then
          ! Insert
          mask = sim%outcome == 1 .and. sim%group == loop_group !Find Ii1
          sample_vector = pack(sim%time, mask)
          if(size(sample_vector) >= 1) then
            new = random_uniform(minval(pack(sim%time, mask)), last(sim%time) - minval(pack(sim%time, mask))) !sample from[Ii1, T]
            sim = Insert_Event(new, step(j), loop_group, sim)
          end if
        else 
          ! Delete
          mask = sim%outcome == step(j) .and. sim%group == loop_group .and. sim%time /= sim%time(1)
          sample_vector = pack(sim%time, mask)
          if(step(j) == -1) then
            if(size(sample_vector) > 1) then
              new = sample(sample_vector)
              sim = Delete_Event(new, step(j), loop_group, sim)
            end if
          else  !若为Thinned事件, 则size需>1，否则 >=1
            if(size(sample_vector) >= 1) then
              new = sample(sample_vector)
              sim = Delete_Event(new, step(j), loop_group, sim)
            end if 
          end if

        end if

        !Update
        mask = sim%outcome == 1 .and. sim%group == loop_group !Find Ii1
        new = random_uniform(minval(pack(sim%time, mask)), last(sim%time) - minval(pack(sim%time, mask))) !sample from[Ii1, T]
        mask = sim%outcome == step(j) .and. sim%group == loop_group .and. sim%time /= sim%time(1)
        sample_vector = pack(sim%time, mask)
        if(size(sample_vector) >= 1) then
          sim = Update_Event(sample(sample_vector), new, step(j), loop_group, sim)
        end if
      end do

      sim = Update_Gaus(sim, loop_group)
      sim = Update_shape_scale(10.0, loop_group, sim)

      call gibbs_sampling(sim)

      beta_sample(i, :) = sim%beta
      pgamma_sample(i, :) = sim%p_gamma
      gamma_sample(i, :) = sim%gamma
      alpha_sample(i, :) = sim%alpha
      tau_sample(i, :) = sim%tau
      omega_sample(i, :) = sim%omega
      eta_sample(i, :) = sim%eta
      ptheta_sample(i, :) = sim%p_theta
      theta_sample(i, :) = sim%theta
      !call sim%output_time
    end do
    print * , "FINALL:", mean(beta_sample(n_iter/2:n_iter, 1))
    call sim%output_parms
    call sim%output_time
    print * , "Shape:", sim%shape_scale
    call write_event_data("data/beta_sample.csv", beta_sample)
    call write_event_data("data/pgamma_sample.csv", pgamma_sample)
    call write_event_data("data/gamma_sample.csv", gamma_sample)
    call write_event_data("data/alpha_sample.csv", alpha_sample)
    call write_event_data("data/tau_sample.csv", tau_sample)
    call write_event_data("data/omega_sample.csv", omega_sample)
    call write_event_data("data/eta_sample.csv", eta_sample)
    call write_event_data("data/ptheta_sample.csv", ptheta_sample)
    call write_event_data("data/theta_sample.csv", theta_sample)
end program main
