program check
    use mod_events, only:event, rand_int, &
                         Delete_Event, Update_Event, Insert_Event, sample, gibbs_sampling, &
                         Update_Gaus, Update_shape_scale
    use iso_fortran_env, only:real32, int32, real64
    use ieee_arithmetic
    use mod_prior, only:prior
    use functional, only: subscript, put => insert, arange, last, sort
    use stdlib_stats_distribution_uniform, only: random_uniform => rvs_uniform
    use stdlib_stats_distribution_exponential, only: rexp => rvs_exp
    use stdlib_stats, only: mean
    use forbear, only : bar_object
    use mod_io, only :write_sim_data


    implicit none
      integer(int32), parameter :: num_groups = 3, N(3)= [10000, 10000, 1000]
  
      real(real32), parameter :: Imin = 0, Time_Last = 50
      real(real32), parameter :: beta(num_groups) = [0.2, 0.1, 0.01]
      real(real32), parameter :: p_gamma(num_groups) = [0.3, 0.2, 0.1]
      real(real32), parameter :: gamma(num_groups) = [0.1, 0.2, 0.3]
      real(real32), parameter :: eta(num_groups) = [0.3, 0.2, 0.1]
      real(real32), parameter :: tau(num_groups) = [0.1, 0.2, 0.3]
      real(real32), parameter :: alpha(num_groups) = [0.1, 0.2, 0.3]
      real(real32), parameter :: omega(num_groups) = [0.003, 0.002, 0.0001]
      real(real32), parameter :: theta(num_groups) = 0.2
      real(real32), parameter :: p_theta(num_groups) = 0.6
      real(real32), parameter :: beta_lambda(num_groups) = [0.003, 0.002, 0.001]

      integer(int32) :: i, current_event(2), vac_group
      integer(int32) :: x(num_groups) = N, y(num_groups) = 0, z(num_groups) = 0, &
                        u(num_groups) = 0, v(num_groups) = 0, w(num_groups) = 0
      
      real(real32) :: t = Imin, time_interval, vaccination_parms(2) = [0.05, 0.01], vaccination_rate(2), vac_del_vector(2)
      real(real32) :: rate_matrix(num_groups, 9), delta_t_matrix(num_groups, 9), delta_t
      real(real32), allocatable :: stored_time(:)
      integer(int32), allocatable :: stored_outcome(:), stored_group(:)

      type(bar_object) :: bar

      stored_time = [real(real32) ::]
      stored_outcome = [integer(int32) ::]
      stored_group = [integer(int32) ::]
      call bar%initialize(filled_char_string='+', prefix_string='progress |', suffix_string='| ', &
                          add_progress_percent=.true., max_value = real(Time_Last, real64),       &
                          min_value = real(Imin, real64), add_date_time = .true., add_progress_speed = .true.)
      call bar%start
      y(1) = 1
      x(1) = x(1) - 1
      do while(t < Time_Last)
        call bar%update(current = real(t, real64))
        !Calculate rate
        rate_matrix(:, 1) = beta_t(t, beta_lambda) * x * sum(y + z + u) !Infection
        rate_matrix(:, 2) = p_theta * theta * y                         !Diagnosis
        rate_matrix(:, 3) = (1 - p_theta) * theta * y                   !Undetected
        rate_matrix(:, 4) = p_gamma * gamma * z                         !Hospitalized
        rate_matrix(:, 5) = (1 - p_gamma) * gamma * z                   !Recovery of Diagnosis
        rate_matrix(:, 9) = omega * z                                   !Death
        rate_matrix(:, 6) = eta * u                                     !Recovery of Undetected
        rate_matrix(:, 7) = tau * w                                     !Suspect of Hospitalized 
        rate_matrix(:, 8) = alpha * v                                   !Suspect of Recovery
        vaccination_rate = vaccination_parms * x(1:2)                        !Rate of Vaccination
        !sample Time by using rate matrix
        delta_t_matrix = rexp_matrix(rate_matrix)
        delta_t = minval(delta_t_matrix)
        current_event = minloc(delta_t_matrix)
        t = t + delta_t
        !Vaccination Event
        vac_del_vector = ieee_value(vac_del_vector, ieee_positive_inf)
        if(vaccination_rate(1) > 0) vac_del_vector(1) = rexp(vaccination_rate(1))
        if(vaccination_rate(2) > 0) vac_del_vector(2) = rexp(vaccination_rate(2))
        vac_group = minloc(vac_del_vector, 1)
        
        if(delta_t > vac_del_vector(vac_group)) then
          x(vac_group) = x(vac_group) - 1
          x(vac_group + 1) = x(vac_group + 1) + 1
          t = t - delta_t + vac_del_vector(vac_group)
          stored_time = [stored_time, t]
          stored_outcome = [stored_outcome, 10]
          stored_group = [stored_group, vac_group]
        else
          stored_time = [stored_time, t]
          stored_outcome = [stored_outcome, current_event(2)]
          stored_group = [stored_group, current_event(1)]
          select case(current_event(2))
            case(1) !Infection
              x(current_event(1)) = x(current_event(1)) - 1
              y(current_event(1)) = y(current_event(1)) + 1
            case(2) !Diagnosis
              y(current_event(1)) = y(current_event(1)) - 1
              z(current_event(1)) = z(current_event(1)) + 1
            case(3) !Undetected
              y(current_event(1)) = y(current_event(1)) - 1
              u(current_event(1)) = u(current_event(1)) + 1
            case(4) !Hospitalized
              w(current_event(1)) = w(current_event(1)) + 1
              z(current_event(1)) = z(current_event(1)) - 1
            case(5) !Recovery of Diagnosis
              v(current_event(1)) = v(current_event(1)) + 1
              z(current_event(1)) = z(current_event(1)) - 1
            case(9) !Death
              z(current_event(1)) = z(current_event(1)) - 1
            case(6) !Recovery of Undetected
              u(current_event(1)) = u(current_event(1)) - 1
              v(current_event(1)) = v(current_event(1)) + 1
            case(7) !Suspect of Hospitalized
              w(current_event(1)) = w(current_event(1)) - 1
              x(current_event(1)) = x(current_event(1)) + 1
            case(8) !Suspect of Recovery
              v(current_event(1)) = v(current_event(1)) - 1
              x(current_event(1)) = x(current_event(1)) + 1
            case default 
              error stop "case out of bound"
          end select
        end if
      end do
      stored_time = [stored_time, Time_Last]
      stored_outcome = [stored_outcome, 11]
      stored_group = [stored_group, 0]
      call write_sim_data("data/simulation.csv", stored_time, stored_outcome, stored_group)

  contains
    pure elemental function beta_t(time, lambda) result(res)
      real(real32), intent(in) :: time
      real(real32), intent(in) :: lambda
      real(real32) :: res
      
      
      res = exp(-time /10) * lambda
      
    end function beta_t
    
    function rexp_matrix(input) result(res)
      real(real32), intent(in) :: input(:, :)
      real(real32) :: res(size(input, 1),size(input, 2))
      integer(int32) :: row, column
    
      do row = 1, size(input, 1)
        do column = 1, size(input, 2)
          if(input(row, column) > 0) then 
            res(row, column) = rexp(input(row, column))
          else 
            res(row, column) = ieee_value(res(row, column), ieee_positive_inf)
          end if
        end do
      end do
    end function rexp_matrix

  end program check