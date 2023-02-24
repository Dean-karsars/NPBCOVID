!------------------------------------------------------------------------------
! China Pharmaceutical University-Department of Biostatistics and Computational Pharmacy
!------------------------------------------------------------------------------
!
! MODULE:  mod_events
!
!> @author
!> Yao Gongzheng}
!
! DESCRIPTION: 
!>  Define the derived types: event
!>  "event" is an array of time points for a sequence of events
!------------------------------------------------------------------------------
module mod_events
  use iso_fortran_env, only:real32, int32, real64
  use ieee_arithmetic
  use functional, only: subscript, put => insert, arange, last
  use mod_prior, only:prior
  use m_mrgrnk, only:sort_index => mrgrnk
  use random_beta_mod, only:random_beta
  use random_gamma_mod, only:random_gamma
  use stdlib_stats_distribution_exponential, only: random_exp => rvs_exp
  use stdlib_stats_distribution_uniform, only: random_uniform => rvs_uniform
  use stdlib_stats_distribution_normal, only: random_uvn => rvs_normal
  use stdlib_array, only: trueloc
  use simple_rnd, only: random_mvn => mnorm_smp
  use stdlib_sorting, only: ord_sort, std_sort => sort
  use stdlib_selection, only:select
  use stdlib_stats_distribution_uniform, only: shuffle
  use mod_matrix, only:inv, det

  implicit none
  
  private
  public :: event, log_likelihood, Gauss_Pred, &
            Insert_Event, rand_int, Delete_Event, &
            Update_Event, sample, gibbs_sampling, Update_Gaus, Update_shape_scale
  
  type :: event
    real(real32), allocatable :: time(:)
    integer(int32), allocatable :: outcome(:) !!"-1:Thinned", "1:Infection", "2:Diagnosis", "3:Undetected", 
                                              !!"4:Hospitalized", "5:Diagnosised to Recovered", "6:Undetected to Recovered"
                                              !!"7:Hospitalsized to Suspect", "8:Recovered to Suspect", "9:Death", "10:Vaccination"
                                              !!"11:Time Last
    integer(int32), allocatable :: group(:)
    integer(int32) :: num_groups, I_min_group
    integer(int32), allocatable :: N(:)
    
    real(real32) :: log_likelihood
    type(prior) :: priors
    real(real32), allocatable :: beta(:), p_gamma(:), gamma(:), eta(:), alpha(:), tau(:), omega(:)
    real(real32) :: p_theta, theta
    real(real32), allocatable :: gaus(:)
    real(real32), allocatable :: shape_scale(:)
    real(real32) :: variance_scale
    integer(int32) :: MPA_Level = 1000
  contains
    procedure, pass(self) :: output_time
    procedure, pass(self) :: output_parms
  end type event

  interface event
    module procedure event_constructor
  end interface

  interface sort
    module procedure event_sort
  end interface

contains
  function event_constructor(time, &
                             outcome, &
                             group, &
                             pop, &
                             I_min_group, &
                             priors, &
                             variance_scale, &
                             shape_scale) result(y)
    
    real(real32), intent(in) :: time(:)
    integer(int32), intent(in) :: outcome(:)
    integer(int32), intent(in) :: group(:)
    integer(int32), intent(in) :: pop(:)
    integer(int32), intent(in) :: I_min_group
    real(real32), intent(in) :: shape_scale(:)
    real(real32), intent(in) :: variance_scale
    type(prior), intent(in) :: priors
    type(event) :: y
    
    !intermediate variables
    integer(int32) :: i
    integer(int32), allocatable :: Infec_loc(:) !Index of Thinned or Infection Event
    integer(int32), allocatable :: idx(:)
    real(real32), allocatable :: cov(:, :)
    real(real32), allocatable :: mean(:)

    real(real32) :: res
    integer :: info
    !automatic reallocate
    y%time = [time] 
    y%outcome = [outcome] 
    y%group = [group] 
    y%gaus = [y%time]
    y%gaus = ieee_value(y%gaus, class = ieee_quiet_nan)

  
    !sort event object according to its time attribute
    call sort(y)
    y%num_groups = size(pop)
    y%N = pop
    


    if(I_min_group < 1 .or. I_min_group > y%num_groups) error stop "I min group must between 1 and num_groups"

    y%I_min_group = I_min_group

    !Gibbs Initialization
    y%priors = priors
    allocate(y%alpha(y%num_groups),   &
             y%beta(y%num_groups),    &
             y%eta(y%num_groups),     &
             y%gamma(y%num_groups),   &
             y%p_gamma(y%num_groups), &
             y%tau(y%num_groups),     &
             y%omega(y%num_groups))

    call gibbs_sampling(y)

    
    !Gaus Initialization
    y%variance_scale = variance_scale
    y%shape_scale = shape_scale
    
    do i = 1, y%num_groups
      Infec_loc = trueloc((y%outcome == 1 .or. y%outcome == -1) .and. y%group == i .and. y%time /= minval(y%time)) !exclude I min
      if(size(Infec_loc) <= y%MPA_Level) then 
        cov = cov_mat(y%time(Infec_loc), y%variance_scale, y%shape_scale(i)) !construct covariance matrix
        !if(i == 1) print *, "group1 covariance matrix size is", size(cov, 1)
        allocate(mean, mold = y%time(Infec_loc)) 
        mean = 0
        !call det(cov, res, info)
        !print * , "group",  i, "det is",  res
        !if(info /= 0) then
          !stop
        !end if
        y%gaus(Infec_loc) = random_mvn(cov, size(Infec_loc), mean) !assign sampled gaus value to corresponding location
        deallocate(mean)
      
      else !MPA
        idx = rand_int(1, size(Infec_loc), y%MPA_Level)
        call std_sort(idx)
        y%gaus(Infec_loc) = MPA_sample(y%time(Infec_loc(idx)), y%time, y%variance_scale, y%shape_scale(i))
      end if

    end do  
    y%log_likelihood = log_likelihood(y)
    print * , y%log_likelihood
  end function event_constructor

  subroutine event_sort(x)
    type(event), intent(inout) :: x

    integer(int32) :: sorted_index(size(x%time))
    call sort_index(x%time, sorted_index)

    x%time = subscript(x%time, sorted_index)
    x%outcome = subscript(x%outcome, sorted_index)
    x%group = subscript(x%group, sorted_index)
    x%gaus = subscript(x%gaus, sorted_index)
  end subroutine event_sort

  pure function log_likelihood(input) result(res)
    !!Evaluate log_likelihood of an event object
    type(event), intent(in) :: input
    real(real32) :: res
    !Intermediate variables
    integer(int32) :: i
    integer(int32) :: x(input%num_groups), y(input%num_groups), z(input%num_groups), &
                      u(input%num_groups), v(input%num_groups), w(input%num_groups)

    logical :: mask(input%num_groups)                  
    !Initialize variables                                
    x = input%N  !易感者
    y = 0        !感染者
    y(input%I_min_group) = 1                         !0号感染者
    x(input%I_min_group) = x(input%I_min_group) - 1  !0号感染者
    z = 0        !确诊者
    u = 0        !未确诊
    v = 0        !住院
    w = 0        !康复
    res = 0
    

    do i = 2, size(input%time) !注意左极限, 先算likelihood, 在更改人数参数
      select case(input%outcome(i))
        case(-1)!Thinned
          res = res + log(input%beta(input%group(i))) + &
                log(real(x(input%group(i)), real32)) + &
                log(real(sum(y) + sum(z) + sum(u), real32)) &
                + log(logit(-1.0 * input%gaus(i)))
          ! Integral Terms
          res = res                                                 -&
                (input%time(i) - input%time(i - 1))                 *& !Time Span
                sum((input%beta * real(x * sum(y + z + u), real32)) +& !Term 1
                    (input%theta * real(y, real32))                 +& !2
                    (input%gamma * real(z, real32))                 +& !3
                    (input%eta * real(u, real32))                   +& !4
                    (input%tau * real(v, real32))                   +& !5
                    (input%alpha * real(w, real32))                 +& !6
                    (input%omega * real(z, real32)))                   !7
                    
          mask = x < 0 .or. y < 0 .or. z < 0 .or. w < 0 .or. u < 0 .or. v < 0
          if(any(mask)) then 
            res =  ieee_value(res, ieee_negative_inf) 
            return
          end if
        case(1) !Infection
          res = res + log(input%beta(input%group(i))) + &
                log(real(x(input%group(i)), real32)) + &
                log(real(sum(y) + sum(z) + sum(u), real32)) &
                + log(logit(input%gaus(i)))
          ! Integral Terms      
          res = res                                                 -&
                (input%time(i) - input%time(i - 1))                 *& !Time Span
                sum((input%beta * real(x * sum(y + z + u), real32)) +& !Term 1
                    (input%theta * real(y, real32))                 +& !2
                    (input%gamma * real(z, real32))                 +& !3
                    (input%eta * real(u, real32))                   +& !4
                    (input%tau * real(v, real32))                   +& !5
                    (input%alpha * real(w, real32))                 +& !6
                    (input%omega * real(z, real32)))                   !7

          x(input%group(i)) = x(input%group(i)) - 1
          y(input%group(i)) = y(input%group(i)) + 1
                              
          mask = x < 0 .or. y < 0 .or. z < 0 .or. w < 0 .or. u < 0 .or. v < 0
          if(any(mask)) then 
            res =  ieee_value(res, ieee_negative_inf) 
            return
          end if
        case(2) !Diagnosis
          res = res + log(input%theta) + log(input%p_theta) + log(real(y(input%group(i)), real32))
          ! Integral Terms      
          res = res                                                 -&
                (input%time(i) - input%time(i - 1))                 *& !Time Span
                sum((input%beta * real(x * sum(y + z + u), real32)) +& !Term 1
                    (input%theta * real(y, real32))                 +& !2
                    (input%gamma * real(z, real32))                 +& !3
                    (input%eta * real(u, real32))                   +& !4
                    (input%tau * real(v, real32))                   +& !5
                    (input%alpha * real(w, real32))                 +& !6
                    (input%omega * real(z, real32)))                   !7

          y(input%group(i)) = y(input%group(i)) - 1
          z(input%group(i)) = z(input%group(i)) + 1

                              
          mask = x < 0 .or. y < 0 .or. z < 0 .or. w < 0 .or. u < 0 .or. v < 0
          if(any(mask)) then 
            res =  ieee_value(res, ieee_negative_inf) 
            return
          end if
        case(3) !Undetected
          res = res + log(input%theta) + log(1 - input%p_theta) + log(real(y(input%group(i)), real32))
          ! Integral Terms      
          res = res                                                 -&
                (input%time(i) - input%time(i - 1))                 *& !Time Span
                sum((input%beta * real(x * sum(y + z + u), real32)) +& !Term 1
                    (input%theta * real(y, real32))                 +& !2
                    (input%gamma * real(z, real32))                 +& !3
                    (input%eta * real(u, real32))                   +& !4
                    (input%tau * real(v, real32))                   +& !5
                    (input%alpha * real(w, real32))                 +& !6
                    (input%omega * real(z, real32)))                   !7
                    
          y(input%group(i)) = y(input%group(i)) - 1
          u(input%group(i)) = u(input%group(i)) + 1
           
          mask = x < 0 .or. y < 0 .or. z < 0 .or. w < 0 .or. u < 0 .or. v < 0
          if(any(mask)) then 
            res =  ieee_value(res, ieee_negative_inf) 
            return
          end if
        case(4) !Hospitalsized
          res = res + log(input%gamma(input%group(i))) + log(input%p_gamma(input%group(i))) + log(real(z(input%group(i)), real32))
          ! Integral Terms      
          res = res                                                 -&
                (input%time(i) - input%time(i - 1))                 *& !Time Span
                sum((input%beta * real(x * sum(y + z + u), real32)) +& !Term 1
                    (input%theta * real(y, real32))                 +& !2
                    (input%gamma * real(z, real32))                 +& !3
                    (input%eta * real(u, real32))                   +& !4
                    (input%tau * real(v, real32))                   +& !5
                    (input%alpha * real(w, real32))                 +& !6
                    (input%omega * real(z, real32)))                   !7

          z(input%group(i)) = z(input%group(i)) - 1
          v(input%group(i)) = v(input%group(i)) + 1

          mask = x < 0 .or. y < 0 .or. z < 0 .or. w < 0 .or. u < 0 .or. v < 0
          if(any(mask)) then 
            res =  ieee_value(res, ieee_negative_inf) 
            return
          end if
        case(5) !Diagnosis to Recovered
          res = res + log(1 - input%gamma(input%group(i))) + log(input%p_gamma(input%group(i))) + &
                log(real(z(input%group(i)), real32))
          ! Integral Terms      
          res = res                                                 -&
                (input%time(i) - input%time(i - 1))                 *& !Time Span
                sum((input%beta * real(x * sum(y + z + u), real32)) +& !Term 1
                    (input%theta * real(y, real32))                 +& !2
                    (input%gamma * real(z, real32))                 +& !3
                    (input%eta * real(u, real32))                   +& !4
                    (input%tau * real(v, real32))                   +& !5
                    (input%alpha * real(w, real32))                 +& !6
                    (input%omega * real(z, real32)))                   !7

          z(input%group(i)) = z(input%group(i)) - 1
          w(input%group(i)) = w(input%group(i)) + 1
                    
          mask = x < 0 .or. y < 0 .or. z < 0 .or. w < 0 .or. u < 0 .or. v < 0
          if(any(mask)) then 
            res =  ieee_value(res, ieee_negative_inf) 
            return
          end if
        case(6) !Undetected to Recovered
          res = res + log(input%eta(input%group(i))) + log(real(u(input%group(i)), real32))
          ! Integral Terms      
          res = res                                                 -&
                (input%time(i) - input%time(i - 1))                 *& !Time Span
                sum((input%beta * real(x * sum(y + z + u), real32)) +& !Term 1
                    (input%theta * real(y, real32))                 +& !2
                    (input%gamma * real(z, real32))                 +& !3
                    (input%eta * real(u, real32))                   +& !4
                    (input%tau * real(v, real32))                   +& !5
                    (input%alpha * real(w, real32))                 +& !6
                    (input%omega * real(z, real32)))                   !7

          u(input%group(i)) = u(input%group(i)) - 1
          w(input%group(i)) = w(input%group(i)) + 1
                              
          mask = x < 0 .or. y < 0 .or. z < 0 .or. w < 0 .or. u < 0 .or. v < 0
          if(any(mask)) then 
            res =  ieee_value(res, ieee_negative_inf) 
            return
          end if
        case(7) !Hospitalsized to Suspect
          res = res + log(input%tau(input%group(i))) + log(real(v(input%group(i)), real32))
          ! Integral Terms      
          res = res                                                 -&
                (input%time(i) - input%time(i - 1))                 *& !Time Span
                sum((input%beta * real(x * sum(y + z + u), real32)) +& !Term 1
                    (input%theta * real(y, real32))                 +& !2
                    (input%gamma * real(z, real32))                 +& !3
                    (input%eta * real(u, real32))                   +& !4
                    (input%tau * real(v, real32))                   +& !5
                    (input%alpha * real(w, real32))                 +& !6
                    (input%omega * real(z, real32)))                   !7

          v(input%group(i)) = v(input%group(i)) - 1
          x(input%group(i)) = x(input%group(i)) + 1
                              
          mask = x < 0 .or. y < 0 .or. z < 0 .or. w < 0 .or. u < 0 .or. v < 0
          if(any(mask)) then 
            res =  ieee_value(res, ieee_negative_inf) 
            return
          end if
        case(8) !Recovered to Suspect
          res = res + log(input%alpha(input%group(i))) + log(real(w(input%group(i)), real32))
          ! Integral Terms      
          res = res                                                 -&
                (input%time(i) - input%time(i - 1))                 *& !Time Span
                sum((input%beta * real(x * sum(y + z + u), real32)) +& !Term 1
                    (input%theta * real(y, real32))                 +& !2
                    (input%gamma * real(z, real32))                 +& !3
                    (input%eta * real(u, real32))                   +& !4
                    (input%tau * real(v, real32))                   +& !5
                    (input%alpha * real(w, real32))                 +& !6
                    (input%omega * real(z, real32)))                   !7

          w(input%group(i)) = w(input%group(i)) - 1
          x(input%group(i)) = x(input%group(i)) + 1
                              
          mask = x < 0 .or. y < 0 .or. z < 0 .or. w < 0 .or. u < 0 .or. v < 0
          if(any(mask)) then 
            res =  ieee_value(res, ieee_negative_inf) 
            return
          end if
        case(9) !Death
          res = res + log(input%omega(input%group(i))) + log(real(z(input%group(i)), real32))
          ! Integral Terms      
          res = res                                                 -&
                (input%time(i) - input%time(i - 1))                 *& !Time Span
                sum((input%beta * real(x * sum(y + z + u), real32)) +& !Term 1
                    (input%theta * real(y, real32))                 +& !2
                    (input%gamma * real(z, real32))                 +& !3
                    (input%eta * real(u, real32))                   +& !4
                    (input%tau * real(v, real32))                   +& !5
                    (input%alpha * real(w, real32))                 +& !6
                    (input%omega * real(z, real32)))                   !7

          z(input%group(i)) = z(input%group(i)) - 1

          mask = x < 0 .or. y < 0 .or. z < 0 .or. w < 0 .or. u < 0 .or. v < 0
          if(any(mask)) then 
            res =  ieee_value(res, ieee_negative_inf) 
            return
          end if
        case(10) !Vaccination
          ! Integral Terms      
          res = res                                                 -&
                (input%time(i) - input%time(i - 1))                 *& !Time Span
                sum((input%beta * real(x * sum(y + z + u), real32)) +& !Term 1
                    (input%theta * real(y, real32))                 +& !2
                    (input%gamma * real(z, real32))                 +& !3
                    (input%eta * real(u, real32))                   +& !4
                    (input%tau * real(v, real32))                   +& !5
                    (input%alpha * real(w, real32))                 +& !6
                    (input%omega * real(z, real32)))                   !7

          x(input%group(i) + 1) = x(input%group(i) + 1) + 1
          x(input%group(i)) = x(input%group(i)) - 1
          !final group does not have vaccination event

          mask = x < 0 .or. y < 0 .or. z < 0 .or. w < 0 .or. u < 0 .or. v < 0
          if(any(mask)) then 
            res =  ieee_value(res, ieee_negative_inf) 
            return
          end if
        case(11) !Time Last 
          ! Integral Terms      
          res = res                                                 -&
                (input%time(i) - input%time(i - 1))                 *& !Time Span
                sum((input%beta * real(x * sum(y + z + u), real32)) +& !Term 1
                    (input%theta * real(y, real32))                 +& !2
                    (input%gamma * real(z, real32))                 +& !3
                    (input%eta * real(u, real32))                   +& !4
                    (input%tau * real(v, real32))                   +& !5
                    (input%alpha * real(w, real32))                 +& !6
                    (input%omega * real(z, real32)))                   !7
                          
          mask = x < 0 .or. y < 0 .or. z < 0 .or. w < 0 .or. u < 0 .or. v < 0
          if(any(mask)) then 
            res =  ieee_value(res, ieee_negative_inf) 
            return
          end if
        case default
          error stop "Log_Likelihood, select case out of bound"
      end select
    end do
  end function log_likelihood
  
  pure function Imin_exclude(x, Imin) result(y)
    !Exclude some eliments from event object
    type(event), intent(in) :: x
    real(real32), intent(in) :: Imin
    type(event) :: y
    logical, allocatable :: mask(:)
    y = x
    mask = x%time >= Imin

    y%time    = pack(array = x%time, mask = mask)
    y%outcome = pack(array = x%outcome, mask = mask)
    y%group   = pack(array = x%group, mask = mask)
    y%gaus    = pack(array = x%gaus, mask = mask)

  end function Imin_exclude
  
  subroutine output_time(self)
    class(event), intent(in) :: self
    integer(int32) :: i
    character(100), allocatable :: outcome_char(:)
    print *, '      Time        |  group | gaus  |      Outcome       '
    print *, '------------------+--------+-------+--------------------'
    allocate(outcome_char(size(self%outcome)))
    do i = 1, size(self%time)
      select case(self%outcome(i))
        case(-1) 
          outcome_char(i) = "Thinned"
        case(1)
          outcome_char(i) = "Infection"
        case(2)
          outcome_char(i) = "Diagnosis"
        case(3)
          outcome_char(i) = "Undetected"
        case(4)
          outcome_char(i) = "Hospitalized"
        case(5)
          outcome_char(i) = "Diagnosis to Recovered" 
        case(6)
          outcome_char(i) = "Undetected to Recovered"
        case(7)
          outcome_char(i) = "Hospitalized to Suspect"
        case(8)
          outcome_char(i) = "Recovered to Suspect"
        case(9)
          outcome_char(i) = "Death"
        case(10)
          outcome_char(i) = "Vaccination"
        case(11)
          outcome_char(i) = "Time Last"
      end select
      write(*,'(f18.5, 5x, i3, x, f9.3, 4x,a)') self%time(i), self%group(i), self%gaus(i), trim(outcome_char(i))
    end do
    deallocate(outcome_char)
  end subroutine output_time

  subroutine output_parms(self)
    class(event), intent(in) :: self
    integer(int32) :: i
    print *, ' Group |   beta    |  p_gamma  |   gamma   |    eta    |   alpha   |    tau    |    omega  |   p_theta  |  theta  '
    print *, '-------+-----------+-----------+-----------+-----------+-----------+-----------+-----------+------------+---------'
    
    do i = 1, self%num_groups
      write(*,'(i3, 5x, 9(f11.3, x))') i,self%beta(i), self%p_gamma(i), self%gamma(i), &
                                      self%eta(i), self%alpha(i), self%tau(i), self%omega(i), &
                                      self%p_theta, self%theta
    end do
    
    !print *, ' p_theta|theta'
    !print *, '--------+-----'
    !write(*,'(2(f7.3, x))') self%p_theta, self%theta
    
  end subroutine output_parms

  pure function sqexp_cov(x, y, alpha, theta) result(res)
    !!squared expontional covariance function
    !!x, y : input values
    !!alpha: variance scale
    !!theta: length scale
    real(real32), intent(in) :: x
    real(real32), intent(in) :: y
    real(real32), intent(in) :: alpha
    real(real32), intent(in) :: theta
    real(real32) :: res

    res = (alpha ** 2) * exp(-1 * ((x - y) ** 2)/(2 * (theta ** 2)))
  end function sqexp_cov

  pure function cov_mat(input, alpha, theta) result(res)
    !!construct covariance matrix for given input array with parms of alpha and theta
    real(real32), intent(in) :: input(:)
    real(real32), intent(in) :: alpha
    real(real32), intent(in) :: theta
    real(real32) :: res(size(input), size(input))
    
    integer(int32) :: i, j

    row: do concurrent(i = 1:size(input))
      column: do concurrent(j = 1:size(input))
        res(i, j) = sqexp_cov(input(i), input(j), alpha, theta)
      end do column
    end do row


    

  end function cov_mat

  function MPA_sample(psedo, true, alpha, theta) result(res)
    !!construct covariance matrix for given input array with parms of alpha and theta
    real(real32), intent(in) :: psedo(:) !!sampled data from the whole data array
    real(real32), intent(in) :: true(:)  !!The whole data array
    real(real32), intent(in) :: alpha    !!sqexp_cov variance parm
    real(real32), intent(in) :: theta    !!sqexp_cov length parm
    
    real(real32) :: res(size(true))
    
    !intermediate virables
    real(real32) :: cov(size(psedo), size(psedo))
    real(real32) :: adj_mat(size(true), size(psedo)), Garb_mat(size(true), size(psedo))
    real(real32) :: mean(size(psedo))
    real(real32) :: psedo_gaus(size(psedo))
    integer(int32) :: info, i, j

    mean = 0
    cov = cov_mat(psedo, alpha, theta)
    psedo_gaus = random_mvn(cov, size(mean), mean)

    row: do concurrent(i = 1:size(true))
      column: do concurrent(j = 1:size(psedo))
        adj_mat(i, j) = sqexp_cov(true(i), psedo(j), alpha, theta)
      end do column 
    end do row
    if(size(psedo) == 1) then
      cov = 1/cov
    end if
    call inv(size(psedo), cov, info)
    if(info /= 0) then
      res = ieee_value(res, ieee_quiet_nan)
      return
    end if

    !! mat[true * psedo] %*% mat[psedo * psedo]
    call sgemm("N",         & !TRANSA
               "N",         & !TRANSB
               size(true),  & !M
               size(psedo), & !N
               size(psedo), & !K
               1.0,         & !ALPHA
               adj_mat,     & !A
               size(true),  & !LDA
               cov,         & !B
               size(psedo), & !LDB
               0,           & !BETA
               Garb_mat,    & !C
               size(true))    !LDC
    
    !! Garb_mat %*% psedo_gaus
    call sgemm("N",         & !TRANSA
               "N",         & !TRANSB
               size(true),  & !M
               1,           & !N
               size(psedo), & !K
               1.0,         & !ALPHA
               Garb_mat,    & !A
               size(true),  & !LDA
               psedo_gaus,  & !B
               size(psedo), & !LDB
               0,           & !BETA
               res,         & !C
               size(true))    !LDC          

  end function MPA_sample


  pure elemental function logit(x) result(res)
    real(real32), intent(in) :: x
    real(real32) :: res

    res = (1 + exp(-1 * x)) ** (-1)
  end function logit

  subroutine gibbs_sampling(input)
    type(event), intent(inout) :: input
    

    integer(int32) :: i
    integer(int32) :: x(input%num_groups), y(input%num_groups), z(input%num_groups), &
                      u(input%num_groups), v(input%num_groups), w(input%num_groups)
    real(real32) :: beta_int(input%num_groups), theta_int, gamma_int(input%num_groups), &
                    eta_int(input%num_groups), alpha_int(input%num_groups), tau_int(input%num_groups), &
                    omega_int(input%num_groups), I2_1, I_submin

    integer(int32) :: beta_shape1(input%num_groups), beta_shape2(input%num_groups), delta(input%num_groups) !shape1:m, shape2:M
    integer(int32) :: theta_parm1(input%num_groups), theta_parm2(input%num_groups) !shape1:q, shape2:Q
    integer(int32) :: pgamma_parm1(input%num_groups), pgamma_parm2(input%num_groups) !shape1:e, shape2:E
    integer(int32) :: f(input%num_groups), a(input%num_groups), b(input%num_groups), g(input%num_groups)

    logical, allocatable :: mask(:)
    !All paramters are initialized as 0
    beta_shape1 = 0
    beta_shape2 = 0
    delta = 0
    theta_parm1 = 0
    theta_parm2 = 0
    pgamma_parm1 = 0
    pgamma_parm2 = 0
    f = 0
    a = 0
    b = 0
    g = 0

    beta_int = 0
    theta_int = 0
    gamma_int = 0
    eta_int = 0
    alpha_int = 0
    tau_int = 0
    omega_int = 0
  
    !Initialize variables                  
    x = input%N  !易感者
    y = 0        !感染者
    y(input%I_min_group) = 1                         !0号感染者
    x(input%I_min_group) = x(input%I_min_group) - 1  !0号感染者
    beta_shape1(input%I_min_group) = 1               !0号感染者
    z = 0        !确诊者
    u = 0        !未确诊
    v = 0        !住院
    w = 0        !康复


    delta(input%I_min_group) = 1
    do i = 2, size(input%time)
      select case(input%outcome(i))
      case(-1)!Thinned
        ! Integral Terms
        beta_int = beta_int                                 +&
                   (input%time(i) - input%time(i - 1))      *& 
                   (real(x * sum(y + z + u), real32)) 
        
        theta_int = theta_int                               +& 
                    sum((input%time(i) - input%time(i - 1)) *&
                    real(y, real32))
              
        gamma_int = gamma_int                               +& 
                    (input%time(i) - input%time(i - 1))     *&
                    real(z, real32)

        eta_int = eta_int                                   +& 
                  (input%time(i) - input%time(i - 1))       *&
                  real(u, real32)                          
      
        alpha_int = alpha_int                               +&
                    (input%time(i) - input%time(i - 1))     *&
                    real(w, real32)                         
         
        tau_int = tau_int                                   +&
                  (input%time(i) - input%time(i - 1))       *&
                  real(v, real32) 
                  
        omega_int = omega_int                               +&
                    (input%time(i) - input%time(i - 1))     *&
                    real(z, real32)
        
        beta_shape2(input%group(i)) = beta_shape2(input%group(i)) + 1  
        
        
      case(1) !Infection
        ! Integral Terms
        beta_int = beta_int                                 +&
                   (input%time(i) - input%time(i - 1))      *& 
                   (real(x * sum(y + z + u), real32)) 
        
        theta_int = theta_int                               +& 
                    sum((input%time(i) - input%time(i - 1)) *&
                    real(y, real32))
              
        gamma_int = gamma_int                               +& 
                    (input%time(i) - input%time(i - 1))     *&
                    real(z, real32)

        eta_int = eta_int                                   +& 
                  (input%time(i) - input%time(i - 1))       *&
                  real(u, real32)                          
      
        alpha_int = alpha_int                               +&
                    (input%time(i) - input%time(i - 1))     *&
                    real(w, real32)                         
         
        tau_int = tau_int                                   +&
                  (input%time(i) - input%time(i - 1))       *&
                  real(v, real32) 
                  
        omega_int = omega_int                               +&
                    (input%time(i) - input%time(i - 1))     *&
                    real(z, real32)

        x(input%group(i)) = x(input%group(i)) - 1
        y(input%group(i)) = y(input%group(i)) + 1
        
        beta_shape1(input%group(i)) = beta_shape1(input%group(i)) + 1

      case(2) !Diagnosis
        ! Integral Terms
        beta_int = beta_int                                 +&
                   (input%time(i) - input%time(i - 1))      *& 
                   (real(x * sum(y + z + u), real32)) 
        
        theta_int = theta_int                               +& 
                    sum((input%time(i) - input%time(i - 1)) *&
                    real(y, real32))
              
        gamma_int = gamma_int                               +& 
                    (input%time(i) - input%time(i - 1))     *&
                    real(z, real32)

        eta_int = eta_int                                   +& 
                  (input%time(i) - input%time(i - 1))       *&
                  real(u, real32)                          
      
        alpha_int = alpha_int                               +&
                    (input%time(i) - input%time(i - 1))     *&
                    real(w, real32)                         
         
        tau_int = tau_int                                   +&
                  (input%time(i) - input%time(i - 1))       *&
                  real(v, real32) 
                  
        omega_int = omega_int                               +&
                    (input%time(i) - input%time(i - 1))     *&
                    real(z, real32)

        y(input%group(i)) = y(input%group(i)) - 1
        z(input%group(i)) = z(input%group(i)) + 1

        theta_parm1(input%group(i)) = theta_parm1(input%group(i)) + 1

      case(3) !Undetected
        ! Integral Terms
        beta_int = beta_int                                 +&
                   (input%time(i) - input%time(i - 1))      *& 
                   (real(x * sum(y + z + u), real32)) 
        
        theta_int = theta_int                               +& 
                    sum((input%time(i) - input%time(i - 1)) *&
                    real(y, real32))
              
        gamma_int = gamma_int                               +& 
                    (input%time(i) - input%time(i - 1))     *&
                    real(z, real32)

        eta_int = eta_int                                   +& 
                  (input%time(i) - input%time(i - 1))       *&
                  real(u, real32)                          
      
        alpha_int = alpha_int                               +&
                    (input%time(i) - input%time(i - 1))     *&
                    real(w, real32)                         
         
        tau_int = tau_int                                   +&
                  (input%time(i) - input%time(i - 1))       *&
                  real(v, real32) 
                  
        omega_int = omega_int                               +&
                    (input%time(i) - input%time(i - 1))     *&
                    real(z, real32)
                  
        y(input%group(i)) = y(input%group(i)) - 1
        u(input%group(i)) = u(input%group(i)) + 1

        theta_parm2(input%group(i)) = theta_parm2(input%group(i)) + 1

      case(4) !Hospitalsized
        ! Integral Terms
        beta_int = beta_int                                 +&
                   (input%time(i) - input%time(i - 1))      *& 
                   (real(x * sum(y + z + u), real32)) 
        
        theta_int = theta_int                               +& 
                    sum((input%time(i) - input%time(i - 1)) *&
                    real(y, real32))
              
        gamma_int = gamma_int                               +& 
                    (input%time(i) - input%time(i - 1))     *&
                    real(z, real32)

        eta_int = eta_int                                   +& 
                  (input%time(i) - input%time(i - 1))       *&
                  real(u, real32)                          
      
        alpha_int = alpha_int                               +&
                    (input%time(i) - input%time(i - 1))     *&
                    real(w, real32)                         
         
        tau_int = tau_int                                   +&
                  (input%time(i) - input%time(i - 1))       *&
                  real(v, real32) 
                  
        omega_int = omega_int                               +&
                    (input%time(i) - input%time(i - 1))     *&
                    real(z, real32)

        z(input%group(i)) = z(input%group(i)) - 1
        v(input%group(i)) = v(input%group(i)) + 1

        pgamma_parm1(input%group(i)) = pgamma_parm1(input%group(i)) + 1

      case(5) !Diagnosis to Recovered
        ! Integral Terms
        beta_int = beta_int                                 +&
                   (input%time(i) - input%time(i - 1))      *& 
                   (real(x * sum(y + z + u), real32)) 
        
        theta_int = theta_int                               +& 
                    sum((input%time(i) - input%time(i - 1)) *&
                    real(y, real32))
              
        gamma_int = gamma_int                               +& 
                    (input%time(i) - input%time(i - 1))     *&
                    real(z, real32)

        eta_int = eta_int                                   +& 
                  (input%time(i) - input%time(i - 1))       *&
                  real(u, real32)                          
      
        alpha_int = alpha_int                               +&
                    (input%time(i) - input%time(i - 1))     *&
                    real(w, real32)                         
         
        tau_int = tau_int                                   +&
                  (input%time(i) - input%time(i - 1))       *&
                  real(v, real32) 
                  
        omega_int = omega_int                               +&
                    (input%time(i) - input%time(i - 1))     *&
                    real(z, real32)

        z(input%group(i)) = z(input%group(i)) - 1
        w(input%group(i)) = w(input%group(i)) + 1

        pgamma_parm2(input%group(i)) = pgamma_parm2(input%group(i)) + 1

      case(6) !Undetected to Recovered
        ! Integral Terms
        beta_int = beta_int                                 +&
                   (input%time(i) - input%time(i - 1))      *& 
                   (real(x * sum(y + z + u), real32)) 
        
        theta_int = theta_int                               +& 
                    sum((input%time(i) - input%time(i - 1)) *&
                    real(y, real32))
              
        gamma_int = gamma_int                               +& 
                    (input%time(i) - input%time(i - 1))     *&
                    real(z, real32)

        eta_int = eta_int                                   +& 
                  (input%time(i) - input%time(i - 1))       *&
                  real(u, real32)                          
      
        alpha_int = alpha_int                               +&
                    (input%time(i) - input%time(i - 1))     *&
                    real(w, real32)                         
         
        tau_int = tau_int                                   +&
                  (input%time(i) - input%time(i - 1))       *&
                  real(v, real32) 
                  
        omega_int = omega_int                               +&
                    (input%time(i) - input%time(i - 1))     *&
                    real(z, real32)

        u(input%group(i)) = u(input%group(i)) - 1
        w(input%group(i)) = w(input%group(i)) + 1

        f(input%group(i)) = f(input%group(i)) + 1


      case(7) !Hospitalsized to Suspect
        ! Integral Terms
        beta_int = beta_int                                 +&
                   (input%time(i) - input%time(i - 1))      *& 
                   (real(x * sum(y + z + u), real32)) 
        
        theta_int = theta_int                               +& 
                    sum((input%time(i) - input%time(i - 1)) *&
                    real(y, real32))
              
        gamma_int = gamma_int                               +& 
                    (input%time(i) - input%time(i - 1))     *&
                    real(z, real32)

        eta_int = eta_int                                   +& 
                  (input%time(i) - input%time(i - 1))       *&
                  real(u, real32)                          
      
        alpha_int = alpha_int                               +&
                    (input%time(i) - input%time(i - 1))     *&
                    real(w, real32)                         
         
        tau_int = tau_int                                   +&
                  (input%time(i) - input%time(i - 1))       *&
                  real(v, real32) 
                  
        omega_int = omega_int                               +&
                    (input%time(i) - input%time(i - 1))     *&
                    real(z, real32)

        v(input%group(i)) = v(input%group(i)) - 1
        x(input%group(i)) = x(input%group(i)) + 1

        a(input%group(i)) = a(input%group(i)) + 1
      case(8) !Recovered to Suspect
        ! Integral Terms
        beta_int = beta_int                                 +&
                   (input%time(i) - input%time(i - 1))      *& 
                   (real(x * sum(y + z + u), real32)) 
        
        theta_int = theta_int                               +& 
                    sum((input%time(i) - input%time(i - 1)) *&
                    real(y, real32))
              
        gamma_int = gamma_int                               +& 
                    (input%time(i) - input%time(i - 1))     *&
                    real(z, real32)

        eta_int = eta_int                                   +& 
                  (input%time(i) - input%time(i - 1))       *&
                  real(u, real32)                          
      
        alpha_int = alpha_int                               +&
                    (input%time(i) - input%time(i - 1))     *&
                    real(w, real32)                         
         
        tau_int = tau_int                                   +&
                  (input%time(i) - input%time(i - 1))       *&
                  real(v, real32) 
                  
        omega_int = omega_int                               +&
                    (input%time(i) - input%time(i - 1))     *&
                    real(z, real32)

        w(input%group(i)) = w(input%group(i)) - 1
        x(input%group(i)) = x(input%group(i)) + 1

        b(input%group(i)) = b(input%group(i)) + 1

      case(9) !Death
        ! Integral Terms
        beta_int = beta_int                                 +&
                   (input%time(i) - input%time(i - 1))      *& 
                   (real(x * sum(y + z + u), real32)) 
        
        theta_int = theta_int                               +& 
                    sum((input%time(i) - input%time(i - 1)) *&
                    real(y, real32))
              
        gamma_int = gamma_int                               +& 
                    (input%time(i) - input%time(i - 1))     *&
                    real(z, real32)

        eta_int = eta_int                                   +& 
                  (input%time(i) - input%time(i - 1))       *&
                  real(u, real32)                          
      
        alpha_int = alpha_int                               +&
                    (input%time(i) - input%time(i - 1))     *&
                    real(w, real32)                         
         
        tau_int = tau_int                                   +&
                  (input%time(i) - input%time(i - 1))       *&
                  real(v, real32) 
                  
        omega_int = omega_int                               +&
                    (input%time(i) - input%time(i - 1))     *&
                    real(z, real32)

        z(input%group(i)) = z(input%group(i)) - 1

        g(input%group(i)) = g(input%group(i)) + 1

      case(10) !Vaccination
        ! Integral Terms
        beta_int = beta_int                                 +&
                   (input%time(i) - input%time(i - 1))      *& 
                   (real(x * sum(y + z + u), real32)) 
        
        theta_int = theta_int                               +& 
                    sum((input%time(i) - input%time(i - 1)) *&
                    real(y, real32))
              
        gamma_int = gamma_int                               +& 
                    (input%time(i) - input%time(i - 1))     *&
                    real(z, real32)

        eta_int = eta_int                                   +& 
                  (input%time(i) - input%time(i - 1))       *&
                  real(u, real32)                          
      
        alpha_int = alpha_int                               +&
                    (input%time(i) - input%time(i - 1))     *&
                    real(w, real32)                         
         
        tau_int = tau_int                                   +&
                  (input%time(i) - input%time(i - 1))       *&
                  real(v, real32) 
                  
        omega_int = omega_int                               +&
                    (input%time(i) - input%time(i - 1))     *&
                    real(z, real32)

        x(input%group(i) + 1) = x(input%group(i) + 1) + 1
        x(input%group(i)) = x(input%group(i)) - 1
        !final group does not have vaccination event
      
      case(11)
        ! Integral Terms
        beta_int = beta_int                                 +&
                   (input%time(i) - input%time(i - 1))      *& 
                   (real(x * sum(y + z + u), real32)) 
        
        theta_int = theta_int                               +& 
                    sum((input%time(i) - input%time(i - 1)) *&
                    real(y, real32))
              
        gamma_int = gamma_int                               +& 
                    (input%time(i) - input%time(i - 1))     *&
                    real(z, real32)

        eta_int = eta_int                                   +& 
                  (input%time(i) - input%time(i - 1))       *&
                  real(u, real32)                          
      
        alpha_int = alpha_int                               +&
                    (input%time(i) - input%time(i - 1))     *&
                    real(w, real32)                         
         
        tau_int = tau_int                                   +&
                  (input%time(i) - input%time(i - 1))       *&
                  real(v, real32) 
                  
        omega_int = omega_int                               +&
                    (input%time(i) - input%time(i - 1))     *&
                    real(z, real32)

      case default
        print *, "error value",input%outcome(i)
        error stop "Gibbs sampling, select case out of bound"
    end select
    end do

    !Sampling
    do i = 1, input%num_groups
      !print *, "beta"
      input%beta(i) = random_gamma(shape = real(beta_shape1(i) + beta_shape2(i) + &
                                                input%priors%beta(i, 1) - delta(i), real32), &
                                   scale = (input%priors%beta(i, 2) + beta_int(i)))

      input%p_gamma(i) = random_beta(input%priors%p_gamma(i, 1) + pgamma_parm1(i), &
                                     input%priors%p_gamma(i, 2) + pgamma_parm2(i))
      !print *, "gamma"
      input%gamma(i) = random_gamma(shape = real(pgamma_parm1(i) + pgamma_parm2(i) + &
                                                 input%priors%gamma(i, 1), real32), &
                                    scale = (input%priors%gamma(i, 2) + gamma_int(i)))
      !print *, "eta"
      !print * , eta_int(1), u(1), input%log_likelihood
      input%eta(i) = random_gamma(shape = real(f(i) + input%priors%eta(i, 1), real32), &
                                  scale = (input%priors%eta(i, 2) + eta_int(i)))   
      !print *, "alpha"  
      input%alpha(i) = random_gamma(shape = real(b(i) + input%priors%alpha(i, 1), real32), &
                                    scale = (input%priors%alpha(i, 2) + alpha_int(i)))
      !print *, "tau"
      input%tau(i) = random_gamma(shape = real(a(i) + input%priors%tau(i, 1), real32), &
                                  scale = (input%priors%tau(i, 2) + tau_int(i)))
      !print *, "omega"
      input%omega(i) = random_gamma(shape = real(g(i) + input%priors%omega(i, 1), real32), &
                                    scale = (input%priors%omega(i, 2) + omega_int(i)))        
    end do

    input%p_theta = random_beta(real(input%priors%p_theta(1) + sum(theta_parm1), real32), &
                                real(input%priors%p_theta(2) + sum(theta_parm2), real32))
    !print *, "theta"
    input%theta = random_gamma(shape = sum(theta_parm1 + theta_parm2) + input%priors%theta(1), &
                               scale = (input%priors%theta(2) + theta_int))   


    I2_1 = random_exp(input%theta + input%beta(input%I_min_group) * (sum(input%N) - 1))    
    
    find_I_submin:do i = 2, size(input%time)
      if(input%outcome(i) == 1) then 
        I_submin = input%time(i)
        exit find_I_submin
      end if
    end do find_I_submin
    
    mask = input%time < (I_submin - I2_1) .and. (input%outcome == 2 .or. input%outcome == 4 .or. input%outcome == 10)
    if(size(pack(mask, mask .eqv. .true.)) >= 1) return !Imin Gibbs 失败
    input%time(1) = I_submin - I2_1
    

    input = Imin_exclude(input, input%time(1))


  end subroutine gibbs_sampling

  function Insert_Event(time, outcome, group, target) result(res)
    real(real32), intent(in) :: time
    integer(int32), intent(in) :: outcome
    integer(int32), intent(in) :: group
    type(event), intent(in) :: target
    type(event) :: res

    integer(int32) :: idx ! new elem's loc
    integer(int32) :: num ! number of old data points
    integer(int32) :: MPA(target%MPA_Level) !MPA sampled idx
    real(real32) :: NaN, unif
    logical, allocatable :: mask(:)
    real(real32), allocatable :: old(:), old_gaus(:), time_interval

    NaN = ieee_value(NaN, class = ieee_quiet_nan)
    res = target
    !Simply insert a time point
    idx = maxval(trueloc(target%time < time)) + 1
    res%time = put(time, idx, res%time)
    res%outcome = put(outcome, idx, res%outcome)
    res%group = put(group, idx, res%group)
    res%gaus = put(NaN, idx, res%gaus)

    num = size(pack(target%outcome, target%outcome == outcome .and. target%group == group))
    
    !Change associated values
    if(outcome == 1 .or. outcome == -1) then !Make Gaussian Predications
      mask = (target%outcome == 1 .or. target%outcome == -1) &
              .and. target%group == group .and. target%time /= target%time(1)
      old = pack(target%time, mask)
      
      if(size(old) < 1) then !不存在old, 无法做出观测
        res = target
        return
      end if
      old_gaus = pack(target%gaus, mask)
      
      if(num < target%MPA_Level) then
        res%gaus(idx) = Gauss_Pred(time, old, old_gaus, target%variance_scale, target%shape_scale(group))
        
      else !MPA: 此时old_gaus数量过大, 故抽取部分样本作为old_gaus 
        
        do while (isnan(res%gaus(idx)))
          MPA = rand_int(1, num, target%MPA_Level)
          call std_sort(MPA)
          res%gaus(idx) = Gauss_Pred(time, subscript(old, MPA), subscript(old_gaus, MPA), &
                                     target%variance_scale, target%shape_scale(group))
        end do
      
      end if
    end if
    
    !Evaluate Likelihood
    res%log_likelihood = log_likelihood(res)
    if(ieee_is_finite(res%log_likelihood)) then
      
      mask = target%outcome == 1 .and. target%group == group !Find Ii1 !Possible bug
      time_interval = last(target%time) - minval(pack(target%time, mask))  !时间区间
      unif = log(real(random_uniform(0.0, 1.0), real32))
      if(log(time_interval)                            +&
         ((res%log_likelihood - target%log_likelihood) -&
         log(real(num + 1, real32))) <= unif) then
        res = target
      end if
    else 
      res = target
    end if
  end function Insert_Event

  function Delete_Event(time, outcome, group, target) result(res)
    real(real32), intent(in) :: time
    integer(int32), intent(in) :: outcome
    integer(int32), intent(in) :: group
    type(event), intent(in) :: target
    type(event) :: res

    integer(int32) :: idx ! delete elem at loc idx
    integer(int32) :: num ! number of old data points
    real(real32) :: unif, time_interval
    logical, allocatable :: mask(:)

    res = target
    !Simply delete a time point
    idx = minval(trueloc(target%time == time .and. target%outcome == outcome .and. target%group == group))
    if(idx <= 1 .or. idx > size(target%time)) then
      res = target
      print *, "idx not found in deletion step, with wrong idx value:", idx
      return
      
    end if
    res%time = [target%time(1:idx - 1), target%time(idx + 1:)]
    res%outcome = [target%outcome(1:idx - 1), target%outcome(idx + 1:)]
    res%group = [target%group(1:idx - 1), target%group(idx + 1:)]
    res%gaus = [target%gaus(1:idx - 1), target%gaus(idx + 1:)]    

    !Evaluate Likelihood
    num = size(pack(target%outcome, target%outcome == outcome .and. target%group == group))
    res%log_likelihood = log_likelihood(res)
    unif = log(real(random_uniform(0.0, 1.0), real32))
    mask = target%outcome == 1 .and. target%group == group !Find Ii1 !Possible bug
    time_interval = last(target%time) - minval(pack(target%time, mask))  !时间区间
    if(ieee_is_finite(res%log_likelihood)) then
      if((log(real(num, real32))                      +&
         (res%log_likelihood - target%log_likelihood) -&
         log(time_interval)) <= unif) then
        res = target
      end if
    else 
      res = target
    end if
  end function Delete_Event

  function Update_Event(old_time, new_time, outcome, group, target) result(res)
    real(real32), intent(in) :: old_time, new_time
    integer(int32), intent(in) :: outcome
    integer(int32), intent(in) :: group
    type(event), intent(in) :: target
    type(event) :: res

    integer(int32) :: old_idx ! delete elem at loc idx
    integer(int32) :: new_idx ! insert elem at loc idx
    integer(int32) :: MPA(target%MPA_Level) !MPA sampled idx
    real(real32) :: NaN, unif
    logical, allocatable :: mask(:)
    real(real32), allocatable :: old(:), old_gaus(:)

    NaN = ieee_value(NaN, class = ieee_quiet_nan)
    res = target
    
    !Simply delete a time point
    old_idx = minval(trueloc(target%time == old_time .and. target%outcome == outcome .and. target%group == group))
    if(old_idx <= 1 .or. old_idx > size(target%time)) then
      res = target
      print *, "idx not found in deletion step, with wrong idx value:", old_idx
      return
      
    end if
    res%time = [target%time(1:old_idx - 1), target%time(old_idx + 1:)]
    res%outcome = [target%outcome(1:old_idx - 1), target%outcome(old_idx + 1:)]
    res%group = [target%group(1:old_idx - 1), target%group(old_idx + 1:)]
    res%gaus = [target%gaus(1:old_idx - 1), target%gaus(old_idx + 1:)]    

    !Simply insert a time point
    new_idx = maxval(trueloc(res%time < new_time)) + 1 !CAUTION：使用res寻找new idx, 而不是target
    if(new_idx <= 1 .or. new_idx > size(target%time)) then
      res = target
      print *, "idx not found in Insert step, with wrong idx value:", new_idx
      return
      
    end if
    res%time = put(new_time, new_idx, res%time)
    res%outcome = put(outcome, new_idx, res%outcome)
    res%group = put(group, new_idx, res%group)
    res%gaus = put(NaN, new_idx, res%gaus)

    
    !Change associated values
    if(outcome == 1 .or. outcome == -1) then !Make Gaussian Predications
      mask = (target%outcome == 1 .or. target%outcome == -1) &
              .and. target%group == group .and. target%time /= target%time(1)
      old = pack(target%time, mask)
      old_gaus = pack(target%gaus, mask)
      
      if(size(old) < target%MPA_Level) then

        res%gaus(new_idx) = Gauss_Pred(new_time, old, old_gaus, target%variance_scale, target%shape_scale(group))

      else !MPA: 此时old_gaus数量过大, 故抽取部分样本作为old_gaus 
        
        do while (isnan(res%gaus(new_idx)))
          MPA = rand_int(1, size(old), target%MPA_Level)
          call std_sort(MPA)
          res%gaus(new_idx) = Gauss_Pred(new_time, subscript(old, MPA), subscript(old_gaus, MPA), &
                                     target%variance_scale, target%shape_scale(group))
        end do
      
      end if
    end if

    !Evaluate Likelihood
    res%log_likelihood = log_likelihood(res)
    if(ieee_is_finite(res%log_likelihood)) then

      unif = log(real(random_uniform(0.0, 1.0), real32))
      if((res%log_likelihood - target%log_likelihood) <= unif) then
        res = target
      end if
    else 
      res = target
    end if

  end function Update_Event

  function Update_Gaus(x, group) result(y)
    type(event), intent(in) :: x
    type(event) :: y
    integer(int32), intent(in) :: group

    real(real32), allocatable :: cov(:, :), mean(:)
    integer(int32), allocatable :: Infec_loc(:), idx(:)
    real(real32), parameter :: nu = 0.5
    real(real32) :: unif

    y = x

    Infec_loc = trueloc((y%outcome == 1 .or. y%outcome == -1) .and. y%group == group .and. y%time /= minval(y%time)) !exclude I min
    if(size(Infec_loc) <= y%MPA_Level) then 
      if(size(Infec_loc) < 1) then
        y = x
        return
      end if
      cov = cov_mat(y%time(Infec_loc), y%variance_scale, y%shape_scale(group)) !construct covariance matrix
        !if(i == 1) print *, "group1 covariance matrix size is", size(cov, 1)
      allocate(mean, mold = y%time(Infec_loc)) 
      mean = 0
      y%gaus(Infec_loc) = sqrt(1 - nu**2) * random_mvn(cov, size(Infec_loc), mean) + nu * y%gaus(Infec_loc)!assign sampled gaus value to corresponding location
      deallocate(mean)
      
    else !MPA
      do while(any(isnan(y%gaus(Infec_loc))))
        idx = rand_int(1, size(Infec_loc), y%MPA_Level)
        call std_sort(idx)
        y%gaus(Infec_loc) = MPA_sample(y%time(Infec_loc(idx)), y%time, y%variance_scale, y%shape_scale(group)) * &
                            sqrt(1 - nu**2) + nu * y%gaus(Infec_loc)
      end do
    end if

  
    !Evaluate Likelihood
    y%log_likelihood = log_likelihood(y)
    if(ieee_is_finite(y%log_likelihood)) then

      unif = log(real(random_uniform(), real32))
      if((y%log_likelihood - x%log_likelihood) <= unif) then
        y = x
      end if
    else 
      y = x
    end if
    
  end function Update_Gaus

  function Update_shape_scale(lambda, group, x) result(y)
    
    real(real32), intent(in) :: lambda
    type(event), intent(in) :: x
    type(event) :: y
    integer(int32), intent(in) :: group

    real(real32), allocatable :: cov(:, :), gaus(:), new_cov(:, :), gaus_trans(:, :), gaus_mat(:, :)
    integer(int32), allocatable :: Infec_loc(:), idx(:)
    real(real32) :: shape, unif, det_old, det_new, up, down
    real(real32), allocatable :: old_garb1(:, :)
    real(real32), allocatable :: old_garb2(:, :)
    real(real32), allocatable :: new_garb1(:, :)
    real(real32), allocatable :: new_garb2(:, :)

    integer(int32) :: info, i
    y = x

    shape = random_uvn(random_exp(lambda), 1.0)
    y%shape_scale(group) = shape

    Infec_loc = trueloc((y%outcome == 1 .or. y%outcome == -1) .and. y%group == group .and. y%time /= minval(y%time)) !exclude I min
    if(size(Infec_loc) <= y%MPA_Level) then 
      if(size(Infec_loc) <= 0) then 
        y = x
        return
      end if
      cov = cov_mat(x%time(Infec_loc), x%variance_scale, x%shape_scale(group)) !construct covariance matrix
      new_cov = cov_mat(y%time(Infec_loc), y%variance_scale, shape) !construct new covariance matrix with new shape value
      gaus = x%gaus(Infec_loc)
    else !MPA
      idx = rand_int(1, size(Infec_loc), x%MPA_Level)
      call std_sort(idx)
      cov = cov_mat(x%time(Infec_loc(idx)), x%variance_scale, x%shape_scale(group))
      new_cov = cov_mat(y%time(Infec_loc(idx)), y%variance_scale, shape)
      gaus = x%gaus(Infec_loc(idx))
    end if

    !Evaluate Likelihood
    if (size(cov, 1) > 1 ) then
      call det(cov, det_old, info)
      if(info /= 0) then
        !print * , "BAD MPA SAMPLING"
        y = x
        return
      end if
      call det(new_cov, det_new, info)
      if(info /= 0) then
        y = x
        return
      end if

      call inv(size(new_cov, 1), new_cov, info)
      if(info /= 0) then
        y = x
        return
      end if
    
      call inv(size(cov, 1), cov, info)  
      if(info /= 0) then
        y = x
        return
      end if
    else 
      det_new = new_cov(1, 1)
      det_old = cov(1, 1)
      cov = 1/cov
      new_cov = 1/new_cov
    end if



    !cov and new_cov are now inversed
    
    allocate(old_garb1(1, size(gaus)), &
             old_garb2(1, 1),          &
             new_garb1(1, size(gaus)), &
             new_garb2(1, 1),          &
             gaus_trans(1, size(gaus)),&
             gaus_mat(size(gaus), 1))
    gaus_trans(1, :) = gaus
    gaus_mat = transpose(gaus_trans)
    !new_garb2 => trans(gaus) %*% inv(new_cov) %*% gaus 
    call sgemm("n", "n", 1, size(gaus), size(gaus), 1.0, gaus_trans, 1, new_cov, size(gaus), 0, new_garb1, 1)
    call sgemm("n", "n", 1, size(gaus), size(gaus), 1.0, new_garb1, 1, gaus_mat, size(gaus), 0, new_garb2, 1)
    !old_garb2 => trans(gaus) %*% inv(cov) %*% gaus 
    call sgemm("n", "n", 1, size(gaus), size(gaus), 1.0, gaus_trans, 1, cov, size(gaus), 0, old_garb1, 1)
    call sgemm("n", "n", 1, size(gaus), size(gaus), 1.0, old_garb1, 1, gaus_mat, size(gaus), 0, old_garb2, 1)

    up = -0.5 * log(det_new) + (-0.5 * new_garb2(1, 1) - lambda * shape)
    down = -0.5 * log(det_old) + (-0.5 * old_garb2(1, 1) - lambda * x%shape_scale(group))

    unif = log(real(random_uniform(), real32))
    if ( (up - down) <= unif ) then
      y = x
      return
    end if
    
  end function Update_shape_scale

  function Gauss_Pred(new, old, gaus, alpha, theta) result(res)
    
    real(real32), intent(in) :: new !!new time point
    real(real32), intent(in) :: old(:) !! old time vector
    real(real32), intent(in) :: gaus(:) !! old gaus vector
    real(real32), intent(in) ::alpha !! variance parm of link function - sqexp_cov
    real(real32), intent(in) :: theta !! scale parm of link function - sqexp_cov

    real(real32) :: res !!resulted gaus prediction

    real(real32) :: K1(1, size(old)), & !K(Iij, X)
                    K2(size(old), size(old)), &!K(X, X)
                    K3(size(old), 1), &!K(X, Iij)
                    Garb(1, size(old)), mean(1), var(1, 1)

    integer(int32) :: i, info

    
    K2 = cov_mat(old, alpha, theta)

    do concurrent(i = 1:size(old))
      K1(1, i) = sqexp_cov(new, old(i), alpha, theta)
    end do

    K3 = transpose(K1)
    if(size(old) == 1) then
      K2 = 1/K2
    end if
    call inv(size(old), K2, info)
    if(info /= 0) then
      res = ieee_value(res, ieee_signaling_nan)
      return
    end if

    call sgemm('n', 'n', 1, size(old), size(old), 1.0, K1, 1, K2, size(old), 0, Garb, 1)
    call sgemm('n', 'n', 1, 1, size(old), 1.0, Garb, 1, gaus, size(old), 0, mean, 1) !mean vector finished
    call sgemm('n', 'n', 1, 1, size(old), 1.0, Garb, 1, K3, size(old), 0, var, 1) !var vector finished

    var = -1.0 * var + sqexp_cov(new, new, alpha, theta)
    if(var(1, 1) <= 0) then
      res = mean(1)
      return
    end if
    res = random_uvn(mean(1), var(1, 1))
    

  end function Gauss_Pred

  function rand_int(a, b, n) result(res)
    !!sample n integers from [a, b]
    
    integer(int32), intent(in) :: a, b, n
    integer(int32) :: res(n)
    res = random_uniform(a, b - a, n)


  end function rand_int

  function sample(x) result(y)
    real(real32), intent(in) :: x(:)
    real(real32) :: y
  
    integer(int32) :: idx
    if(size(x) == 1) then
      y = x(1)
      return
    end if
    idx = maxval(rand_int(1, size(x), 1))
    y = x(idx)
  end function sample
end module mod_events