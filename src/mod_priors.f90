module mod_prior
    use iso_fortran_env, only:real32, int32
    implicit none
    
    private
    public :: prior

    type :: prior
        real(real32), allocatable :: beta(:, :), p_gamma(:, :), gamma(:, :), eta(:, :), alpha(:, :), tau(:, :), omega(:, :)
        real(real32), allocatable :: p_theta(:), theta(:)
    contains
        procedure, pass(self) :: output
    end type prior

    interface prior
        module procedure :: prior_constructor
    end interface

contains
    pure function prior_constructor(beta, p_gamma, gamma, eta, &
                                    alpha, tau, omega, p_theta, theta) result(res)
        
        real(real32), intent(in):: beta(:, :), p_gamma(:, :), gamma(:, :), &
                                   eta(:, :), alpha(:, :), tau(:, :), omega(:, :)
        real(real32), intent(in) :: p_theta(:), theta(:)
        type(prior) :: res
        
        res%alpha = alpha
        res%beta = beta
        res%tau = tau
        res%omega = omega
        res%gamma = gamma
        res%eta = eta
        res%p_theta = p_theta
        res%theta = theta
        res%p_gamma = p_gamma

    end function prior_constructor

    subroutine output(self)
        class(prior), intent(in) :: self
        integer(int32) :: i
        print *, ' Group | beta  |p_gamma| gamma |  eta  | alpha |  tau  |  omega'
        print *, '-------+-------+-------+-------+-------+-------+-------+-------'
        
        do i = 1, size(self%beta)
          write(*,'(i3, 5x, 7(f7.3, x, f7.3))') i,self%beta(i, :), self%p_gamma(i, :), self%gamma(i, :), &
                                         self%eta(i, :), self%alpha(i, :), self%tau(i, :), self%omega(i,:)
        end do
        print *, ' p_theta|theta'
        print *, '--------+-----'
        write(*,'(2(f7.3, x, f7.3))') self%p_theta, self%theta
    end subroutine output

end module mod_prior