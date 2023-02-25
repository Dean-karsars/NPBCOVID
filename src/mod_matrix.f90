module mod_matrix
    use iso_fortran_env, only:real32
    use ieee_arithmetic
    implicit none
    
contains
    subroutine inv(ndim, Amat, info)
        ! ndim is dimension of Amat,the input Amat will be written which is result
        implicit none
        integer,parameter::dp = 4
        integer::i
        integer,intent(out)::info
        integer,intent(in)::ndim   ! in 代表这个值只能是个输入值,不会被改变
    !    IPIV   : INTEGER. Array,DIMENSION at least max(1,n). The pivot indices that define
    !    the permutation matrix P; row i of the matrix was interchanged with
    !    row ipiv(i). Corresponds to the single precision factorization (if info=
    !    0 and iter ≥ 0) or the double precision factorization (if info= 0 and
    !    iter < 0).
        integer,allocatable::ipiv(:)
        real(dp),parameter::zone = (1.00,0.00)
    !    Amat  :
    !    Overwritten by the factors L and U from the factorization of A = P*L*U;
    !    the unit diagonal elements of L are not stored.
    !    If iterative refinement has been successfully used (info= 0 and iter≥
    !    0),then A is unchanged.
    !    If double precision factorization has been used (info= 0 and iter <
    !    0),then the array A contains the factors L and U from the factorization
    !    A = P*L*U; the unit diagonal elements of L are not stored.
        real(dp),intent(inout):: Amat(ndim,ndim)
    !    Bmat  :
    !    Overwritten by the solution matrix X for dgesv,sgesv,zgesv,zgesv.
        real(dp),allocatable::Bmat(:,:)
        allocate(ipiv(ndim))
        allocate(Bmat(ndim,ndim))
        ipiv=0
        ! unit matrix
        Bmat= (0.0,0.0)
        do i=1,ndim
            Bmat(i,i)= zone
        end do
        call sgesv(ndim,ndim,Amat,ndim,ipiv,Bmat,ndim,info)
        !if(info.ne.0) print *,'Function INV Failed: sgesv error', "ndim = ", ndim
        Amat=Bmat
        return
    end subroutine inv

    ! performs matrix-matrix multiply
! C := alpha*op( A )*op( B ) + beta*C
    function mat_mul(A,B)  result(C)
        ! nmatdim is dimension of matrix
        implicit none
        integer,parameter::dp = real32 ! double precision(精度控制)
        real(dp)::ALPHA
        real(dp)::BETA 
        real(dp),intent(in) ::A(:, :)
        real(dp),intent(in) ::B(:, :)
        !complex(dp)::mat_mul(nmatdim,nmatdim)
        real(dp) , allocatable :: C(:, :)
        integer :: m, n, k
        alpha = 1.0
        beta = 0.0
        allocate(C(m,n))
        C(:,:) = (0.00,0.00)
        m = size(A, 1) !rows of A
        n = size(B, 2) !columns of B
        k = size(A, 2) !columns of A
        

        call SGEMM('N','N', m, n, k, ALPHA, A, m, B, k, BETA, C, m)
        
    end function mat_mul

    subroutine det(mat, res, info)
        real(4), intent(in) :: mat(:, :)
        integer, intent(in) :: info
        real(4), intent(inout) :: res
        integer :: i,j,n = rank(mat)
        integer :: ipiv(rank(mat))


        call SGETRF(n,n,mat,n,ipiv,info)
        
        ! Check for errors
        if (info /= 0) then
            !print *, 'Error: DGETRF failed'
            res = ieee_value(res, ieee_quiet_nan)
            return
        end if

        ! Compute the determinant as the product of diagonal elements of U 
        ! and the sign change due to row interchanges
        res = mat(1,1)

        ! Loop over diagonal elements of U  
        do i =2,n 
            res = res * mat(i,i)
        end do

        ! Loop over pivot indices and count sign changes  
        do i=1,n 
            if (ipiv(i) /= i) then 
                res = -res 
            end if 
        end do
       
    end subroutine det

    subroutine cholesky_det(mat, res, info)
        real(real32), intent(in) :: mat(:, :)
        integer, intent(out) :: info
        real(real32), intent(inout) :: res
        
        real(real32) :: det_a ! 行列式det(A)
        real(real32) :: det_l ! 行列式det(L)
        integer :: i, j! 循环变量 
        integer :: n ! 矩阵的维度
        n = size(mat, 1) 

        call spotrf('L', n, mat, n, info)
        ! 检查返回值是否正常
        if (info /= 0) then
            !print *, 'dpotrf failed with info = ', info
            !stop
            res = ieee_value(res, ieee_quiet_nan)
            return
        end if

        
        ! 计算det(L)为主对角线元素的乘积
        det_l = 1.0
        do i = 1,n
            det_l = det_l * mat(i,i)
        end do

        ! 计算det(A)为(det(L))^2
        det_a = det_l**2
        res = det_a

    end subroutine cholesky_det
end module mod_matrix