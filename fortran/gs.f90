MODULE gs

  implicit none

CONTAINS

  subroutine serialgs(n,M,b,x,r,iter,eps,maxit)

    ! SCALAR ARGUMENTS
    integer :: n,iter,maxit
    real(8) :: eps,r

    ! ARRAY ARGUMENTS
    real(8) :: M(n,n),b(n),x(n)
       
    intent(in   ) :: M,b,n,eps,maxit
    intent(out  ) :: r,iter
    intent(inout) :: x

    ! LOCAL SCALARS
    integer :: i
    real(8) :: prev(n)

    r = 1.0D+20

    iter = 0

    do while (r .gt. eps .and. iter .lt. maxit)

       ! UPDATING X

       do i = 1,n
          prev(i) = x(i)
       end do

       call updatex(1,n,n,M,b,x)

       ! STOPPING CRITERIUM

       r = 0.0D0
       do i = 1,n
          r = max(r,abs(prev(i) - x(i)))
          if (r .gt. eps) then
             exit
          end if
       end do

       iter = iter + 1

    end do

  end subroutine serialgs

  subroutine updatex(start,end,n,M,b,x)

    ! This subroutine applies the Gauss-Seidel step in a subset of the
    ! rows and variables of the system.

    ! ARGUMENTS
    integer :: start,end,n
    real(8) :: M(n,n),b(n),x(n)

    intent(in   ) :: start,end,n,M,b
    intent(inout) :: x

    ! LOCAL SCALARS
    integer :: i,j
    real(8) :: s

    do i = start,end
       s = b(i)
       do j = 1,n
          if (i .ne. j) then
             s = s - M(i,j) * x(j)
          end if
       end do
       x(i) = s / M(i,i)
    end do

  end subroutine updatex

  integer function ceil(x)

    ! SCALAR ARGUMENTS
    real(8), intent(in) :: x

    ceil = INT(x)

    if (ceil .ne. x) then
       return
    end if

    return
    
  end function ceil

  subroutine sharedgs(n,M,b,x,r,iter,eps,maxit)

    ! Shared memory Gauss-Seidel method, using OpenMP routines.

    ! SCALAR ARGUMENTS
    integer :: n,iter,maxit
    real(8) :: eps

    ! ARRAY ARGUMENTS
    real(8) :: M(n,n),b(n),x(n),r
       
    intent(in   ) :: M,b,n,eps,maxit
    intent(out  ) :: r,iter
    intent(inout) :: x

    ! LOCAL SCALARS
    integer :: i,np,tn,start,end,it

    ! LOCAL ARRAYS
    real(8) :: res(n)

    ! FUNCTIONS
    integer omp_get_num_threads,omp_get_thread_num

    do i = 1,n
       res(i) = 1.0D+20
    end do

    iter = 0

    !$OMP PARALLEL PRIVATE(np,tn,start,end,it)

    np = omp_get_num_threads()

    tn = omp_get_thread_num()

    start = 1 + ceil(tn * n / (1.0D0 * np))

    if (tn .eq. np - 1) then
       end = n
    else
       end = ceil((tn + 1) * n / (1.0D0 * np))
    end if

    call partialgs(n,M,b,x,res,start,end,it,eps,maxit)

    !$OMP ATOMIC
    iter = iter + it

    !$OMP END PARALLEL

    r = maxval(res)

  end subroutine sharedgs

  subroutine partialgs(n,M,b,x,r,start,end,iter,eps,maxit)

    ! SCALAR ARGUMENTS
    integer :: n,start,end,iter,maxit
    real(8) :: eps

    ! ARRAY ARGUMENTS
    real(8) :: M(n,n),b(n),x(n),r(n)
       
    intent(in   ) :: M,b,n,eps,start,end,maxit
    intent(out  ) :: r,iter
    intent(inout) :: x
    
    ! LOCAL SCALARS
    integer :: i
    real(8) :: res

    ! LOCAL ARRAYS
    real(8) :: prev(start:end)

    iter = 0

    res = maxval(r)

    do while (res .gt. eps .and. iter .lt. maxit)

       do i = start,end
          prev(i) = x(i)
       end do

       call updatex(start,end,n,M,b,x)
       
       ! STOPPING CRITERIUM
       
       res = 0.0D0
       do i = start,end
          r(i) = abs(prev(i) - x(i))
          res = max(res,r(i))
       end do
       
       if (res .lt. eps) then
          res = maxval(r)
       end if

       iter = iter + 1
       
    end do
    
  end subroutine partialgs

END MODULE
