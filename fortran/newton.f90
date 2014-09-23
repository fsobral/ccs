module newton

  implicit none

contains

  subroutine tnmin_serial(n,x,func,grad,hess,epslin,epsopt,maxit)
    
    use gs

    implicit none

    ! This subroutine uses the truncated Newton's method for
    ! minimizing a function without constraints. This serial version
    ! uses serial Gauss-Seidel for solving the linear systems.

    ! SCALAR ARGUMENTS
    integer :: n,maxit
    real(8) :: epslin,epsopt

    ! ARRAY ARGUMENTS
    real(8) :: x(n)

    ! FUNCTION ARGUMENTS
    INTERFACE
       subroutine func(n,x,f)
         integer, intent( in) :: n
         real(8), intent(out) :: f
         real(8), intent( in) :: x(n)
       end subroutine func
    end INTERFACE

    INTERFACE
       subroutine grad(n,x,g)
         integer, intent( in) :: n
         real(8), intent(out) :: g(n)
         real(8), intent( in) :: x(n)
       end subroutine grad
    end INTERFACE
    
    INTERFACE
       subroutine hess(n,x,h)
         integer, intent( in) :: n
         real(8), intent(out) :: h(n,n)
         real(8), intent( in) :: x(n)
       end subroutine hess
    end INTERFACE

    ! LOCAL SCALARS
    integer :: i,iter,initer,inmaxit
    real(8) :: f,gnorm,r

    ! LOCAL ARRAYS
    real(8), allocatable :: d(:),g(:),h(:,:)
    
    
    allocate(d(n),g(n),h(n,n)) ! Pegar o erro aqui!

    call func(n,x,f)
    call grad(n,x,g)

    gnorm = 0.0D0
    do i = 1,n
       g(i) = - g(i)
       gnorm = max(gnorm,abs(g(i)))
    end do

    iter = 0

    do while (gnorm .gt. epsopt .and. iter .le. maxit) 

       write(*,FMT=9000)iter,f,gnorm

       call hess(n,x,h)

       ! Calculates the step

       call serialgs(n,h,g,d,r,initer,epsopt,inmaxit)

       ! Updates the point

       do i = 1,n
          x(i) = x(i) + d(i)
       end do

       call func(n,x,f)
       call grad(n,x,g)

       gnorm = 0.0D0
       do i = 1,n
          g(i) = - g(i)
          gnorm = max(gnorm,abs(g(i)))
       end do

       iter = iter + 1

    end do

    deallocate(d,g,h)

9000 FORMAT('Iter:',1X,I010,1X,'F=',1X,E10.4,1X,'NORMG=',1X,E10.4)

  end subroutine tnmin_serial

end module newton
