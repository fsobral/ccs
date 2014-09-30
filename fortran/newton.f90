module newton

  implicit none

contains

  subroutine tnmin_serial(n,x,func,grad,hess,epslin,epsopt,maxit,flag)
    
    use gs

    implicit none

    ! This subroutine uses the truncated Newton's method for
    ! minimizing a function without constraints. This serial version
    ! uses serial Gauss-Seidel for solving the linear systems.

    ! SCALAR ARGUMENTS
    integer :: n,maxit,flag
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
    integer :: i,iter,initer,inmaxit,status
    real(8) :: f,gnorm,r,inepsilon

    ! LOCAL ARRAYS
    real(8), allocatable :: d(:),g(:),h(:,:)
    
    
    allocate(d(n),g(n),h(n,n),STAT=status)

    if (status .ne. 0) then
       write(*,*) 'Memory problems'
       flag = -1

       stop
    end if

    ! Initial calculations

    call func(n,x,f)
    call grad(n,x,g)

    gnorm = 0.0D0
    do i = 1,n
       g(i) = - g(i)
       gnorm = max(gnorm,abs(g(i)))
    end do

    iter = 1

    inepsilon = 1.0D0;

    inmaxit = maxit / 10

    ! Main loop

    do while (gnorm .gt. epsopt .and. iter .le. maxit) 

       write(*,FMT=9000)iter,f,gnorm

       ! Calculates the step

       call hess(n,x,h)

       inepsilon = max(epsopt,inepsilon / 2.0D0)

       call serialgs(n,h,g,d,r,initer,inepsilon,inmaxit)

       write(*,FMT=9010)initer,r

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

    write(*,FMT=9000)iter,f,gnorm

    deallocate(d,g,h)

9000 FORMAT('Iter:',1X,I0.5,1X,'F=',1X,E10.4,1X,'NORMG=',1X,E10.4)
9010 FORMAT('Inner GS took',1X,I0.5,1X,'iterations.',1X,'RESIDUAL=',1X,E10.4)

  end subroutine tnmin_serial

  ! ---------------------------------------------------------------- !
  ! ---------------------------------------------------------------- !

  subroutine tnmin_async(n,x,d,func,grad,hess,epslin,epsopt,maxit,flag)

    use gs ! function 'ceil'

    implicit none

    ! SCALAR ARGUMENTS
    integer :: n,maxit,flag
    real(8) :: epslin,epsopt

    ! ARRAY ARGUMENTS
    real(8) :: x(n),d(n)

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
    integer :: np,tn,start,end

    ! FUNCTIONS
    integer omp_get_num_threads,omp_get_thread_num

    !$OMP PARALLEL PRIVATE(np,tn,start,end)

    np = omp_get_num_threads()

    tn = omp_get_thread_num()

    start = 1 + ceil(tn * n / (1.0D0 * np))

    if (tn .eq. np - 1) then
       end = n
    else
       end = ceil((tn + 1) * n / (1.0D0 * np))
    end if

    call tnmin_partial(n,x,d,func,grad,hess,start,end,epslin,epsopt,&
         maxit,flag)

    !$OMP END PARALLEL

  end subroutine tnmin_async
  

  ! ---------------------------------------------------------------- !
  ! ---------------------------------------------------------------- !

  subroutine tnmin_partial(n,x,d,func,grad,hess,start,end,epslin,epsopt,&
       maxit,flag)
    
    use gs

    implicit none

    ! SCALAR ARGUMENTS
    integer :: start,n,maxit,flag,end
    real(8) :: epslin,epsopt

    ! ARRAY ARGUMENTS
    real(8) :: x(n),d(n)

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
    integer :: i,iter,initer,inmaxit,status,proc
    real(8) :: f,gnorm,inepsilon

    ! LOCAL ARRAYS
    real(8), allocatable :: g(:),h(:,:)
    real(8) :: r(start:end)

    ! FUNCTIONS
    integer omp_get_thread_num

    allocate(g(n),h(n,n),STAT=status)

    if (status .ne. 0) then
       write(*,*) 'Memory problems'
       flag = -1

       stop
    end if

    ! Initial calculations

    call func(n,x,f)
    call grad(n,x,g)

    gnorm = 0.0D0
    do i = 1,n
       g(i) = - g(i)
       gnorm = max(gnorm,abs(g(i)))
       d(i) = 0.0D0
    end do

    iter = 1

    inepsilon = 5.0D-1;

    inmaxit = 10

    proc = omp_get_thread_num()

    ! Main loop

    do while (gnorm .gt. epsopt .and. iter .le. maxit) 

       write(*,FMT=8000) proc,iter,f,gnorm

       ! Calculates the step

       call hess(n,x,h)

       inepsilon = min(gnorm,sqrt(inepsilon))

       ! Teste do Emerson: apenas 1 iteracao
       !call updatex(start,end,n,h,g,d)
       ! Meu teste: criterio decrescente
       do i = start,end
          r(i) = 1.0D+20
       end do
       call partialgs(n,h,g,d,r,start,end,initer,inepsilon,inmaxit)

       write(*,FMT=8010) proc,initer,maxval(r)

       ! Updates the point

       do i = start,end
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

    write(*,FMT=8000) proc,iter,f,gnorm

    deallocate(g,h)

8000 FORMAT('Proc',1X,I2,1X,'Iter:',1X,I0.10,1X,'F=',1X,E10.4,1X,&
          'NORMG=',1X,E10.4)
8010 FORMAT('Proc',1X,I2,1X,'Inner GS took',1X,I0.5,1X,&
          'iterations.',1X,'RESIDUAL=',1X,E10.4)

  end subroutine tnmin_partial
  

end module newton
