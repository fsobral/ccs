program test_newton

  use newton

  implicit none

  integer :: n,maxit,flag,i
  real(8),allocatable :: x(:)
  real(8) :: epslin,epsopt

  epslin = 1.0D-8
  epsopt = 1.0D-8
  maxit = 1000

  write(*,*) 'Digite n'
  read(*,*) n

  allocate(x(n))

  do i = 1,n
     x(i) = 5.0D0
  end do

  call tnmin_serial(n,x,func,grad,hess,epslin,epsopt,maxit,flag)

  deallocate(x)

contains

  subroutine func(n,x,f)
    integer, intent( in) :: n
    real(8), intent(out) :: f
    real(8), intent( in) :: x(n)

    integer :: i

    f = 0.0D0
    do i = 1,n
       f = f + x(i) ** 2.0D0
    end do

  end subroutine func
  
  subroutine grad(n,x,g)
    integer, intent( in) :: n
    real(8), intent(out) :: g(n)
    real(8), intent( in) :: x(n)

    integer :: i

    do i = 1,n
       g(i) = 2.0D0 * x(i)
    end do

  end subroutine grad
  
  subroutine hess(n,x,h)
    integer, intent( in) :: n
    real(8), intent(out) :: h(n,n)
    real(8), intent( in) :: x(n)

    integer :: i,j

    do j = 1,n
       do i = 1,j-1
          h(i,j) = 0.0D0
          h(j,i) = 0.0D0
       end do
       h(j,j) = 2.0D0
    end do

  end subroutine hess

end program test_newton
