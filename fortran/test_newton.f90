program test_newton

  use newton

  implicit none

  integer :: n,maxit,flag,i
  real(8),allocatable :: x(:),d(:)
  real(8) :: epslin,epsopt

  epslin = 1.0D-8
  epsopt = 1.0D-8
  maxit = 100000000

  write(*,*) 'Digite n'
  read(*,*) n

  allocate(x(n),d(n))

  do i = 1,n
     x(i) = -1.0D0
     d(i) =  0.0D0
  end do

  call tnmin_async(n,x,d,func,grad,hess,epslin,epsopt,maxit,flag)

  deallocate(x,d)

contains

  subroutine func(n,x,f)
    integer, intent( in) :: n
    real(8), intent(out) :: f
    real(8), intent( in) :: x(n)

    integer :: i

    f = 0.0D0
    do i = 1,n - 1
       f = f + (1.0D0 - x(i)) ** 2.0D0 &
             + 100.0D0 * (x(i + 1) - x(i) ** 2.0D0) ** 2.0D0
    end do

  end subroutine func
  
  subroutine grad(n,x,g)
    integer, intent( in) :: n
    real(8), intent(out) :: g(n)
    real(8), intent( in) :: x(n)

    integer :: i

    do i = 1,n - 1
       g(i) = - 2.0D0 * (1.0D0 - x(i)) &
              - 4.0D+2 * (x(i + 1) - x(i) ** 2.0D0) * x(i)
    end do
    g(n) = 0.0D0
    do i = 2,n
       g(i) = g(i) + 2.0D+2 * (x(i) - x(i - 1) ** 2.0D0)
    end do

  end subroutine grad
  
  subroutine hess(n,x,h)
    integer, intent( in) :: n
    real(8), intent(out) :: h(n,n)
    real(8), intent( in) :: x(n)

    integer :: i,j

    do j = 2,n - 1
       h(j,j + 1) = - 4.0D+2 * x(j)
       h(j + 1,j) = h(j,j + 1)
       do i = j + 2,n
          h(i,j) = 0.0D0
          h(j,i) = 0.0D0
       end do
       h(j,j) = 2.02D+2 + 1.2D+3 * x(j) ** 2.0D0 - 4.0D+2 * x(j + 1)
    end do
    h(1,1) = 2.0D0 + 1.2D+3 * x(1) ** 2.0D0 - 4.0D+2 * x(2)
    h(n,n) = 2.0D+2

  end subroutine hess

end program test_newton
