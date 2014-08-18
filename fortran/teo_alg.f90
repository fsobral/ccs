PROGRAM teoric_algorithm

  use gs

  implicit none

  ! PARAMETERS
  integer, parameter :: MAXIT = 1000

  ! LOCAL ARRAYS
  real(8), allocatable :: M(:,:),b(:),x(:),rv(:),sol(:)
  integer, allocatable :: seed(:)

  ! LOCAL SCALARS
  real(8) :: r,mr,msdiff,mt,miter
  integer :: i,k,n,status,rsize,t1,t2,clock_rate,clock_max,iter
  
  write(*,*) "Type dimension..."
  read(*,*) n

  allocate(M(n,n),b(n),x(n),rv(n),sol(n),STAT=status)

  if (status .ne. 0) then
     write(*,*) "Memory problems"
     stop
  end if

  CALL RANDOM_SEED(SIZE=rsize)

  allocate(seed(rsize))
  do i = 1,rsize
     seed(i) = 123456789
  end do

  CALL RANDOM_SEED(PUT=seed)

  write(*,*) "Starting serial CCS..."

  call genDDRandSystem(n,M,b,sol)

  do i = 1,n
     x(i) = 0.0D0
  end do

  CALL SYSTEM_CLOCK(t1,clock_rate,clock_max)

  CALL serialccs(5*n,n,M,b,x,r,iter,1.0D-5,MAXIT)
     
  CALL SYSTEM_CLOCK(t2,clock_rate,clock_max)

  write(*,9000) r,(t2 - t1) / (1.0 * clock_rate),iter,maxval(abs(x - sol))

9000 FORMAT('RESIDUAL:',1X,D12.5,1X,'TIME (in secs):',1X,F7.4,&
          1X,'ITERATIONS:',1X,I8,1X,'DIFF',1X,D12.5)

contains

  subroutine genDDRandSystem(n,A,b,sol)
    
    ! Generates a random diagonal dominant linear system

    ! ARGUMENTS
    integer :: n
    real(8) :: A(n,n),b(n),sol(n)

    intent(out) :: A,b,sol
    intent(in ) :: n

    ! LOCAL SCALARS
    integer :: i,j
    real(8) :: rnumber,s

    do i = 1,n
       s = 0.0D0
       do j = 2,n
          CALL RANDOM_NUMBER(rnumber)
          M(i,j) = - n + 2 * n * rnumber
          s = s + abs(M(i,j))
       end do
       M(i,i) = s
    end do

    do i = 1,n
       CALL RANDOM_NUMBER(rnumber)
       sol(i) = - n + 2 * n * rnumber
       b(i) = 0.0D0
    end do

    do j = 1,n
       do i = 1,n
          b(i) = b(i) + A(i,j) * sol(j)
       end do
    end do       

  end subroutine genDDRandSystem

END PROGRAM teoric_algorithm
