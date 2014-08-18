PROGRAM gs_comp

  use gs

  implicit none

  ! PARAMETERS
  integer, parameter :: MAXIT = 1000
  real(8), parameter :: TOL = 1.0D-5

  ! LOCAL ARRAYS
  real(8), allocatable :: M(:,:),b(:),x(:),rv(:),sol(:)
  integer, allocatable :: seed(:)

  ! LOCAL SCALARS
  real(8) :: r,mr,msdiff,mt,miter
  integer :: i,k,n,status,rsize,t1,t2,clock_rate,clock_max,iter, &
       trials
  
  write(*,*) "Type dimension and number of trials..."
  read(*,*) n,trials

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

  write(*,*) "Starting serial..."

  CALL RANDOM_SEED(PUT=seed)

  mt = 0.0D0
  mr = 0.0D0
  miter = 0.0D0
  msdiff = 0.0D0

  do k = 1,trials

     call genDDRandSystem(n,M,b,sol)

     do i = 1,n
        x(i) = 0.0D0
     end do

     CALL SYSTEM_CLOCK(t1,clock_rate,clock_max)

     CALL serialgs(n,M,b,x,r,iter,TOL,MAXIT)

     CALL SYSTEM_CLOCK(t2,clock_rate,clock_max)
     
     miter = miter + iter / (1.0D0 * trials)
     mt = mt + (t2 - t1) / (1.0D0 * clock_rate * trials)
     mr = mr + r / (1.0D0 * trials)
     msdiff = msdiff + maxval(abs(x - sol)) / (1.0D0 * trials)

  end do

  write(*,9000) mr,mt,INT(miter),msdiff
     
  write(*,*) "Starting parallel..."

  CALL RANDOM_SEED(PUT=seed)

  mt = 0.0D0
  mr = 0.0D0
  miter = 0.0D0
  msdiff = 0.0D0

  do k = 1,trials

     call genDDRandSystem(n,M,b,sol)

     do i = 1,n
        x(i) = 0.0D0
     end do

     CALL SYSTEM_CLOCK(t1,clock_rate,clock_max)

     CALL sharedgs(n,M,b,x,r,iter,TOL,MAXIT)
 
     CALL SYSTEM_CLOCK(t2,clock_rate,clock_max)

     miter = miter + iter / (1.0D0 * trials)
     mt = mt + (t2 - t1) / (1.0D0 * clock_rate * trials)
     mr = mr + r / (1.0D0 * trials)
     msdiff = msdiff + maxval(abs(x - sol)) / (1.0D0 * trials)

  end do

  write(*,9000) mr,mt,INT(miter),msdiff
!  write(*,2000) x
!  write(*,2000) sol

2000 FORMAT(/,'Point:',/,3(1X,E20.10))
9000 FORMAT('RESIDUAL:',1X,D12.5,1X,'TIME (in secs):',1X,F10.4,&
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
       M(i,i) = s + 1.0D0
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

END PROGRAM gs_comp
