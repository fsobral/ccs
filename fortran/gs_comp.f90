PROGRAM gs_comp

  use gs
  use gentests

  implicit none

  ! PARAMETERS
  integer, parameter :: MAXIT = 100000
  real(8), parameter :: TOL = 1.0D-5

  ! LOCAL ARRAYS
  real(8), allocatable :: M(:,:),b(:),x(:),rv(:),sol(:)
  integer, allocatable :: seed(:)

  ! LOCAL SCALARS
  real(8) :: r,mr,msdiff,mt,miter
  integer :: i,k,n,status,rsize,t1,t2,clock_rate,clock_max,iter, &
       trials,strials,nfailed,type
  
  write(*,*) 'Type dimension and number of trials...'
  read(*,*) n,trials
  write(*,*) '1- Diagonal Dominant 2- Stochastic Least Squares'
  write(*,*) '3- Symmetric Definite Positive'
  write(*,*) '4- Symmetric Semi-Definite Positive'
  write(*,*) '5- Symmetric Indefinite'
  read(*,*) type

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

  nfailed = 0
  mt = 0.0D0
  mr = 0.0D0
  miter = 0.0D0
  msdiff = 0.0D0

  do k = 1,trials

     if (type .eq. 1) then
        call genDDRandSystem(n,M,b,sol)
     else if (type .eq. 2) then
        call genRandStochSystem(n,M,b,sol)
     else if (type .eq. 3) then
        call genRandSDPSystem(n,M,b,sol,0.0D0,1.0D0 * n)
     else if (type .eq. 4) then
        call genRandSDPSystem(n,M,b,sol,4.0D-1,1.0D0 * n)
     else if (type .eq. 5) then
        call genRandSDPSystem(n,M,b,sol,2.0D-1,- 1.0D0 * n)
     else
        write(*,*) 'Invalid type.'
        exit
     end if

     do i = 1,n
        x(i) = 0.0D0
     end do

     CALL SYSTEM_CLOCK(t1,clock_rate,clock_max)

     CALL serialgs(n,M,b,x,r,iter,TOL,MAXIT)

     CALL SYSTEM_CLOCK(t2,clock_rate,clock_max)
     
     if (r .le. TOL) then
        miter = miter + iter
        mt = mt + (t2 - t1)
        mr = mr + r
        msdiff = msdiff + maxval(abs(x - sol))
     else
        nfailed = nfailed + 1
     end if

  end do

  strials = trials - nfailed

  write(*,9000) mr / (1.0D0 * strials),mt / (1.0D0 * clock_rate * strials), &
       INT(miter / (1.0D0 * strials)),msdiff / (1.0D0 * strials), &
       (1.0D+2 * nfailed) / trials
     
  write(*,*) "Starting parallel..."

  CALL RANDOM_SEED(PUT=seed)

  nfailed = 0
  mt = 0.0D0
  mr = 0.0D0
  miter = 0.0D0
  msdiff = 0.0D0

  do k = 1,trials

     if (type .eq. 1) then
        call genDDRandSystem(n,M,b,sol)
     else if (type .eq. 2) then
        call genRandStochSystem(n,M,b,sol)
     else if (type .eq. 3) then
        call genRandSDPSystem(n,M,b,sol,0.0D0,1.0D0 * n)
     else if (type .eq. 4) then
        call genRandSDPSystem(n,M,b,sol,4.0D-1,1.0D0 * n)
     else if (type .eq. 5) then
        call genRandSDPSystem(n,M,b,sol,2.0D-1,- 1.0D0 * n)
     else
        write(*,*) 'Invalid type.'
        exit
     end if

     do i = 1,n
        x(i) = 0.0D0
     end do

     CALL SYSTEM_CLOCK(t1,clock_rate,clock_max)

     CALL sharedgs(n,M,b,x,r,iter,TOL,MAXIT)
 
     CALL SYSTEM_CLOCK(t2,clock_rate,clock_max)

     if (r .le. TOL) then
        miter = miter + iter
        mt = mt + (t2 - t1)
        mr = mr + r
        msdiff = msdiff + maxval(abs(x - sol))
     else
        nfailed = nfailed + 1
     end if

  end do

  strials = trials - nfailed

  write(*,9000) mr / (1.0D0 * strials),mt / (1.0D0 * clock_rate * strials), &
       INT(miter / (1.0D0 * strials)),msdiff / (1.0D0 * strials), &
       (1.0D+2 * nfailed) / trials
     
!  write(*,2000) x
!  write(*,2000) sol

!2000 FORMAT(/,'Point:',/,3(1X,E20.10))
9000 FORMAT('RESIDUAL:',1X,D12.5,1X,'TIME (in secs):',1X,F10.4,&
          1X,'ITERATIONS:',1X,I8,1X,'DIFF',1X,D12.5,1X,&
          'FAILED:',1X,F6.2,'%')

END PROGRAM gs_comp
