module gentests

  implicit none

contains

! ---------------------------------------------------------- !
! ---------------------------------------------------------- !

  subroutine genDDRandSystem(n,A,b,sol)

    implicit none
    
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
          A(i,j) = - n + 2 * n * rnumber
          s = s + abs(A(i,j))
       end do
       A(i,i) = s + 1.0D0
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

! ---------------------------------------------------------- !
! ---------------------------------------------------------- !

  subroutine genRandStochSystem(n,A,b,sol)
    
    implicit none

    ! Generates a linear system associated with the gradient of a
    ! least squares problem created with a stochastic matrix (See
    ! Nesterov's paper).

    ! PARAMETERS
    real(8), parameter :: density = 0.3D0
    real(8), parameter :: rho = 2.0D0

    ! ARGUMENTS
    integer :: n
    real(8) :: A(n,n),b(n),sol(n)

    intent(out) :: A,b,sol
    intent(in ) :: n

    ! LOCAL SCALARS
    integer :: i,j,k,nlinks,start,pos
    real(8) :: rnumber,s

    ! LOCAL ARRAYS
    real(8) :: Eb(n,n)
    integer :: E(n,n),links(n * n)
    
    start = 1

    do i = 1,n*n
       links(i) = i
    end do

    do j = 1,n
       do i = 1,n
          E(i,j) = 0
       end do
    end do

    do i = 1,n
       links((i - 1) * n + i) = links(start)
       start = start + 1
    end do

    nlinks = 0

    do while (nlinks / (1.0D0 * n ** 2) .lt. density)

       CALL RANDOM_NUMBER(rnumber)
       pos = INT(rnumber * (n ** 2 - start)) + start
       
       i = 1 + INT((links(pos) - 1) / n)
       j = links(pos) - (i - 1) * n

       E(i,j) = 1
       links(pos) = links(start)
       start = start + 1

       nlinks = nlinks + 1

    end do

    do j = 1,n
       s = 0.0D0
       do i = 1,n
          s = s + E(i,j)
       end do
       if (s .eq. 0.0D0) then
          do i = 1,n
             Eb(i,j) = 1.0D0 / (1.0D0 * n)
          end do
       else
          do i = 1,n
             Eb(i,j) = E(i,j) / s
          end do
       end if
    end do

    do i = 1,n
       do j = 1,n          
          s = - Eb(j,i) - Eb(i,j) + 2.0D0 * rho
          if (i .eq. j) then
             s = s + 1.0D0
          end if
          do k = 1,n
             s = s + Eb(k,i) * Eb(k,j)
          end do
          A(i,j) = s
       end do
    end do

    do i = 1,n
       b(i) = 2.0D0 * rho
    end do

  end subroutine genRandStochSystem

! ---------------------------------------------------------- !
! ---------------------------------------------------------- !

  subroutine genRandSDPSystem(n,A,b,sol,semiprob,magnitude)

    ! This subroutine creates a symmetric linear system.
    !
    ! On entry:
    !
    ! n: dimension of the problem
    !
    ! semiprob: probability that an eigenvalue of the matrix A will be
    ! "near" zero
    !
    ! magnitude: range of the random values. If negative, the matrix A
    ! may be indefinite
    !
    ! On exit
    !
    ! The (n x n) system Ax = b with solution 'sol'


    ! ARGUMENTS
    integer :: n
    real(8) :: magnitude,semiprob
    real(8) :: A(n,n),b(n),sol(n)

    intent(out) :: A,b,sol
    intent(in ) :: n,magnitude,semiprob

    ! LOCAL SCALARS
    integer :: i,j
    real(8) :: rnumber

    ! Creates a symmetric matrix
    call genRandHouseDP(n,A,semiprob,magnitude)

    do i = 1,n
       CALL RANDOM_NUMBER(rnumber)
       sol(i) = - abs(magnitude) + 2 * abs(magnitude) * rnumber
       b(i) = 0.0D0
    end do

    do j = 1,n
       do i = 1,n
          b(i) = b(i) + A(i,j) * sol(j)
       end do
    end do       

  end subroutine genRandSDPSystem

! ---------------------------------------------------------- !
! ---------------------------------------------------------- !

  subroutine genRandCholL(n,L,prob)
    
    ! This subroutine does not work!!!!

    ! PARAMETERS
    real(8),parameter :: epsdiag = 1.0D-5

    ! ARRAY ARGUMENTS
    real(8) :: L(n,n)

    ! SCALAR ARGUMENTS
    integer :: n
    real(8) :: prob

    intent( in) :: prob,n
    intent(out) :: L

    ! LOCAL SCALARS
    integer :: i,j
    real(8) :: rnumber

    do j = 1,n
       do i = 1,j - 1 
          L(i,j) = 0.0D0
       end do
       ! The diagonal element must be non-negative
       CALL RANDOM_NUMBER(rnumber)
       if (rnumber .lt. prob) then
          CALL RANDOM_NUMBER(rnumber)
          L(j,j) = (1.0D0 + rnumber) * epsdiag
       else
          CALL RANDOM_NUMBER(rnumber)
          L(j,j) = 1.0D0 + n * rnumber
       end if
       do i = j+1,n
          CALL RANDOM_NUMBER(rnumber)
          L(i,j) = - n + 2.0D0 * n * rnumber
       end do
    end do

  end subroutine genRandCholL

! ---------------------------------------------------------- !
! ---------------------------------------------------------- !

  subroutine genHouseholderOrt(n,Q,magnitude)
    ! This subroutine generates an randomly generated orthogonal
    ! Householder matrix:
    !
    ! I - 2 * (u * u^T) / (u^T * u)
    !
    ! and returns such matrix in Q

    ! SCALAR ARGUMENTS
    integer :: n
    real(8) :: magnitude

    ! ARRAY ARGUMENTS
    real(8) :: Q(n,n)

    intent( in) :: n,magnitude
    intent(out) :: Q

    ! LOCAL SCALARS
    integer :: i,j
    real(8) :: rnumber,dp
    
    ! LOCAL ARRAYS
    real(8) :: u(n)

    dp = 0.0D0
    do i = 1,n
       CALL RANDOM_NUMBER(rnumber)
       u(i) = - abs(magnitude) + 2.0D0 * abs(magnitude) * rnumber
       dp = dp + u(i) ** 2.0D0
    end do
    
    do j = 1,n
       do i = 1,n
          Q(i,j) = - 2.0D0 * u(i) * u(j) / dp
       end do
       Q(j,j) = Q(j,j) + 1.0D0
    end do

  end subroutine genHouseholderOrt

! ---------------------------------------------------------- !
! ---------------------------------------------------------- !

  subroutine genRandDP(n,A,prob)

    ! ARRAY ARGUMENTS
    real(8) :: A(n,n)

    ! SCALAR ARGUMENTS
    integer :: n
    real(8) :: prob

    intent( in) :: prob,n
    intent(out) :: A

    ! LOCAL SCALARS
    integer :: i,j,k

    ! LOCAL ARRAYS
    real(8) :: L(n,n)

    call genRandCholL(n,L,prob)

    do j = 1,n
       do i = 1,n
          A(i,j) = 0.0D0
       end do
    end do

    do j = 1,n
       do i = 1,j
          do k = 1,n
             A(k,j) = A(k,j) + L(k,i) * L(j,i)
          end do
       end do
    end do

  end subroutine genRandDP
  
! ---------------------------------------------------------- !
! ---------------------------------------------------------- !

  subroutine genRandHouseDP(n,A,prob,magnitude)
    ! This subroutine generates a random symmetric matrix A with
    ! dimension n.
    !
    ! prob: is the probability of a eigenvalue to be very near to zero
    !
    ! magnitude: if is positive, defines the range of random
    ! values. If is negative, also defines this range, but the matrix
    ! has chance to be indefinite


    ! PARAMETERS
    real(8),parameter :: epsdiag = 1.0D-5

    ! ARRAY ARGUMENTS
    real(8) :: A(n,n)

    ! SCALAR ARGUMENTS
    integer :: n
    real(8) :: prob,magnitude

    intent( in) :: prob,n,magnitude
    intent(out) :: A

    ! LOCAL SCALARS
    integer :: i,j,k
    real(8) :: rnumber

    ! LOCAL ARRAYS
    real(8) :: Q(n,n),d(n)

    call genHouseholderOrt(n,Q,magnitude)

    ! Generates a random diagonal matrix D and calculates Q * D * Q^t

    do i = 1,n
       CALL RANDOM_NUMBER(rnumber)
       if (rnumber .lt. prob) then
          ! Semi-definite positiveness probability
          CALL RANDOM_NUMBER(rnumber)
          d(i) = (1.0D0 + rnumber) * epsdiag
       else if (magnitude .gt. 0.0D0) then
          ! Matrix is positive definite
          CALL RANDOM_NUMBER(rnumber)
          d(i) = epsdiag + magnitude * rnumber
       else
          ! Matrix is indefinite
          CALL RANDOM_NUMBER(rnumber)
          d(i) = - abs(magnitude) + 2.0D0 * abs(magnitude) * rnumber          
       end if
    end do
    
    do j = 1,n
       do i = 1,n
          A(i,j) = 0.0D0
       end do
    end do

    do j = 1,n
       do i = 1,n
          do k = 1,n
             A(k,j) = A(k,j) + Q(k,i) * Q(j,i) * d(i)
          end do
       end do
    end do

  end subroutine genRandHouseDP

end module gentests
