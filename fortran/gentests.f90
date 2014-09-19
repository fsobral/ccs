module gentests

  implicit none

contains

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

  subroutine genRandPositiveMatrix(n,A,b,sol)

    ! PARAMETERS
    real(8) :: semiprob = 8.0D-1

    ! ARGUMENTS
    integer :: n
    real(8) :: A(n,n),b(n),sol(n)

    intent(out) :: A,b,sol
    intent(in ) :: n

    ! LOCAL ARRAYS
    real(8) :: L(n,n)
    
    ! LOCAL SCALARS
    integer :: i,j,k
    real(8) :: rnumber

    do j = 1,n
       ! The diagonal element must be non-negative
       CALL RANDOM_NUMBER(rnumber)
       L(j,j) = 1.0D0 + n * rnumber
       CALL RANDOM_NUMBER(rnumber)
       if (rnumber .le. semiprob) then
          L(j,j) = 1.0D-1 + rnumber * 1.0D-1
       end if
       do i = j+1,n
          CALL RANDOM_NUMBER(rnumber)
          L(i,j) = - n + 2.0D0 * n * rnumber
       end do
    end do

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

  end subroutine genRandPositiveMatrix

end module gentests
