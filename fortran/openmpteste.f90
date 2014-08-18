program teste

  integer i
  integer omp_get_thread_num,omp_get_num_threads

  real(8) :: x(15)
  integer :: start,end

  !$OMP PARALLEL PRIVATE(i,start,end)
  i = omp_get_thread_num()
  write(*,*) 'hello world', i

  start = 1 + i * INT(15 / omp_get_num_threads())
  if (i .eq. omp_get_num_threads() - 1) then
     end = 15
  else
     end = (i + 1) * INT(15 / omp_get_num_threads())
  end if

  call alteravetor(end - start + 1,x(start:end))
  !$OMP END PARALLEL

  write(*,2000) x

2000 FORMAT(/,'Point:',/,3(1X,E20.10))

contains

  subroutine alteravetor(n,x)

    integer :: n
    real(8) :: x(n)

    integer :: i

    do i = 1,n
       x(i) = omp_get_thread_num()
    end do

  end subroutine alteravetor

end program teste
