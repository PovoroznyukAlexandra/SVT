! Solving a linear system with ILU preconditioning and PGMRES method
program linear_solver
    implicit none
    real*8, dimension(:), allocatable :: a, b, alu, w, sol, vv
    integer, dimension(:), allocatable :: ia, ja, jlu, ju, levs, jw
    integer :: iout, n, nnz, ierr, iwk, maxits, im, lfil
    integer :: i
    integer :: nmax = 264751, nzmax = 3541545, iounit = 21
    real*8 :: beg, fin, eps
    
    allocate(ia(nmax+1), ja(nzmax), a(nzmax))
    open(iounit, file = "/Users/deti/Desktop/3_курс/6_сем/СВТ/Task_1/4.dat")
    call readsk(nmax,nzmax,n,nnz,a,ja,ia,iounit,ierr)
    
    ! b_i = sin(i)
    allocate(b(n))
    do i = 1, n
        b(i) = sin(real(i))
    end do
    
    ! ILU
    lfil = 1
    iwk = nnz * 100
    allocate(jw(3*n), w(n), alu(iwk), jlu(iwk), levs(iwk), ju(n))
    call cpu_time(beg)
    call iluk(n,a,ja,ia,0,alu,jlu,ju,levs,iwk,w,jw,ierr)
    call cpu_time(fin)
    print *, "ILU time: ", fin-beg
    
    ! PGMRES
    im = 15
    eps = 1e-6
    maxits = n * n
    iout = 22
    open(iout, file = "/Users/deti/Desktop/3_курс/6_сем/СВТ/Task_1/logs.txt")
    allocate(sol(n), vv(n * (im + 1)))
    call cpu_time(beg)
    call pgmres(n, im, b, sol, vv, eps, maxits, iout, a, ja, ia, alu, jlu, ju, ierr)
    call cpu_time(fin)
    print *, "PGMRES time: ", fin-beg
    
    ! Free memory
    deallocate(ia, ja, a, jw, w, alu, jlu, levs, ju, sol, vv, b)
end program linear_solver