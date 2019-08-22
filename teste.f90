program teste
    use mpi
   
    ! call system('mkdir temp')
    ! open(10,file='temp/test.txt',status="replace")
    ! write(10,*) 'lalalala'
    ! print*,'a'
    
    implicit none
    real :: sen(2), res(8)
    integer ( kind = 4 ) status(MPI_STATUS_SIZE)
    integer ( kind = 4 ) ierr
    integer ( kind = 4 ) np
    integer ( kind = 4 ) tag
    integer ( kind = 4 ) id

    call MPI_Init ( ierr )
    call MPI_Comm_rank ( MPI_COMM_WORLD, id, ierr )                            
    call MPI_Comm_size ( MPI_COMM_WORLD, np, ierr )

    res = [0,0,0,0,0,0,0,0]
    sen = [1.0,1.0]
    sen = sen*id +1

    call MPI_ALLGATHER(sen, 2, MPI_REAL, res, 2, MPI_REAL , MPI_COMM_WORLD, ierr)

    print '("id ", f3.1, " res ", f3.1," ", f3.1," ", f3.1," ", f3.1," ", f3.1," ", f3.1," ", f3.1," ", f3.1)', &
        id, res(1),res(2),res(3), res(4) ,res(5) ,res(6), res(7), res(8) 

    
    
    call MPI_Finalize ( ierr )
end program teste
