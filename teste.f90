

program teste

    ! use mpi
   
    ! call system('mkdir temp')
    ! open(10,file='temp/test.txt',status="replace")
    ! write(10,*) 'lalalala'
    ! print*,'a'
    
    implicit none
    real :: a(10) 
    real, parameter :: PI = 3.1415, kb = 1.38064852E-23
    integer :: b(10)

    b = [0,1,0,1,0,0,1,0,1,0]

    a = [0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0]
    print*, "a = ", a*180.0
    a = a*pi
    print*, "rem a por 2*pi", modulo(a,2*pi)*180/pi

    include "teste_include.f90"
    

end program teste