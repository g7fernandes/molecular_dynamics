program teste
    ! use mod1
   
    ! call system('mkdir temp')
    ! open(10,file='temp/test.txt',status="replace")
    ! write(10,*) 'lalalala'
    ! print*,'a'
    
    implicit none
    real :: a,b,c,d
    character(LEN=*), parameter :: fmt = '(f16.0, ", ",f32.16, ", ",f32.16, ", ",f32.16, ", ",f32.16 )'
    character(170) :: out

    a = 12.0
    b = 1.2
    c = sqrt(2.0)
    d = -sqrt(3.0)
    write(out,fmt) a,b,c,d,d
    out = trim(out)
    write(*,*) out
    
    
    
end program teste
