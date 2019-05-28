! Este módulo imprime a matriz
module matprint
    contains
    subroutine printm(m,r)
        use mod1
        real(dp), dimension(:,:), intent(in) :: m
        integer,dimension(2) :: s       
        integer :: i,j  
        character(*), optional, intent(in) :: r 

        s = shape(m)
        if(present(r)) then
            print*, r,' ='
        else
            print*, ' ='
        end if
        
        do i = 1,s(1)
            do j = 1,s(2)
                write(*,'(EN14.3)',advance='no') m(i,j)
            end do
            write(*,*)
        end do
        write(*,*) ' '
    end subroutine printm
end module matprint    
