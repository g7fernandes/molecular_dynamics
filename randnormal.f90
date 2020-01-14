module randnormal
    contains
    
    function GaussDeviate() 
    ! BOX MULLER 
        real :: GaussDeviate
        real :: a1, a2, s, r, b1
        logical, save :: iset = .false.
        real, save :: b2
        
        r = 3.0
        
        if (.not. iset) then
            do while (r >= 1)
            ! Generate random numbers between -1 e 1
                call random_seed()
                call random_number(a1)
                call random_number(a2)
                a1 = 2*a1-1            
                a2 = 2*a2-1
                r = a1**2 + a2**2
            end do
            s = sqrt(-2*log(r)/r) ! polar transformation
            b1 = a1*s
            b2 = a2*s
            iset = .true.
            GaussDeviate = b1
        else
            iset = .false.
            GaussDeviate = b2
        end if 
    
    end function

end module randnormal
