module mod1
!     use linkedlist
    integer, parameter:: dp=kind(0.d0)                   ! double precision
    ! vetor de strings
    type string
        character(len=:), allocatable :: str
        integer :: quanty                   
    end type string
end module mod1
