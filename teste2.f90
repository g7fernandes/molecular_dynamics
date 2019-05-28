! A derived type for storing data.
module mod1
    integer, parameter:: dp=kind(0.d0)                   ! double precision
    ! vetor de strings
    type string
        character(len=:), allocatable :: str
        integer :: quanty
    end type string
end module mod1    

module data
    use linkedlist
    implicit none

    private
    public :: data_t
    public :: data_ptr

    ! Data is stored in data_t
    type :: data_t
    real :: x
    end type data_t

    ! A trick to allow us to store pointers in the list
    type :: data_ptr
        type(data_t), pointer :: p
    end type data_ptr
    
end module data

module mod3
    use linkedlist  
    type :: container
        type(list_t), pointer :: list
    end type container
end module mod3
    
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

module mod2
    contains
    
    subroutine testes(m,ptr,ptrn)
        use linkedlist
        use mod3
        type(container),intent(in), allocatable,dimension(:,:) :: m
        integer,pointer :: ptr,ptrn
        print*,'linha 72'
        ptrn = transfer(list_get(m(2,1)%list), ptr)
        print *, 'Head node data:', ptrn
    end subroutine testes
    
    subroutine testchange(node2,m)
        use linkedlist
        use mod3
        type(container),intent(in), allocatable,dimension(:,:) :: m
        type(list_t), pointer :: node2
        call list_change(node2,m(2,1)%list)
    end subroutine testchange

    subroutine sub2(y)
        integer, intent(inout) :: y
        
        y  = y*2
    end subroutine sub2
    
    subroutine sub1(x,y)
        integer, intent(in) :: x
        integer, intent(out) :: y
        
        y = x
        call sub2(y)
        
        
    end subroutine sub1
    
    
end module mod2
! 
! ! A simple generic linked list test program
! program list_test
!   use linkedlist
!   use data
!   implicit none
! 
!   type(list_t), pointer :: ll => null()
!   type(data_t), target :: dat_a
!   type(data_t), target :: dat_b
!   type(data_ptr) :: ptr
! 
!   ! Initialize two data objects
!   dat_a%x = 17.5
!   dat_b%x = 3.0
! 
!   ! Initialize the list with dat_a
!   ptr%p => dat_a
!   call list_init(ll, DATA=transfer(ptr, list_data))
!   print *, 'Initializing list with data:', ptr%p
! 
!   ! Insert dat_b into the list
!   ptr%p => dat_b
!   call list_insert(ll, DATA=transfer(ptr, list_data))
!   print *, 'Inserting node with data:', ptr%p
! 
!   ! Get the head node
!   ptr = transfer(list_get(ll), ptr)
!   print *, 'Head node data:', ptr%p
! 
!   ! Get the next node
!   ptr = transfer(list_get(list_next(ll)), ptr)
!   print *, 'Second node data:', ptr%p
! 
!   ! Free the list
!   call list_free(ll)
! end program list_test 





program teste2
    use linkedlist
    use data
    use matprint
    use mod3
    use mod1
    use mod2
    use mpi 

    implicit none
!     type(container), allocatable,dimension(:,:) :: m
    integer,target :: d1,d2,d3,d4,d5,d6,d7,d8
    integer,pointer :: ptr
    character(5) :: c
    integer :: d
    logical :: laux
    real(dp) :: a(4), b(16)
    real(dp), allocatable, dimension(:) :: r
    integer :: i,j,jcell(10)
    type(container), allocatable,dimension(:,:) :: m,n
    integer ( kind = 4 ) id, root, dest
    integer ( kind = 4 ) status(MPI_STATUS_SIZE)
    integer ( kind = 4 ) ierr
    integer ( kind = 4 ) p
    integer ( kind = 4 ) tag
    
    allocate(m(2,1))
    allocate(n(2,1))
    d1 = 10
    ptr => d1
    call list_init(m(1,1)%list, DATA=transfer(ptr, list_data))

    d2 =20 
    ptr => d2
    call list_insert(m(1,1)%list, DATA=transfer(ptr, list_data))

    d3 = 30
    ptr => d3
    call list_insert(m(1,1)%list, DATA=transfer(ptr, list_data))
    
    d4 = 40
    ptr => d4
    call list_insert(m(1,1)%list, DATA=transfer(ptr, list_data))

    d4 = 50
    ptr => d4
    call list_insert(m(1,1)%list, DATA=transfer(ptr, list_data)) 
     
     
    d5 = 1
     ptr => d5
    call list_init(m(2,1)%list, DATA=transfer(ptr, list_data))

    d5 =2
     ptr => d5
    call list_insert(m(2,1)%list, DATA=transfer(ptr, list_data))

    d5 = 3
     ptr => d5
    call list_insert(m(2,1)%list, DATA=transfer(ptr, list_data))
    
    d5 = 4
     ptr => d5
    call list_insert(m(2,1)%list, DATA=transfer(ptr, list_data))  
    
    CALL RANDOM_SEED()
    CALL RANDOM_NUMBER(r)
    print*, r

    call MPI_Init ( ierr )
    call MPI_Comm_rank ( MPI_COMM_WORLD, id, ierr )
    call MPI_Comm_size ( MPI_COMM_WORLD, p, ierr )
    
    if (id == 0 .or. id == 1) then
        a = [1,2,3,4]
        ! CALL MPI_SEND(BUF, COUNT, DATATYPE, DEST, TAG, COMM, IERROR)
        tag = 2
        dest = 2
        call MPI_send(a,4,MPI_DOUBLE_PRECISION,dest,tag,MPI_COMM_WORLD, ierr)
        a = a*10
        tag = 3
        dest = 3
        call MPI_send(a,4,MPI_DOUBLE_PRECISION,dest,tag,MPI_COMM_WORLD, ierr)

    end if

    ! call MPI_Barrier(MPI_COMM_WORLD,ierr)
    
    ! if (id == 1 .or. id == 0) then
    !     n(1,1) = m(1,1)
    ! else 
    !     n(:,1) = m(:,1)
    ! end if 
    ! deallocate(m)

    ! print*, ' ' 
    ! print*, '\nProcess',id
    ! print*, ' ' 
    ! print*, 'ONE'
    ! ptr = transfer(list_get(n(1,1)%list), ptr)
    ! print *, 'Head node data:', ptr
    
    ! ! Get the next node
    ! ptr = transfer(list_get(list_next(n(1,1)%list)), ptr)
    ! print *, 'Second node data:', ptr  

    ! print*, 'Associated n(1,2)?',associated(list_next(n(1,1)%list))
    ! ptr = transfer(list_get(list_next(n(1,1)%list)), ptr)
    ! print *, ' data:', ptr

    call MPI_Finalize ( ierr )
    ! print*, 'OTHER'
    ! ptr = transfer(list_get(m(2,1)%list), ptr)
    ! print *, 'Head node data:', ptr
    
    ! ! Get the next node
    ! ptr = transfer(list_get(list_next(m(2,1)%list)), ptr)
    ! print *, 'Second node data:', ptr  

    
!     ptr => d5
!     print*, associated(m(2,1)%list)
!     do d5 = 1,4
!         if (associated(m(2,1)%list)) then
!             print*,d5
! !             ptr => d5
!             print*, transfer(ptr, list_data)
!             call list_insert(m(2,1)%list, DATA=transfer(ptr, list_data))
!             print*,d5
!         else
!             print*,d5
! !             ptr => d5
!             print*, transfer(ptr, list_data)
!             call list_init(m(2,1)%list, DATA=transfer(ptr, list_data))
!             print*,d5
!         end if
!     end do
! 
!       
!     print*, "AAA"
!     
!     
! !     d = 3
! !     ptr => d
! !     call list_init(m(2)%list, DATA=transfer(ptr, list_data))
! !   
! !     d = 4
! !     ptr => d
! !     call list_insert(m(2)%list, DATA=transfer(ptr, list_data))
!     
! 
! 
! !     ptr = transfer(list_get(m(2)%list), ptr)
! !     print *, 'Head node data:', ptr
! !     
! !     ! Get the next node
! !     ptr = transfer(list_get(list_next(m(2)%list)), ptr)
! !     print *, 'Second node data:', ptr  
! ! 
!     print *, "OUTRA"
!     
! !     ptr = transfer(list_get(m(2,1)%list), ptr)
! !     print *, 'Head node data:', ptr
!     call testes(m,ptr,ptr)
!         
!     ptr = transfer(list_get(list_next(m(2,1)%list)), ptr)
!     !print *, 'Second node data:', ptr  
!     
!     node => list_next(m(2,1)%list)
! !     node2 => node
!     ptr = transfer(list_get(node), ptr)
!     print *, 'Second node data:', ptr  
!     
!     node => list_next(node)
!     node3 => list_next(node)
!     
! !     print*,associated(node)
!     ptr = transfer(list_get(node), ptr)
!     print *, 'Third node data:', ptr  
! 
!     node => list_next(node)
!     ptr = transfer(list_get(node), ptr)
!     print *, 'Fourth node data:', ptr  
!     
!     node => list_next(node)
!     if (associated(node)) then
!     print*,associated(node)
!     ptr = transfer(list_get(node), ptr)
!     print *, 'Fifth node data:', ptr    
!     end if
! 
!     print *, 'UMA'
!     
!     ptr = transfer(list_get(m(1,1)%list), ptr)
!     print *, 'Head node data:', ptr
!         
!  !   ptr = transfer(list_get(list_next(m(1,1)%list)), ptr)
!  !   print *, 'Second node data:', ptr  
!     
!     node => list_next(m(1,1)%list)
!     
!     ptr = transfer(list_get(node), ptr)
!     print *, 'Second node data:', ptr  
!     
! 
!     node2 => list_next(node)    
!     node => list_next(node)
! 
!    ! print*,associated(node)
!     ptr = transfer(list_get(node), ptr)
!     print *, 'Third node data:', ptr  
! 
!     node => list_next(node)
!      
!     ptr = transfer(list_get(node), ptr)
!     print *, 'Fourth node data:', ptr  
!     
!     node => list_next(node)
!     print*,associated(node)
!     if (associated(node)) then
!   !  print*,associated(node)
!     ptr = transfer(list_get(node), ptr)
!     print *, 'Fifth node data:', ptr    
!     
!         node => list_next(node)
!         print*,associated(node)
!         if (associated(node)) then
!     !  print*,associated(node)
!         ptr = transfer(list_get(node), ptr)
!         print *, 'Sixth node data:', ptr  
!         
!             node => list_next(node)
!             print*,associated(node)
!             if (associated(node)) then
!         !  print*,associated(node)
!             ptr = transfer(list_get(node), ptr)
!             print *, 'Seventh node data:', ptr    
!             end if
!         
!         end if
!         
!     
!     end if
!     
! 
!     
!     print *, "mudando o terceiro elemento"
!     
!    call testchange(node2,m)
! 
!     ptr = transfer(list_get(m(2,1)%list), ptr)
!     print *, 'Head node data:', ptr
!         
!     ptr = transfer(list_get(list_next(m(2,1)%list)), ptr)
!     !print *, 'Second node data:', ptr  
!     
!     node => list_next(m(2,1)%list)
!     ptr = transfer(list_get(node), ptr)
!     print *, 'Second node data:', ptr  
!     
!     node => list_next(node)
!     
! !     print*,associated(node)
!     ptr = transfer(list_get(node), ptr)
!     print *, 'Third node data:', ptr  
! 
!     node => list_next(node)
!     ptr = transfer(list_get(node), ptr)
!     print *, 'Fourth node data:', ptr  
!     
!     node => list_next(node)
!     if (associated(node)) then
!     print*,associated(node)
!     ptr = transfer(list_get(node), ptr)
!     print *, 'Fifth node data:', ptr    
!     end if
!     
!     print *, '\n original'
!     
!         ptr = transfer(list_get(m(1,1)%list), ptr)
!     print *, 'Head node data:', ptr
!         
!     ptr = transfer(list_get(list_next(m(1,1)%list)), ptr)
!     print *, 'Second node data:', ptr  
!     
!     node => list_next(m(1,1)%list)
!     ptr = transfer(list_get(node), ptr)
!    ! print *, 'Second node data:', ptr  
!     
!     node => list_next(node)
!     
!    ! print*,associated(node)
!     ptr = transfer(list_get(node), ptr)
!     print *, 'Third node data:', ptr  
!     
!     print*,'TESTE DOIDO'
!     print*, ptr
!     b = ptr
!     print*, 'b=',b
!     node => list_next(node)
!     ptr = transfer(list_get(node), ptr)
!     print *, 'Fourth node data:', ptr  
!     
!     node => list_next(node)
!     print*,associated(node),'forth'
!     if (associated(node)) then
!   !  print*,associated(node)
!     ptr = transfer(list_get(node), ptr)
!     print *, 'Fifth node data:', ptr    
!     
!         node => list_next(node)
!         print*,associated(node)
!         if (associated(node)) then
!     !  print*,associated(node)
!         ptr = transfer(list_get(node), ptr)
!         print *, 'Sixth node data:', ptr  
!         
!             node => list_next(node)
!             print*,associated(node)
!             if (associated(node)) then
!         !  print*,associated(node)
!             ptr = transfer(list_get(node), ptr)
!             print *, 'Seventh node data:', ptr    
!             end if
!         
!         end if
!         
!     
!     end if
!     

    !ptr = transfer(list_get(node), ptr)
    !print *, 'Fifth node data:', ptr 

!     open(10,file='grupo1.csv',status="old")
!     do i = 1,4
!     read(10,*) a(i,1),a(i,2)
!     end do
    
    
    
!     call printm(a,'a')
!     j = 2
!     jcell = (/((i*j),i = 1,10)/)    
!     
!     laux = .true.
!     do while (laux)
!         print*, 'abc'
!     end do
!    PRINT *, "TESTE", (.NOT. 0)
    c = 'abcde'
    if ('abcde' == c) then
        print*, "ola"
    else
        print*, "oi"
    end if 

    
end program teste2        

