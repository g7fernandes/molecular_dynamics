! A simple generic linked list test program

module data
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

program list_test
  use linkedlist
  use data
  implicit none

  integer :: cont
  type(list_t), pointer :: ll => null(), node
  type(data_t), target :: dat_a
  type(data_t), target :: dat_b
  type(data_ptr) :: ptr


  print*, "enter"
  read(*,*)

  do cont = 1, 100000
    ! Initialize two data objects
    dat_a%x = 17.5
    dat_b%x = 3.0

    ! Initialize the list with dat_a
    allocate(ptr%p)
    ptr%p%x = 17.5 ! dat_a%x
    call list_init(ll, DATA=transfer(ptr, list_data))
    print *, 'Initializing list with data:', ptr%p
    ! deallocate(ptr%p)

    allocate(ptr%p)
    ! Insert dat_b into the list
    ptr%p%x = 3.0 ! dat_b%x
    call list_insert(ll, DATA=transfer(ptr, list_data))
    print *, 'Inserting node with data:', ptr%p
    ! deallocate(ptr%p)

    allocate(ptr%p)
    ! Insert dat_b into the list
    ptr%p%x = 40.0 ! dat_b%x
    call list_insert(ll, DATA=transfer(ptr, list_data))
    print *, 'Inserting node with data:', ptr%p
    ! deallocate(ptr%p)
    
    dat_a%x = 0
    dat_b%x = 0

    ! Get the head node
    ptr = transfer(list_get(ll), ptr)
    print *, 'Head node data:', ptr%p
    deallocate(ptr%p)
    ! Get the next node
    node => list_next(ll)
    ptr = transfer(list_get(list_next(ll)), ptr)
    print *, 'Second node data:', ptr%p
    deallocate(ptr%p)
    ! Get the next node
    node => list_next(node)
    ptr = transfer(list_get(node), ptr)
    print *, 'Third node data:', ptr%p
    deallocate(ptr%p)
    ! Free the list

    ! print*, "blablabla"
    !     ! Get the head node
    ! ptr = transfer(list_get(ll), ptr)
    ! print *, 'Head node data:', ptr%p
    ! deallocate(ptr%p)
    ! ! Get the next node
    ! node => list_next(ll)
    ! ptr = transfer(list_get(list_next(ll)), ptr)
    ! print *, 'Second node data:', ptr%p
    ! deallocate(ptr%p)
    ! ! Get the next node
    ! node => list_next(node)
    ! ptr = transfer(list_get(node), ptr)
    ! print *, 'Third node data:', ptr%p
    ! deallocate(ptr%p)
    ! ! Free the list

    call list_free(ll)
  end do 
  print*, "enter"
  read(*,*)

!   print*, "enter"
!   read(*,*)
!   do cont = 1, 100000
!     ! Initialize two data objects
!     dat_a%x = 17.5
!     dat_b%x = 3.0

!     ! Initialize the list with dat_a
!     ptr%p => dat_a
!     call list_init(ll, DATA=transfer(ptr, list_data))
!     print *, 'Initializing list with data:', ptr%p

!     ! Insert dat_b into the list
!     ptr%p => dat_b
!     call list_insert(ll, DATA=transfer(ptr, list_data))
!     print *, 'Inserting node with data:', ptr%p

!     ! Get the head node
!     ptr = transfer(list_get(ll), ptr)
!     print *, 'Head node data:', ptr%p

!     ! Get the next node
!     ptr = transfer(list_get(list_next(ll)), ptr)
!     print *, 'Second node data:', ptr%p

!     ! Free the list
!     call list_free(ll)
!     end do 



end program list_test