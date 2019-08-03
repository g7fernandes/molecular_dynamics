module data
    use mod1
    use linkedlist  

    private
    public :: data_t
    public :: data_ptr

    ! Data is stored in data_t
    type :: data_t
        real(dp),dimension(2) :: x !posição
        real(dp),dimension(2) :: v !velocidade
        real(dp),dimension(2) :: F !força nela
        integer               :: grupo !grupo que a partícula pertence
        integer               :: n !numero da partícula
        ! integer, dimension(2) :: mic !minimal image convenson (ajuda a contar quantas vezes a partcula atravessou a borda numa direção)
        logical               :: flag ! bandeira auxiliar
    end type data_t

    ! A trick to allow us to store pointers in the list
    type :: data_ptr
        type(data_t), pointer :: p
    end type data_ptr


    
    ! pra fazer vetor de ponteiro
!     type :: container
!         type(list_t), pointer :: list => null()
!     end type container
  
end module data
