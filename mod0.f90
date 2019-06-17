module mod0
    use mod1
    use linkedlist
        type :: container
            type(list_t), pointer :: list => null()
        end type container

        ! Propriedades de cada grupo de partícula (indicado pelo índice do vetor com este tipo)
        ! nos permite mais fases e economiza espaço
        type :: prop_grupo
            real(dp) :: m ! mass
            real(dp) :: epsilon 
            real(dp) :: sigma 
            real(dp) :: rs !raio sólido 
            real(dp)  :: x_lockdelay    ! só vai poder mudar de posição a partir de t = x_lockdelay
            real(dp) :: fric_term
        end type prop_grupo

            ! LSTR  lista de transferência pro MPI
        type :: lstr
            real(dp),allocatable :: lstrdb_N(:)
            real(dp),allocatable :: lstrdb_S(:)
            real(dp),allocatable :: lstrdb_E(:)
            real(dp),allocatable :: lstrdb_W(:)
            real(dp),allocatable :: lstrdb_D(:) ! Diagonal. Para transferir nas direções NE, NS, SE, SW no encontro de 4 processos

            integer,allocatable :: lstrint_N(:)
            integer,allocatable :: lstrint_S(:)
            integer,allocatable :: lstrint_E(:)
            integer,allocatable :: lstrint_W(:)
            integer,allocatable :: lstrint_D(:)
        end type lstr
end module mod0