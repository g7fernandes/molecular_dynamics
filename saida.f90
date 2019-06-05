module saida
    contains
    subroutine vec2csv(v,n,d,prop,step,t,nimpre,start)
        use mod1
        !N é o número de partículas e d é a dimensão da propriedade
        ! cada linha é uma partícula 
        implicit none
        integer, intent(in) :: n,d
        real(dp), intent(in) :: v(n,d), t, start
        integer :: i,step, nimpre
        character(*) :: prop
        character(4) :: extensao = '.csv', passo
        real(dp) :: time
        real(dp), save :: timep, etc, dtimepp

        write(passo,'(i0)') step
        open(10,file='temp/'//prop//extensao//'.'//trim(passo),status="replace")
        if (d == 1) then
            do i = 1,n
                write(10,*) v(i,1)
            end do
        else if (d == 2) then
            do i = 1,n
                write(10,*) v(i,1),',',v(i,2)
            end do
        else
            do i = 1,n
                write(10,*) v(i,1),',',v(i,2),',',v(i,3)
            end do
        end if
        close(10)

        if (prop == "position") then 
            if (step == 0) timep = 0
            call cpu_time(time)    
            if (step == 1) then 
                etc = ((time - start)/step+ (time-timep))*0.5*nimpre - (time - start)
                dtimepp = (time-timep)
                print '("Salvo arquivo ", A, "  t = ", f10.3, "  ETC: ", f10.3, "s" )',prop//extensao//'.'//trim(passo),t,etc                                    
            else if (step > 1) then
                etc = ((time - start)*2/step + ((time-timep)*4 + dtimepp*4))*(nimpre-step)/10 
                print '("Salvo arquivo ", A, "  t = ", f10.3, "  ETC: ", f10.3, "s" )',prop//extensao//'.'//trim(passo),t,etc                                    
                dtimepp = (time-timep)
            else 
                print '("Salvo arquivo ", A, "  t = ", f10.3, "  ETC: ", "unknown" )',prop//extensao//'.'//trim(passo),t
            end if

            timep = time 
        else 
            if (step > 0) then 
                print '("Salvo arquivo ", A, "  t = ", f10.3, "  ETC: ", f10.3, "s" )',prop//extensao//'.'//trim(passo),t,etc                                    
            else 
                print '("Salvo arquivo ", A, "  t = ", f10.3, "  ETC: ", "unknown" )',prop//extensao//'.'//trim(passo),t
            end if
            timep = time
        end if     
        ! print*, 'Salvo arquivo ',prop//extensao//'.'//trim(passo), "t =", t, "ETC: ", etc 
        
    end subroutine vec2csv 
    
    subroutine linked2vec(malha,domx,domy,nxv,aux1)
       
        use linkedlist
        use mod1
        use data
        use mod0
        
        type(container), allocatable,dimension(:,:),intent(in) :: malha
        integer :: i,j,aux1
        integer, intent(in) :: domx(2), domy(2)
        type(list_t), pointer :: node
        real(dp), intent(out) :: nxv(:)
        type(data_ptr) :: ptr

        aux1 = 1
        ! real(dp),intent(inout) :: celula(:,:)
        
        do i = domy(1), domy(2)
            do j = domx(1), domx(2)
               ! print *, 'posição', i, ',', j
                node => list_next(malha(i,j)%list)
                do while (associated(node))
                    ptr = transfer(list_get(node), ptr)
                    nxv(aux1:aux1+4) = [ real(ptr%p%n, kind(0.d0)), ptr%p%x(1),ptr%p%x(2), &
                    ptr%p%v(1), ptr%p%v(2)]
                    aux1 = aux1 + 5
                    node => list_next(node)
                end do
            end do
        end do
        aux1 = aux1 -1
        
    end subroutine linked2vec    
    
    subroutine nancheck(malha,mesh,t) 
        use linkedlist
        use mod1
        use data
        use mod0

        real(dp), intent(in) :: t
        integer :: i 
        integer, intent(in) :: mesh(:)
        type(container), allocatable,dimension(:,:),intent(in) :: malha
        type(data_ptr) :: ptr
        type(list_t), pointer :: node
        do i = 1,mesh(2)+2 
            do j = 1,mesh(1)+2
                node => list_next(malha(i,j)%list)
                do while (associated(node))
                    ptr = transfer(list_get(node), ptr)
                    !calcula a energia cinética atual
                    node => list_next(node)
                    if (isnan(ptr%p%x(1))) then
                        print*, 'NaN em part',ptr%p%n,'x, t=',t
                    end if
                    if (isnan(ptr%p%x(2))) then
                        print*, 'NaN em part',ptr%p%n,'y, t=',t
                    end if
                end do
            end do
        end do        
    end subroutine nancheck
end module saida 