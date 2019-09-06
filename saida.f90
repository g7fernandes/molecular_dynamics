module saida
    contains
    subroutine vec2csv(v,n,d,prop,step,t,nimpre,start,concat0)
        use mod1
        !N é o número de partículas e d é a dimensão da propriedade
        ! cada linha é uma partícula 
        implicit none
        logical, optional :: concat0
        logical :: laux
        integer, intent(in) :: n,d
        real(dp), intent(in) :: v(n,d), t, start
        integer :: i,step, nimpre, ic1, cpu_countrate, horas, min 
        real(dp) :: sec
        character(*) :: prop
        character(4) :: extensao = '.csv', passo
        real(dp) :: time
        real(dp), save :: timep, etc, dtimepp
        character(LEN=*),parameter :: fmt5 = '(f32.16, ", ",f32.16, ", ",f32.16, ", ",f32.16, ", ",f32.16 )'
        character(LEN=*),parameter :: fmt6 = '(f5.0, ", ",f32.16, ", ",f32.16, ", ",f32.16, ", ",f32.16, ", ",f32.16 )'
        character(LEN=*),parameter :: fmt7 = '(f5.0, ", ",f5.0, ", ",f32.16, ", ",f32.16, ", ",f32.16, ", ",f32.16, ", ",f32.16 )'

        laux = .false.
        if (present(concat0)) laux = concat0
        write(passo,'(i0)') step
        if (d <= 3) then
            open(10,file='temp/'//prop//extensao//'.'//trim(passo),status="replace")
        else
            open(10,file='temp2/'//prop//extensao//'.'//trim(passo),status="replace")
        end if

        if (.not. laux) then
            if (d == 1) then
                do i = 1,n
                    write(10,*) v(i,1)
                end do
            else if (d == 2) then
                do i = 1,n
                    write(10,*) v(i,1),',',v(i,2)
                end do
            else if (d == 3) then
                do i = 1,n
                    write(10,*) v(i,1),',',v(i,2),',',v(i,3)
                end do
            else if (d == 4) then
                do i = 1,n
                    write(10,*) v(i,1),',',v(i,2),',',v(i,3),',',v(i,4)
                end do
            else if (d == 5) then 
                do i = 1,n
                    write(10,fmt5) v(i,1),v(i,2),v(i,3),v(i,4),v(i,5)
                end do          
            else if (d == 6) then 
                do i = 1,n
                    write(10,fmt6) v(i,1),v(i,2),v(i,3),v(i,4),v(i,5),v(i,6)
                end do          
            else if (d == 7) then 
                do i = 1,n
                    write(10,fmt7) v(i,1),v(i,2),v(i,3),v(i,4),v(i,5),v(i,6),v(i,7)
                end do          
            end if
        else 
            if (d == 1) then
                do i = 1,n
                    write(10,*) v(i,1),',',0
                end do
            else if (d == 2) then
                do i = 1,n
                    write(10,*) v(i,1),',',v(i,2),',',0
                end do
            else if (d == 3) then
                do i = 1,n
                    write(10,*) v(i,1),',',v(i,2),',',v(i,3),',',0
                end do
            else if (d == 4) then 
                do i = 1,n
                    write(10,fmt5) v(i,1),v(i,2),v(i,3),v(i,4),0
                end do          
            else if (d == 5) then 
                do i = 1,n
                    write(10,fmt6) v(i,1),v(i,2),v(i,3),v(i,4),v(i,5),0
                end do          
            else if (d == 6) then 
                do i = 1,n
                    write(10,fmt7) v(i,1),v(i,2),v(i,3),v(i,4),v(i,5),v(i,6),0
                end do          
            end if
        end if 
        
        close(10)

        if (prop == "position") then 
            if (step == 0) timep = 0
            ! call cpu_time(time)    
                call system_clock(ic1,cpu_countrate)
                time = real(ic1,kind(0.d0))/real(cpu_countrate,kind(0.d0))
            if (step == 2) then 
                call system_clock(ic1,cpu_countrate)
                time = real(ic1,kind(0.d0))/real(cpu_countrate,kind(0.d0))
                etc = ((time - start)/real(step,kind(0.d0)) + (time-timep))*0.5*real(nimpre,kind(0.d0)) - (time - start)
                dtimepp = (time-timep)
                print '("Salvo arquivo ", A, "  t = ", f10.3, "  ETC: ", f10.3, "s" )',prop//extensao//'.'//trim(passo),t,etc                                    
            else if (step > 2) then
                call system_clock(ic1,cpu_countrate)
                time = real(ic1,kind(0.d0))/real(cpu_countrate,kind(0.d0))
                etc = ((time - start)*6/real(step,kind(0.d0)) + ((time-timep)*2 + dtimepp*2))* & 
                    (real(nimpre,kind(0.d0))-real(step,kind(0.d0)))/10 
                horas = int(etc/3600)
                min = (int(etc) - horas*3600)/60 
                sec = etc - real(horas*3600 + min*60,kind(0.d0))
                print '("Salvo arquivo ", A, "  t = ", f10.3, "  ETC: ", i3,":",i2,":",f4.1 )' &
                    ,prop//extensao//'.'//trim(passo),t,horas, min, sec                                    
                dtimepp = (time-timep)  
            else 
                print '("Salvo arquivo ", A, "  t = ", f10.3, "  ETC: ", "unknown" )',prop//extensao//'.'//trim(passo),t
            end if

            timep = time 
        else 
            if (step > 2) then 
                print '("Salvo arquivo ", A, "  t = ", f10.3, "  ETC: ", f10.3, "s" )',prop//extensao//'.'//trim(passo),t,etc                                    
            else 
                print '("Salvo arquivo ", A, "  t = ", f10.3, "  ETC: ", "unknown" )',prop//extensao//'.'//trim(passo),t
            end if
           ! time = real(ic1,kind(0.d0))/real(cpu_countrate,kind(0.d0))
           ! timep = time
        end if     
        ! print*, 'Salvo arquivo ',prop//extensao//'.'//trim(passo), "t =", t, "ETC: ", etc 
        
    end subroutine vec2csv 
    
    subroutine linked2vec(malha,mesh,domx,domy,nxv,aux1)
       
        use linkedlist
        use mod1
        use data
        use mod0
        
        type(container), allocatable,dimension(:,:),intent(in) :: malha
        integer :: i,j,aux1
        integer, intent(in) :: domx(2), domy(2), mesh(2)
        type(list_t), pointer :: node
        real(dp), intent(out) :: nxv(:)
        type(data_ptr) :: ptr

        aux1 = 1
        ! real(dp),intent(inout) :: celula(:,:)
        
        do i = domy(1), domy(2)
            do j = domx(1), domx(2)
                if(i > 1 .and. i < mesh(2)+2 .and. j > 1 .and. j < mesh(1)+2) then                  
                    node => list_next(malha(i,j)%list)
                    do while (associated(node))
                        ptr = transfer(list_get(node), ptr)
                        nxv(aux1:aux1+4) = [ real(ptr%p%n, kind(0.d0)), ptr%p%x(1),ptr%p%x(2), &
                        ptr%p%v(1), ptr%p%v(2)]
                        aux1 = aux1 + 5
                        node => list_next(node)
                    end do
                end if
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