module fisica
    contains
    !Calcula velocidades distribuidas de acordo com MaxwellBolzmann
    ! factor = sqrt(<vd²>) = sqrt(kb*T/m) por componente
    subroutine MaxwellBoltzmann(malha,mesh,factor)
        use linkedlist
        use mod1
        use randnormal
        use mod0
        use data
        
        real(dp),intent(in) :: factor(2)
        integer :: i,j
        type(container), allocatable,dimension(:,:),intent(in) :: malha
        integer, intent(in) :: mesh(2)
        type(list_t), pointer :: node
        type(data_ptr) :: ptr

        do i = 2,mesh(2)+1
            do j = 2,mesh(1)+1
                node => list_next(malha(i,j)%list)
                do while (associated(node))
                    ptr = transfer(list_get(node), ptr)
                    ptr%p%v(1) = ptr%p%v(1) + factor(1)*GaussDeviate()
                    ptr%p%v(2) = ptr%p%v(2) + factor(2)*GaussDeviate()
                    ptr%p%F = [0,0]
                    node => list_next(node)
                end do
            end do
        end do

    end subroutine MaxwellBoltzmann
    
    !computa o potencial total 
    subroutine comp_pot(mesh,malha,propriedade,r_cut,domx,domy, ids, id, t, nRfu) 
        !Computa e imprime o potencial e forças em cada partícula, intermoleculares

        use linkedlist
        use mod1
        use data
        use mod0
        use mpi
        
        type(prop_grupo), allocatable,dimension(:),intent(in) :: propriedade
        type(container), allocatable,dimension(:,:),intent(in) :: malha
        real(dp), intent(in) :: t,r_cut
        real(dp) :: sigma, epsil, sigma_a, epsil_a,sigma_b, epsil_b, rcut, fric_term !fric_term = força de ficção
        real(dp) :: x1(2),v1(2),x2(2),p1(2), rs1, rs2, coss, sine, u   
        real(dp), intent(inout), allocatable, dimension(:) :: nRfu
        integer :: i,j,ct = 0, n1, n2 !,ptr, ptrn
        integer, intent(in) :: mesh(:),domx(2),domy(2), id
        real(dp) :: Fi(2)=0,r, aux2(2),fR(2), fric_term1, fric_term2, dox(2), doy(2)
        type(list_t), pointer :: node, next_node
        type(data_ptr) :: ptr,ptrn
        integer ( kind = 4 ), intent(in) :: ids(8)
        
 
        !Lennard Jones
        dox = domx 
        doy = domy 
        if (sum(ids(1:4)) > -4) then
            !caso paralelo
            if (domx(1) > 1) dox(1) = domx(1) - 1
            if (domy(1) > 1) doy(1) = domy(1) - 1
            if (domx(2) < mesh(1)+2) dox(2) = domx(2) + 1
            if (domy(2) < mesh(2)+2) doy(2) = domy(2) + 1
        else
            dox = domx 
            doy = domy 
        end if 
        
        do i = doy(1),doy(2) ! i é linha
            do j = dox(1),dox(2)
                node => list_next(malha(i,j)%list)
                if (associated(node)) then
                    ptr = transfer(list_get(node), ptr)
                end if
                do while (associated(node))
                    ! print*, "Encontrada", ct, "celulas, id=", id
                    ptr = transfer(list_get(node), ptr) !particula selecionada
                    ptr%p%flag = .true. ! indica ao comp_x que a partícula precisa ser calculada
                    x1 = ptr%p%x
                    v1 = ptr%p%v
                    n1 = ptr%p%n
                    m1 = propriedade(ptr%p%grupo)%m
                    p1 = [v1(1), v1(2)]*m1 
                    rs1 = propriedade(ptr%p%grupo)%rs !raio sólido 
                    fric_term1 = propriedade(ptr%p%grupo)%fric_term
                    sigma_a = propriedade(ptr%p%grupo)%sigma
                    ! rcut = r_cut*sigma
                    epsil_a = propriedade(ptr%p%grupo)%epsilon 

                    ! print '("x1  =", f10.6, " ", f10.6, " n ", i2, " cell ",i2," ",i2)',x1(1),x1(2), ptr%p%n,i,j
                    ! if (id == 0) read(*,*)
                    !calcular a força desta com todas as outras partículas
                    next_node => list_next(node) ! próxima partícula da célula
                    node => list_next(node) ! a ser computado com a particula selecionada

                    if (i /= doy(1) .and. i /= doy(2) .and. j /= dox(1) .and. j /= dox(2)) then
                        nRfu(n1*6-5) = real(ptr%p%n,kind(0.d0))
                    end if
                    nRfu(n1*6-1:n1*6) = [p1(1), p1(2)] 
                    do while (associated(node))
                        ! NA PRÓPRIA CELULA 
                        ptrn = transfer(list_get(node), ptrn) !outra particula selecionada
                        ! print*, "L 108"
                        sigma_b = propriedade(ptrn%p%grupo)%sigma
                        epsil_b = propriedade(ptrn%p%grupo)%epsilon 
                        rs2 = propriedade(ptrn%p%grupo)%rs 
                        ! Lorenz-Betherlot rule for mixing epsilon sigma 
                        if (sigma_a > sigma_b) then
                            rcut = r_cut*sigma_a + rs1 + rs2
                            sigma = 0.5*(sigma_a + sigma_b)
                            epsil = sqrt(epsil_a *epsil_b )
                        else if (sigma_a < sigma_b) then 
                            rcut = r_cut*sigma_b + rs1 + rs2
                            sigma = 0.5*(sigma_a + sigma_b)
                            epsil = sqrt(epsil_a *epsil_b )
                        else 
                            rcut = r_cut*sigma_a + rs1 + rs2
                            sigma = sigma_a
                            epsil = epsil_a
                        end if 
                        x2 = ptrn%p%x
                        n2 = ptrn%p%n
                        r = sqrt((x1(1)-x2(1))**2 + (x1(2)-x2(2))**2)
                        coss = (x1(1)-x2(1))/r 
                        sine = (x1(2)-x2(2))/r 
                        r = r - rs1 - rs2 !raio
                        ! print*, "L 389 r", r, "id",id
                        ! print*, "L 129"
                        if (r <= rcut) then
                            aux2 = -(1/r**2)*(sigma/r)**6*(1-2*(sigma/r)**6)*24*epsil* & 
                            [(x1(1)-x2(1)) - (rs1+rs2)*coss, (x1(2)-x2(2)) - (rs1+rs2)*sine]  
                            u = 4*epsil*((sigma/r)**12 - (sigma/r)**6)  
                            nRfu(n1*6-4:n1*6-2) = [aux2(2)*r*coss, aux2(1)*r*sine, u] &
                            + nRfu(n1*6-4:n1*6-2)
                            nRfu(n2*6-4:n2*6-2) = [aux2(2)*r*coss, aux2(1)*r*sine, u] & 
                            + nRfu(n2*6-4:n2*6-2)
                            if (i /= doy(1) .and. i /= doy(2) .and. j /= dox(1) .and. j /= dox(2)) then
                                nRfu(n2*6-5) = real(n2,kind(0.d0))
                            end if
                            ! print*, "nrfu2", nRfu(n1*6-5:n1*6)
                            !read(*,*)
                        end if
                        node => list_next(node) ! próxima partícula da célula
                    end do

                    !Células ao redor  !i é linha, j é coluna
                    if (i /= mesh(2)+2) then ! se não for a última linha
                        if (j == mesh(1)+2) then ! se for a última coluna
                            node => list_next(malha(i+1,j)%list) !interagirá com a próxima linha apenas
                            do while (associated(node))
                                ptrn = transfer(list_get(node), ptrn) !outra particula selecionada
                                sigma_b = propriedade(ptrn%p%grupo)%sigma
                                epsil_b = propriedade(ptrn%p%grupo)%epsilon 
                                rs2 = propriedade(ptrn%p%grupo)%rs 
                                fric_term2 = propriedade(ptrn%p%grupo)%fric_term
                                ! Lorenz-Betherlot rule for mixing epsilon sigma 
                                if (sigma_a > sigma_b) then
                                    rcut = r_cut*sigma_a + (rs1 + rs2)
                                    sigma = 0.5*(sigma_a + sigma_b)
                                    epsil = sqrt(epsil_a*epsil_b)
                                else if (sigma_a < sigma_b) then 
                                    rcut = r_cut*sigma_b + (rs1 + rs2)
                                    sigma = 0.5*(sigma_a + sigma_b)
                                    epsil = sqrt(epsil_a*epsil_b)
                                else 
                                    rcut = r_cut*sigma_a + (rs1 + rs2)
                                    sigma = sigma_a
                                    epsil = epsil_a
                                end if 
                                x2 = ptrn%p%x
                                n2 = ptrn%p%n
                                r = sqrt((x1(1)-x2(1))**2 + (x1(2)-x2(2))**2)
                                coss = (x1(1)-x2(1))/r 
                                sine = (x1(2)-x2(2))/r 
                                r = r - rs1 - rs2 !raio
                                ! print*, "L 666 r", r, "id",id
                                if (r <= rcut) then
                                    aux2 = -(1/r**2)*(sigma/r)**6*(1-2*(sigma/r)**6)*24*epsil* & 
                                        [(x1(1)-x2(1)) - (rs1+rs2)*coss, (x1(2)-x2(2)) - (rs1+rs2)*sine]
                                    nRfu(n1*6-4:n1*6-2) = [aux2(2)*r*coss, aux2(1)*r*sine, u] &
                                        + nRfu(n1*6-4:n1*6-2)
                                    nRfu(n2*6-4:n2*6-2) = [aux2(2)*r*coss, aux2(1)*r*sine, u] & 
                                        + nRfu(n2*6-4:n2*6-2)
                                    if (i /= doy(1) .and. i /= doy(2) .and. j /= dox(1) .and. j /= dox(2)) then
                                        nRfu(n2*6-5) = real(n2,kind(0.d0))
                                    end if
                                end if
                                node => list_next(node) ! próxima partícula da célula
                            end do
                            
                            if (j /= 1) then !se não for a primeira coluna 
                                node => list_next(malha(i+1,j-1)%list) !interagirá com a próxima linha apenas
                                do while (associated(node))
                                    ptrn = transfer(list_get(node), ptrn) !outra particula selecionada
                                    sigma_b = propriedade(ptrn%p%grupo)%sigma
                                    epsil_b = propriedade(ptrn%p%grupo)%epsilon 
                                    rs2 = propriedade(ptrn%p%grupo)%rs 
                                    fric_term2 = propriedade(ptrn%p%grupo)%fric_term
                                    ! Lorenz-Betherlot rule for mixing epsilon sigma 
                                    if (sigma_a > sigma_b) then
                                        rcut = r_cut*sigma_a + (rs1 + rs2)
                                        sigma = 0.5*(sigma_a + sigma_b)
                                        epsil = sqrt(epsil_a*epsil_b)
                                    else if (sigma_a < sigma_b) then 
                                        rcut = r_cut*sigma_b + (rs1 + rs2)
                                        sigma = 0.5*(sigma_a + sigma_b)
                                        epsil = sqrt(epsil_a*epsil_b)
                                    else 
                                        rcut = r_cut*sigma_a + (rs1 + rs2)
                                        sigma = sigma_a
                                        epsil = epsil_a
                                    end if 

                                    x2 = ptrn%p%x
                                    n2 = ptrn%p%n
                                    r = sqrt((x1(1)-x2(1))**2 + (x1(2)-x2(2))**2)
                                    coss = (x1(1)-x2(1))/r 
                                    sine = (x1(2)-x2(2))/r 
                                    r = r - rs1 - rs2 !raio
                                    ! print*, "L 710 r", r, "id",id
                                    if (r <= rcut) then
                                        aux2 = -(1/r**2)*(sigma/r)**6*(1-2*(sigma/r)**6)*24*epsil* & 
                                            [(x1(1)-x2(1)) - (rs1+rs2)*coss, (x1(2)-x2(2)) - (rs1+rs2)*sine]

                                        nRfu(n1*6-4:n1*6-2) = [aux2(2)*r*coss, aux2(1)*r*sine, u] &
                                        + nRfu(n1*6-4:n1*6-2)
                                        nRfu(n2*6-4:n2*6-2) = [aux2(2)*r*coss, aux2(1)*r*sine, u] & 
                                        + nRfu(n2*6-4:n2*6-2)
                                        if (i /= doy(1) .and. i /= doy(2) .and. j /= dox(1) .and. j /= dox(2)) then
                                            nRfu(n2*6-5) = real(n2,kind(0.d0))
                                        end if
                                    end if
                                    node => list_next(node) ! próxima partícula da célula
                                end do                            
                            end if
                        else
                             !interagirá com a próxima linha e coluna, e na diagonal
                            node => list_next(malha(i,j+1)%list) 
                            do while (associated(node))
                                ptrn = transfer(list_get(node), ptrn) !outra particula selecionada
                                sigma_b = propriedade(ptrn%p%grupo)%sigma
                                epsil_b = propriedade(ptrn%p%grupo)%epsilon 
                                rs2 = propriedade(ptrn%p%grupo)%rs
                                fric_term2 = propriedade(ptrn%p%grupo)%fric_term 
                                ! Lorenz-Betherlot rule for mixing epsilon sigma 
                                if (sigma_a > sigma_b) then
                                    rcut = r_cut*sigma_a + (rs1 + rs2)
                                    sigma = 0.5*(sigma_a + sigma_b)
                                    epsil = sqrt(epsil_a*epsil_b)
                                else if (sigma_a < sigma_b) then 
                                    rcut = r_cut*sigma_b + (rs1 + rs2)
                                    sigma = 0.5*(sigma_a + sigma_b)
                                    epsil = sqrt(epsil_a*epsil_b)
                                else 
                                    rcut = r_cut*sigma_a + (rs1 + rs2)
                                    sigma = sigma_a
                                    epsil = epsil_a
                                end if 

                                x2 = ptrn%p%x
                                n2 = ptrn%p%n
                                r = sqrt((x1(1)-x2(1))**2 + (x1(2)-x2(2))**2)
                                coss = (x1(1)-x2(1))/r 
                                sine = (x1(2)-x2(2))/r 
                                r = r - rs1 - rs2 !raio
                                ! print*, "L 757 r", r, "id",id
                                if (r <= rcut) then
                                    aux2 = -(1/r**2)*(sigma/r)**6*(1-2*(sigma/r)**6)*24*epsil* & 
                                        [(x1(1)-x2(1)) - (rs1+rs2)*coss, (x1(2)-x2(2)) - (rs1+rs2)*sine] 

                                    nRfu(n1*6-4:n1*6-2) = [aux2(2)*r*coss, aux2(1)*r*sine, u] &
                                    + nRfu(n1*6-4:n1*6-2)
                                    nRfu(n2*6-4:n2*6-2) = [aux2(2)*r*coss, aux2(1)*r*sine, u] & 
                                    + nRfu(n2*6-4:n2*6-2)
                                    if (i /= doy(1) .and. i /= doy(2) .and. j /= dox(1) .and. j /= dox(2)) then
                                        nRfu(n2*6-5) = real(n2,kind(0.d0))
                                    end if
                                end if
                                node => list_next(node) ! próxima partícula da célula                                
                            end do

                            node => list_next(malha(i+1,j)%list) !interagirá com a próxima linha 
                            do while (associated(node))
                                ptrn = transfer(list_get(node), ptrn) !outra particula selecionada
                                sigma_b = propriedade(ptrn%p%grupo)%sigma
                                epsil_b = propriedade(ptrn%p%grupo)%epsilon 
                                rs2 = propriedade(ptrn%p%grupo)%rs 
                                fric_term2 = propriedade(ptrn%p%grupo)%fric_term
                                ! Lorenz-Betherlot rule for mixing epsilon sigma 
                                if (sigma_a > sigma_b) then
                                    rcut = r_cut*sigma_a + (rs1 + rs2)
                                    sigma = 0.5*(sigma_a + sigma_b)
                                    epsil = sqrt(epsil_a*epsil_b)
                                else if (sigma_a < sigma_b) then 
                                    rcut = r_cut*sigma_b + (rs1 + rs2)
                                    sigma = 0.5*(sigma_a + sigma_b)
                                    epsil = sqrt(epsil_a*epsil_b)
                                else 
                                    rcut = r_cut*sigma_a + (rs1 + rs2)
                                    sigma = sigma_a
                                    epsil = epsil_a
                                end if 
                                x2 = ptrn%p%x
                                n2 = ptrn%p%n
                                r = sqrt((x1(1)-x2(1))**2 + (x1(2)-x2(2))**2)
                                coss = (x1(1)-x2(1))/r 
                                sine = (x1(2)-x2(2))/r 
                                r = r - rs1 - rs2 !raio
                                ! print*, "L 798 r", r, "id",id
                                if (r <= rcut) then
                                    aux2 = -(1/r**2)*(sigma/r)**6*(1-2*(sigma/r)**6)*24*epsil* & 
                                        [(x1(1)-x2(1)) - (rs1+rs2)*coss, (x1(2)-x2(2)) - (rs1+rs2)*sine] 
                                    nRfu(n1*6-4:n1*6-2) = [aux2(2)*r*coss, aux2(1)*r*sine, u] &
                                        + nRfu(n1*6-4:n1*6-2)
                                    nRfu(n2*6-4:n2*6-2) = [aux2(2)*r*coss, aux2(1)*r*sine, u] & 
                                        + nRfu(n2*6-4:n2*6-2)
                                    if (i /= doy(1) .and. i /= doy(2) .and. j /= dox(1) .and. j /= dox(2)) then
                                        nRfu(n2*6-5) = real(n2,kind(0.d0))
                                    end if
                                end if
                                node => list_next(node) ! próxima partícula da célula                                
                            end do
                            
                            node => list_next(malha(i+1,j+1)%list) !interagirá com a próxima linha e coluna
                            do while (associated(node))
                                ptrn = transfer(list_get(node), ptrn) !outra particula selecionada
                                sigma_b = propriedade(ptrn%p%grupo)%sigma
                                epsil_b = propriedade(ptrn%p%grupo)%epsilon 
                                rs2 = propriedade(ptrn%p%grupo)%rs 
                                fric_term2 = propriedade(ptrn%p%grupo)%fric_term
                                ! Lorenz-Betherlot rule for mixing epsilon sigma 
                                if (sigma_a > sigma_b) then
                                    rcut = r_cut*sigma_a + (rs1 + rs2)
                                    sigma = 0.5*(sigma_a + sigma_b)
                                    epsil = sqrt(epsil_a*epsil_b)
                                else if (sigma_a < sigma_b) then 
                                    rcut = r_cut*sigma_b + (rs1 + rs2)
                                    sigma = 0.5*(sigma_a + sigma_b)
                                    epsil = sqrt(epsil_a*epsil_b)
                                else 
                                    rcut = r_cut*sigma_a + (rs1 + rs2)
                                    sigma = sigma_a
                                    epsil = epsil_a
                                end if 

                                x2 = ptrn%p%x
                                n2 = ptrn%p%n
                                r = sqrt((x1(1)-x2(1))**2 + (x1(2)-x2(2))**2)
                                coss = (x1(1)-x2(1))/r 
                                sine = (x1(2)-x2(2))/r 
                                r = r - rs1 - rs2 !raio
                                ! print*, "L 841 r", r, "id",id
                                if (r <= rcut) then
                                    aux2 = -(1/r**2)*(sigma/r)**6*(1-2*(sigma/r)**6)*24*epsil* & 
                                        [(x1(1)-x2(1)) - (rs1+rs2)*coss, (x1(2)-x2(2)) - (rs1+rs2)*sine] 
                                    nRfu(n1*6-4:n1*6-2) = [real(ptr%p%n,kind(0.d0)), aux2(2)*r*coss, aux2(1)*r*sine, u] &
                                        + nRfu(n1*6-4:n1*6-2)
                                    nRfu(n2*6-4:n2*6-2) = [real(ptr%p%n,kind(0.d0)), aux2(2)*r*coss, aux2(1)*r*sine, u] & 
                                        + nRfu(n2*6-4:n2*6-2)
                                    if (i /= doy(1) .and. i /= doy(2) .and. j /= dox(1) .and. j /= dox(2)) then
                                        nRfu(n2*6-5) = real(n2,kind(0.d0))
                                    end if
                                end if
                                node => list_next(node) ! próxima partícula da célula                                
                            end do
                            
                            if (j /= 1) then !se não for a primeira coluna 
                                node => list_next(malha(i+1,j-1)%list) !interagirá com a próxima linha apenas
                                do while (associated(node))
                                    ptrn = transfer(list_get(node), ptrn) !outra particula selecionada
                                    sigma_b = propriedade(ptrn%p%grupo)%sigma
                                    epsil_b = propriedade(ptrn%p%grupo)%epsilon 
                                    rs2 = propriedade(ptrn%p%grupo)%rs 
                                    fric_term2 = propriedade(ptrn%p%grupo)%fric_term
                                    ! Lorenz-Betherlot rule for mixing epsilon sigma 
                                    if (sigma_a > sigma_b) then
                                        rcut = r_cut*sigma_a + (rs1 + rs2)
                                        sigma = 0.5*(sigma_a + sigma_b)
                                        epsil = sqrt(epsil_a*epsil_b)
                                    else if (sigma_a < sigma_b) then 
                                        rcut = r_cut*sigma_b + (rs1 + rs2)
                                        sigma = 0.5*(sigma_a + sigma_b)
                                        epsil = sqrt(epsil_a*epsil_b)
                                    else 
                                        rcut = r_cut*sigma_a + (rs1 + rs2)
                                        sigma = sigma_a
                                        epsil = epsil_a
                                    end if 

                                    x2 = ptrn%p%x
                                    n2 = ptrn%p%n
                                    r = sqrt((x1(1)-x2(1))**2 + (x1(2)-x2(2))**2)
                                    coss = (x1(1)-x2(1))/r 
                                    sine = (x1(2)-x2(2))/r 
                                    r = r - rs1 - rs2 !raio
                                    ! print*, "L 885 r", r, "id",id
                                    if (r <= rcut) then
                                        aux2 = -(1/r**2)*(sigma/r)**6*(1-2*(sigma/r)**6)*24*epsil* & 
                                            [(x1(1)-x2(1)) - (rs1+rs2)*coss, (x1(2)-x2(2)) - (rs1+rs2)*sine] 
                                        nRfu(n1*6-4:n1*6-2) = [aux2(2)*r*coss, aux2(1)*r*sine, u] &
                                            + nRfu(n1*6-4:n1*6-2)
                                        nRfu(n2*6-4:n2*6-2) = [aux2(2)*r*coss, aux2(1)*r*sine, u] & 
                                            + nRfu(n2*6-4:n2*6-2)
                                        if (i /= doy(1) .and. i /= doy(2) .and. j /= dox(1) .and. j /= dox(2)) then
                                            nRfu(n2*6-5) = real(n2,kind(0.d0))
                                        end if
                                    end if
                                    node => list_next(node) ! próxima partícula da célula
                                end do                            
                            end if
                        end if
                        
                    else ! se for a última lina, só interage com a celua ao lado 
                        node => list_next(malha(i,j+1)%list) !interagirá com a próxima linha e coluna
                        do while (associated(node))
                            ptrn = transfer(list_get(node), ptrn) !outra particula selecionada
                            sigma_b = propriedade(ptrn%p%grupo)%sigma
                            epsil_b = propriedade(ptrn%p%grupo)%epsilon 
                            rs2 = propriedade(ptrn%p%grupo)%rs 
                            fric_term2 = propriedade(ptrn%p%grupo)%fric_term
                            ! Lorenz-Betherlot rule for mixing epsilon sigma 
                            if (sigma_a > sigma_b) then
                                rcut = r_cut*sigma_a + (rs1 + rs2)
                                sigma = 0.5*(sigma_a + sigma_b)
                                epsil = sqrt(epsil_a*epsil_b)
                            else if (sigma_a < sigma_b) then 
                                rcut = r_cut*sigma_b + (rs1 + rs2)
                                sigma = 0.5*(sigma_a + sigma_b)
                                epsil = sqrt(epsil_a*epsil_b)
                            else 
                                rcut = r_cut*sigma_a + (rs1 + rs2)
                                sigma = sigma_a
                                epsil = epsil__a
                            end if 

                            x2 = ptrn%p%x
                            n2 = ptrn%p%n
                            r = sqrt((x1(1)-x2(1))**2 + (x1(2)-x2(2))**2)
                            coss = (x1(1)-x2(1))/r 
                            sine = (x1(2)-x2(2))/r 
                            r = r - rs1 - rs2 !raio
                            ! print*, "L 931 r", r, "id",id
                            if (r <= rcut) then
                                aux2 = -(1/r**2)*(sigma/r)**6*(1-2*(sigma/r)**6)*24*epsil* & 
                                    [(x1(1)-x2(1)) - (rs1+rs2)*coss, (x1(2)-x2(2)) - (rs1+rs2)*sine] 
                                nRfu(n1*6-4:n1*6-2) = [aux2(2)*r*coss, aux2(1)*r*sine, u] &
                                    + nRfu(n1*6-4:n1*6-2)
                                nRfu(n2*6-4:n2*6-2) = [-aux2(2)*r*coss, -aux2(1)*r*sine, u] & 
                                    + nRfu(n2*6-4:n2*6-2)
                                if (i /= doy(1) .and. i /= doy(2) .and. j /= dox(1) .and. j /= dox(2)) then
                                    nRfu(n2*6-5) = real(n2,kind(0.d0))
                                end if
                            end if
                            node => list_next(node) ! próxima partícula da célula                                
                        end do
                    end if
                    node => next_node
                end do
            end do
        end do
        aux1 = aux1 - 1
        ! print*, "ID, AUX1", id, aux1
        ! print*, "FIM"
    end subroutine comp_pot    
 
    function comp_Kglobal(malha,domx,domy,propriedade,np,id,t) result(K)
        use linkedlist
        use mod1
        use data
        use mod0
        use mpi

        type(prop_grupo), allocatable,dimension(:),intent(in) :: propriedade
        integer, intent(in) :: domx(2),domy(2)
        type(container), allocatable,dimension(:,:),intent(in) :: malha
        type(data_ptr) :: ptr
        real(dp) :: kb = 1.38064852E-23,K, auxres(8), nump, aux(2)
        type(list_t), pointer :: node
        integer ( kind = 4 ) :: np, ierr, id
        real(dp), intent(in) :: t
        
        K = 0; nump = 0;
        
        do i = domy(1),domy(2)
            do j = domx(1),domx(2)
               ! print *, 'posição', i, ',', j
                node => list_next(malha(i,j)%list)
                do while (associated(node))
                    ptr = transfer(list_get(node), ptr)
                    !calcula a energia cinética atual
                    if (propriedade(ptr%p%grupo)%x_lockdelay <= t) then
                        K = 0.5*propriedade(ptr%p%grupo)%m*(ptr%p%v(1)**2+ptr%p%v(2)**2) + K
                        nump = nump +1
                    end if
                    node => list_next(node)
                end do
            end do
        end do        
        
        ! então o processo é paralelo, vamos juntar tudo 
        aux = [K, nump]
        
        if (np > 1) then
            call MPI_GATHER(aux, 2, MPI_DOUBLE_PRECISION, auxres, 2, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
            if (id == 0) then
                K = auxres(1)+auxres(3)+auxres(5)+auxres(7)
                nump = auxres(2)+auxres(4)+auxres(6)+auxres(8)
                auxres = [K,nump,K,nump,K,nump,K,nump]
            end if
            call MPI_SCATTER(auxres, 2, MPI_DOUBLE_PRECISION, aux, 2, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
        end if
        K = aux(1)/aux(2)
   
    end function comp_Kglobal
        
    !Corrige energia cinética/configura temepratura global
    subroutine corr_Kglobal(malha,domx,domy,Td,propriedade, np,id,tt)
        use linkedlist
        use mod1
        use data
        use mod0
        use mpi
 
        type(prop_grupo), allocatable,dimension(:),intent(in) :: propriedade
        real(dp),intent(in) :: Td!energia cinética
        real(dp) :: T,K,kb = 1.38064852E-23, beta
        integer, intent(in) :: domx(2),domy(2)
        type(container), allocatable,dimension(:,:),intent(in) :: malha
        type(data_ptr) :: ptr
        type(list_t), pointer :: node
        integer ( kind = 4 ) :: np, id
        real(dp), intent(in) :: tt
        !calcula temperatura atual
        ! T = (2/(Nf*kb))*Ekin com Nf = número de graus de liberadade 
        T = comp_Kglobal(malha,domx,domy,propriedade,np,id,tt)

        beta = sqrt(Td/T) ! aqui o T é o Beta
        do i = domy(1),domy(2)
            do j = domx(1),domx(2)
               ! print *, 'posição', i, ',', j
                node => list_next(malha(i,j)%list)
                do while (associated(node))
                    ptr = transfer(list_get(node), ptr)
                    if (propriedade(ptr%p%grupo)%x_lockdelay <= t) then
                        !calcula a energia cinética atual
                        ptr%p%v = beta*ptr%p%v ! aqui o T é o beta 
                    end if
                    node => list_next(node)
                end do
            end do
        end do        
    end subroutine corr_Kglobal
    
    function comp_K(malha,domx,domy,s_cells,propriedade,np,id,t) result(K)
        use linkedlist
        use mod1
        use data
        use mod0
        use mpi

        type(prop_grupo), allocatable,dimension(:),intent(in) :: propriedade
        integer, intent(in) :: domx(2),domy(2), s_cells(4) !selected cells
        type(container), allocatable,dimension(:,:),intent(in) :: malha
        type(data_ptr) :: ptr
        real(dp) :: kb = 1.38064852E-23,K, auxres(8), nump, aux(2)
        type(list_t), pointer :: node
        integer ( kind = 4 ) :: np, ierr, id
        real(dp), intent(in) :: t
        
        K = 0; nump = 0;
        
        do i = domy(1),domy(2)
            do j = domx(1),domx(2)
                if (j >= s_cells(1) .and. j <= s_cells(2) .and. i >= s_cells(3) .and. i <= s_cells(4)) then
                    node => list_next(malha(i,j)%list)
                    do while (associated(node))
                        ptr = transfer(list_get(node), ptr)
                        !calcula a energia cinética atual
                        if (propriedade(ptr%p%grupo)%x_lockdelay <= t) then
                            K = 0.5*propriedade(ptr%p%grupo)%m*(ptr%p%v(1)**2+ptr%p%v(2)**2) + K
                            nump = nump +1
                        end if
                        node => list_next(node)
                    end do
                end if
            end do
        end do        
        
        ! então o processo é paralelo, vamos juntar tudo 
        aux = [K, nump]
        
        if (np > 1) then
            call MPI_GATHER(aux, 2, MPI_DOUBLE_PRECISION, auxres, 2, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
            if (id == 0) then
                K = auxres(1)+auxres(3)+auxres(5)+auxres(7)
                nump = auxres(2)+auxres(4)+auxres(6)+auxres(8)
                auxres = [K,nump,K,nump,K,nump,K,nump]
            end if
            call MPI_SCATTER(auxres, 2, MPI_DOUBLE_PRECISION, aux, 2, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
        end if
        K = aux(1)/aux(2)
   
    end function comp_K 

    subroutine corr_K(malha,domx,domy,cold_cells,hot_cells,Td_hot,Td_cold,propriedade, np,id,tt)
        use linkedlist
        use mod1
        use data
        use mod0
        use mpi
 
        type(prop_grupo), allocatable,dimension(:),intent(in) :: propriedade
        real(dp),intent(in) :: Td_hot, Td_cold!energia cinética
        real(dp) :: T_hot, T_cold,K,kb = 1.38064852E-23, betac, betah
        integer, intent(in) :: domx(2),domy(2), cold_cells(4), hot_cells(4)
        type(container), allocatable,dimension(:,:),intent(in) :: malha
        type(data_ptr) :: ptr
        type(list_t), pointer :: node
        integer ( kind = 4 ) :: np, id
        real(dp), intent(in) :: tt
        !calcula temperatura atual
        ! T = (2/(Nf*kb))*Ekin com Nf = número de graus de liberadade 
        T_hot =  comp_K(malha,domx,domy,hot_cells,propriedade,np,id,tt)
        T_cold = comp_K(malha,domx,domy,cold_cells,propriedade,np,id,tt)
        ! print*, "T =", T, "Td =", Td
        ! if (id  == 0) read(*,*)    
        ! call MPI_barrier(MPI_COMM_WORLD, ierr)
        betac = sqrt(Td_cold/T_cold) 
        betah = sqrt(Td_hot/T_hot)
        do i = domy(1),domy(2)
            do j = domx(1),domx(2)
                if (j >= cold_cells(1) .and. j <= cold_cells(2) .and. i >= cold_cells(3) .and. i <= cold_cells(4)) then
                    node => list_next(malha(i,j)%list)
                    do while (associated(node))
                        ptr = transfer(list_get(node), ptr)
                        if (propriedade(ptr%p%grupo)%x_lockdelay <= t) then
                            !calcula a energia cinética atual
                            ptr%p%v = betac*ptr%p%v ! aqui o T é o beta 
                        end if
                        node => list_next(node)
                    end do
                else if (j >= hot_cells(1) .and. j <= hot_cells(2) .and. i >= hot_cells(3) .and. i <= hot_cells(4)) then
                    node => list_next(malha(i,j)%list)
                    do while (associated(node))
                        ptr = transfer(list_get(node), ptr)
                        if (propriedade(ptr%p%grupo)%x_lockdelay <= t) then
                            !calcula a energia cinética atual
                            ptr%p%v = betah*ptr%p%v ! aqui o T é o beta 
                        end if
                        node => list_next(node)
                    end do
                end if
            end do
        end do        
    end subroutine corr_K

    !força de ficção 
    function comp_fric(r,v,fric_term) result(f)
        use mod1
        real(dp), intent(in) :: r(2),v(2),fric_term
        real(dp) :: f(2)
        
        f = (fric_term/(r(1)**2 + r(2)**2))*(v(1)*r(1)+v(2)*r(2))*[r(1),r(2)]
    end function comp_fric
        
    
    ! atualiza forças
    subroutine comp_F(GField, mesh,malha,propriedade,r_cut,domx,domy, ids, id, t)
        use linkedlist
        use mod1
        use data
        use mod0
        use mpi
        
        type(prop_grupo), allocatable,dimension(:),intent(in) :: propriedade
        type(container), allocatable,dimension(:,:),intent(in) :: malha
        real(dp), intent(in) :: GField(2), t
        real(dp) :: sigma, epsil, sigma_a, epsil_a,sigma_b, epsil_b, rcut,r_cut, fric_term !fric_term = força de ficção
        real(dp) :: x1(2),v1(2),x2(2),v2(2), rs1, rs2, coss, sine   
        integer :: i,j,ct = 0 !,ptr, ptrn
        integer, intent(in) :: mesh(:),domx(2),domy(2), id
        real(dp) :: Fi(2)=0,r, aux2(2),fR(2), fric_term1, fric_term2, dox(2), doy(2)
        type(list_t), pointer :: node, next_node
        type(data_ptr) :: ptr,ptrn
        integer ( kind = 4 ), intent(in) :: ids(8)
 
        !Lennard Jones
        fR = [0,0]
        dox = domx 
        doy = domy 
        if (sum(ids(1:4)) > -4) then
            !caso paralelo
            if (domx(1) > 1) dox(1) = domx(1) - 1
            if (domy(1) > 1) doy(1) = domy(1) - 1
            if (domx(2) < mesh(1)+2) dox(2) = domx(2) + 1
            if (domy(2) < mesh(2)+2) doy(2) = domy(2) + 1
        else
            dox = domx 
            doy = domy 
        end if 
        
        do i = doy(1),doy(2) ! i é linha
            do j = dox(1),dox(2)
                node => list_next(malha(i,j)%list)
                ! if (associated(node)) then
                !     ptr = transfer(list_get(node), ptr)
                !     ! print*, "L 253", ptr%p%n
                ! end if
                do while (associated(node))
                    ! print*, "Encontrada", ct, "celulas, id=", id
                    ptr = transfer(list_get(node), ptr) !particula selecionada
                    ptr%p%flag = .true. ! indica ao comp_x que a partícula precisa ser calculada
                    x1 = ptr%p%x
                    v1 = ptr%p%v
                    m1 = propriedade(ptr%p%grupo)%m
                    rs1 = propriedade(ptr%p%grupo)%rs !raio sólido 
                    fric_term1 = propriedade(ptr%p%grupo)%fric_term
                    sigma_a = propriedade(ptr%p%grupo)%sigma
                    ! rcut = r_cut*sigma
                    epsil_a = propriedade(ptr%p%grupo)%epsilon 
                    if (propriedade(ptr%p%grupo)%x_lockdelay <= t) then
                        ptr%p%F = ptr%p%F + GField*m1
                    end if
                    ! print '("x1  =", f10.6, " ", f10.6, " n ", i2, " cell ",i2," ",i2)',x1(1),x1(2), ptr%p%n,i,j
                    ! if (id == 0) read(*,*)
                    !calcular a força desta com todas as outras partículas
                    next_node => list_next(node) ! próxima partícula da célula
                    node => list_next(node) ! a ser computado com a particula selecionada
             
                    ! NA PRÓPRIA CELULA 
                    do while (associated(node))
                        ptrn = transfer(list_get(node), ptrn) !outra particula selecionada
                        sigma_b = propriedade(ptrn%p%grupo)%sigma
                        epsil_b = propriedade(ptrn%p%grupo)%epsilon 
                        rs2 = propriedade(ptrn%p%grupo)%rs 
                        fric_term2 = propriedade(ptrn%p%grupo)%fric_term
                        ! Lorenz-Betherlot rule for mixing epsilon sigma 
                        if (sigma_a > sigma_b) then
                            rcut = r_cut*sigma_a + rs1 + rs2
                            sigma = 0.5*(sigma_a + sigma_b)
                            epsil = sqrt(epsil_a *epsil_b )
                        else if (sigma_a < sigma_b) then 
                            rcut = r_cut*sigma_b + rs1 + rs2
                            sigma = 0.5*(sigma_a + sigma_b)
                            epsil = sqrt(epsil_a *epsil_b )
                        else 
                            rcut = r_cut*sigma_a + rs1 + rs2
                            sigma = sigma_a
                            epsil = epsil_a
                        end if 
                        x2 = ptrn%p%x
                        v2 = ptrn%p%v
                        r = sqrt((x1(1)-x2(1))**2 + (x1(2)-x2(2))**2)
                        coss = (x1(1)-x2(1))/r 
                        sine = (x1(2)-x2(2))/r 
                        r = r - rs1 - rs2 !raio
                        ! print*, "L 389 r", r, "id",id
                        if (r <= rcut) then
                            aux2 = -(1/r**2)*(sigma/r)**6*(1-2*(sigma/r)**6)*24*epsil* & 
                            [(x1(1)-x2(1)) - (rs1+rs2)*coss, (x1(2)-x2(2)) - (rs1+rs2)*sine]  
                            ! print*, "L 395 r", r, "id",id
  !print '("(x1(1)-x2(1)) ", f7.3 ," (rs1+rs2)*coss ", f7.3, " (x1(2)-x2(2)) ", f7.3, " (rs1+rs2)*sine ", f7.3 )', (x1(1)-x2(1)), &
  !(rs1+rs2)*coss, (x1(2)-x2(2)), (rs1+rs2)*sine
                            fric_term = (fric_term1+fric_term2)/2
                            if (fric_term > 0) then
                                fR = comp_fric([-(x1(1)-x2(1))+ (rs1+rs2)*coss,-(x1(2)-x2(2))+(rs1+rs2)*sine], &
                                [-(v1(1)-v2(1)),-(v1(2)-v2(2))],fric_term)
                            end if
                            ptr%p%F = aux2 + ptr%p%F +fR
                            ptrn%p%F = -aux2 + ptrn%p%F - fR
                         end if
                        node => list_next(node) ! próxima partícula da célula
                    end do

                    !Células ao redor  !i é linha, j é coluna
                    if (i /= mesh(2)+2) then ! se não for a última linha
                        if (j == mesh(1)+2) then ! se for a última coluna
                            node => list_next(malha(i+1,j)%list) !interagirá com a próxima linha apenas
                            do while (associated(node))
                                ptrn = transfer(list_get(node), ptrn) !outra particula selecionada
                                sigma_b = propriedade(ptrn%p%grupo)%sigma
                                epsil_b = propriedade(ptrn%p%grupo)%epsilon 
                                rs2 = propriedade(ptrn%p%grupo)%rs 
                                fric_term2 = propriedade(ptrn%p%grupo)%fric_term
                                ! Lorenz-Betherlot rule for mixing epsilon sigma 
                                if (sigma_a > sigma_b) then
                                    rcut = r_cut*sigma_a + (rs1 + rs2)
                                    sigma = 0.5*(sigma_a + sigma_b)
                                    epsil = sqrt(epsil_a*epsil_b)
                                else if (sigma_a < sigma_b) then 
                                    rcut = r_cut*sigma_b + (rs1 + rs2)
                                    sigma = 0.5*(sigma_a + sigma_b)
                                    epsil = sqrt(epsil_a*epsil_b)
                                else 
                                    rcut = r_cut*sigma_a + (rs1 + rs2)
                                    sigma = sigma_a
                                    epsil = epsil_a
                                end if 

                                x2 = ptrn%p%x
                                v2 = ptrn%p%v
                                r = sqrt((x1(1)-x2(1))**2 + (x1(2)-x2(2))**2)
                                coss = (x1(1)-x2(1))/r 
                                sine = (x1(2)-x2(2))/r 
                                r = r - rs1 - rs2 !raio
                                ! print*, "L 666 r", r, "id",id
                                if (r <= rcut) then
                                    aux2 = -(1/r**2)*(sigma/r)**6*(1-2*(sigma/r)**6)*24*epsil* & 
                                        [(x1(1)-x2(1)) - (rs1+rs2)*coss, (x1(2)-x2(2)) - (rs1+rs2)*sine]
                                    fric_term = (fric_term1+fric_term2)/2 
                                    if (fric_term > 0) then
                                        fR = comp_fric([-(x1(1)-x2(1))+ (rs1+rs2)*coss,-(x1(2)-x2(2))+(rs1+rs2)*sine], &
                                [-(v1(1)-v2(1)),-(v1(2)-v2(2))],fric_term)
                                    end if
                                    ptr%p%F = aux2 + ptr%p%F +fR
                                    ptrn%p%F = -aux2 + ptrn%p%F - fR
                                end if
                                node => list_next(node) ! próxima partícula da célula
                            end do
                            
                            if (j /= 1) then !se não for a primeira coluna 
                                node => list_next(malha(i+1,j-1)%list) !interagirá com a próxima linha apenas
                                do while (associated(node))
                                    ptrn = transfer(list_get(node), ptrn) !outra particula selecionada
                                    sigma_b = propriedade(ptrn%p%grupo)%sigma
                                    epsil_b = propriedade(ptrn%p%grupo)%epsilon 
                                    rs2 = propriedade(ptrn%p%grupo)%rs 
                                    fric_term2 = propriedade(ptrn%p%grupo)%fric_term
                                    ! Lorenz-Betherlot rule for mixing epsilon sigma 
                                    if (sigma_a > sigma_b) then
                                        rcut = r_cut*sigma_a + (rs1 + rs2)
                                        sigma = 0.5*(sigma_a + sigma_b)
                                        epsil = sqrt(epsil_a*epsil_b)
                                    else if (sigma_a < sigma_b) then 
                                        rcut = r_cut*sigma_b + (rs1 + rs2)
                                        sigma = 0.5*(sigma_a + sigma_b)
                                        epsil = sqrt(epsil_a*epsil_b)
                                    else 
                                        rcut = r_cut*sigma_a + (rs1 + rs2)
                                        sigma = sigma_a
                                        epsil = epsil_a
                                    end if 

                                    x2 = ptrn%p%x
                                    v2 = ptrn%p%v
                                    r = sqrt((x1(1)-x2(1))**2 + (x1(2)-x2(2))**2)
                                    coss = (x1(1)-x2(1))/r 
                                    sine = (x1(2)-x2(2))/r 
                                    r = r - rs1 - rs2 !raio
                                    ! print*, "L 710 r", r, "id",id
                                    if (r <= rcut) then
                                        aux2 = -(1/r**2)*(sigma/r)**6*(1-2*(sigma/r)**6)*24*epsil* & 
                                            [(x1(1)-x2(1)) - (rs1+rs2)*coss, (x1(2)-x2(2)) - (rs1+rs2)*sine]
                                        fric_term = (fric_term1+fric_term2)/2
                                        if (fric_term > 0) then
                                            fR = comp_fric([-(x1(1)-x2(1))+ (rs1+rs2)*coss,-(x1(2)-x2(2))+(rs1+rs2)*sine], &
                                [-(v1(1)-v2(1)),-(v1(2)-v2(2))],fric_term)
                                        end if
                                        ptr%p%F = aux2 + ptr%p%F +fR
                                        ptrn%p%F = -aux2 + ptrn%p%F - fR
                                    end if
                                    node => list_next(node) ! próxima partícula da célula
                                end do                            
                            end if
                        else
                             !interagirá com a próxima linha e coluna, e na diagonal
                            node => list_next(malha(i,j+1)%list) 
                            
                            do while (associated(node))
                                ptrn = transfer(list_get(node), ptrn) !outra particula selecionada
                                sigma_b = propriedade(ptrn%p%grupo)%sigma
                                epsil_b = propriedade(ptrn%p%grupo)%epsilon 
                                rs2 = propriedade(ptrn%p%grupo)%rs
                                fric_term2 = propriedade(ptrn%p%grupo)%fric_term 
                                ! Lorenz-Betherlot rule for mixing epsilon sigma 
                                if (sigma_a > sigma_b) then
                                    rcut = r_cut*sigma_a + (rs1 + rs2)
                                    sigma = 0.5*(sigma_a + sigma_b)
                                    epsil = sqrt(epsil_a*epsil_b)
                                else if (sigma_a < sigma_b) then 
                                    rcut = r_cut*sigma_b + (rs1 + rs2)
                                    sigma = 0.5*(sigma_a + sigma_b)
                                    epsil = sqrt(epsil_a*epsil_b)
                                else 
                                    rcut = r_cut*sigma_a + (rs1 + rs2)
                                    sigma = sigma_a
                                    epsil = epsil_a
                                end if 

                                x2 = ptrn%p%x
                                v2 = ptrn%p%v
                                
                                r = sqrt((x1(1)-x2(1))**2 + (x1(2)-x2(2))**2)
                                coss = (x1(1)-x2(1))/r 
                                sine = (x1(2)-x2(2))/r 
                                r = r - rs1 - rs2 !raio
                                ! print*, "L 757 r", r, "id",id
                                if (r <= rcut) then
                                    aux2 = -(1/r**2)*(sigma/r)**6*(1-2*(sigma/r)**6)*24*epsil* & 
                                        [(x1(1)-x2(1)) - (rs1+rs2)*coss, (x1(2)-x2(2)) - (rs1+rs2)*sine] 
                                    fric_term = (fric_term1+fric_term2)/2
                                    if (fric_term > 0) then
                                        fR = comp_fric([-(x1(1)-x2(1))+ (rs1+rs2)*coss,-(x1(2)-x2(2))+(rs1+rs2)*sine], &
                                [-(v1(1)-v2(1)),-(v1(2)-v2(2))],fric_term)
                                    end if
                                    ptr%p%F = aux2 + ptr%p%F +fR
                                    ptrn%p%F = -aux2 + ptrn%p%F - fR
                                end if
                                node => list_next(node) ! próxima partícula da célula                                
                            end do

                            node => list_next(malha(i+1,j)%list) !interagirá com a próxima linha 
                            do while (associated(node))
                                ptrn = transfer(list_get(node), ptrn) !outra particula selecionada
                                sigma_b = propriedade(ptrn%p%grupo)%sigma
                                epsil_b = propriedade(ptrn%p%grupo)%epsilon 
                                rs2 = propriedade(ptrn%p%grupo)%rs 
                                fric_term2 = propriedade(ptrn%p%grupo)%fric_term
                                ! Lorenz-Betherlot rule for mixing epsilon sigma 
                                if (sigma_a > sigma_b) then
                                    rcut = r_cut*sigma_a + (rs1 + rs2)
                                    sigma = 0.5*(sigma_a + sigma_b)
                                    epsil = sqrt(epsil_a*epsil_b)
                                else if (sigma_a < sigma_b) then 
                                    rcut = r_cut*sigma_b + (rs1 + rs2)
                                    sigma = 0.5*(sigma_a + sigma_b)
                                    epsil = sqrt(epsil_a*epsil_b)
                                else 
                                    rcut = r_cut*sigma_a + (rs1 + rs2)
                                    sigma = sigma_a
                                    epsil = epsil_a
                                end if 
                                x2 = ptrn%p%x
                                v2 = ptrn%p%v
                                r = sqrt((x1(1)-x2(1))**2 + (x1(2)-x2(2))**2)
                                coss = (x1(1)-x2(1))/r 
                                sine = (x1(2)-x2(2))/r 
                                r = r - rs1 - rs2 !raio
                                ! print*, "L 798 r", r, "id",id
                                if (r <= rcut) then
                                    aux2 = -(1/r**2)*(sigma/r)**6*(1-2*(sigma/r)**6)*24*epsil* & 
                                        [(x1(1)-x2(1)) - (rs1+rs2)*coss, (x1(2)-x2(2)) - (rs1+rs2)*sine] 
                                    fric_term = (fric_term1+fric_term2)/2
                                    if (fric_term > 0) then
                                        fR = comp_fric([-(x1(1)-x2(1))+ (rs1+rs2)*coss,-(x1(2)-x2(2))+(rs1+rs2)*sine], &
                                [-(v1(1)-v2(1)),-(v1(2)-v2(2))],fric_term)
                                    end if
                                    ptr%p%F = aux2 + ptr%p%F +fR
                                    ptrn%p%F = -aux2 + ptrn%p%F - fR
                                end if
                                node => list_next(node) ! próxima partícula da célula                                
                            end do
                            
                            node => list_next(malha(i+1,j+1)%list) !interagirá com a próxima linha e coluna
                            do while (associated(node))
                                ptrn = transfer(list_get(node), ptrn) !outra particula selecionada
                                sigma_b = propriedade(ptrn%p%grupo)%sigma
                                epsil_b = propriedade(ptrn%p%grupo)%epsilon 
                                rs2 = propriedade(ptrn%p%grupo)%rs 
                                fric_term2 = propriedade(ptrn%p%grupo)%fric_term
                                ! Lorenz-Betherlot rule for mixing epsilon sigma 
                                if (sigma_a > sigma_b) then
                                    rcut = r_cut*sigma_a + (rs1 + rs2)
                                    sigma = 0.5*(sigma_a + sigma_b)
                                    epsil = sqrt(epsil_a*epsil_b)
                                else if (sigma_a < sigma_b) then 
                                    rcut = r_cut*sigma_b + (rs1 + rs2)
                                    sigma = 0.5*(sigma_a + sigma_b)
                                    epsil = sqrt(epsil_a*epsil_b)
                                else 
                                    rcut = r_cut*sigma_a + (rs1 + rs2)
                                    sigma = sigma_a
                                    epsil = epsil_a
                                end if 

                                x2 = ptrn%p%x
                                v2 = ptrn%p%v
                                r = sqrt((x1(1)-x2(1))**2 + (x1(2)-x2(2))**2)
                                coss = (x1(1)-x2(1))/r 
                                sine = (x1(2)-x2(2))/r 
                                r = r - rs1 - rs2 !raio
                                ! print*, "L 841 r", r, "id",id
                                if (r <= rcut) then
                                    aux2 = -(1/r**2)*(sigma/r)**6*(1-2*(sigma/r)**6)*24*epsil* & 
                                        [(x1(1)-x2(1)) - (rs1+rs2)*coss, (x1(2)-x2(2)) - (rs1+rs2)*sine] 
                                    fric_term = (fric_term1+fric_term2)/2
                                    if (fric_term > 0) then
                                        fR = comp_fric([-(x1(1)-x2(1))+ (rs1+rs2)*coss,-(x1(2)-x2(2))+(rs1+rs2)*sine], &
                                [-(v1(1)-v2(1)),-(v1(2)-v2(2))],fric_term)
                                    end if
                                    ptr%p%F = aux2 + ptr%p%F +fR
                                    ptrn%p%F = -aux2 + ptrn%p%F - fR
                                end if
                                node => list_next(node) ! próxima partícula da célula                                
                            end do
                            
                            if (j /= 1) then !se não for a primeira coluna 
                                node => list_next(malha(i+1,j-1)%list) !interagirá com a próxima linha apenas
                                do while (associated(node))
                                    ptrn = transfer(list_get(node), ptrn) !outra particula selecionada
                                    sigma_b = propriedade(ptrn%p%grupo)%sigma
                                    epsil_b = propriedade(ptrn%p%grupo)%epsilon 
                                    rs2 = propriedade(ptrn%p%grupo)%rs 
                                    fric_term2 = propriedade(ptrn%p%grupo)%fric_term
                                    ! Lorenz-Betherlot rule for mixing epsilon sigma 
                                    if (sigma_a > sigma_b) then
                                        rcut = r_cut*sigma_a + (rs1 + rs2)
                                        sigma = 0.5*(sigma_a + sigma_b)
                                        epsil = sqrt(epsil_a*epsil_b)
                                    else if (sigma_a < sigma_b) then 
                                        rcut = r_cut*sigma_b + (rs1 + rs2)
                                        sigma = 0.5*(sigma_a + sigma_b)
                                        epsil = sqrt(epsil_a*epsil_b)
                                    else 
                                        rcut = r_cut*sigma_a + (rs1 + rs2)
                                        sigma = sigma_a
                                        epsil = epsil_a
                                    end if 

                                    x2 = ptrn%p%x
                                    v2 = ptrn%p%v
                                    r = sqrt((x1(1)-x2(1))**2 + (x1(2)-x2(2))**2)
                                    coss = (x1(1)-x2(1))/r 
                                    sine = (x1(2)-x2(2))/r 
                                    r = r - rs1 - rs2 !raio
                                    ! print*, "L 885 r", r, "id",id
                                    if (r <= rcut) then
                                        aux2 = -(1/r**2)*(sigma/r)**6*(1-2*(sigma/r)**6)*24*epsil* & 
                                            [(x1(1)-x2(1)) - (rs1+rs2)*coss, (x1(2)-x2(2)) - (rs1+rs2)*sine] 
                                        fric_term = (fric_term1+fric_term2)/2
                                        if (fric_term > 0) then
                                            fR = comp_fric([-(x1(1)-x2(1))+ (rs1+rs2)*coss,-(x1(2)-x2(2))+(rs1+rs2)*sine], &
                                [-(v1(1)-v2(1)),-(v1(2)-v2(2))],fric_term)
                                        end if
                                        ptr%p%F = aux2 + ptr%p%F +fR
                                        ptrn%p%F = -aux2 + ptrn%p%F - fR
                                    end if
                                    node => list_next(node) ! próxima partícula da célula
                                end do                            
                            end if
                        end if
                        
                    else ! se for a última lina, só interage com a celua ao lado 
                        node => list_next(malha(i,j+1)%list) !interagirá com a próxima linha e coluna
                        do while (associated(node))
                            ptrn = transfer(list_get(node), ptrn) !outra particula selecionada
                            sigma_b = propriedade(ptrn%p%grupo)%sigma
                            epsil_b = propriedade(ptrn%p%grupo)%epsilon 
                            rs2 = propriedade(ptrn%p%grupo)%rs 
                            fric_term2 = propriedade(ptrn%p%grupo)%fric_term
                            ! Lorenz-Betherlot rule for mixing epsilon sigma 
                            if (sigma_a > sigma_b) then
                                rcut = r_cut*sigma_a + (rs1 + rs2)
                                sigma = 0.5*(sigma_a + sigma_b)
                                epsil = sqrt(epsil_a*epsil_b)
                            else if (sigma_a < sigma_b) then 
                                rcut = r_cut*sigma_b + (rs1 + rs2)
                                sigma = 0.5*(sigma_a + sigma_b)
                                epsil = sqrt(epsil_a*epsil_b)
                            else 
                                rcut = r_cut*sigma_a + (rs1 + rs2)
                                sigma = sigma_a
                                epsil = epsil__a
                            end if 

                            x2 = ptrn%p%x
                            v2 = ptrn%p%v
                            r = sqrt((x1(1)-x2(1))**2 + (x1(2)-x2(2))**2)
                            coss = (x1(1)-x2(1))/r 
                            sine = (x1(2)-x2(2))/r 
                            r = r - rs1 - rs2 !raio
                            ! print*, "L 931 r", r, "id",id
                            if (r <= rcut) then
                                aux2 = -(1/r**2)*(sigma/r)**6*(1-2*(sigma/r)**6)*24*epsil* & 
                                    [(x1(1)-x2(1)) - (rs1+rs2)*coss, (x1(2)-x2(2)) - (rs1+rs2)*sine] 
                                fric_term = (fric_term1+fric_term2)/2
                                if (fric_term > 0) then
                                    fR = comp_fric([-(x1(1)-x2(1))+ (rs1+rs2)*coss,-(x1(2)-x2(2))+(rs1+rs2)*sine], &
                                [-(v1(1)-v2(1)),-(v1(2)-v2(2))],fric_term)
                                end if
                                ptr%p%F = aux2 + ptr%p%F +fR
                                ptrn%p%F = -aux2 + ptrn%p%F - fR
                            end if
                            node => list_next(node) ! próxima partícula da célula                                
                        end do
                    end if
                    node => next_node
                end do
            end do
        end do
         !print*, 'Linha 143'
        ! if (id == 0) read(*,*)
        ! call MPI_barrier(MPI_COMM_WORLD, ierr) 
    end subroutine comp_F
    
    ! atualiza posições
    subroutine comp_x(icell,jcell,malha,N,mesh,propriedade, dx_max,t,dt,ids,LT,domx,domy,wall,id, np)
        use linkedlist
        use mod1
        use data
        use mod0
        use mpi
        use matprint

        character(1) :: north, south, east, west
        type(prop_grupo), allocatable,dimension(:),intent(in) :: propriedade
        type(container), allocatable,dimension(:,:),intent(in) :: malha
        ! integer, intent(inout) :: mic(:)
        integer :: i,j,k, cell(2), status(MPI_STATUS_SIZE), count, destD,dx1,dx2,dy1,dy2
        integer, intent(in) :: N,mesh(2),domx(2),domy(2)
        type(list_t), pointer :: node, previous_node
        real(dp), intent(in) :: dt, t, dx_max
        character(4), intent(in) :: wall
        type(data_ptr) :: ptr
        real(dp) :: x(2),m, dx(2)
        real(dp), intent(in) :: icell(:), jcell(:)
        integer ( kind = 4 ), intent(in) :: ids(8), id, np
        integer ( kind = 4 ) :: tag,cont_db(8),cont_int(8)
        type(lstr) :: LT
        !  real(dp),intent(inout) :: celula(:,:)
         ! IDS correspondem às celulas [N,S,E,W]
        cont_db = 0
        cont_int = 0

        north = wall(1:1)        
        south = wall(2:2)
        east = wall(3:3)
        west = wall(4:4)
        !  if (id == 0) read(*,*) 
      
         ! Esvazia as celulas emprestadas 
        if (domy(1) > 1) then
            i = domy(1)-1
            do j = domx(1),domx(2)
                previous_node => malha(i,j)%list
                node => list_next(malha(i,j)%list) 
                ! print*, associated(node), "A", id
                do while (associated(node)) 
                    ptr = transfer(list_get(node), ptr)
                    deallocate(ptr%p)
                    call list_remove(previous_node)
                    node => list_next(previous_node)
                    ! node => list_next(node)
                end do
            end do
        else if (domy(1) == 1) then !esvazia celulas fantasma
            i = domy(1)
            do j = domx(1),domx(2)
                previous_node => malha(i,j)%list
                node => list_next(malha(i,j)%list) 
                ! print*, associated(node), "A", id
                do while (associated(node)) 
                    ptr = transfer(list_get(node), ptr)
                    deallocate(ptr%p)
                    call list_remove(previous_node)
                    node => list_next(previous_node)
                    ! node => list_next(node)
                end do
            end do            
        end if 
        if (domy(2) < mesh(2)+2) then
            i = domy(2)+1
            do j = domx(1),domx(2)
                previous_node => malha(i,j)%list
                node => list_next(malha(i,j)%list)
                ! print*, associated(node), "B", id
                do while (associated(node)) 
                    ptr = transfer(list_get(node), ptr)
                    deallocate(ptr%p)
                    call list_remove(previous_node)
                    node => list_next(previous_node)
                    ! node => list_next(node)
                end do 
            end do
        else if  (domy(2) == mesh(2)+2) then 
            i = domy(2)
            do j = domx(1),domx(2)
                previous_node => malha(i,j)%list
                node => list_next(malha(i,j)%list)
                ! print*, associated(node), "B", id
                do while (associated(node)) 
                    ptr = transfer(list_get(node), ptr)
                    deallocate(ptr%p)
                    call list_remove(previous_node)
                    node => list_next(previous_node)
                    ! node => list_next(node)
                end do 
            end do
        end if 
        if (domx(2) < mesh(1)+2) then
            j = domx(2)+1
            do i = domy(1),domy(2)
                previous_node => malha(i,j)%list
                node => list_next(malha(i,j)%list) 
                ! print*, associated(node), "D", id
                do while (associated(node))
                    ptr = transfer(list_get(node), ptr)
                    deallocate(ptr%p)
                    call list_remove(previous_node)
                    node => list_next(previous_node)
                    ! node => list_next(node)
                end do 
            end do
        else if (domx(2) == mesh(1)+2) then 
            j = domx(2)
            do i = domy(1),domy(2)
                previous_node => malha(i,j)%list
                node => list_next(malha(i,j)%list) 
                ! print*, associated(node), "D", id
                do while (associated(node))
                    ptr = transfer(list_get(node), ptr)
                    deallocate(ptr%p)
                    call list_remove(previous_node)
                    node => list_next(previous_node)
                    ! node => list_next(node)
                end do 
            end do
        end if  
        if (domx(1) > 1) then
            j = domx(1)-1
            do i = domy(1),domy(2)
                node => list_next(malha(i,j)%list)
                previous_node => malha(i,j)%list
                do while (associated(node)) 
                    ptr = transfer(list_get(node), ptr)
                    deallocate(ptr%p)
                    call list_remove(previous_node)
                    node => list_next(previous_node)
                    ! node => list_next(node)
                end do
            end do
        else if (domx(1) == 1) then
            j = domx(1)
            do i = domy(1),domy(2)
                node => list_next(malha(i,j)%list)
                previous_node => malha(i,j)%list
                do while (associated(node)) 
                    ptr = transfer(list_get(node), ptr)
                    deallocate(ptr%p)
                    call list_remove(previous_node)
                    node => list_next(previous_node)
                    ! node => list_next(node)
                end do
            end do
        end if 
        ! PARA 4 processos. Limpa o que está na diagonal adjacente
        if (np > 1 .and. sum(ids(5:8)) > -4) then  
            if (domy(1) /= 1) i = domy(1) -1
            if (domy(2) /= mesh(2)+2) i = domy(2) +1
            if (domx(1) /= 1) j = domx(1) -1
            if (domx(2) /= mesh(1) + 2) j = domx(2) +1 
            previous_node => malha(i,j)%list
            node => list_next(malha(i,j)%list)
            do while (associated(node)) 
                ptr = transfer(list_get(node), ptr)
                deallocate(ptr%p)
                call list_remove(previous_node)
                node => list_next(previous_node)
                ! node => list_next(node)
            end do
        end if 
        ! print*, "cont_int", cont_int
        ! if (id == 0) read(*,*) 
        ! call MPI_barrier(MPI_COMM_WORLD, ierr)
        dx1 = domx(1)
        dx2 = domx(2)
        dy1 = domy(1)
        dy2 = domy(2)
        if (dy1 == 1) dy1 = 2
        if (dx1 == 1) dx1 = 2
        if (dy2 == mesh(2)+2) dy2 = mesh(2)+1
        if (dx2 == mesh(1)+2) dx2 = mesh(1)+1

        do i = dy1, dy2
            do j = dx1,dx2
                previous_node => malha(i,j)%list
                node => list_next(malha(i,j)%list)
                 ! ptr%p%flag = true significa que ainda não foi computado para evitar loop infinito
                
                if (associated(node))  then
                    ptr = transfer(list_get(node), ptr)
                end if
                
                do while (associated(node) .and. ptr%p%flag)
                    ! print*, "L PART", i,j, "id", id   
                    m = propriedade(ptr%p%grupo)%m
                    ! aqui vemos se já passou o tempo que as posições ficaram travadas. 
                    if (propriedade(ptr%p%grupo)%x_lockdelay > t) ptr%p%flag = .false. 
                         
                    if (ptr%p%flag) then 
                        dx(1) = dt*ptr%p%v(1) + ptr%p%F(1)*dt**2/(2*m)
                        dx(2) = dt*ptr%p%v(2) + ptr%p%F(2)*dt**2/(2*m)
                    else 
                        dx = [0,0]
                    end if 
                    ptr%p%flag = .false.
                    if ((dx(1)**2 + dx(2)**2) >= dx_max) then
                        print '("Particulas rápidas demais! dx =", f18.5, " ", f18.5, " | n ", i4)', dx(1), dx(2), ptr%p%n 
                        print*, "F",ptr%p%F, "v", ptr%p%v
                        print*, "x", ptr%p%x(1), ptr%p%x(2)
                        call system('killall lennard.out')
                        dx = [dx(1)/dx(1),dx(1)/dx(1)]*dx_max
                        read(*,*)
                    end if
                    if (isnan(dx(1)) .or. isnan(dx(2))) then
                        print '("NaN NaN Nan Batman! dx =", f18.5, " ", f18.5, " | n ", i4)', dx(1), dx(2), ptr%p%n 
                        print*, "F",ptr%p%F, "v", ptr%p%v
                        print*, "x", ptr%p%x(1), ptr%p%x(2)
                        call system('killall lennard.out')
                    end if
                    
                    ptr%p%x(1) = ptr%p%x(1) + dx(1) !dt*ptr%p%v(1) + ptr%p%F(1)*dt**2/(2*m)
                    ptr%p%x(2) = ptr%p%x(2) + dx(2) !dt*ptr%p%v(2) + ptr%p%F(2)*dt**2/(2*m)
                    ! print*,'F_0  = ',ptr%p%F(1),ptr%p%F(2), "n", ptr%p%n, "id", id 
                    ! print '("x_0  = [",f10.3,", ",f10.3, "] ", "i, j =", i2, " ", i2, " n ",i2 )',ptr%p%x(1),ptr%p%x(2),i,j, ptr%p%n!"id", id, "n", ptr%p%n
                    ! print '("v_0  = [",f10.3,", ",f10.3, "] ", "i, j =", i2, " ", i2, " n ",i2 )',ptr%p%v(1),ptr%p%v(2),i,j, ptr%p%n!"id", id, "n", ptr%p%n
                    ! Ordena as celulas
                    x = ptr%p%x
                    ! teste se a partícula pode pular celulas
                    ! Teste se escapou para a esquerda
                    if (x(1) <= jcell(j-1)) then
                        cell = [i,j-1]
                        !para esquerda para cima
                        if (x(2) <= icell(i-1)) then
                            cell = [i-1,j-1]
                        !para esquerda e para baixo
                        else if (x(2) > icell(i)) then
                            cell = [i+1,j-1]
                        end if
                        
                    ! Teste se escapou para a direita
                    else if (x(1) > jcell(j)) then
                        cell = [i,j+1]
                        !para direita para cima
                        if (x(2) <= icell(i-1)) then
                            cell = [i-1,j+1]
                        !para direita e para baixo
                        else if (x(2) > icell(i)) then
                            cell = [i+1,j+1]
                            ! print*, "icell(i)", icell(i), "jcell(j)", jcell(j), "x=", x
                            ! ! print*, "L PUTAMERDA"
                        end if       
                     
                    else if (x(2) <=  icell(i-1)) then
                        ! print*,'!Teste se escapou para cima'
                        ! print*,x(2) ,'<',  icell(i-1)
                        cell = [i-1, j]
                    !Teste se escapou para baixo
                    else if (x(2) >  icell(i)) then
                        ! print*,'!Teste se escapou para baixo'
                        ! print*,x(2), '>',  icell(i)
                        cell = [i+1,j]
                    else 
                        cell = [i,j]
                    end if 

                    ! print*, "Contint(1)", cont_int(1)
                    
                    ! VERIFICAÇÃO SE MUDOUD DE DOMÍNIO
                    if (cell(1) /= i .or. cell(2) /= j)  then !Mudou de celula
                        ! print*, "MUDOU"
                        ! if (id == 0) read(*,*)
                        ! Se a partícula chegar na fronteira do domínio que um processador
                        ! cuida, então esta partícula é colocada na lista linkada referente
                        ! a este célula com list_change. Além disso ela precisa ser 
                        ! transferida para o outro processo vizinho para que ele possa 
                        ! calcular sua influência no domínio.
                        ! Se a partícula ultrapassar o domínio, então ela é transferida
                        ! para o próximo processo e dealocada da lista com list_remove 
                        ! ! ! print*, "L FORÇA", ptr%p%F, id
                        if (cell(2) >= domx(2) .and. cell(2) < mesh(1)+2 ) then
                            ! print*, "L 538",  cell(1), cell(2),"part", ptr%p%n,"i,j", i, j
                            ! 6 elementos
                            LT%lstrdb_E(cont_db(3)+1:cont_db(3)+6) = [x(1),x(2), &
                            ptr%p%v(1),ptr%p%v(2), ptr%p%F(1),ptr%p%F(2)]
                            ! 4 elementos
                            ! print*, "id",id,"transferindo para o leste",  [cell(1),cell(2),ptr%p%n,ptr%p%grupo]
                            LT%lstrint_E(cont_int(3)+1:cont_int(3)+4) = [cell(1),cell(2),ptr%p%n,ptr%p%grupo]
                            ! print*, "id",id,"transferindo para o leste",   LT%lstrint_E, "cont", cont_int(3)
                            ! if (cell(2) > domx(2)) then
                            !     ! print*, 'affs'
                            !     ! read(*,*)
                            !     print*, "Removido", ptr%p%n, "de malha", i,j
                            !     print*, "Irá para", LT%lstrint_E

                            cont_db(3) = cont_db(3) + 6
                            cont_int(3) = cont_int(3) + 4
                        end if 
                        if (cell(2) <= domx(1) .and. cell(2) > 1) then
                            ! print*, "L 554",  cell(1), cell(2), "part", ptr%p%n, "i,j", i, j
                            ! print*, "domx", domx
                            LT%lstrdb_W(cont_db(4)+1:cont_db(4)+6) = [x(1), x(2), ptr%p%v(1),ptr%p%v(2), ptr%p%F(1),ptr%p%F(2)]
                            LT%lstrint_W(cont_int(4)+1:cont_int(4)+4) = [cell(1),cell(2),ptr%p%n,ptr%p%grupo]
                            ! print*, "id",id,"transferindo para o oeste",  [cell(1),cell(2),ptr%p%n,ptr%p%grupo]
                            ! if (cell(2) < domx(1)) then
                            !     ! print*, 'affs'
                            !     ! read(*,*)
                            !     ! print*, "Removido", ptr%p%n, "de malha", i,j
                            !     print*, "Removido", ptr%p%n, "de malha", i,j
                            !     print*, "Irá para", LT%lstrint_W
        
                            cont_db(4) = cont_db(4) + 6
                            cont_int(4) = cont_int(4) + 4
                        end if 
                        if (cell(1) >= domy(2) .and. cell(1) < mesh(2)+2) then
                            ! print*, "L 567",  cell(1), cell(2),  "part", ptr%p%n, "i,j", i, j
                            ! print*, "id",id,"transferindo para o norte",  [cell(1),cell(2),ptr%p%n,ptr%p%grupo]
                            LT%lstrdb_N(cont_db(1)+1:cont_db(1)+6) = [x(1),x(2), ptr%p%v(1),ptr%p%v(2), ptr%p%F(1),ptr%p%F(2)]                
                            LT%lstrint_N(cont_int(1)+1:cont_int(1)+4) = [cell(1),cell(2),ptr%p%n,ptr%p%grupo]
                            ! if (cell(1) > domy(1)) then
                            !     ! read(*,*)
                            !     ! print*, 'affs'
                            !     ! call list_remove(previous_node)
                            !     print*, "Removido", ptr%p%n, "de malha", i,j
                            !     print*, "Irá para", LT%lstrint_N(cont_int(1)+1:cont_int(1)+4), "cont", cont_int(1)
                            ! !     node => list_next(previous_node)
                            ! ! else 
                            ! !     call list_change(previous_node,malha(cell(1),cell(2))%list)
                            ! !     node => list_next(previous_node)
                            ! end if 
                            cont_db(1) = cont_db(1) + 6
                            cont_int(1) = cont_int(1) + 4
                        end if 
                        if (cell(1) <= domy(1) .and. cell(1) > 1) then
                            ! print*, "L 580", cell(1), cell(2), "part", ptr%p%n, "i,j", i, j
                            LT%lstrdb_S(cont_db(2)+1:cont_db(2)+6) = [x(1),x(2), ptr%p%v(1),ptr%p%v(2), ptr%p%F(1),ptr%p%F(2)]
                            LT%lstrint_S(cont_int(2)+1:cont_int(2)+4) = [cell(1),cell(2),ptr%p%n,ptr%p%grupo]
                            ! ! print*, "L 583"
                            ! print*, "id",id,"transferindo para o sul",  [cell(1),cell(2),ptr%p%n,ptr%p%grupo]
                            ! if (cell(1) < domy(1)) then
                            !     ! read(*,*)
                            !     ! print*, 'affs'
                            !     ! call list_remove(previous_node)
                            !     print*, "Removido", ptr%p%n, "de malha", i,j
                            !     print*, "Irá para", LT%lstrint_S(cont_int(2)+1:cont_int(2)+4), "cont", cont_int(1)
                            ! !     node => list_next(previous_node)
                            ! ! else 
                            ! !     !!!! ! print*, "L 590"
                            ! !     call list_change(previous_node,malha(cell(1),cell(2))%list)
                            ! !     node => list_next(previous_node)
                            ! !     !!!! ! print*, "L 624"
                            ! end if 
                            cont_db(2) = cont_db(2) + 6
                            cont_int(2) = cont_int(2) + 4
                        end if

                        ! CASO PERIODICO

                        if (east == 'p' .and. cell(2) == mesh(1)+2 ) then
                            ! print*, "L 538",  cell(1), cell(2),"part", ptr%p%n,"i,j", i, j
                            ! 6 elementos
                            LT%lstrdb_W(cont_db(4)+1:cont_db(4)+6) = [x(1)- jcell(mesh(1)+1),x(2), &
                            ptr%p%v(1),ptr%p%v(2), ptr%p%F(1),ptr%p%F(2)]
                            ! 4 elementos
                            ! print*, "id",id,"transferindo para o leste",  [cell(1),cell(2),ptr%p%n,ptr%p%grupo]
                            LT%lstrint_W(cont_int(4)+1:cont_int(4)+4) = [cell(1), 2,ptr%p%n,ptr%p%grupo]

                            cont_db(4) = cont_db(4) + 6
                            cont_int(4) = cont_int(4) + 4
                        end if 
                        if (west == 'p' .and. cell(2) == 1) then
                            ! print*, "L 554",  cell(1), cell(2), "part", ptr%p%n, "i,j", i, j
                            ! print*, "domx", domx
                            LT%lstrdb_E(cont_db(3)+1:cont_db(3)+6) = [x(1)+ jcell(mesh(1)+1), x(2), ptr%p%v(1),ptr%p%v(2), ptr%p%F(1),ptr%p%F(2)]
                            LT%lstrint_E(cont_int(3)+1:cont_int(3)+4) = [cell(1),mesh(1)+1,ptr%p%n,ptr%p%grupo]

                            cont_db(3) = cont_db(3) + 6
                            cont_int(3) = cont_int(3) + 4
                        end if 
                        if (north == 'p' .and. cell(1) == mesh(2)+2) then
                            ! print*, "L 567",  cell(1), cell(2),  "part", ptr%p%n, "i,j", i, j
                            ! print*, "id",id,"transferindo para o norte",  [cell(1),cell(2),ptr%p%n,ptr%p%grupo]
                            LT%lstrdb_S(cont_db(2)+1:cont_db(2)+6) = [x(1),x(2) - icell(mesh(2)+1), ptr%p%v(1),ptr%p%v(2), ptr%p%F(1),ptr%p%F(2)]                
                            LT%lstrint_S(cont_int(2)+1:cont_int(2)+4) = [2,cell(2),ptr%p%n,ptr%p%grupo]

                            cont_db(2) = cont_db(2) + 6
                            cont_int(2) = cont_int(2) + 4
                        end if 
                        if (south == 'p' .and. cell(1) == 1) then
                            ! print*, "L 580", cell(1), cell(2), "part", ptr%p%n, "i,j", i, j
                            LT%lstrdb_N(cont_db(1)+1:cont_db(1)+6) = [x(1),icell(mesh(2)+1) + x(2), ptr%p%v(1),ptr%p%v(2), ptr%p%F(1),ptr%p%F(2)]
                            LT%lstrint_N(cont_int(1)+1:cont_int(1)+4) = [mesh(2)+1,cell(2),ptr%p%n,ptr%p%grupo]
 
                            cont_db(1) = cont_db(1) + 6
                            cont_int(1) = cont_int(1) + 4
                        end if                        

                        !!! FIM DA VERIFICAÇÃO SE MUDOU DE DOMÍNIO !!!

                        ! Antes aqui tinha um else para o caso de não mudar de processo. 
                        ! Removi todos os call list_remove e deixei que a partícula fosse para
                        ! a área emprestada que será limpa na próxima execução de call comp_x. Isto pois a célula precisa
                        ! sentir a presença da partícula que ela perdeu mas ficou ao lado 
                        ! print*, "mudou2"
                        ! if (id == 0) read(*,*)

                        call list_change(previous_node,malha(cell(1),cell(2))%list)
                        node => list_next(previous_node)
 
                    else !Não mudou de celula
                        ! Mas tem que avisar quais ainda estão na
                        ! vizinhança
                        
                        ! print*, "NAO mudou cell",cell, "domx", domx, "domy", domy
                        if (cell(2) == domx(2)) then
                            !!! ! print*, "L 538",  cell(1), cell(2), domy(1), domy(2), "part", ptr%p%n
                            ! 6 elementos
                            LT%lstrdb_E(cont_db(3)+1:cont_db(3)+6) = [x(1),x(2), &
                                ptr%p%v(1),ptr%p%v(2), ptr%p%F(1),ptr%p%F(2)]
                            ! 4 elementos
                            LT%lstrint_E(cont_int(3)+1:cont_int(3)+4) = [cell(1),cell(2),ptr%p%n,ptr%p%grupo]
                            ! print*, "L id",id,"transferindo para o leste",  [cell(1),cell(2),ptr%p%n,ptr%p%grupo]
                            cont_db(3) = cont_db(3) + 6
                            cont_int(3) = cont_int(3) + 4

                        else if (cell(2) == domx(1)) then
                            !!! ! print*, "L 554",  cell(1), cell(2), domy(1), domy(2), "part", ptr%p%n
                            LT%lstrdb_W(cont_db(4)+1:cont_db(4)+6) = [x(1),x(2), ptr%p%v(1),ptr%p%v(2), ptr%p%F(1),ptr%p%F(2)]
                            LT%lstrint_W(cont_int(4)+1:cont_int(4)+4) = [cell(1),cell(2),ptr%p%n,ptr%p%grupo]
                            ! print*, "L id",id,"transferindo para o oeste",  [cell(1),cell(2),ptr%p%n,ptr%p%grupo]
                            cont_db(4) = cont_db(4) + 6
                            cont_int(4) = cont_int(4) + 4
                        end if
                        
                        if (cell(1) == domy(1)) then
                            !!! ! print*, "L 567",  cell(1), cell(2), domy(1), domy(2), "part", ptr%p%n
                            LT%lstrdb_S(cont_db(2)+1:cont_db(2)+6) = [x(1),x(2), ptr%p%v(1),ptr%p%v(2), ptr%p%F(1),ptr%p%F(2)]                
                            LT%lstrint_S(cont_int(2)+1:cont_int(2)+4) = [cell(1),cell(2),ptr%p%n,ptr%p%grupo]
                            ! print*, "L id",id,"transferindo para o sul",  [cell(1),cell(2),ptr%p%n,ptr%p%grupo]
                            cont_db(2) = cont_db(2) + 6
                            cont_int(2) = cont_int(2) + 4

                        elseif (cell(1) == domy(2)) then
                            !!! ! print*, "L 580", cell(1), cell(2), domy(1), domy(2), "part", ptr%p%n
                            LT%lstrdb_N(cont_db(1)+1:cont_db(1)+6) = [x(1),x(2), ptr%p%v(1),ptr%p%v(2), ptr%p%F(1),ptr%p%F(2)]
                            LT%lstrint_N(cont_int(1)+1:cont_int(1)+4) = [cell(1),cell(2),ptr%p%n,ptr%p%grupo]
                        !    print*, "L id",id,"transferindo para o norte",  [cell(1),cell(2),ptr%p%n,ptr%p%grupo]
                            cont_db(1) = cont_db(1) + 6
                            cont_int(1) = cont_int(1) + 4
                        end if

                        ! CASO PERIODICO 

                        if (east == 'p' .and. cell(2) == mesh(1)+1) then
                            !!! ! print*, "L 538",  cell(1), cell(2), domy(1), domy(2), "part", ptr%p%n
                            ! 6 elementos
                            LT%lstrdb_W(cont_db(4)+1:cont_db(4)+6) = [x(1)- jcell(mesh(1)+1),x(2),x(2), &
                                ptr%p%v(1),ptr%p%v(2), ptr%p%F(1),ptr%p%F(2)]
                            ! 4 elementos
                            LT%lstrint_W(cont_int(4)+1:cont_int(4)+4) = [cell(1),1,ptr%p%n,ptr%p%grupo]
                            ! print*, "L id",id,"transferindo para o oeste",  [cell(1),cell(2),ptr%p%n,ptr%p%grupo]
                            cont_db(4) = cont_db(4) + 6
                            cont_int(4) = cont_int(4) + 4

                        else if (west == 'p' .and. cell(2) == 2) then
                            !!! ! print*, "L 554",  cell(1), cell(2), domy(1), domy(2), "part", ptr%p%n
                            LT%lstrdb_E(cont_db(3)+1:cont_db(3)+6) = [x(1)+ jcell(mesh(1)+1),x(2), ptr%p%v(1),ptr%p%v(2), ptr%p%F(1),ptr%p%F(2)]
                            LT%lstrint_E(cont_int(3)+1:cont_int(3)+4) = [cell(1),mesh(1)+2,ptr%p%n,ptr%p%grupo]
                            ! print*, "L id",id,"transferindo para o leste",  [cell(1),cell(2),ptr%p%n,ptr%p%grupo]
                            cont_db(3) = cont_db(3) + 6
                            cont_int(3) = cont_int(3) + 4
                        end if
                        
                        if (south == 'p' .and. cell(1) == mesh(2)+1) then
                            !!! ! print*, "L 567",  cell(1), cell(2), domy(1), domy(2), "part", ptr%p%n
                            LT%lstrdb_N(cont_db(2)+1:cont_db(2)+6) = [x(1),icell(mesh(2)+1) + x(2), ptr%p%v(1),ptr%p%v(2), ptr%p%F(1),ptr%p%F(2)]                
                            LT%lstrint_N(cont_int(2)+1:cont_int(2)+4) = [1,cell(2),ptr%p%n,ptr%p%grupo]
                            ! print*, "L id",id,"transferindo para o norte",  [cell(1),cell(2),ptr%p%n,ptr%p%grupo]
                            cont_db(1) = cont_db(1) + 6
                            cont_int(1) = cont_int(1) + 4

                        elseif (north == 'p' .and. cell(1) == 2) then
                            !!! ! print*, "L 580", cell(1), cell(2), domy(1), domy(2), "part", ptr%p%n
                            LT%lstrdb_S(cont_db(1)+1:cont_db(1)+6) = [x(1),x(2) - icell(mesh(2)+1), ptr%p%v(1),ptr%p%v(2), ptr%p%F(1),ptr%p%F(2)]
                            LT%lstrint_S(cont_int(1)+1:cont_int(1)+4) = [mesh(2)+2,cell(2),ptr%p%n,ptr%p%grupo]
                        !    print*, "L id",id,"transferindo para o sul",  [cell(1),cell(2),ptr%p%n,ptr%p%grupo]
                            cont_db(2) = cont_db(2) + 6
                            cont_int(2) = cont_int(2) + 4
                        end if
                       
                        previous_node => node
                        node => list_next(node)
                        ! print*, "L740 FORÇA", ptr%p%F, id
                    end if
                    if (associated(node)) ptr = transfer(list_get(node), ptr)
                end do      
            end do
        end do
        
        ! transferir os termos que estão no encontro de 4 processos 
        ! print*, "L 794", id
        ! if (id == 0) read(*,*) 
        ! call MPI_barrier(MPI_COMM_WORLD, ierr)
        if (np > 1 .and. sum(ids(5:8)) > -4) then
            do i = 5,8
                !identifica qual o processo na diagonal
                if (ids(i) > -1) then
                    destD = i 
                end if
            end do
            ! ! print*, "L 787", ids, "D", destD
            ! Aqui está com suporte a 4 processadores por enquanto
            ! Primeiro transfere para as celulas emprestadas
            
            if (domy(1) /= 1) i = domy(1)
            if (domy(2) /= mesh(2)+2) i = domy(2)
            if (domx(1) /= 1) j = domx(1)
            if (domx(2) /= mesh(1) + 2) j = domx(2)

            previous_node => malha(i,j)%list
            node => list_next(malha(i,j)%list)

            ! if (associated(node)) ptr = transfer(list_get(node), ptr) !para que isto?
            ! if (associated(node))  print*, "transferencia diagonal", i,j, "ID", id
            do while(associated(node))
                ptr = transfer(list_get(node), ptr)
                x = ptr%p%x ! não sei pq essas linhas não estavam aqui antes. Atentar
                cell = [i, j]
                ! ! print*, "L 806", cell(1), cell(2), domy(1), domy(2), "part", ptr%p%n
                LT%lstrdb_D(cont_db(destD)+1:cont_db(destD)+6) = [x(1),x(2), ptr%p%v(1),ptr%p%v(2), ptr%p%F(1),ptr%p%F(2)]
                LT%lstrint_D(cont_int(destD)+1:cont_int(destD)+4) = [cell(1),cell(2),ptr%p%n,ptr%p%grupo]
                ! print*, "L id",id,"transferindo diagonal",  [cell(1),cell(2),ptr%p%n], "p/ id", ids(destD)
                cont_db(destD) = cont_db(destD) + 6
                cont_int(destD) = cont_int(destD) + 4 
                previous_node => node
                node => list_next(node)
            end do
            
            ! Segundo transfere para as celulas do domínio caso a partícula tenha mudado de processo
            if (domy(1) /= 1) i = domy(1) -1
            if (domy(2) /= mesh(2)+2) i = domy(2) +1
            if (domx(1) /= 1) j = domx(1) -1
            if (domx(2) /= mesh(1) + 2) j = domx(2) +1 
        
            previous_node => malha(i,j)%list
            node => list_next(malha(i,j)%list)
            ! if (associated(node)) ptr = transfer(list_get(node), ptr)
            ! if (associated(node))  print*, "mudança diagonal", i,j, "ID", id
            do while(associated(node))
                ptr = transfer(list_get(node), ptr)
                x = ptr%p%x
                cell = [i, j]
                ! ! print*, "L 806", cell(1), cell(2), domy(1), domy(2), "part", ptr%p%n
                LT%lstrdb_D(cont_db(destD)+1:cont_db(destD)+6) = [x(1),x(2), ptr%p%v(1),ptr%p%v(2), ptr%p%F(1),ptr%p%F(2)]
                LT%lstrint_D(cont_int(destD)+1:cont_int(destD)+4) = [cell(1),cell(2),ptr%p%n,ptr%p%grupo]
                ! print*, "L id",id,"mudando diagonal",  [cell(1),cell(2),ptr%p%n], "p/ id", ids(destD)
                cont_db(destD) = cont_db(destD) + 6
                cont_int(destD) = cont_int(destD) + 4 
                previous_node => node
                node => list_next(node)
            end do

            !Caso Periódico. Só será necessário se os domínios dos processos formarem região 2x2 
            if (wall == 'pppp') then
                ! Primeiro transfere para as celulas emprestadas (periodico)
                ! O x aqui tem outra função
                !extremidades arestas norte e sul
                if (domy(1) == 1) then
                    i = domy(1) + 1
                    x = [0.0_dp, icell(mesh(2)+1)]
                    cell = [mesh(2)+2, 1]
                end if
                if (domy(2) == mesh(2)+2) then 
                    i = domy(2) - 1
                    x = [0.0_dp, -icell(mesh(2)+1)]
                    cell = [1, 1]
                end if
                if (domx(1) /= 1) j = domx(1)  
                if (domx(2) /= mesh(1) + 2) j = domx(2)
                cell(2) = j

                previous_node => malha(i,j)%list
                node => list_next(malha(i,j)%list)
    
                ! if (associated(node)) ptr = transfer(list_get(node), ptr)
                ! if (associated(node))  print*, "transferencia diagonal", i,j, "ID", id
                do while(associated(node))
                    ptr = transfer(list_get(node), ptr)
                    LT%lstrdb_D(cont_db(destD)+1:cont_db(destD)+6) = [x(1)+ptr%p%x(1),x(2)+ptr%p%x(2), ptr%p%v(1),ptr%p%v(2), & 
                    ptr%p%F(1),ptr%p%F(2)]
                    LT%lstrint_D(cont_int(destD)+1:cont_int(destD)+4) = [cell(1),cell(2),ptr%p%n,ptr%p%grupo]
                    ! print*, "L id",id,"transferindo diagonal",  [cell(1),cell(2),ptr%p%n], "p/ id", ids(destD)
                    cont_db(destD) = cont_db(destD) + 6
                    cont_int(destD) = cont_int(destD) + 4 
                    previous_node => node
                    node => list_next(node)
                end do

                !extremidades leste oeste
                if (domy(1) /= 1) i = domy(1)
                if (domy(2) /= mesh(2)+2) i = domy(2)
                if (domx(1) == 1) then 
                    j = 2
                    x = [jcell(mesh(1)+2), 0.0_dp]
                    cell = [i, mesh(1) + 2]
                end if
                if (domx(2) == mesh(1) + 2) then 
                    j = domx(2) - 1
                    x = [-jcell(mesh(1)+1), 0.0_dp]
                    cell = [i, 1]
                end if

                previous_node => malha(i,j)%list
                node => list_next(malha(i,j)%list)
    
                ! if (associated(node)) ptr = transfer(list_get(node), ptr)
                ! if (associated(node))  print*, "transferencia diagonal", i,j, "ID", id
                do while(associated(node))
                    ptr = transfer(list_get(node), ptr)
                    x = ptr%p%x ! não sei pq essas linhas não estavam aqui antes. Atentar
                    ! ! print*, "L 806", cell(1), cell(2), domy(1), domy(2), "part", ptr%p%n
                    LT%lstrdb_D(cont_db(destD)+1:cont_db(destD)+6) = [x(1)+ptr%p%x(1),x(2)+ptr%p%x(2), ptr%p%v(1),ptr%p%v(2), &
                        ptr%p%F(1),ptr%p%F(2)]
                    LT%lstrint_D(cont_int(destD)+1:cont_int(destD)+4) = [cell(1),cell(2),ptr%p%n,ptr%p%grupo]
                    ! print*, "L id",id,"transferindo diagonal",  [cell(1),cell(2),ptr%p%n], "p/ id", ids(destD)
                    cont_db(destD) = cont_db(destD) + 6
                    cont_int(destD) = cont_int(destD) + 4 
                    previous_node => node
                    node => list_next(node)
                end do
                
                !extremidades nas pontas
                if (domy(1) == 1) then
                    i = domy(1) +1
                    cell(1) = mesh(2) + 2
                    x(2)  =  icell(mesh(2)+1)
                end if
                if (domy(2) == mesh(2)+2) then
                    i = domy(2) - 1
                    cell(1) = 1
                    x(2)  =  - icell(mesh(2)+1)
                end if
                if (domx(1) == 1) then 
                    j = domx(1)+1
                    cell(2) = mesh(1) + 2
                    x(1) = jcell(mesh(1)+1)
                end if
                if (domx(2) == mesh(1) + 2) then 
                    j = domx(2)-1
                    cell(2) = 1
                    x(1) = - jcell(mesh(1)+2)
                end if
                previous_node => malha(i,j)%list
                node => list_next(malha(i,j)%list)
    
                ! if (associated(node)) ptr = transfer(list_get(node), ptr)
                ! if (associated(node))  print*, "transferencia diagonal", i,j, "ID", id
                do while(associated(node))
                    ptr = transfer(list_get(node), ptr)
                    x = ptr%p%x ! não sei pq essas linhas não estavam aqui antes. Atentar
                    cell = [i, j]
                    ! ! print*, "L 806", cell(1), cell(2), domy(1), domy(2), "part", ptr%p%n
                    LT%lstrdb_D(cont_db(destD)+1:cont_db(destD)+6) = [x(1)+ptr%p%x(1),x(2)+ptr%p%x(2), ptr%p%v(1),ptr%p%v(2), &
                         ptr%p%F(1),ptr%p%F(2)]
                    LT%lstrint_D(cont_int(destD)+1:cont_int(destD)+4) = [cell(1),cell(2),ptr%p%n,ptr%p%grupo]
                    ! print*, "L id",id,"transferindo diagonal",  [cell(1),cell(2),ptr%p%n], "p/ id", ids(destD)
                    cont_db(destD) = cont_db(destD) + 6
                    cont_int(destD) = cont_int(destD) + 4 
                    previous_node => node
                    node => list_next(node)
                end do

                ! Segundo transfere para as celulas do domínio caso a partícula tenha mudado de processo
                ! (periodico)

                ! norte e sul
 
                if (domy(1) == 1) then
                    i = domy(1)
                    x = [0.0_dp, icell(mesh(2)+1)]
                    cell = [mesh(2)+1, 1]
                end if
                if (domy(2) == mesh(2)+2) then     
                    i = domy(2)
                    x = [0.0_dp, -icell(mesh(2)+1)]
                    cell = [1, 0]
                end if
                if (domx(1) /= 1) j = domx(1)   
                if (domx(2) /= mesh(1) + 2) j = domx(2)  
                cell(2) = j  
            
                previous_node => malha(i,j)%list
                node => list_next(malha(i,j)%list)
                ! if (associated(node)) ptr = transfer(list_get(node), ptr)
                ! if (associated(node))  print*, "mudança diagonal", i,j, "ID", id
                do while(associated(node))
                    ptr = transfer(list_get(node), ptr)
                    x = ptr%p%x ! não sei pq essas linhas não estavam aqui antes. Atentar
                    ! ! print*, "L 806", cell(1), cell(2), domy(1), domy(2), "part", ptr%p%n
                    LT%lstrdb_D(cont_db(destD)+1:cont_db(destD)+6) = [x(1)+ptr%p%x(1),x(2)+ptr%p%x(2), ptr%p%v(1),ptr%p%v(2), &
                         ptr%p%F(1),ptr%p%F(2)]
                    LT%lstrint_D(cont_int(destD)+1:cont_int(destD)+4) = [cell(1),cell(2),ptr%p%n,ptr%p%grupo]
                    ! print*, "L id",id,"mudando diagonal",  [cell(1),cell(2),ptr%p%n], "p/ id", ids(destD)
                    cont_db(destD) = cont_db(destD) + 6
                    cont_int(destD) = cont_int(destD) + 4 
                    previous_node => node
                    node => list_next(node)
                end do

                !extremidades leste oeste
                if (domy(1) /= 1) i = domy(1)
                if (domy(2) /= mesh(2)+2) i = domy(2)
                if (domx(1) == 1) then 
                    j = 1
                    x = [jcell(mesh(1)+1), 0.0_dp]
                    cell = [i, mesh(1) + 2]
                end if
                if (domx(2) == mesh(1) + 2) then 
                    j = domx(2) 
                    x = [-jcell(mesh(1)+1), 0.0_dp]
                    cell = [i, 2]
                end if
            
                previous_node => malha(i,j)%list
                node => list_next(malha(i,j)%list)
                ! if (associated(node)) ptr = transfer(list_get(node), ptr)
                ! if (associated(node))  print*, "mudança diagonal", i,j, "ID", id
                do while(associated(node))
                    ptr = transfer(list_get(node), ptr)
                    x = ptr%p%x ! não sei pq essas linhas não estavam aqui antes. Atentar
                    cell = [i, j]
                    ! ! print*, "L 806", cell(1), cell(2), domy(1), domy(2), "part", ptr%p%n
                    LT%lstrdb_D(cont_db(destD)+1:cont_db(destD)+6) = [x(1)+ptr%p%x(1),x(2)+ptr%p%x(2), ptr%p%v(1),ptr%p%v(2), &
                        ptr%p%F(1),ptr%p%F(2)]
                    LT%lstrint_D(cont_int(destD)+1:cont_int(destD)+4) = [cell(1),cell(2),ptr%p%n,ptr%p%grupo]
                    ! print*, "L id",id,"mudando diagonal",  [cell(1),cell(2),ptr%p%n], "p/ id", ids(destD)
                    cont_db(destD) = cont_db(destD) + 6
                    cont_int(destD) = cont_int(destD) + 4 
                    previous_node => node
                    node => list_next(node)
                end do

                !extremidades nas pontas
                if (domy(1) == 1) then
                    i = domy(1)
                    cell(1) = mesh(2) + 2
                    x(2)  =  icell(mesh(2)+1)
                end if 
                if (domy(2) == mesh(2)+2) i = domy(2)
                    i = domy(2)
                    cell(1) = 1
                    x(2)  =  - icell(mesh(2)+1)
                end if
                if (domx(1) == 1) then 
                    j = domx(1)
                    cell(2) = mesh(1) + 2
                    x(1) = jcell(mesh(1)+1)
                end if
                if (domx(2) == mesh(1) + 2) then 
                    j = domx(2)
                    cell(2) = 1
                    x(1) = - jcell(mesh(1)+2)
                end if
            
                previous_node => malha(i,j)%list
                node => list_next(malha(i,j)%list)
                ! if (associated(node)) ptr = transfer(list_get(node), ptr)
                ! if (associated(node))  print*, "mudança diagonal", i,j, "ID", id
                do while(associated(node))
                    ptr = transfer(list_get(node), ptr)
                    x = ptr%p%x ! não sei pq essas linhas não estavam aqui antes. Atentar
                    cell = [i, j]
                    ! ! print*, "L 806", cell(1), cell(2), domy(1), domy(2), "part", ptr%p%n
                    LT%lstrdb_D(cont_db(destD)+1:cont_db(destD)+6) = [x(1)+ptr%p%x(1),x(2)+ptr%p%x(2), ptr%p%v(1),ptr%p%v(2), & 
                        ptr%p%F(1),ptr%p%F(2)]
                    LT%lstrint_D(cont_int(destD)+1:cont_int(destD)+4) = [cell(1),cell(2),ptr%p%n,ptr%p%grupo]
                    ! print*, "L id",id,"mudando diagonal",  [cell(1),cell(2),ptr%p%n], "p/ id", ids(destD)
                    cont_db(destD) = cont_db(destD) + 6
                    cont_int(destD) = cont_int(destD) + 4 
                    previous_node => node
                    node => list_next(node)
                end do

            ! Fim da transfência para diagonal
        end if
        
        ! print*, "L 854", id
        ! if (id == 0) read(*,*)
        ! call MPI_barrier(MPI_COMM_WORLD, ierr)
        ! print*, "L SOMA", (ids(1) + ids(2) + ids(3) +ids(4)), "id", id  
        ! print*, "enviando"
        if (np > 1) then !se for paralelo 
            ! print*, 'L 862 >> id', id
            ! call MPI_SEND(BUF, COUNT, DATATYPE, DEST, TAG, COMM, IERROR)
            tag = 1
            if (ids(1) >= 0) call MPI_SEND(LT%lstrdb_N, cont_db(1), MPI_DOUBLE_PRECISION, ids(1), tag,MPI_COMM_WORLD, ierr)
            if (ids(2) >= 0) call MPI_SEND(LT%lstrdb_S, cont_db(2), MPI_DOUBLE_PRECISION, ids(2), tag,MPI_COMM_WORLD, ierr)   
            if (ids(3) >= 0) call MPI_SEND(LT%lstrdb_E, cont_db(3), MPI_DOUBLE_PRECISION, ids(3), tag,MPI_COMM_WORLD, ierr)   
            if (ids(4) >= 0) call MPI_SEND(LT%lstrdb_W, cont_db(4), MPI_DOUBLE_PRECISION, ids(4), tag,MPI_COMM_WORLD, ierr)   
            if (sum(ids(5:8)) > -4) then 
                call MPI_SEND(LT%lstrdb_D, cont_db(destD), MPI_DOUBLE_PRECISION, ids(destD), tag,MPI_COMM_WORLD, ierr)
            end if 
            tag = 2 
            if (ids(1) >= 0) call MPI_SEND(LT%lstrint_N, cont_int(1), MPI_integer, ids(1), tag,MPI_COMM_WORLD, ierr)
            if (ids(2) >= 0) call MPI_SEND(LT%lstrint_S, cont_int(2), MPI_integer, ids(2), tag,MPI_COMM_WORLD, ierr)
            if (ids(3) >= 0) call MPI_SEND(LT%lstrint_E, cont_int(3), MPI_integer, ids(3), tag,MPI_COMM_WORLD, ierr)
            if (ids(4) >= 0) call MPI_SEND(LT%lstrint_W, cont_int(4), MPI_integer, ids(4), tag,MPI_COMM_WORLD, ierr)
            if (sum(ids(5:8)) > -4) then 
                call MPI_SEND(LT%lstrint_D, cont_int(destD), MPI_integer, ids(destD), tag,MPI_COMM_WORLD, ierr)
            end if

            cont_db = 0
            cont_int = 0
            ! print*, "L 935 TUDO ENVIADO", id
            tag = 1
            if (ids(1) >= 0) then  
                call MPI_probe(ids(1),tag, MPI_COMM_WORLD, status, ierr)
                CALL MPI_Get_count ( status , MPI_DOUBLE_PRECISION , cont_db(1) , ierr )
                call MPI_RECV (LT%lstrdb_N, cont_db(1), MPI_DOUBLE_PRECISION,ids(1), tag, MPI_COMM_WORLD, status, ierr)
            end if 
            
            if (ids(2) >= 0) then 
                call MPI_probe(ids(2),tag, MPI_COMM_WORLD, status, ierr)
                CALL MPI_Get_count ( status , MPI_DOUBLE_PRECISION ,  cont_db(2) , ierr )
                call MPI_RECV (LT%lstrdb_S, cont_db(2), MPI_DOUBLE_PRECISION,ids(2), tag, MPI_COMM_WORLD, status, ierr)
            
            end if 
            if (ids(3) >= 0) then
                call MPI_probe(ids(3),tag, MPI_COMM_WORLD, status, ierr)
                CALL MPI_Get_count ( status , MPI_DOUBLE_PRECISION ,  cont_db(3) , ierr )
                call MPI_RECV (LT%lstrdb_E, cont_db(3), MPI_DOUBLE_PRECISION,ids(3), tag, MPI_COMM_WORLD, status, ierr)
            end if 

            if (ids(4) >= 0) then
                call MPI_probe(ids(4),tag, MPI_COMM_WORLD, status, ierr)
                CALL MPI_Get_count ( status , MPI_DOUBLE_PRECISION ,  cont_db(4) , ierr )
                call MPI_RECV (LT%lstrdb_W, cont_db(4), MPI_DOUBLE_PRECISION,ids(4), tag, MPI_COMM_WORLD, status, ierr)
            end if 

            tag = 2
            if (ids(1) >= 0) then
                call MPI_probe(ids(1),tag, MPI_COMM_WORLD, status, ierr)
                CALL MPI_Get_count ( status , MPI_integer ,  cont_int(1) , ierr )
                call MPI_RECV (LT%lstrint_N, cont_int(1), MPI_integer,ids(1), tag, MPI_COMM_WORLD, status, ierr)
            end if 
            if (ids(2) >= 0) then
                call MPI_probe(ids(2),tag, MPI_COMM_WORLD, status, ierr)
                CALL MPI_Get_count ( status , MPI_integer ,  cont_int(2) , ierr )
                call MPI_RECV (LT%lstrint_S, cont_int(2), MPI_integer,ids(2), tag, MPI_COMM_WORLD, status, ierr)
            end if 
            if (ids(3) >= 0) then
                call MPI_probe(ids(3),tag, MPI_COMM_WORLD, status, ierr)
                CALL MPI_Get_count ( status , MPI_integer ,  cont_int(3) , ierr )
                call MPI_RECV (LT%lstrint_E, cont_int(3), MPI_integer,ids(3), tag, MPI_COMM_WORLD, status, ierr)
            end if
            if (ids(4) >= 0) then
                call MPI_probe(ids(4),tag, MPI_COMM_WORLD, status, ierr)
                CALL MPI_Get_count ( status , MPI_integer ,  cont_int(4) , ierr )
                call MPI_RECV (LT%lstrint_W, cont_int(4), MPI_integer,ids(4), tag, MPI_COMM_WORLD, status, ierr)
            end if
            if (sum(ids(5:8)) > -4) then
                tag = 1
                ! RECEBE DA DIAGONAL (suporte para 4 processos)
                ! call MPI_probe(,tag, MPI_COMM_WORLD, status, ierr)
                call MPI_probe(ids(destD),tag, MPI_COMM_WORLD, status, ierr)
                CALL MPI_Get_count ( status , MPI_DOUBLE_PRECISION ,  cont_db(destD) , ierr )
                call MPI_RECV (LT%lstrdb_D, cont_db(destD), MPI_DOUBLE_PRECISION,ids(destD), tag, MPI_COMM_WORLD, status, ierr)
                tag = 2
                ! Recebe INT da DIAGONAL. Aqui está com suporte para 4 processadores
                call MPI_probe(ids(destD),tag, MPI_COMM_WORLD, status, ierr)
                CALL MPI_Get_count ( status , MPI_DOUBLE_PRECISION ,  cont_int(destD) , ierr )
                call MPI_RECV (LT%lstrint_D, cont_int(destD), MPI_DOUBLE_PRECISION,ids(destD), tag, MPI_COMM_WORLD, status, ierr)

                ! DIAGONAL, tem que alterar para o caso com >4 processos
                j = destD
                do i = 0, cont_db(j)/6
                    if (cont_db(j) > 0) then
                        ! print*, "LA 1108"
                        ! print*, "D al", id,  LT%lstrdb_D(i*6+1:i*6+2)
                        allocate(ptr%p)
                        ptr%p%x = LT%lstrdb_D(i*6+1:i*6+2)
                        ! print*, ptr%p%x
                        ptr%p%v = LT%lstrdb_D(i*6+3:i*6+4)
                        ptr%p%grupo = LT%lstrint_D(i*4+4) 
                        ptr%p%F = LT%lstrdb_D(i*6+5:i*6+6)
                        ! ! print*, "L FORÇÆ", ptr%p%F
                        ptr%p%n = LT%lstrint_D(i*4+3) !identidade da partícula importante para imprimir
                        call list_insert(malha(LT%lstrint_D(i*4+1), &
                            LT%lstrint_D(i*4+2))%list, data=transfer(ptr, list_data)) 
                        ! print*, ptr%p%n, "D allocado em", LT%lstrint_D(i*4+1), LT%lstrint_D(i*4+2), "id =", id, "posição",ptr%p%x
                        ! print*, "W allocado F = ",ptr%p%F
                        cont_db(j) = cont_db(j) - 6
                    end if
                end do
            end if

            ! print*, "cont db", cont_db

            ! if (id == 0)  read(*,*)
            ! call MPI_barrier(MPI_COMM_WORLD, ierr)
            ! print*, "Tudo recebido", id
            !  call MPI_barrier(MPI_COMM_WORLD, ierr)
            !  print*, 'ID APOS BARREIRA', id
            ! print*, "L 1033"
            do j = 1,4
                do i = 0, cont_db(j)/6
                    !! ! print*, "L 840" , i, j, "id", id, "cont_db/6", cont_db(j)/6
                    if (j == 1 .and. cont_db(1) > 0 ) then
                        ! print*, "LA 1133" 
                        allocate(ptr%p)
                        
                        ptr%p%x = LT%lstrdb_N(i*6+1:i*6+2)
                        ptr%p%v = LT%lstrdb_N(i*6+3:i*6+4)
                        ptr%p%grupo = LT%lstrint_N(i*4+4) 
                        ptr%p%F = LT%lstrdb_N(i*6+5:i*6+6)
                        ! ! print*, "L FORÇÆ", ptr%p%F
                        ptr%p%n = LT%lstrint_N(i*4+3) !identidade da partícula importante para imprimir
                        ! !print*, 'LLLL', LT%lstrint_N(i*4+1), LT%lstrint_N(i*4+2)
                        call list_insert(malha(LT%lstrint_N(i*4+1), &
                            LT%lstrint_N(i*4+2))%list, data=transfer(ptr, list_data)) 
                        ! print*, ptr%p%n, "N allocado em", LT%lstrint_N(i*4+1), LT%lstrint_N(i*4+2), "id =", id
                        ! print*, "N allocado F = ",ptr%p%F
                        cont_db(j) = cont_db(j) - 6
                        
                    end if
                    if (j == 2  .and. cont_db(2) > 0) then
                        ! print*, "LA 1151", id
                        allocate(ptr%p)
                        ! print*, "L 859 <<<", id
                        ptr%p%x = LT%lstrdb_S(i*6+1:i*6+2)
                        ptr%p%v = LT%lstrdb_S(i*6+3:i*6+4)
                        ptr%p%grupo = LT%lstrint_S(i*4+4) 
                        ptr%p%F = LT%lstrdb_S(i*6+5:i*6+6)
                        ! ! print*, "L FORÇÆ", ptr%p%F
                        ptr%p%n = LT%lstrint_S(i*4+3) !identidade da partícula importante para imprimir
                        ! print*, "L 865 <<<", LT%lstrint_S(i*6+1:i*6+6)
                        call list_insert(malha(LT%lstrint_S(i*4+1), &
                            LT%lstrint_S(i*4+2))%list, data=transfer(ptr, list_data)) 
                        ! print*, ptr%p%n, "S allocado em", LT%lstrint_S(i*4+1), LT%lstrint_S(i*4+2), "id =", id
                        ! print*, "S allocado F = ",ptr%p%F
                        !! ! print*, "L 869 <<<"
                        cont_db(j) = cont_db(j) - 6
                    end if
                    if (j == 3  .and. cont_db(3) > 0) then
                        ! print*, "LA 1170"
                        allocate(ptr%p)
                        ptr%p%x = LT%lstrdb_E(i*6+1:i*6+2)
                        ptr%p%v = LT%lstrdb_E(i*6+3:i*6+4)
                        ptr%p%grupo = LT%lstrint_E(i*4+4) 
                        ptr%p%F = LT%lstrdb_E(i*6+5:i*6+6)
                        ! ! print*, "L FORÇÆ", ptr%p%F
                        ptr%p%n = LT%lstrint_E(i*4+3) !identidade da partícula importante para imprimir
                        call list_insert(malha(LT%lstrint_E(i*4+1), &
                            LT%lstrint_E(i*4+2))%list, data=transfer(ptr, list_data)) 
                        ! print*, ptr%p%n, "E allocado em", LT%lstrint_E(i*4+1), LT%lstrint_E(i*4+2), "x=",ptr%p%x
                        ! print*, "E allocado", LT%lstrint_E
                        cont_db(j) = cont_db(j) - 6
                    end if
                    if (j == 4 .and. cont_db(4) > 0 ) then
                        ! print*, "LA 1185"
                        allocate(ptr%p)
                        ptr%p%x = LT%lstrdb_W(i*6+1:i*6+2)
                        ! print*, ptr%p%x
                        ptr%p%v = LT%lstrdb_W(i*6+3:i*6+4)
                        ptr%p%grupo = LT%lstrint_W(i*4+4) 
                        ptr%p%F = LT%lstrdb_W(i*6+5:i*6+6)
                        ! ! print*, "L FORÇÆ", ptr%p%F
                        ptr%p%n = LT%lstrint_W(i*4+3) !identidade da partícula importante para imprimir
                        call list_insert(malha(LT%lstrint_W(i*4+1), &
                            LT%lstrint_W(i*4+2))%list, data=transfer(ptr, list_data)) 
                        ! print*, ptr%p%n, "W allocado em", LT%lstrint_W(i*4+1), LT%lstrint_W(i*4+2), "id =", id
                        ! print*, "W allocado F = ",ptr%p%F
                        cont_db(j) = cont_db(j) - 6
                    end if  
                end do
            end do   
            ! print*, "Tudo Alocado"
              ! if (id == 0) read(*,*) 
            ! call MPI_barrier(MPI_COMM_WORLD, ierr)
        end if

        ! print*, "L 1115 fim comp_x", id
        ! if (id == 0)  read(*,*)
        ! call MPI_barrier(MPI_COMM_WORLD, ierr)
    end subroutine comp_x
    
    ! atualiza velocidades
    subroutine comp_v(malha,mesh,dt,t,propriedade,domx,domy)
        use linkedlist
        use mod1
        use data
        use mod0

        type(prop_grupo), allocatable,dimension(:),intent(in) :: propriedade
        type(container), allocatable,dimension(:,:),intent(in) :: malha
        integer :: i,j, bsh(4)
        integer, intent(in) :: mesh(2),domx(2),domy(2)
        type(list_t), pointer :: node 
        real(dp), intent(in) :: dt, t
        type(data_ptr) :: ptr

        bsh = [0,0,0,0]
        if (domx(1) == 1) bsh(1) = 1
        if (domx(2) == mesh(1)+2) bsh(2) = -1
        if (domy(1) == 1) bsh(3) = 1
        if (domy(2) == mesh(2)+2) bsh(4) = -1

        do i = (domy(1)+bsh(3)), (domy(2)+bsh(4))
            do j = (domx(1)+bsh(1)), (domx(2)+bsh(2))
                node => list_next(malha(i,j)%list)
                do while (associated(node))
                    ptr = transfer(list_get(node), ptr)
                    ! print*, "n",  ptr%p%n, "F =",ptr%p%F, "v =",ptr%p%v
                    if (propriedade(ptr%p%grupo)%x_lockdelay <= t) then
                        m = propriedade(ptr%p%grupo)%m
                        ptr%p%v(1) = ptr%p%v(1) + ptr%p%F(1)*dt/(2*propriedade(ptr%p%grupo)%m)
                        ptr%p%v(2) = ptr%p%v(2) + ptr%p%F(2)*dt/(2*propriedade(ptr%p%grupo)%m)    
                        ! print '("{{ n ", i2, " F ", f8.3, " ", f8.3)',ptr%p%n,ptr%p%F(1),ptr%p%F(2)
                        ! print '("{{ n ", i2, " v ", f8.3, " ", f8.3)',ptr%p%n,ptr%p%v(1),ptr%p%v(2)
                    end if
                    ptr%p%F = [0,0]
                    node => list_next(node)
                end do
            end do
        end do
    end subroutine comp_v
    
    subroutine walls(icell,jcell,mesh,malha,domx,domy,wall,subx,suby,np,id)
        use linkedlist
        use mod1
        use data
        use mod0
        use, intrinsic :: ieee_arithmetic, only: IEEE_Value, IEEE_QUIET_NAN
        use, intrinsic :: iso_fortran_env, only: real32

        real(real32) :: nan
        type(container), allocatable,dimension(:,:),intent(in) :: malha
        real(dp),intent(in) :: icell(:), jcell(:)
        character(1) :: north, south, east, west
        character(4), intent(in) :: wall
        type(list_t), pointer :: node, previous_node
        integer :: i,j,cell(2)
        integer, intent(in) :: mesh(:), np, subx, suby
        type(data_ptr) :: ptr,ptrn
        integer, intent(in) :: domx(2), domy(2), id
        nan = IEEE_VALUE(nan, IEEE_QUIET_NAN)

        north = wall(1:1)        
        south = wall(2:2)
        east = wall(3:3)
        west = wall(4:4)
        
        if (domy(2) == mesh(2)+2) then 
           !! ! print*, "L 980", id
            if (north == 'e') then !elastic
                i = mesh(2) +2
                do j = domx(1),domx(2)
                    previous_node => malha(i,j)%list
                    node => list_next(malha(i,j)%list)
                    
                    do while (associated(node))
                        ! print*, "cell",i,j, "north", id
                        ptr = transfer(list_get(node), ptr)
                        ptr%p%x(2) = 2*icell(mesh(2)+1) - ptr%p%x(2) 
                        ptr%p%v(2) = -ptr%p%v(2)
                        call list_change(previous_node,malha(i-1,j)%list)
                        node => list_next(previous_node)    
                    end do
                end do
            else if (north == 'o') then
                i = mesh(2) +2
                do j = domx(1), domx(2)
                    previous_node => malha(i,j)%list               
                    node => list_next(malha(i,j)%list)
                    do while (associated(node))
                        ! print*, "cell",i,j
                        ptr = transfer(list_get(node), ptr)
                        ptr%p%x(2) = nan
                        ptr%p%v(2) = nan
                        call list_remove(previous_node)
                        node => list_next(previous_node)                    
                    end do
                end do           
            else !periodic 
                ! Para o caso periódico, as partículas que estão nas celulas fantasmas
                ! vão ter que ir para o outro lado da região de cálculo para a região não- 
                ! fantasma. Mas não deverá ser removida de região fantasma já que elas in-
                ! fluenciam a região ao lado. Isto deverá durar uma iteração e serão apagadas
                ! da região fantasma no comp_x. 
                ! Para o caso paralelo quando mudar de domínio deverão ser transferidas como
                ! no comp_x
                if (suby == 1) then ! caso não paralelo ou o processo calcula ao longo de todo y
                    i = mesh(2) +2
                    do j = domx(1),domx(2)
                        previous_node => malha(i,j)%list
                        node => list_next(malha(i,j)%list)
                        do while (associated(node))
                            ptr = transfer(list_get(node), ptr)
                            if (.not. ptr%p%flag) then 
                                ptr%p%x(2) = -icell(mesh(2)+1) + ptr%p%x(2) 
                                call list_change(previous_node,malha(2,j)%list)
                                node => list_next(previous_node)  
                            else 
                                previous_node => node 
                                node => list_next(node)
                            end if 
                            
                        end do
                    end do
                    i = mesh(2) + 1 ! copiada pra celula fantasma do lado oposto 
                    do j = domx(1),domx(2)
                        previous_node => malha(i,j)%list
                        node => list_next(malha(i,j)%list)
                        do while (associated(node))
                            ! print*, "cell",i,j, "north", id
                            ptr = transfer(list_get(node), ptr)
                            allocate(ptrn%p)
                            ptrn%p = ptr%p 
                            ptrn%p%x(2) = -icell(mesh(2)+1) + ptr%p%x(2)
                            !precisamos de um marcador para particula não voltar  
                            ptrn%p%flag = .true. 
                            call list_insert(malha(1,j)%list, data=transfer(ptrn,list_data))
                            node => list_next(node)    
                        end do
                    end do                   
                end if
                ! print*, "oi"
            end if
        end if
        
       !! print*, "Linha 1018",id
        if (domy(1) == 1) then
           !! ! print*, "L 1019", id
            if (south == 'e') then
                i = 1
                do j = domx(1), domx(2)
                    previous_node => malha(i,j)%list
                    node => list_next(malha(i,j)%list)
                    do while (associated(node))
                        ! print*, "cell",i,j, "south", id
                        ptr = transfer(list_get(node), ptr)
                        ptr%p%x(2) = -ptr%p%x(2)  
                        ptr%p%v(2) = -ptr%p%v(2)
                        call list_change(previous_node,malha(2,j)%list)
                        node => list_next(previous_node)                    
                    end do
                end do
            else if (south == 'o') then
                i = 1
                do j = domx(1), domx(2)
                    previous_node => malha(i,j)%list              
                    node => list_next(malha(i,j)%list)
                    do while (associated(node))
                        ptr = transfer(list_get(node), ptr)
                        ptr%p%x(2) = nan
                        ptr%p%v(2) = nan
                        call list_remove(previous_node)
                        node => list_next(previous_node)                    
                    end do
                end do           
            else !periodic 
                if (suby == 1) then ! caso não paralelo ou o processo calcula ao longo de todo y
                    i = 1
                    do j = domx(1),domx(2)
                        previous_node => malha(i,j)%list
                        node => list_next(malha(i,j)%list)
                
                        do while (associated(node))
                            ! print*, "BB"
                            ! print*, "cell",i,j, "north", id
                            ptr = transfer(list_get(node), ptr)
                            if (.not. ptr%p%flag) then 
                                ptr%p%x(2) = icell(mesh(2)+1) + ptr%p%x(2) 
                                ! print*,ptr%p%n 
                                call list_change(previous_node,malha(mesh(2)+1,j)%list)
                                node => list_next(previous_node)    
                            else 
                                previous_node => node 
                                node => list_next(node)
                            end if 
                        end do
                    end do
                    i = 2 ! copiada pra celula fantasma do lado oposto 
                    do j = domx(1),domx(2)
                        previous_node => malha(i,j)%list
                        node => list_next(malha(i,j)%list)
                        do while (associated(node))
                            ! print*, "cell",i,j, "north", id
                            ptr = transfer(list_get(node), ptr)
                            allocate(ptrn%p)
                            ptrn%p = ptr%p 
                            ptrn%p%x(2) = icell(mesh(2)+1) + ptr%p%x(2)
                            !precisamos de um marcador para particula não voltar  
                            ptrn%p%flag = .true. 
                            call list_insert(malha(mesh(2)+2,j)%list, data=transfer(ptrn,list_data))
                            node => list_next(node)    
                        end do
                    end do                   
                end if
            end if 
        end if 
       !! ! print*, "L 1057", id
        if (domx(1) == 1) then
            if (west == 'e') then
                j = 1
                do i = domy(1), domy(2)
                    previous_node => malha(i,j)%list
                    node => list_next(malha(i,j)%list)
                    do while (associated(node))
                        ! print*, "cell",i,j, "west", id
                        ptr = transfer(list_get(node), ptr)
                        ! print*, "n", ptr%p%x 
                        ptr%p%x(1) = -ptr%p%x(1)
                        ptr%p%v(1) = -ptr%p%v(1)
                        call list_change(previous_node,malha(i,2)%list)
                        node => list_next(previous_node)    
                    end do
                end do    
            else if (west == 'o') then
                j = 1
                do i = domy(1), domy(2)
                    previous_node => malha(i,j)%list
                    
                    node => list_next(malha(i,j)%list)
                    !  ! !print*, 'Linha 236', associated(previous_node)
                    do while (associated(node))
                        !  print*,'N'!,x(ptr,2)
                        !read(*,*)
                        ptr = transfer(list_get(node), ptr)
                        ptr%p%x(2) = nan 
                        ptr%p%v(2) = nan
                        call list_remove(previous_node)
                        node => list_next(previous_node)                    
                    end do
                end do           
            else !periodic 
                if (subx == 1) then
                    j = 1
                    do i = domy(1), domy(2)
                        previous_node => malha(i,j)%list
                        node => list_next(malha(i,j)%list)
                        do while (associated(node))
                            ! print*, "cell",i,j, "west", id
                            ptr = transfer(list_get(node), ptr)
                            ! print*, "n", ptr%p%x 
                            if (.not. ptr%p%flag) then 
                                ptr%p%x(1) = ptr%p%x(1) + jcell(mesh(1)+1)
                                call list_change(previous_node,malha(i,mesh(1)+1)%list)
                                node => list_next(previous_node)    
                            else 
                                previous_node => node 
                                node => list_next(node)
                            end if
                        end do
                    end do  

                    j = 2
                    do i = domy(1), domy(2)
                        previous_node => malha(i,j)%list
                        node => list_next(malha(i,j)%list)
                        do while (associated(node))
                            ! print*, "cell",i,j, "west", id
                            ptr = transfer(list_get(node), ptr)
                            ! print*, "n", ptr%p%x 
                            allocate(ptrn%p)
                            ptrn%p = ptr%p 
                            ptrn%p%x(1) = ptr%p%x(1) + jcell(mesh(1)+1)
                            ptr%p%flag = .true. 
                            call list_insert(malha(i,mesh(1)+2)%list, data=transfer(ptr,list_data))
                            node => list_next(node)    
                        end do
                    end do  

                end if
            end if
        end if
        if (domx(2) == mesh(1) +2) then
            if (east == 'e') then
                j = mesh(1)+2
                do i = domy(1), domy(2)
                    previous_node => malha(i,j)%list
                    node => list_next(malha(i,j)%list)
                    do while (associated(node))
                        ! print*, "cell",i,j, "east", id
                        ptr = transfer(list_get(node), ptr)
                        ptr%p%x(1) = 2*jcell(mesh(1)+1) - ptr%p%x(1) 
                        ptr%p%v(1) = -ptr%p%v(1)
                        call list_change(previous_node,malha(i,j-1)%list)
                        node => list_next(previous_node)   
                    end do
                end do
            else if (east == 'o') then
                j = mesh(1)+2
                do i = domy(1), domy(2)
                    previous_node => malha(i,j)%list
                    node => list_next(malha(i,j)%list)
                    do while (associated(node))
                        ptr = transfer(list_get(node), ptr)
                        ptr%p%x(1) = nan 
                        ptr%p%v(1) = nan
                        call list_remove(previous_node)
                        node => list_next(previous_node)                    
                    end do
                end do           
            else !periodic 
                if (subx == 1) then
                    j = mesh(1)+2
                    do i = domy(1), domy(2)
                        previous_node => malha(i,j)%list
                        node => list_next(malha(i,j)%list)
                        do while (associated(node))
                            ! print*, "cell",i,j, "west", id
                            ptr = transfer(list_get(node), ptr)
                            ! print*, "n", ptr%p%x 
                            if (.not. ptr%p%flag) then 
                                ptr%p%x(1) = ptr%p%x(1) - jcell(mesh(1)+1)
                                call list_change(previous_node,malha(i,2)%list)
                                node => list_next(previous_node)    
                            else 
                                previous_node => node 
                                node => list_next(node)
                            end if
                        end do
                    end do  
                    j = mesh(1)+1
                    do i = domy(1), domy(2)
                        previous_node => malha(i,j)%list
                        node => list_next(malha(i,j)%list)
                        do while (associated(node))
                            ! print*, "cell",i,j, "west", id
                            ptr = transfer(list_get(node), ptr)
                            ! print*, "n", ptr%p%x 
                            allocate(ptrn%p)
                            ptrn%p = ptr%p 
                            ptrn%p%x(1) = ptr%p%x(1) - jcell(mesh(1)+1)
                            ptrn%p%flag = .true. 
                            call list_insert(malha(i,1)%list, data=transfer(ptr,list_data))
                            node => list_next(node)    
                        end do
                    end do  
                end if
            end if
        end if
        ! print*,"linha 1146 ok"
        ! call MPI_barrier(MPI_COMM_WORLD, ierr) 
     !    !read(*,*)
    end subroutine walls
    
    subroutine distribuicao_inicial(malha,N,mesh)
        use linkedlist
        use mod1
        use data
        use mod0
        
        type(container), allocatable,dimension(:,:),intent(in) :: malha
        integer :: i,j,aux1
        integer, intent(in) :: N,mesh(2)
        type(list_t), pointer :: node
        type(data_ptr) :: ptr
        
        
        do i = 1,mesh(2)+2
            do j = 1,mesh(1)+2
                node => list_next(malha(i,j)%list)
                aux1 = 0
                do while (associated(node))
                    aux1 = aux1 +1
                    ptr = transfer(list_get(node), ptr)
                    ! celula(ptr,:) = [j,i]
                    print *, 'posição', i, ',', j
                    print *, ptr%p%x
                    node => list_next(node)
                    if (aux1 > N) then
                        print*, 'deu ruim'
                        exit
                    end if
                end do
            end do
        end do
    end subroutine distribuicao_inicial

    subroutine clean_mesh(malha, mesh, domx, domy,id, tudo)
        !  This subroutine cleans the unaccessible cells of each process bounded by domx, domy
        use linkedlist
        use mod1
        use data
        use mod0

        type(container), allocatable,dimension(:,:) :: malha
        integer :: i,j, k = 0
        integer, intent(in) :: mesh(2),domx(2),domy(2)
        type(list_t), pointer :: node, previous_node
        type(data_ptr) :: ptr
        logical, intent(in) :: tudo 
        ! print*, "AAA", id
        if (tudo) k = 1
        do i = 2-k,mesh(1)+1+k
            do j = 2-k, mesh(2)+1+k
                if (((i < domy(1)-1 .or. i > domy(2)+1) .and. (j < domx(1)-1 .or. j > domx(2)+1)) .or. tudo) then 
                    node => list_next(malha(i,j)%list)
                    ! print*, "i,j,id", i,j, id
                    do while (associated(node)) 
                        ptr = transfer(list_get(node), ptr)
                        deallocate(ptr%p)
                        node => list_next(node)
                    end do
                    node => list_next(malha(i,j)%list) 
                    if (associated(node)) then
                        call list_free(malha(i,j)%list)
                        call list_init(malha(i,j)%list)
                        ! if (.not. tudo) then
                        !     print*, "BBB"
                        ! else
                        !     print*, "CCC"
                        ! end if 
                    end if
                end if
            end do
        end do

        if (tudo) then
            do i = 2-k,mesh(1)+1+k
                do j = 2-k, mesh(2)+1+k
                    deallocate(malha(i,j)%list)
                end do
            end do
            ! dellocate( )
        end if

    end subroutine clean_mesh

end module fisica


program main
    use mod1
    use linkedlist
    use data
    use mod0
    use matprint
    use saida
    use m_config
    use fisica
    use randnormal
    use mpi
    use, intrinsic :: ieee_arithmetic, only: IEEE_Value, IEEE_QUIET_NAN
    use, intrinsic :: iso_fortran_env, only: real32

    implicit none
!    Variáveis
    integer :: N,Ntype,i=1, nimpre,j = 1, ii, quant = 0,mesh(2), cont = 1, aux1 = 0,cont2 = 1,domx(2), domy(2), aux3
    integer :: subx, suby, NMPT, j2, cold_cells(4), hot_cells(4), nimpre_init, print_TC
    integer, target :: k
    real(dp), dimension(:,:), allocatable :: v, x, celula, rFUp !força n e n+1, velocidades e posições, [r*F, potencial e momento]
    real(dp), dimension(:), allocatable :: icell,jcell, nxv, nxv_send, nRfu, nRfu_send !dimensões das celulas e vetor de resultado pra imprimir
    real(dp) :: t=0,t_fim,dt,printstep, sigma, epsil, rcut,aux2,start = 0,finish,Td,kb = 1.38064852E-23,vd(2)
    real(dp) :: GField(2), temp_Td(3), dimX, dimY, dx_max, Td_hot, Td_cold
    integer,allocatable :: interv(:), interv_Td(:), grupo(:), rcounts(:), displs(:) !, mic(:,:), mic_rcv(:) ! mic são vetores para contar quantas vezes atravessou a borda periodica
    type(container), allocatable,dimension(:,:) :: malha
    type(prop_grupo),allocatable,dimension(:) :: propriedade
    type(list_t), pointer :: node
 !   Variáveis do arquivo de configuração 
    type(CFG_t) :: my_cfg
    character(5) :: particle
    character(4) :: wall
    character(20) :: nome,arquivox, arquivov = '%' !nome de partícula deve ter até 20 caracteres
    type(string) :: part_nomes(10) ! vetor de strings
    character(1) :: optio
    type(data_ptr) :: ptr !integer,pointer :: ptr,ptrn
    logical :: laux, termostato
    character(len=32) :: arg
    character :: arg1(1)
    type(lstr) :: LT
    integer ( kind = 4 ) status(MPI_STATUS_SIZE)
    integer ( kind = 4 ) ierr
    integer ( kind = 4 ) np
    integer ( kind = 4 ) tag
    integer ( kind = 4 ) id, ids(8) 
    real(real32) :: nan

    nan = IEEE_VALUE(nan, IEEE_QUIET_NAN)
    

    ! depois ver se não da pra fazer tudo serialmente
    call MPI_Init ( ierr )
    call MPI_Comm_rank ( MPI_COMM_WORLD, id, ierr )                            
    call MPI_Comm_size ( MPI_COMM_WORLD, np, ierr )

    ! argumento para checar NAN no nancheck
    if(id == 0) call system('python verify_settings.py')
    call MPI_barrier(MPI_COMM_WORLD, ierr)
    call get_command_argument(1, arg)
    arg1 = trim(arg)
    
!    CONDIÇÕES INICIAIS/CORTORNO E PROPRIEDADES
    ! Propriedades que cobrem todo um grupo de partículas
    !  serão armazenadas na estrutura propriedade(grupo)
!    Inicia o método de configuração
    call CFG_add(my_cfg, "global%N",1,&
      "number of particles")
    call CFG_add(my_cfg, "global%Ntype",1,&
      "number of types")
    call CFG_add(my_cfg, "global%dt",0.0_dp,&
      "time step")
    call CFG_add(my_cfg, "global%t_fim",0.0_dp,&
      "simulation time") 
    call CFG_add(my_cfg, "global%nimpre_init",0,&
      "start nimpre") 
    call CFG_add(my_cfg, "global%nimpre",1,&
      "number of result files")
    call CFG_add(my_cfg, "global%dimX",1.0_dp,&
      "y-dimension of the region of calculus")  
    call CFG_add(my_cfg, "global%dimY",1.0_dp,&
      "y-dimension of the region of calculus")    
    call CFG_add(my_cfg, "global%mesh",(/0, 0/),&
      "mesh size")  !celulas
    call CFG_add(my_cfg, "global%sigma",1.0_dp,&
      "Leonard-Jones sigma factor")        
    call CFG_add(my_cfg, "global%epsilon",1.0_dp,&
      "Leonard-Jones epsilon factor")   
    call CFG_add(my_cfg, "global%rcut",1.0_dp,&
      "Cut radius")
    call CFG_add(my_cfg,"global%wall",'eeee', &
        "wall's periodic vs elastic") 
    call CFG_add(my_cfg,"global%Td",-1.0_dp, &
        "Thermostat")
    call CFG_add(my_cfg,"global%temp_Td",(/0.0_dp, 0.0_dp, 0.0_dp/), &
        "Thermostat iteractions time")
    call CFG_add(my_cfg,"global%Td_hot",-1.0_dp, &
        "Thermostat hot wall")
    call CFG_add(my_cfg,"global%Td_cold",-1.0_dp, &
        "Thermostat cold wall")
    call CFG_add(my_cfg,"global%cold_cells",(/0, 0, 0, 0/), &
        "Cold Cells")     
    call CFG_add(my_cfg,"global%hot_cells",(/0, 0, 0, 0/), &
        "Hot Cells")
    call CFG_add(my_cfg,"global%vd",(/0.0_dp, 0.0_dp/), &
        "Mean vd")                 
    call CFG_add(my_cfg,"global%NMPT",1, &
        "Max Number of particles that will change process per iteraction")
    call CFG_add(my_cfg,"global%GField",(/0.0_dp, 0.0_dp/), &
        "Uniform Gravitational Field")
    call CFG_add(my_cfg,"global%print_TC",1, &
        "Print transport coefficient data")
   
    call CFG_read_file(my_cfg, "settings.ini")
    call CFG_get(my_cfg,"global%N",N)
    call CFG_get(my_cfg,"global%Ntype",Ntype)
    call CFG_get(my_cfg,"global%t_fim",t_fim)
    call CFG_get(my_cfg,"global%nimpre_init",nimpre_init)
    call CFG_get(my_cfg,"global%dt",dt)
    call CFG_get(my_cfg,"global%nimpre",nimpre)
    call CFG_get(my_cfg,"global%dimX", dimX)
    call CFG_get(my_cfg,"global%dimY", dimY)
    call CFG_get(my_cfg,"global%mesh",mesh)
    call CFG_get(my_cfg,"global%rcut",rcut)
    call CFG_get(my_cfg,"global%wall",wall)
    call CFG_get(my_cfg,"global%Td",Td)
    call CFG_get(my_cfg,"global%temp_Td",temp_Td)
    call CFG_get(my_cfg,"global%cold_cells",cold_cells)
    call CFG_get(my_cfg,"global%hot_cells",hot_cells)
    call CFG_get(my_cfg,"global%Td_cold",Td_cold)
    call CFG_get(my_cfg,"global%Td_hot",Td_hot)
    call CFG_get(my_cfg,"global%vd",vd)
    call CFG_get(my_cfg,"global%NMPT",NMPT)
    call CFG_get(my_cfg,"global%GField",GField)
    call CFG_get(my_cfg,"global%print_TC",print_TC)
    

    if (NMPT < 1) NMPT = N
    printstep = ((temp_Td(2)-temp_Td(1))/dt)/temp_Td(3) 
    if (printstep < 1) printstep = 1
    !aloca as Variáveis com base no número de partículas
    allocate(v(N,2),x(N,2),grupo(N),interv(nimpre),interv_Td(nint(printstep)), & 
        jcell(mesh(2)+1),icell(mesh(1)+1),celula(N,2),propriedade(Ntype))
    allocate(LT%lstrdb_N(NMPT*6),LT%lstrdb_S(NMPT*6),LT%lstrdb_E(NMPT*6),LT%lstrdb_W(NMPT*6),LT%lstrdb_D(NMPT*6), &
    LT%lstrint_N(NMPT*4),LT%lstrint_S(NMPT*4),LT%lstrint_E(NMPT*4),LT%lstrint_W(NMPT*4), LT%lstrint_D(NMPT*4)) 

    !if (wall == 'pppp') then 
    !    allocate(mic(N))
    !    mic = 0
    !    if (id == 0) allocate(mic_rcv(N))
    !end if

    ! allocate(nxv(N*5),nxv_send(N*5))

    aux2 = dimX/mesh(1)
    jcell = (/((i*aux2),i = 0,mesh(1))/) ! direção x
    
    aux2 = dimY/mesh(2)
    icell = (/((i*aux2),i = 0,mesh(2))/) ! direção Y
   
   
    do i = 0,(Ntype-1)     
        if (id == 0) write(*,'(a,i0)') 'par_',i
        write(particle,'(a,i0)') 'par_',i
        call CFG_add(my_cfg, particle//"%m",1.1_dp,&
           "mass of "//particle)
        call CFG_add(my_cfg, particle//"%v",(/0.0_dp, 0.0_dp/), &
            "global velocity of "//particle)         
        call CFG_add(my_cfg, particle//"%x","arquivo.csv", &
            "position of "//particle)   
        call CFG_add(my_cfg, particle//"%v_file","arquivo.csv", &
            "velecities of "//particle)     
        call CFG_add(my_cfg, particle//"%nome",'abcde', &
            "name of "//particle)  
        call CFG_add(my_cfg, particle//"%quantidade",1, &
            "quantity of "//particle)   
        call CFG_add(my_cfg, particle//"%epsilon",1.1_dp, &
            "sigma "//particle)   
        call CFG_add(my_cfg, particle//"%sigma",1.1_dp, &
            "epsilon"//particle)  
        call CFG_add(my_cfg, particle//"%x_lockdelay",1.1_dp, &
            "change in position delay "//particle)       
        call CFG_add(my_cfg, particle//"%rs",1.1_dp, &
            "solid radius "//particle) 
        call CFG_add(my_cfg,particle//"%fric_term",1.0_dp, &
            "Friction term")    
    end do
    
    dx_max = 10*dimx !pra definir critério de estabilidade no uso de malha
    do i = 0,(Ntype-1)
        write(particle,'(a,i0)') 'par_',i
    
        call CFG_get(my_cfg, particle//"%quantidade", quant)
        call CFG_get(my_cfg, particle//"%nome", nome)
        part_nomes(i+1)%str = nome
        call CFG_get(my_cfg, particle//"%x", arquivox)
        call CFG_get(my_cfg, particle//"%v_file", arquivov)
        call CFG_get(my_cfg, particle//"%m", propriedade(i+1)%m)
        call CFG_get(my_cfg, particle//"%epsilon", propriedade(i+1)%epsilon)
        call CFG_get(my_cfg, particle//"%sigma", propriedade(i+1)%sigma)
        call CFG_get(my_cfg, particle//"%x_lockdelay", propriedade(i+1)%x_lockdelay)
        call CFG_get(my_cfg, particle//"%rs", propriedade(i+1)%rs)
        call CFG_get(my_cfg, particle//"%fric_term",propriedade(i+1)%fric_term)

        ! le o arquivo com posições
        if (dx_max < propriedade(i+1)%sigma*rcut/2) then
            dx_max = propriedade(i+1)%sigma*rcut/2
        end if 
        open(20,file=arquivox,status='old')
        if (arquivov(1:1) /= '%') then ! usa arquivo com velocidades
            open(30,file=arquivov,status='old')
            do j = 1,quant
                grupo(cont) = i+1
                read(20,*) x(cont,:) 
                read(30,*) v(cont,:)
                cont = cont+1
            end do   
            close(30)     
        else 
            do j = 1,quant
                call CFG_get(my_cfg, particle//"%v", v(cont,:))    
                grupo(cont) = i+1
                read(20,*) x(cont,:) 
                cont = cont+1
            end do
        end if
        close(20)
    end do
    dx_max = dx_max**2
   ! mostra informações lidas
    if (id == 0) then
        call CFG_write(my_cfg, "settings.txt") 
        print*, "Nomes das partículas"
        do i = 1,Ntype
            print*, part_nomes(i)%str
        end do
    end if     
    
    !--------------------------------------------------!
    ! CRIA LISTAS LIGADAS

    !cada elemento da matriz representa uma celula da malha
    allocate(malha(mesh(2)+2,mesh(1)+2)) ! o primeiro e últimos índices correspondem a 
                                         ! celulas fantasma
    !inicialização
    k = 0
    do i = 1,mesh(2)+2 !linha
        do j = 1,mesh(1)+2
            call list_init(malha(i,j)%list) !, DATA=transfer(ptr, list_data))
        end do
    end do
    ! Vetor que indica quando imprimir
    printstep = nint(t_fim/dt)/nimpre
    interv = (/(printstep*i,i=0,nimpre) /)
    ! Vetor que indica quando fazer velocity scaling
    printstep = nint(temp_Td(3))
    interv_Td = (/(printstep*i,i=0,nint( ((temp_Td(2)-temp_Td(1))/dt)/temp_Td(3) )) /)
    interv_Td = interv_Td + nint(temp_Td(1)/dt)
    i = 0
    
    if (Td == 0) then
        Td = (2/(2*N*kb))*sum(0.5*propriedade(ptr%p%grupo)%m*(v(:,1)**2+v(:,2)**2))
    else if (Td < 0 .and. (Td_hot*Td_hot) <= 0) then
        interv_Td = -1
    end if
    !Aloca as partículas nas celulas
    do k = 1,cont-1
        laux = .true.
        j = 2
        do while (laux)
            if (x(k,1) > jcell(j)) then
                j = j+1
            else
                laux = .false.       
            end if 
        end do
        
        laux = .true.
        i = 2
        do while (laux)
            if (x(k,2) > icell(i)) then
                i = i+1
            else
                laux = .false. 
            end if 
        end do        
        !i e j contém agora a posição da partícula na matriz malha
        allocate(ptr%p)
        ptr%p%x = x(k,:)
        ptr%p%v = v(k,:)
        ptr%p%grupo = grupo(k)
        ptr%p%F = [0,0] !passo 1 do SV
        ptr%p%n = k !identidade da partícula importante para imprimir
        ! ptr%p%mic = [0,0]
        call list_insert(malha(i,j)%list, data=transfer(ptr, list_data)) 
    end do
    deallocate(grupo)

    ! adiciona às particulas velocidade inicial de acordo com distribuição maxwell boltzmann
    if (vd(1) /= 0 .and. vd(2) /= 0) call MaxwellBoltzmann(malha,mesh, vd)
    if (vd(1) == 0 .and. vd(2) == 0 .and. Td > 0) then
        vd = sqrt([Td, Td])
        call MaxwellBoltzmann(malha,mesh,vd)
    end if 
    ! adicionamos uma celula para corrigir o fato de estarmos usando fantasmas
    cold_cells = cold_cells + 1
    hot_cells = hot_cells + 1
    !-------------------------------------------------!
    ! ITERAÇÕES
    
    ! imprime condições iniciais

    if (id == 0) then  
        j = 0 
        call system('mkdir temp') !pasta temporária para armazenar os resultados
        ! call linked2vec(malha,domx,domy,nxv,aux1)
        call vec2csv(x,N,2,'position',j,t,nimpre,start)
        call vec2csv(v,N,2,'velocity',j,t,nimpre,start)
        ! call vec2csv(celula,N,2,'cell',j)
        j = j + 1

        ! print*,'V_tot0 = ',comp_pot(mesh,malha,propriedade,rcut)
    end if         

    !! AQUI COMEÇA O PARALELISMO                                           

    ! só vai funcionar com serial e np par
    ! Descobir a melhor forma de dividir a malha
    
    if (dimX > 3*dimY .and. np > 1) then
        subx = 4
        suby = 1
    else if (dimY > 3*dimX .and. np > 1) then
        subx = 1
        suby = 4
    else if (np > 1) then
    !if (np > 1) then
        subx = 2
        suby = 2
    else 
        subx = 1
        suby = 1
    end if  

    print*, 'subx', subx, 'suby', suby
    np = subx*suby
    ids = [-1,-1,-1,-1,-1,-1,-1,-1]

    if (np > 2) then
        L_O: do j = 0, (suby-1)
            do i = 0, (subx-1)
                if (id == (i+j*suby)) then 
                    ! ids = [N,S,E,W, NE, NW, SE, SW]
                    ! Depois tentar para caso geral 
                    ! ids = [(j-1)*subx +i, (j+1)*subx +i, &
                    !         j*subx + i +1, -(i>0)*(j*subx + i)-1]
                    ! if (ids(1) < 0) ids(1) = -1
                    ! if (ids(4) < 0) ids(4) = -1 
                    ! if (j == suby-1) ids(2) = -1
                    ! if (i == subx -1) ids(3) = -1 
                    if (subx == 2) then
                        if (id == 0) then 
                            ids = [2, -1, 1,-1, 3, -1, -1, -1]
                        else if (id == 1) then
                            ids = [3, -1, -1, 0, -1, 2, -1, -1]
                        else if (id == 2) then
                            ids = [-1, 0, 3, -1, -1, -1, 1, -1]
                        else 
                            ids = [-1, 1, -1, 2, -1, -1, -1, 0]
                        end if
                    else if (subx == 4) then
                        if (id == 0) then 
                            ids = [-1, -1, 1,-1, -1, -1, -1, -1]
                        else if (id == 1) then
                            ids = [-1, -1, 2, 0, -1, -1, -1, -1]
                        else if (id == 2) then
                            ids = [-1, -1, 3, 1, -1, -1, -1, -1]
                        else 
                            ids = [-1, -1, -1, 2, -1, -1, -1, -1]
                        end if
                    else if (suby == 4) then
                        if (id == 0) then 
                            ids = [1, -1, -1, -1, -1, -1, -1, -1]
                        else if (id == 1) then
                            ids = [2, 0, -1, -1, -1, -1, -1, -1]
                        else if (id == 2) then
                            ids = [3, 1, -1, -1, -1, -1, -1, -1]
                        else 
                            ids = [-1, 2, -1, -1, -1, -1, -1, -1]
                        end if
                    end if 
                    domx = [int(i*(mesh(1)+2)/subx)+1, int((i+1)*(mesh(1)+2)/(subx))]
                    domy = [int(j*(mesh(2)+2)/suby)+1, int((j+1)*(mesh(2)+2)/(suby))]
                    exit L_O
                end if
            end do
        end do L_O
    else 
        domx = [1,mesh(1)+2]
        domy = [1,mesh(2)+2]
    end if

    print '("id = ",I2, "; domx = [",I3," ",I3, "]; domy = [",I3," ",I3, "]; '& 
    'icell = [",F10.3," ",F10.3,"] jcell = [",F10.3," ",F10.3,"]")' &
    , id, domx(1),domx(2), domy(1), domy(2), icell(domy(1)),  icell(domy(2)-1), jcell(domx(1)),  jcell(domx(2)-1)
    print '("id = ", I2, " ids = [", I2, " ", I2, " ", I2, " ", I2,"]")', id, ids(1), ids(2), ids(3), ids(4)
    ! libera a memória utilizada pelos processos não-raiz
    ! pois eles não vão imprimir x e v pra csv
    if (id /= 0) then
        deallocate(v,x)
    end if
    ! aloca alguns vetores utilizados para o mpi_gatherv
    allocate(rcounts(np), displs(np))

    !print '("Mesh divided in ", I2, "x", I2, " subregions."  )', subx, suby
    ! print*, 'ID ', id, 'iterv ', interv
    call clean_mesh(malha, mesh, domx, domy,id,.false.)
    ! print*, "bbb", id
  
    ! if (id == 0) then
    !     print*,'Press Return to continue'
    !    read(*,*)
    ! end if
    call MPI_barrier(MPI_COMM_WORLD, ierr)
    
    ! id = 0 é o processo raiz 
    if (id == 0)  call cpu_time(start)

    ! Agora dividimos para cada id uma região. 
    i = 0
    j = 1
    j2 = 1
   
    if (nimpre_init > 0) then
        j = nimpre_init
        t = t_fim*nimpre_init/nimpre
        i = interv(j)
    end if

    if (print_TC == 1) then
        allocate(nRfu(N*5),nRfu_send(N*6),rFUp(N,6))
        call system('mkdir temp2')
    end if

    do while (t_fim > t)
        ! print*, "L 1990"
        call comp_F(GField, mesh,malha,propriedade,rcut,domx,domy,ids,id,t)  !altera Força
        ! IDS são os ids das regiões vizinhas, vetor(4) int posições, 
        ! DOMX e DOMY são vetores(2) int com o domínio (quat. de celulas)
        ! que o processo vai cuidar (subdivisões)
        ! print*, "L 1702", id
        ! call MPI_barrier(MPI_COMM_WORLD, ierr)
        ! print*, "L 1999", id
        !  if (id == 0) print*, "TEMPO", t, '<<<<<<<<<<', id, i
        call comp_x(icell,jcell,malha,N,mesh,propriedade, dx_max,t,dt,ids,LT,domx,domy,wall,id, np) ! altera posição
        ! print*, "L 2002", id
        ! if (id == 0) read(*,*)
        ! call MPI_barrier(MPI_COMM_WORLD, ierr)
        
        call walls(icell,jcell,mesh,malha,domx,domy,wall,subx,suby,np,id) ! altera posição e malha
        
        ! print*, "L 1714", id
        ! if (id == 0) read(*,*)
        ! call MPI_barrier(MPI_COMM_WORLD, ierr)
        call comp_F(GField, mesh,malha,propriedade,rcut,domx,domy,ids,id,t) !altera força
        ! print*, "L 1717"
        call comp_v(malha,mesh,dt,t,propriedade,domx,domy) !altera velocidade
        ! print*, "1721"

        t = t + dt
        cont2 = cont2+1

        if (i == interv_Td(j2)) then
            if (Td > 0) then
                call corr_Kglobal(malha,domx,domy,Td,propriedade, np,id,t)
            end if
            if (Td_hot > 0 .or. Td_cold > 0) then
                call corr_K(malha,domx,domy,cold_cells,hot_cells,Td_hot,Td_cold,propriedade, np,id,t)
            end if
            j2 = j2 + 1
        end if 
        
        if (i == interv(j)) then
            deallocate(LT%lstrdb_N, LT%lstrdb_S, LT%lstrdb_E, LT%lstrdb_W, LT%lstrdb_D, &
            LT%lstrint_N, LT%lstrint_S, LT%lstrint_E, LT%lstrint_W, LT%lstrint_D)
            allocate(nxv(N*5),nxv_send(N*5))
            ! call MPI_barrier(MPI_COMM_WORLD, ierr) 
            ! print*, 'MANDANDO PRINTAR', id, i
            
            call linked2vec(malha,mesh,domx,domy,nxv_send,aux1)
            ! ! ! print*, "L nxv_send", nxv_send  
            !! print*, "LINKED 2 VEC, id", id, "aux1", aux1
            if (np > 1) then ! paralelo
                ! call MPI_GATHER(sbuf, scount, MPI_integer, rbuf, rcount, MPI_integer, root, MPI_COMM_WORLD, ierr)
                call MPI_GATHER(aux1,     1,    MPI_integer, rcounts, 1,   MPI_integer,  0,   MPI_COMM_WORLD, ierr)

                ! print*, "Aux1 =", aux1, "id =", id
                ! if (id == 0) print*, "rcounts =", rcounts
                ! !print*, 'L IMPRE >', id
                displs = 0

                if (id == 0) then
                    do ii = 2,np
                        displs(ii) = displs(ii-1) + rcounts(ii-1)
                    end do
                    ! print*, "displs", displs
                end if
                ! print*,'E_tot0 = ',(comp_pot(mesh,malha,propriedade,rcut) +comp_K(malha,mesh,propriedade))
                ! MPI_GATHERV    (sbuf,   scount,  stype, rbuf, rcounts,  displs,   rtype,  root,  comm,  ierr)
                call MPI_GATHERV(nxv_send, aux1, MPI_DOUBLE_PRECISION, nxv, rcounts, &
                 displs, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
            else 
                nxv = nxv_send
            end if 
            
            if (id == 0) then
                do ii = 0, N-1
                    x(int(nxv(ii*5+1)),:) = [nxv(ii*5+2),nxv(ii*5+3)]
                    v(int(nxv(ii*5+1)),:) = [nxv(ii*5+4),nxv(ii*5+5)]
                end do
                
                call vec2csv(x,N,2,'position',j,t,nimpre,start)
                call vec2csv(v,N,2,'velocity',j,t,nimpre,start)
            
            end if
            ! if (id == 0) read(*,*)
            ! call MPI_barrier(MPI_COMM_WORLD, ierr) 
            
            deallocate(nxv,nxv_send)

            if (print_TC == 1) then 
                nRfu = 0.0_dp

                call comp_pot(mesh,malha,propriedade,rcut,domx,domy, ids, id, t, nRfu) 

                !! print*, "LINKED 2 VEC, id", id, "aux1", aux1
                if (np > 1) then ! paralelo
                    aux3 = 1
                    do ii = 1,N
                        if (nRfu(ii*6-5) /= 0) then
                            nRfu_send(aux3:aux3+5) = nRfu(ii*6-5:ii*6)
                            aux3 = aux3 + 6
                        end if
                    end do
                    aux1 = aux3 -1

                    ! call MPI_GATHER(sbuf, scount, MPI_integer, rbuf, rcount, MPI_integer, root, MPI_COMM_WORLD, ierr)
                    call MPI_GATHER(aux1,     1,    MPI_integer, rcounts, 1,   MPI_integer,  0,   MPI_COMM_WORLD, ierr)
    
                    displs = 0
    
                    if (id == 0) then
                        do ii = 2,np
                            displs(ii) = displs(ii-1) + rcounts(ii-1)
                        end do
                        ! print*, "displs", displs
                    end if
                    ! if (id == 0) then 
                    !     print*, "L 3099"
                    !     read(*,*)
                    ! end if
                    ! call MPI_barrier(MPI_COMM_WORLD, ierr) 

                    ! MPI_GATHERV    (sbuf,   scount,  stype, rbuf, rcounts,  displs,   rtype,  root,  comm,  ierr)
                    call MPI_GATHERV(nRfu_send, aux1, MPI_DOUBLE_PRECISION, nRfu, rcounts, &
                     displs, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
                end if 
                ! print*, "AA"          

                if (id == 0) then
                    do ii = 0, N-1
                        ! print*, [nRfu(ii*6+2),nRfu(ii*6+3),nRfu(ii*6+4),nRfu(ii*6+5),nRfu(ii*6+6)]
                        rFUp(int(nRfu(ii*6+1)),:) = [nRfu(ii*6+2),nRfu(ii*6+3),nRfu(ii*6+4),nRfu(ii*6+5),nRfu(ii*6+6)]
                    end do
                    ! print*, "AAA"
                    call vec2csv(rFUp,N,5,'rF_u_P',j,t,nimpre,start)
                    
                end if 
                ! print*, "B"
                ! read(*,*)
                ! deallocate(nRfu,nRfu_send,rFUp)
                ! print*, "C"
            end if 
            allocate(LT%lstrdb_N(NMPT*6),LT%lstrdb_S(NMPT*6),LT%lstrdb_E(NMPT*6),LT%lstrdb_W(NMPT*6),LT%lstrdb_D(NMPT*6), &
            LT%lstrint_N(NMPT*4),LT%lstrint_S(NMPT*4),LT%lstrint_E(NMPT*4),LT%lstrint_W(NMPT*4), LT%lstrint_D(NMPT*4))
            j = j+1
        end if
        i = i+1
        ! print*, "L 1754 >", id
        ! if (id  == 0) read(*,*)    
        ! call MPI_barrier(MPI_COMM_WORLD, ierr)
        ! print*, t
    end do

    call clean_mesh(malha, mesh, domx, domy,id,.true.)
    deallocate(LT%lstrdb_N, LT%lstrdb_S, LT%lstrdb_E, LT%lstrdb_W, LT%lstrdb_D, &
            LT%lstrint_N, LT%lstrint_S, LT%lstrint_E, LT%lstrint_W, LT%lstrint_D)
    ! print*, 'L 3104'
    ! if (print_TC == 1) deallocate(nRfu,nRfu_send,rFUp)
    
    if (id == 0)  then
        call cpu_time(finish)
        print '("Time = ",f10.3," seconds.")',(finish-start)
        open(unit=22,file='settings.txt',status="old", position="append", action="write")
        write(22,*) "#:Execution time = ",(finish-start)," seconds."
        close(22)
        call system('python csv2vtk_particles.py')
    end if
    
    ! print*, 'L 1808'
    call MPI_Finalize ( ierr )
end program main                