subroutine comp_FT(GField,Hfield,theta,Tor,pr,mesh,malha,propriedade,r_cut,domx,domy, ids, id, t,dt,partlst)
    use linkedlist
    use mod1
    use data
    use mod0
    use mpi
    
    type(prop_grupo), allocatable,dimension(:),intent(in) :: propriedade
    type(container), allocatable,dimension(:,:),intent(in) :: malha
    real(dp), intent(in) :: dt, GField(2),hfield(5), pr(5), t ! outras propriedades da particula que gira [A, B, alpha, beta, phase]
    real(dp), intent(inout), dimension(:) :: Tor ! torque
    real(dp), intent(in), dimension(:) :: theta ! ângulo das partículas que giram, tangencial
    real(dp) :: sigma, epsil, sigma_a, epsil_a,sigma_b, epsil_b, rcut,r_cut, fric_term !fric_term = força de ficção
    real(dp) :: x1(2),v1(2),x2(2),v2(2), rs1, rs2, coss, sine, A, B, alpha,beta,ph   
    integer :: i,j, ct = 0, dox(2), doy(2) !,ptr, ptrn
    integer, intent(in) :: mesh(:),domx(2),domy(2), id
    real(dp) :: Fi(2)=0,r, aux2(2),aux3,fR(2), fric_term1, fric_term2
    integer, save :: pa
    type(list_t), pointer :: node, next_node
    type(data_ptr) :: ptr,ptrn
    integer, intent(out), dimension(:) :: partlst
    integer ( kind = 4 ), intent(in) :: ids(8)

    partlst = 0
    A =     pr(1)
    B =     pr(2)
    alpha = pr(3)
    beta =  pr(4)
    ph =    pr(5)

    ! quanta a quantidade de partículas antes das que giram
    ! print*, "t,2*dt", t, 2*dt
    if (t < 2*dt) then
        rs = 0
        i = 0
        pa = 0
        do while (rs == 0)
            i = i+1
            rs = propriedade(i)%rs
            pa = pa + propriedade(i)%quant
        end do
        pa = pa - propriedade(i)%quant
        ! print*, "pa calculado!",id,pa
    end if 
    ! print*, "pa = ", pa
    

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
                ! PRINT*, "L 759"
                ptr%p%flag = .true. ! indica ao comp_x que a partícula precisa ser calculada
                ! PRINT*, "L 764"
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

                if (rs1 > 0) then
                    partlst(ptr%p%n - pa) = 1
                    Tor(-pa +ptr%p%n) = Tor(-pa +ptr%p%n) + magtorque(theta(-pa+ptr%p%n),hfield,t)
                end if
                ! print '("x1  =", f10.6, " ", f10.6, " n ", i2, " cell ",i2," ",i2)',x1(1),x1(2), ptr%p%n,i,j
                ! if (id == 0) read(*,*)
                !calcular a força desta com todas as outras partículas
                next_node => list_next(node) ! próxima partícula da célula
                node => list_next(node) ! a ser computado com a particula selecionada
         
                ! PRINT*, "L 780"
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
                       if (rs1 == 0 .and. rs2 == 0) then
                            aux2 = -(1/r**2)*(sigma/r)**6*(1-2*(sigma/r)**6)*24*epsil* & 
                            [(x1(1)-x2(1)) - (rs1+rs2)*coss, (x1(2)-x2(2)) - (rs1+rs2)*sine]  
                            ! print*, "L 395 r", r, "id",id

                            fric_term = (fric_term1+fric_term2)/2
                            if (fric_term > 0) then
                                fR = comp_fric([-(x1(1)-x2(1))+ (rs1+rs2)*coss,-(x1(2)-x2(2))+(rs1+rs2)*sine], &
                                [-(v1(1)-v2(1)),-(v1(2)-v2(2))],fric_term)
                            end if
                            ptr%p%F = aux2 + ptr%p%F +fR
                            ptrn%p%F = -aux2 + ptrn%p%F - fR
                       else if (rs1 > 0 .and. rs2 == 0) then
                            gamma = theta(-pa +ptr%p%n) + atan2(coss,sine)

                            aux2 = - ( (6 * epsil * sigma**6 * (A * sin(alpha * gamma) + 1) * ((B * sin(beta * gamma) + r)**6 - 2 * sigma**6))/(B * sin(beta * gamma) + r)**13 ) * [(x1(1)-x2(1)) - (rs1+rs2)*coss, (x1(2)-x2(2)) - (rs1+rs2)*sine]  

                            aux3 =  (1/r) * (epsil * (A * sin(alpha * gamma) + 1)* ((6 * beta * B * sigma**6 * cos(beta * gamma))/(B * sin(beta * gamma) + r)**7 - (12 * beta * B * sigma**12 * cos(beta * gamma))/(B * sin(beta * gamma) + r)**13) + alpha * A * epsil * cos(alpha * gamma) * (sigma**12/(B * sin(beta * gamma) + r)**12 - sigma**6/(B * sin(beta * gamma) + r)**6))

                            Tor(-pa +ptr%p%n) =  -aux3*(r+rs1+rs2) + Tor(-pa +ptr%p%n)
                            ptr%p%F = aux2 + ptr%p%F 
                            ptrn%p%F = -aux2 + ptrn%p%F - aux3*[sine,-coss] 
                            ! print*, "r1 >0"
                            ! print*, "ptr%p%F" , ptr%p%F, "particula"
                            ! print*, "ptrn%p%F" , ptrn%p%F
                            ! print*, "aux2", aux2, "aux3", aux3
                            ! print*, "-aux3*[sine,-cos]", -aux3*[sine,-coss] 

                            ! read(*,*)
                            ! partlst(ptr%p%n - pa) = 1
                                
                        else if (rs2 > 0 .and. rs1 == 0) then 
                            gamma = theta(-pa +ptrn%p%n) + atan2(coss,sine)

                            aux2 = - ( (6 * epsil * sigma**6 * (A * sin(alpha * gamma) + 1) * ((B * sin(beta * gamma) + r)**6 - 2 * sigma**6))/(B * sin(beta * gamma) + r)**13 ) * [(x1(1)-x2(1)) - (rs1+rs2)*coss, (x1(2)-x2(2)) - (rs1+rs2)*sine]  

                            aux3 =  (1/r) * (epsil * (A * sin(alpha * gamma) + 1)* ((6 * beta * B * sigma**6 * cos(beta * gamma))/(B * sin(beta * gamma) + r)**7 - (12 * beta * B * sigma**12 * cos(beta * gamma))/(B * sin(beta * gamma) + r)**13) + alpha * A * epsil * cos(alpha * gamma) * (sigma**12/(B * sin(beta * gamma) + r)**12 - sigma**6/(B * sin(beta * gamma) + r)**6))

                            
                            Tor(-pa +ptrn%p%n) =  -aux3*(r+rs1+rs2) + Tor(-pa +ptrn%p%n)
                            ptr%p%F = aux2 + ptr%p%F + aux3*[sine,-coss] 
                            ptrn%p%F = -aux2 + ptrn%p%F 
                            ! print*, "r2 > 0"
                            ! print*, "ptr%p%F" , ptr%p%F
                            ! print*, "ptrn%p%F" , ptrn%p%F, "particula"
                            ! print*, "aux2", aux2, "aux3", aux3
                            ! print*, "aux3*[sine,-cos]", aux3*[sine,-coss] 
                            ! read(*,*)
                            partlst(ptrn%p%n - pa) = 1
                        else
                            ! Não vou estudar aqui a interação entre duas partículas, elas vão se repelir como átomos
                            aux2 = -(1/r**2)*(sigma*(1+B)/r)**6*(1-2*(sigma*(1+B)/r)**6)*24*epsil*(1+A)* & 
                            [(x1(1)-x2(1)) - (rs1+rs2)*coss, (x1(2)-x2(2)) - (rs1+rs2)*sine]  
                            ! print*, "L 395 r", r, "id",id
                            ! implementar 24*ep0*(1/(r-rs)**2) * ((sig0/(r-rs))**6 *(1-2* (sig0/(r-rs))**6 )) * (r-rs) + (Np/(2*deltad))*((np.log((deltad*2 + sig0)/(r-rs) ) -1) + (deltad*2+sig0))*heaviside(r-(rs+sig0) , 1)*(1-heaviside(r-(deltad*2+rs+sig0) , 1))
                            ! talvez o fricterm seja interessante em algum momento, vou deixar aqui
                            fric_term = (fric_term1+fric_term2)/2
                            if (fric_term > 0) then
                                fR = comp_fric([-(x1(1)-x2(1))+ (rs1+rs2)*coss,-(x1(2)-x2(2))+(rs1+rs2)*sine], &
                                [-(v1(1)-v2(1)),-(v1(2)-v2(2))],fric_term)
                            end if

                            ptr%p%F = aux2 + ptr%p%F +fR
                            ptrn%p%F = -aux2 + ptrn%p%F - fR
                        end if
                    end if
                    node => list_next(node) ! próxima partícula da célula
                end do
                ! print*, "L 871"
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
                                if (rs1 == 0 .and. rs2 == 0) then
                                    aux2 = -(1/r**2)*(sigma/r)**6*(1-2*(sigma/r)**6)*24*epsil* & 
                                    [(x1(1)-x2(1)) - (rs1+rs2)*coss, (x1(2)-x2(2)) - (rs1+rs2)*sine]  
                                    ! print*, "L 395 r", r, "id",id
    
                                    fric_term = (fric_term1+fric_term2)/2
                                    if (fric_term > 0) then
                                        fR = comp_fric([-(x1(1)-x2(1))+ (rs1+rs2)*coss,-(x1(2)-x2(2))+(rs1+rs2)*sine], &
                                        [-(v1(1)-v2(1)),-(v1(2)-v2(2))],fric_term)
                                    end if
                                    ptr%p%F = aux2 + ptr%p%F +fR
                                    ptrn%p%F = -aux2 + ptrn%p%F - fR
                                else if (rs1 > 0 .and. rs2 == 0) then
    
                                    gamma = theta(-pa +ptr%p%n) + atan2(coss,sine)

                                    aux2 = - ( (6 * epsil * sigma**6 * (A * sin(alpha * gamma) + 1) * ((B * sin(beta * gamma) + r)**6 - 2 * sigma**6))/(B * sin(beta * gamma) + r)**13 ) * [(x1(1)-x2(1)) - (rs1+rs2)*coss, (x1(2)-x2(2)) - (rs1+rs2)*sine]  

                                    aux3 =  (1/r) * (epsil * (A * sin(alpha * gamma) + 1)* ((6 * beta * B * sigma**6 * cos(beta * gamma))/(B * sin(beta * gamma) + r)**7 - (12 * beta * B * sigma**12 * cos(beta * gamma))/(B * sin(beta * gamma) + r)**13) + alpha * A * epsil * cos(alpha * gamma) * (sigma**12/(B * sin(beta * gamma) + r)**12 - sigma**6/(B * sin(beta * gamma) + r)**6))
    
                                    Tor(-pa +ptr%p%n) =  -aux3*(r+rs1+rs2) + Tor(-pa +ptr%p%n)
                                    ptr%p%F = aux2 + ptr%p%F 
                                    ptrn%p%F = -aux2 + ptrn%p%F - aux3*[sine,-coss] 
                                    ! partlst(ptr%p%n - pa) = 1
                                else if (rs2 > 0 .and. rs1 == 0) then 
                                    gamma = theta(-pa +ptrn%p%n) + atan2(coss,sine)

                                    aux2 = - ( (6 * epsil * sigma**6 * (A * sin(alpha * gamma) + 1) * ((B * sin(beta * gamma) + r)**6 - 2 * sigma**6))/(B * sin(beta * gamma) + r)**13 ) * [(x1(1)-x2(1)) - (rs1+rs2)*coss, (x1(2)-x2(2)) - (rs1+rs2)*sine]  

                                    aux3 =  (1/r) * (epsil * (A * sin(alpha * gamma) + 1)* ((6 * beta * B * sigma**6 * cos(beta * gamma))/(B * sin(beta * gamma) + r)**7 - (12 * beta * B * sigma**12 * cos(beta * gamma))/(B * sin(beta * gamma) + r)**13) + alpha * A * epsil * cos(alpha * gamma) * (sigma**12/(B * sin(beta * gamma) + r)**12 - sigma**6/(B * sin(beta * gamma) + r)**6))

                                    Tor(-pa +ptrn%p%n) =  -aux3*(r+rs1+rs2) + Tor(-pa +ptrn%p%n)
                                    ptr%p%F = aux2 + ptr%p%F + aux3*[sine,-coss] 
                                    ptrn%p%F = -aux2 + ptrn%p%F 
                                    partlst(ptrn%p%n - pa) = 1
                                else
                                    ! Não vou estudar aqui a interação entre duas partículas, elas vão se repelir como átomos
                                    aux2 = -(1/r**2)*(sigma*(1+B)/r)**6*(1-2*(sigma*(1+B)/r)**6)*24*epsil*(1+A)* & 
                                    [(x1(1)-x2(1)) - (rs1+rs2)*coss, (x1(2)-x2(2)) - (rs1+rs2)*sine]  
                                    ! print*, "L 395 r", r, "id",id
    
                                    ! talvez o fricterm seja interessante em algum momento, vou deixar aqui
                                    fric_term = (fric_term1+fric_term2)/2
                                    if (fric_term > 0) then
                                        fR = comp_fric([-(x1(1)-x2(1))+ (rs1+rs2)*coss,-(x1(2)-x2(2))+(rs1+rs2)*sine], &
                                        [-(v1(1)-v2(1)),-(v1(2)-v2(2))],fric_term)
                                    end if
    
                                    ptr%p%F = aux2 + ptr%p%F +fR
                                    ptrn%p%F = -aux2 + ptrn%p%F - fR
                                    
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
                                v2 = ptrn%p%v
                                r = sqrt((x1(1)-x2(1))**2 + (x1(2)-x2(2))**2)
                                coss = (x1(1)-x2(1))/r 
                                sine = (x1(2)-x2(2))/r 
                                r = r - rs1 - rs2 !raio
                                ! print*, "L 710 r", r, "id",id
                                if (r <= rcut) then
                                    if (rs1 == 0 .and. rs2 == 0) then
                                        aux2 = -(1/r**2)*(sigma/r)**6*(1-2*(sigma/r)**6)*24*epsil* & 
                                        [(x1(1)-x2(1)) - (rs1+rs2)*coss, (x1(2)-x2(2)) - (rs1+rs2)*sine]  
                                        ! print*, "L 395 r", r, "id",id
        
                                        fric_term = (fric_term1+fric_term2)/2
                                        if (fric_term > 0) then
                                            fR = comp_fric([-(x1(1)-x2(1))+ (rs1+rs2)*coss,-(x1(2)-x2(2))+(rs1+rs2)*sine], &
                                            [-(v1(1)-v2(1)),-(v1(2)-v2(2))],fric_term)
                                        end if
                                        ptr%p%F = aux2 + ptr%p%F +fR
                                        ptrn%p%F = -aux2 + ptrn%p%F - fR
                                else if (rs1 > 0 .and. rs2 == 0) then
        
                                        gamma = theta(-pa +ptr%p%n) + atan2(coss,sine)
                                        aux2 = - ( (6 * epsil * sigma**6 * (A * sin(alpha * gamma) + 1) * ((B * sin(beta * gamma) + r)**6 - 2 * sigma**6))/(B * sin(beta * gamma) + r)**13 ) * [(x1(1)-x2(1)) - (rs1+rs2)*coss, (x1(2)-x2(2)) - (rs1+rs2)*sine]  

                                        aux3 =  (1/r) * (epsil * (A * sin(alpha * gamma) + 1)* ((6 * beta * B * sigma**6 * cos(beta * gamma))/(B * sin(beta * gamma) + r)**7 - (12 * beta * B * sigma**12 * cos(beta * gamma))/(B * sin(beta * gamma) + r)**13) + alpha * A * epsil * cos(alpha * gamma) * (sigma**12/(B * sin(beta * gamma) + r)**12 - sigma**6/(B * sin(beta * gamma) + r)**6))
        
                                        Tor(-pa +ptr%p%n) =  -aux3*(r+rs1+rs2) + Tor(-pa +ptr%p%n)
                                        ptr%p%F = aux2 + ptr%p%F 
                                        ptrn%p%F = -aux2 + ptrn%p%F - aux3*[sine,-coss] 
                                        ! partlst(ptr%p%n - pa) = 1
                                    else if (rs2 > 0 .and. rs1 == 0) then 
                                        gamma = theta(-pa +ptrn%p%n) + atan2(coss,sine)
                                        aux2 = - ( (6 * epsil * sigma**6 * (A * sin(alpha * gamma) + 1) * ((B * sin(beta * gamma) + r)**6 - 2 * sigma**6))/(B * sin(beta * gamma) + r)**13 ) * [(x1(1)-x2(1)) - (rs1+rs2)*coss, (x1(2)-x2(2)) - (rs1+rs2)*sine]  

                                        aux3 =  (1/r) * (epsil * (A * sin(alpha * gamma) + 1)* ((6 * beta * B * sigma**6 * cos(beta * gamma))/(B * sin(beta * gamma) + r)**7 - (12 * beta * B * sigma**12 * cos(beta * gamma))/(B * sin(beta * gamma) + r)**13) + alpha * A * epsil * cos(alpha * gamma) * (sigma**12/(B * sin(beta * gamma) + r)**12 - sigma**6/(B * sin(beta * gamma) + r)**6))
        
                                        Tor(-pa +ptrn%p%n) =  -aux3*(r+rs1+rs2) + Tor(-pa +ptrn%p%n)
                                        ptr%p%F = aux2 + ptr%p%F + aux3*[sine,-coss] 
                                        ptrn%p%F = -aux2 + ptrn%p%F 
                                        partlst(ptrn%p%n - pa) = 1
                                    else
                                        ! Não vou estudar aqui a interação entre duas partículas, elas vão se repelir como átomos
                                        aux2 = -(1/r**2)*(sigma*(1+B)/r)**6*(1-2*(sigma*(1+B)/r)**6)*24*epsil*(1+A)* & 
                                        [(x1(1)-x2(1)) - (rs1+rs2)*coss, (x1(2)-x2(2)) - (rs1+rs2)*sine]  
                                        ! print*, "L 395 r", r, "id",id
        
                                        ! talvez o fricterm seja interessante em algum momento, vou deixar aqui
                                        fric_term = (fric_term1+fric_term2)/2
                                        if (fric_term > 0) then
                                            fR = comp_fric([-(x1(1)-x2(1))+ (rs1+rs2)*coss,-(x1(2)-x2(2))+(rs1+rs2)*sine], &
                                            [-(v1(1)-v2(1)),-(v1(2)-v2(2))],fric_term)
                                        end if
        
                                        ptr%p%F = aux2 + ptr%p%F +fR
                                        ptrn%p%F = -aux2 + ptrn%p%F - fR
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
                            v2 = ptrn%p%v
                            
                            r = sqrt((x1(1)-x2(1))**2 + (x1(2)-x2(2))**2)
                            coss = (x1(1)-x2(1))/r 
                            sine = (x1(2)-x2(2))/r 
                            r = r - rs1 - rs2 !raio
                            ! print*, "L 757 r", r, "id",id
                            if (r <= rcut) then
                                if (rs1 == 0 .and. rs2 == 0) then
                                    aux2 = -(1/r**2)*(sigma/r)**6*(1-2*(sigma/r)**6)*24*epsil* & 
                                    [(x1(1)-x2(1)) - (rs1+rs2)*coss, (x1(2)-x2(2)) - (rs1+rs2)*sine]  
                                    ! print*, "L 395 r", r, "id",id
    
                                    fric_term = (fric_term1+fric_term2)/2
                                    if (fric_term > 0) then
                                        fR = comp_fric([-(x1(1)-x2(1))+ (rs1+rs2)*coss,-(x1(2)-x2(2))+(rs1+rs2)*sine], &
                                        [-(v1(1)-v2(1)),-(v1(2)-v2(2))],fric_term)
                                    end if
                                    ptr%p%F = aux2 + ptr%p%F +fR
                                    ptrn%p%F = -aux2 + ptrn%p%F - fR
                            else if (rs1 > 0 .and. rs2 == 0) then
    
                                    gamma = theta(-pa +ptr%p%n) + atan2(coss,sine)
                                    aux2 = - ( (6 * epsil * sigma**6 * (A * sin(alpha * gamma) + 1) * ((B * sin(beta * gamma) + r)**6 - 2 * sigma**6))/(B * sin(beta * gamma) + r)**13 ) * [(x1(1)-x2(1)) - (rs1+rs2)*coss, (x1(2)-x2(2)) - (rs1+rs2)*sine]  

                                    aux3 =  (1/r) * (epsil * (A * sin(alpha * gamma) + 1)* ((6 * beta * B * sigma**6 * cos(beta * gamma))/(B * sin(beta * gamma) + r)**7 - (12 * beta * B * sigma**12 * cos(beta * gamma))/(B * sin(beta * gamma) + r)**13) + alpha * A * epsil * cos(alpha * gamma) * (sigma**12/(B * sin(beta * gamma) + r)**12 - sigma**6/(B * sin(beta * gamma) + r)**6))
    
                                    Tor(-pa +ptr%p%n) =  -aux3*(r+rs1+rs2) + Tor(-pa +ptr%p%n)
                                    ptr%p%F = aux2 + ptr%p%F 
                                    ptrn%p%F = -aux2 + ptrn%p%F - aux3*[sine,-coss] 
                                    ! partlst(ptr%p%n - pa) = 1
                                else if (rs2 > 0 .and. rs1 == 0) then 
                                    gamma = theta(-pa +ptrn%p%n) + atan2(coss,sine)
                                    aux2 = - ( (6 * epsil * sigma**6 * (A * sin(alpha * gamma) + 1) * ((B * sin(beta * gamma) + r)**6 - 2 * sigma**6))/(B * sin(beta * gamma) + r)**13 ) * [(x1(1)-x2(1)) - (rs1+rs2)*coss, (x1(2)-x2(2)) - (rs1+rs2)*sine]  

                                    aux3 =  (1/r) * (epsil * (A * sin(alpha * gamma) + 1)* ((6 * beta * B * sigma**6 * cos(beta * gamma))/(B * sin(beta * gamma) + r)**7 - (12 * beta * B * sigma**12 * cos(beta * gamma))/(B * sin(beta * gamma) + r)**13) + alpha * A * epsil * cos(alpha * gamma) * (sigma**12/(B * sin(beta * gamma) + r)**12 - sigma**6/(B * sin(beta * gamma) + r)**6))
    
                                    Tor(-pa +ptrn%p%n) =  -aux3*(r+rs1+rs2) + Tor(-pa +ptrn%p%n)
                                    ptr%p%F = aux2 + ptr%p%F + aux3*[sine,-coss] 
                                    ptrn%p%F = -aux2 + ptrn%p%F 
                                    partlst(ptrn%p%n - pa) = 1
                                else
                                    ! Não vou estudar aqui a interação entre duas partículas, elas vão se repelir como átomos
                                    aux2 = -(1/r**2)*(sigma*(1+B)/r)**6*(1-2*(sigma*(1+B)/r)**6)*24*epsil*(1+A)* & 
                                    [(x1(1)-x2(1)) - (rs1+rs2)*coss, (x1(2)-x2(2)) - (rs1+rs2)*sine]  
                                    ! print*, "L 395 r", r, "id",id
    
                                    ! talvez o fricterm seja interessante em algum momento, vou deixar aqui
                                    fric_term = (fric_term1+fric_term2)/2
                                    if (fric_term > 0) then
                                        fR = comp_fric([-(x1(1)-x2(1))+ (rs1+rs2)*coss,-(x1(2)-x2(2))+(rs1+rs2)*sine], &
                                        [-(v1(1)-v2(1)),-(v1(2)-v2(2))],fric_term)
                                    end if
    
                                    ptr%p%F = aux2 + ptr%p%F +fR
                                    ptrn%p%F = -aux2 + ptrn%p%F - fR
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
                            v2 = ptrn%p%v
                            r = sqrt((x1(1)-x2(1))**2 + (x1(2)-x2(2))**2)
                            coss = (x1(1)-x2(1))/r 
                            sine = (x1(2)-x2(2))/r 
                            r = r - rs1 - rs2 !raio
                            ! print*, "L 798 r", r, "id",id
                            if (r <= rcut) then
                                if (rs1 == 0 .and. rs2 == 0) then
                                    aux2 = -(1/r**2)*(sigma/r)**6*(1-2*(sigma/r)**6)*24*epsil* & 
                                    [(x1(1)-x2(1)) - (rs1+rs2)*coss, (x1(2)-x2(2)) - (rs1+rs2)*sine]  
                                    ! print*, "L 395 r", r, "id",id
    
                                    fric_term = (fric_term1+fric_term2)/2
                                    if (fric_term > 0) then
                                        fR = comp_fric([-(x1(1)-x2(1))+ (rs1+rs2)*coss,-(x1(2)-x2(2))+(rs1+rs2)*sine], &
                                        [-(v1(1)-v2(1)),-(v1(2)-v2(2))],fric_term)
                                    end if
                                    ptr%p%F = aux2 + ptr%p%F +fR
                                    ptrn%p%F = -aux2 + ptrn%p%F - fR
                                    
                            else if (rs1 > 0 .and. rs2 == 0) then
    
                                    gamma = theta(-pa +ptr%p%n) + atan2(coss,sine)
                            
                                    aux2 = - ( (6 * epsil * sigma**6 * (A * sin(alpha * gamma) + 1) * ((B * sin(beta * gamma) + r)**6 - 2 * sigma**6))/(B * sin(beta * gamma) + r)**13 ) * [(x1(1)-x2(1)) - (rs1+rs2)*coss, (x1(2)-x2(2)) - (rs1+rs2)*sine]  

                                    aux3 =  (1/r) * (epsil * (A * sin(alpha * gamma) + 1)* ((6 * beta * B * sigma**6 * cos(beta * gamma))/(B * sin(beta * gamma) + r)**7 - (12 * beta * B * sigma**12 * cos(beta * gamma))/(B * sin(beta * gamma) + r)**13) + alpha * A * epsil * cos(alpha * gamma) * (sigma**12/(B * sin(beta * gamma) + r)**12 - sigma**6/(B * sin(beta * gamma) + r)**6))
                            
                                    Tor(-pa +ptr%p%n) =  -aux3*(r+rs1+rs2) + Tor(-pa +ptr%p%n)
                                    ptr%p%F = aux2 + ptr%p%F 
                                    ptrn%p%F = -aux2 + ptrn%p%F - aux3*[sine,-coss] 
                                    ! partlst(ptr%p%n - pa) = 1
                                else if (rs2 > 0 .and. rs1 == 0) then 
                                    gamma = theta(-pa +ptrn%p%n) + atan2(coss,sine)
                                    aux2 = - ( (6 * epsil * sigma**6 * (A * sin(alpha * gamma) + 1) * ((B * sin(beta * gamma) + r)**6 - 2 * sigma**6))/(B * sin(beta * gamma) + r)**13 ) * [(x1(1)-x2(1)) - (rs1+rs2)*coss, (x1(2)-x2(2)) - (rs1+rs2)*sine]  

                                    aux3 =  (1/r) * (epsil * (A * sin(alpha * gamma) + 1)* ((6 * beta * B * sigma**6 * cos(beta * gamma))/(B * sin(beta * gamma) + r)**7 - (12 * beta * B * sigma**12 * cos(beta * gamma))/(B * sin(beta * gamma) + r)**13) + alpha * A * epsil * cos(alpha * gamma) * (sigma**12/(B * sin(beta * gamma) + r)**12 - sigma**6/(B * sin(beta * gamma) + r)**6))
    
                                    Tor(-pa +ptrn%p%n) =  -aux3*(r+rs1+rs2) + Tor(-pa +ptrn%p%n)
                                    ptr%p%F = aux2 + ptr%p%F + aux3*[sine,-coss] 
                                    ptrn%p%F = -aux2 + ptrn%p%F 
                                    partlst(ptrn%p%n - pa) = 1
                                else
                                    ! Não vou estudar aqui a interação entre duas partículas, elas vão se repelir como átomos
                                    aux2 = -(1/r**2)*(sigma*(1+B)/r)**6*(1-2*(sigma*(1+B)/r)**6)*24*epsil*(1+A)* & 
                                    [(x1(1)-x2(1)) - (rs1+rs2)*coss, (x1(2)-x2(2)) - (rs1+rs2)*sine]  
                                    ! print*, "L 395 r", r, "id",id
    
                                    ! talvez o fricterm seja interessante em algum momento, vou deixar aqui
                                    fric_term = (fric_term1+fric_term2)/2
                                    if (fric_term > 0) then
                                        fR = comp_fric([-(x1(1)-x2(1))+ (rs1+rs2)*coss,-(x1(2)-x2(2))+(rs1+rs2)*sine], &
                                        [-(v1(1)-v2(1)),-(v1(2)-v2(2))],fric_term)
                                    end if
    
                                    ptr%p%F = aux2 + ptr%p%F +fR
                                    ptrn%p%F = -aux2 + ptrn%p%F - fR
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
                            v2 = ptrn%p%v
                            r = sqrt((x1(1)-x2(1))**2 + (x1(2)-x2(2))**2)
                            coss = (x1(1)-x2(1))/r 
                            sine = (x1(2)-x2(2))/r 
                            r = r - rs1 - rs2 !raio
                            ! print*, "L 841 r", r, "id",id
                            if (r <= rcut) then
                                if (rs1 == 0 .and. rs2 == 0) then
                                    aux2 = -(1/r**2)*(sigma/r)**6*(1-2*(sigma/r)**6)*24*epsil* & 
                                    [(x1(1)-x2(1)) - (rs1+rs2)*coss, (x1(2)-x2(2)) - (rs1+rs2)*sine]  
                                    ! print*, "L 395 r", r, "id",id
    
                                    fric_term = (fric_term1+fric_term2)/2
                                    if (fric_term > 0) then
                                        fR = comp_fric([-(x1(1)-x2(1))+ (rs1+rs2)*coss,-(x1(2)-x2(2))+(rs1+rs2)*sine], &
                                        [-(v1(1)-v2(1)),-(v1(2)-v2(2))],fric_term)
                                    end if
                                    ptr%p%F = aux2 + ptr%p%F +fR
                                    ptrn%p%F = -aux2 + ptrn%p%F - fR
                                else if (rs1 > 0 .and. rs2 == 0) then
    
                                    gamma = theta(-pa +ptr%p%n) + atan2(coss,sine)
                                    aux2 = - ( (6 * epsil * sigma**6 * (A * sin(alpha * gamma) + 1) * ((B * sin(beta * gamma) + r)**6 - 2 * sigma**6))/(B * sin(beta * gamma) + r)**13 ) * [(x1(1)-x2(1)) - (rs1+rs2)*coss, (x1(2)-x2(2)) - (rs1+rs2)*sine]  

                                    aux3 =  (1/r) * (epsil * (A * sin(alpha * gamma) + 1)* ((6 * beta * B * sigma**6 * cos(beta * gamma))/(B * sin(beta * gamma) + r)**7 - (12 * beta * B * sigma**12 * cos(beta * gamma))/(B * sin(beta * gamma) + r)**13) + alpha * A * epsil * cos(alpha * gamma) * (sigma**12/(B * sin(beta * gamma) + r)**12 - sigma**6/(B * sin(beta * gamma) + r)**6))
    
                                    Tor(-pa +ptr%p%n) =  -aux3*(r+rs1+rs2) + Tor(-pa +ptr%p%n)
                                    ptr%p%F = aux2 + ptr%p%F 
                                    ptrn%p%F = -aux2 + ptrn%p%F - aux3*[sine,-coss] 
                                    ! partlst(ptr%p%n - pa) = 1
                                else if (rs2 > 0 .and. rs1 == 0) then 
                                    gamma = theta(-pa +ptrn%p%n) + atan2(coss,sine)
                                    
                                    aux2 = - ( (6 * epsil * sigma**6 * (A * sin(alpha * gamma) + 1) * ((B * sin(beta * gamma) + r)**6 - 2 * sigma**6))/(B * sin(beta * gamma) + r)**13 ) * [(x1(1)-x2(1)) - (rs1+rs2)*coss, (x1(2)-x2(2)) - (rs1+rs2)*sine]  

                                    aux3 =  (1/r) * (epsil * (A * sin(alpha * gamma) + 1)* ((6 * beta * B * sigma**6 * cos(beta * gamma))/(B * sin(beta * gamma) + r)**7 - (12 * beta * B * sigma**12 * cos(beta * gamma))/(B * sin(beta * gamma) + r)**13) + alpha * A * epsil * cos(alpha * gamma) * (sigma**12/(B * sin(beta * gamma) + r)**12 - sigma**6/(B * sin(beta * gamma) + r)**6))
    
                                    Tor(-pa +ptrn%p%n) =  -aux3*(r+rs1+rs2) + Tor(-pa +ptrn%p%n)
                                    ptr%p%F = aux2 + ptr%p%F + aux3*[sine,-coss] 
                                    ptrn%p%F = -aux2 + ptrn%p%F 
                                    partlst(ptrn%p%n - pa) = 1
                                else
                                    ! Não vou estudar aqui a interação entre duas partículas, elas vão se repelir como átomos
                                    aux2 = -(1/r**2)*(sigma*(1+B)/r)**6*(1-2*(sigma*(1+B)/r)**6)*24*epsil*(1+A)* & 
                                    [(x1(1)-x2(1)) - (rs1+rs2)*coss, (x1(2)-x2(2)) - (rs1+rs2)*sine]  
                                    ! print*, "L 395 r", r, "id",id
    
                                    ! talvez o fricterm seja interessante em algum momento, vou deixar aqui
                                    fric_term = (fric_term1+fric_term2)/2
                                    if (fric_term > 0) then
                                        fR = comp_fric([-(x1(1)-x2(1))+ (rs1+rs2)*coss,-(x1(2)-x2(2))+(rs1+rs2)*sine], &
                                        [-(v1(1)-v2(1)),-(v1(2)-v2(2))],fric_term)
                                    end if
    
                                    ptr%p%F = aux2 + ptr%p%F +fR
                                    ptrn%p%F = -aux2 + ptrn%p%F - fR
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
                                v2 = ptrn%p%v
                                r = sqrt((x1(1)-x2(1))**2 + (x1(2)-x2(2))**2)
                                coss = (x1(1)-x2(1))/r 
                                sine = (x1(2)-x2(2))/r 
                                r = r - rs1 - rs2 !raio
                                ! print*, "L 885 r", r, "id",id
                                if (r <= rcut) then
                                    if (rs1 == 0 .and. rs2 == 0) then
                                        aux2 = -(1/r**2)*(sigma/r)**6*(1-2*(sigma/r)**6)*24*epsil* & 
                                        [(x1(1)-x2(1)) - (rs1+rs2)*coss, (x1(2)-x2(2)) - (rs1+rs2)*sine]  
                                        ! print*, "L 395 r", r, "id",id
        
                                        fric_term = (fric_term1+fric_term2)/2
                                        if (fric_term > 0) then
                                            fR = comp_fric([-(x1(1)-x2(1))+ (rs1+rs2)*coss,-(x1(2)-x2(2))+(rs1+rs2)*sine], &
                                            [-(v1(1)-v2(1)),-(v1(2)-v2(2))],fric_term)
                                        end if
                                        ptr%p%F = aux2 + ptr%p%F +fR
                                        ptrn%p%F = -aux2 + ptrn%p%F - fR
                                else if (rs1 > 0 .and. rs2 == 0) then
        
                                        gamma = theta(-pa +ptr%p%n) + atan2(coss,sine)
                                        
                                        aux2 = - ( (6 * epsil * sigma**6 * (A * sin(alpha * gamma) + 1) * ((B * sin(beta * gamma) + r)**6 - 2 * sigma**6))/(B * sin(beta * gamma) + r)**13 ) * [(x1(1)-x2(1)) - (rs1+rs2)*coss, (x1(2)-x2(2)) - (rs1+rs2)*sine]  

                                        aux3 =  (1/r) * (epsil * (A * sin(alpha * gamma) + 1)* ((6 * beta * B * sigma**6 * cos(beta * gamma))/(B * sin(beta * gamma) + r)**7 - (12 * beta * B * sigma**12 * cos(beta * gamma))/(B * sin(beta * gamma) + r)**13) + alpha * A * epsil * cos(alpha * gamma) * (sigma**12/(B * sin(beta * gamma) + r)**12 - sigma**6/(B * sin(beta * gamma) + r)**6))
        
                                        Tor(-pa +ptr%p%n) =  -aux3*(r+rs1+rs2) + Tor(-pa +ptr%p%n)
                                        ptr%p%F = aux2 + ptr%p%F 
                                        ptrn%p%F = -aux2 + ptrn%p%F - aux3*[sine,-coss] 
                                        ! partlst(ptr%p%n - pa) = 1
                                    else if (rs2 > 0 .and. rs1 == 0) then 
                                        gamma = theta(-pa +ptrn%p%n) + atan2(coss,sine)

                                        aux2 = - ( (6 * epsil * sigma**6 * (A * sin(alpha * gamma) + 1) * ((B * sin(beta * gamma) + r)**6 - 2 * sigma**6))/(B * sin(beta * gamma) + r)**13 ) * [(x1(1)-x2(1)) - (rs1+rs2)*coss, (x1(2)-x2(2)) - (rs1+rs2)*sine]  

                                        aux3 =  (1/r) * (epsil * (A * sin(alpha * gamma) + 1)* ((6 * beta * B * sigma**6 * cos(beta * gamma))/(B * sin(beta * gamma) + r)**7 - (12 * beta * B * sigma**12 * cos(beta * gamma))/(B * sin(beta * gamma) + r)**13) + alpha * A * epsil * cos(alpha * gamma) * (sigma**12/(B * sin(beta * gamma) + r)**12 - sigma**6/(B * sin(beta * gamma) + r)**6))
        
                                        Tor(-pa +ptrn%p%n) =  -aux3*(r+rs1+rs2) + Tor(-pa +ptrn%p%n)
                                        ptr%p%F = aux2 + ptr%p%F + aux3*[sine,-coss] 
                                        ptrn%p%F = -aux2 + ptrn%p%F 
                                        partlst(ptrn%p%n - pa) = 1
                                    else
                                        ! Não vou estudar aqui a interação entre duas partículas, elas vão se repelir como átomos
                                        aux2 = -(1/r**2)*(sigma*(1+B)/r)**6*(1-2*(sigma*(1+B)/r)**6)*24*epsil*(1+A)* & 
                                        [(x1(1)-x2(1)) - (rs1+rs2)*coss, (x1(2)-x2(2)) - (rs1+rs2)*sine]  
                                        ! print*, "L 395 r", r, "id",id
        
                                        ! talvez o fricterm seja interessante em algum momento, vou deixar aqui
                                        fric_term = (fric_term1+fric_term2)/2
                                        if (fric_term > 0) then
                                            fR = comp_fric([-(x1(1)-x2(1))+ (rs1+rs2)*coss,-(x1(2)-x2(2))+(rs1+rs2)*sine], &
                                            [-(v1(1)-v2(1)),-(v1(2)-v2(2))],fric_term)
                                        end if
        
                                        ptr%p%F = aux2 + ptr%p%F +fR
                                        ptrn%p%F = -aux2 + ptrn%p%F - fR
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
                        v2 = ptrn%p%v
                        r = sqrt((x1(1)-x2(1))**2 + (x1(2)-x2(2))**2)
                        coss = (x1(1)-x2(1))/r 
                        sine = (x1(2)-x2(2))/r 
                        r = r - rs1 - rs2 !raio
                        ! print*, "L 931 r", r, "id",id
                        if (r <= rcut) then
                            if (rs1 == 0 .and. rs2 == 0) then
                                 aux2 = -(1/r**2)*(sigma/r)**6*(1-2*(sigma/r)**6)*24*epsil* & 
                                 [(x1(1)-x2(1)) - (rs1+rs2)*coss, (x1(2)-x2(2)) - (rs1+rs2)*sine]  
                                 ! print*, "L 395 r", r, "id",id
 
                                 fric_term = (fric_term1+fric_term2)/2
                                 if (fric_term > 0) then
                                     fR = comp_fric([-(x1(1)-x2(1))+ (rs1+rs2)*coss,-(x1(2)-x2(2))+(rs1+rs2)*sine], &
                                     [-(v1(1)-v2(1)),-(v1(2)-v2(2))],fric_term)
                                 end if
                                 ptr%p%F = aux2 + ptr%p%F +fR
                                 ptrn%p%F = -aux2 + ptrn%p%F - fR
                            else if (rs1 > 0 .and. rs2 == 0) then
                                 gamma = theta(-pa +ptr%p%n) + atan2(coss,sine)
 
                                aux2 = - ( (6 * epsil * sigma**6 * (A * sin(alpha * gamma) + 1) * ((B * sin(beta * gamma) + r)**6 - 2 * sigma**6))/(B * sin(beta * gamma) + r)**13 ) * [(x1(1)-x2(1)) - (rs1+rs2)*coss, (x1(2)-x2(2)) - (rs1+rs2)*sine]  

                                aux3 =  (1/r) * (epsil * (A * sin(alpha * gamma) + 1)* ((6 * beta * B * sigma**6 * cos(beta * gamma))/(B * sin(beta * gamma) + r)**7 - (12 * beta * B * sigma**12 * cos(beta * gamma))/(B * sin(beta * gamma) + r)**13) + alpha * A * epsil * cos(alpha * gamma) * (sigma**12/(B * sin(beta * gamma) + r)**12 - sigma**6/(B * sin(beta * gamma) + r)**6))
 
                                 ! print*, "r1 > 0 gamma", gamma
                                 ! print*, "aux2, aux3",aux2, aux3
                                 ! print*, "x1 =", x1
                                 ! print*, "x2 =", x2
                                 ! read(*,*)
                                Tor(-pa +ptr%p%n) =  -aux3*(r+rs1+rs2) + Tor(-pa +ptr%p%n)
                                ptr%p%F = aux2 + ptr%p%F 
                                ptrn%p%F = -aux2 + ptrn%p%F - aux3*[sine,-coss] 
                                 ! partlst(ptr%p%n - pa) = 1
 
                            else if (rs2 > 0 .and. rs1 == 0) then 
                                gamma = theta(-pa +ptrn%p%n) + atan2(coss,sine)

                                aux2 = - ( (6 * epsil * sigma**6 * (A * sin(alpha * gamma) + 1) * ((B * sin(beta * gamma) + r)**6 - 2 * sigma**6))/(B * sin(beta * gamma) + r)**13 ) * [(x1(1)-x2(1)) - (rs1+rs2)*coss, (x1(2)-x2(2)) - (rs1+rs2)*sine]  

                                aux3 =  (1/r) * (epsil * (A * sin(alpha * gamma) + 1)* ((6 * beta * B * sigma**6 * cos(beta * gamma))/(B * sin(beta * gamma) + r)**7 - (12 * beta * B * sigma**12 * cos(beta * gamma))/(B * sin(beta * gamma) + r)**13) + alpha * A * epsil * cos(alpha * gamma) * (sigma**12/(B * sin(beta * gamma) + r)**12 - sigma**6/(B * sin(beta * gamma) + r)**6))

                                Tor(-pa +ptrn%p%n) =  -aux3*(r+rs1+rs2) + Tor(-pa +ptrn%p%n)
                                ptr%p%F = aux2 + ptr%p%F + aux3*[sine,-coss] 
                                ptrn%p%F = -aux2 + ptrn%p%F 
                                partlst(ptrn%p%n - pa) = 1
                            else
                                ! Não vou estudar aqui a interação entre duas partículas, elas vão se repelir como átomos
                                aux2 = -(1/r**2)*(sigma*(1+B)/r)**6*(1-2*(sigma*(1+B)/r)**6)*24*epsil*(1+A)* & 
                                [(x1(1)-x2(1)) - (rs1+rs2)*coss, (x1(2)-x2(2)) - (rs1+rs2)*sine]  
                                ! print*, "L 395 r", r, "id",id

                                ! talvez o fricterm seja interessante em algum momento, vou deixar aqui
                                fric_term = (fric_term1+fric_term2)/2
                                if (fric_term > 0) then
                                    fR = comp_fric([-(x1(1)-x2(1))+ (rs1+rs2)*coss,-(x1(2)-x2(2))+(rs1+rs2)*sine], &
                                    [-(v1(1)-v2(1)),-(v1(2)-v2(2))],fric_term)
                                end if

                                ptr%p%F = aux2 + ptr%p%F +fR
                                ptrn%p%F = -aux2 + ptrn%p%F - fR
                            end if
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
end subroutine comp_FT    


! subroutine comp_FT(GField,Hfield,theta,Tor,pr,mesh,malha,propriedade,r_cut,domx,domy, ids, id, t,dt,partlst)
!     use linkedlist
!     use mod1
!     use data
!     use mod0
!     use mpi
    
!     type(prop_grupo), allocatable,dimension(:),intent(in) :: propriedade
!     type(container), allocatable,dimension(:,:),intent(in) :: malha
!     real(dp), intent(in) :: dt, GField(2),hfield(5), pr(5), t ! outras propriedades da particula que gira [A, B, alpha, beta, phase]
!     real(dp), intent(inout), dimension(:) :: Tor ! torque
!     real(dp), intent(in), dimension(:) :: theta ! ângulo das partículas que giram, tangencial
!     real(dp) :: sigma, epsil, sigma_a, epsil_a,sigma_b, epsil_b, rcut,r_cut, fric_term !fric_term = força de ficção
!     real(dp) :: x1(2),v1(2),x2(2),v2(2), rs1, rs2, coss, sine, A, B, alpha,beta,ph   
!     integer :: i,j, ct = 0, dox(2), doy(2) !,ptr, ptrn
!     integer, intent(in) :: mesh(:),domx(2),domy(2), id
!     real(dp) :: Fi(2)=0,r, aux2(2),aux3,fR(2), fric_term1, fric_term2
!     integer, save :: pa
!     type(list_t), pointer :: node, next_node
!     type(data_ptr) :: ptr,ptrn
!     integer, intent(out), dimension(:) :: partlst
!     integer ( kind = 4 ), intent(in) :: ids(8)

!     partlst = 0
!     A =     pr(1)
!     B =     pr(2)
!     alpha = pr(3)
!     beta =  pr(4)
!     ph =    pr(5)

!     ! quanta a quantidade de partículas antes das que giram
!     ! print*, "t,2*dt", t, 2*dt
!     if (t < 2*dt) then
!         rs = 0
!         i = 0
!         pa = 0
!         do while (rs == 0)
!             i = i+1
!             rs = propriedade(i)%rs
!             pa = pa + propriedade(i)%quant
!         end do
!         pa = pa - propriedade(i)%quant
!         ! print*, "pa calculado!",id,pa
!     end if 
!     ! print*, "pa = ", pa
    

!     !Lennard Jones
!     fR = [0,0]
!     dox = domx 
!     doy = domy 
!     if (sum(ids(1:4)) > -4) then
!         !caso paralelo
!         if (domx(1) > 1) dox(1) = domx(1) - 1
!         if (domy(1) > 1) doy(1) = domy(1) - 1
!         if (domx(2) < mesh(1)+2) dox(2) = domx(2) + 1
!         if (domy(2) < mesh(2)+2) doy(2) = domy(2) + 1
!     else
!         dox = domx 
!         doy = domy 
!     end if 
    
!     do i = doy(1),doy(2) ! i é linha
!         do j = dox(1),dox(2)
!             node => list_next(malha(i,j)%list)
!             ! if (associated(node)) then
!             !     ptr = transfer(list_get(node), ptr)
!             !     ! print*, "L 253", ptr%p%n
!             ! end if
!             do while (associated(node))
!                 ! print*, "Encontrada", ct, "celulas, id=", id
!                 ptr = transfer(list_get(node), ptr) !particula selecionada
!                 ! PRINT*, "L 759"
!                 ptr%p%flag = .true. ! indica ao comp_x que a partícula precisa ser calculada
!                 ! PRINT*, "L 764"
!                 x1 = ptr%p%x
!                 v1 = ptr%p%v
!                 m1 = propriedade(ptr%p%grupo)%m
!                 rs1 = propriedade(ptr%p%grupo)%rs !raio sólido 
!                 fric_term1 = propriedade(ptr%p%grupo)%fric_term
!                 sigma_a = propriedade(ptr%p%grupo)%sigma
!                 ! rcut = r_cut*sigma
!                 epsil_a = propriedade(ptr%p%grupo)%epsilon 
!                 if (propriedade(ptr%p%grupo)%x_lockdelay <= t) then
!                     ptr%p%F = ptr%p%F + GField*m1
!                 end if

!                 if (rs1 > 0) then
!                     partlst(ptr%p%n - pa) = 1
!                     Tor(-pa +ptr%p%n) = Tor(-pa +ptr%p%n) + magtorque(theta(-pa+ptr%p%n),hfield,t)
!                 end if
!                 ! print '("x1  =", f10.6, " ", f10.6, " n ", i2, " cell ",i2," ",i2)',x1(1),x1(2), ptr%p%n,i,j
!                 ! if (id == 0) read(*,*)
!                 !calcular a força desta com todas as outras partículas
!                 next_node => list_next(node) ! próxima partícula da célula
!                 node => list_next(node) ! a ser computado com a particula selecionada
         
!                 ! PRINT*, "L 780"
!                 ! NA PRÓPRIA CELULA 
!                 do while (associated(node))
!                     ptrn = transfer(list_get(node), ptrn) !outra particula selecionada
!                     sigma_b = propriedade(ptrn%p%grupo)%sigma
!                     epsil_b = propriedade(ptrn%p%grupo)%epsilon 
!                     rs2 = propriedade(ptrn%p%grupo)%rs 
!                     fric_term2 = propriedade(ptrn%p%grupo)%fric_term
!                     ! Lorenz-Betherlot rule for mixing epsilon sigma 
!                     if (sigma_a > sigma_b) then
!                         rcut = r_cut*sigma_a + rs1 + rs2
!                         sigma = 0.5*(sigma_a + sigma_b)
!                         epsil = sqrt(epsil_a *epsil_b )
!                     else if (sigma_a < sigma_b) then 
!                         rcut = r_cut*sigma_b + rs1 + rs2
!                         sigma = 0.5*(sigma_a + sigma_b)
!                         epsil = sqrt(epsil_a *epsil_b )
!                     else 
!                         rcut = r_cut*sigma_a + rs1 + rs2
!                         sigma = sigma_a
!                         epsil = epsil_a
!                     end if 
!                     x2 = ptrn%p%x
!                     v2 = ptrn%p%v
!                     r = sqrt((x1(1)-x2(1))**2 + (x1(2)-x2(2))**2)
!                     coss = (x1(1)-x2(1))/r 
!                     sine = (x1(2)-x2(2))/r 
!                     r = r - rs1 - rs2 !raio
!                     ! print*, "L 389 r", r, "id",id
!                     if (r <= rcut) then
!                        if (rs1 == 0 .and. rs2 == 0) then
!                             aux2 = -(1/r**2)*(sigma/r)**6*(1-2*(sigma/r)**6)*24*epsil* & 
!                             [(x1(1)-x2(1)) - (rs1+rs2)*coss, (x1(2)-x2(2)) - (rs1+rs2)*sine]  
!                             ! print*, "L 395 r", r, "id",id

!                             fric_term = (fric_term1+fric_term2)/2
!                             if (fric_term > 0) then
!                                 fR = comp_fric([-(x1(1)-x2(1))+ (rs1+rs2)*coss,-(x1(2)-x2(2))+(rs1+rs2)*sine], &
!                                 [-(v1(1)-v2(1)),-(v1(2)-v2(2))],fric_term)
!                             end if
!                             ptr%p%F = aux2 + ptr%p%F +fR
!                             ptrn%p%F = -aux2 + ptrn%p%F - fR
!                        else if (rs1 > 0 .and. rs2 == 0) then
!                             gamma = theta(-pa +ptr%p%n) + atan2(coss,sine)

!                             aux2 = -24*epsil*(1+A*sin(alpha*gamma+ph))*(1/r**2) * & 
!                                 ( (sigma*(1+B*sin(beta*gamma))/r)**6 * (1-2* (sigma*(1+B*sin(beta*gamma))/r)**6 )) * &
!                                 [(x1(1)-x2(1)) - (rs1+rs2)*coss, (x1(2)-x2(2)) - (rs1+rs2)*sine]  
!                             aux3 = (1/r)* ( (4*epsil*(A*sin(alpha*gamma+ph) + 1 ))* (sigma/r)**6*(B*sin(beta*gamma)+1)**5 * &
!                                 ((sigma/r)**6 * 12*beta*B*cos(beta*gamma)*(B*sin(beta*gamma)+1)**6 - 6*beta*B*cos(beta*gamma)) + &
!                                 4*alpha*A*epsil*cos(alpha*gamma+ph)*(sigma/r)**6*(B*sin(beta*gamma)+1)**6 * &
!                                 ((sigma/r)**6 * (B*sin(beta*gamma)+1)**6 - 1)) 


!                             Tor(-pa +ptr%p%n) =  -aux3*(r+rs1+rs2) + Tor(-pa +ptr%p%n)
!                             ptr%p%F = aux2 + ptr%p%F 
!                             ptrn%p%F = -aux2 + ptrn%p%F - aux3*[sine,-coss] 
!                             ! print*, "r1 >0"
!                             ! print*, "ptr%p%F" , ptr%p%F, "particula"
!                             ! print*, "ptrn%p%F" , ptrn%p%F
!                             ! print*, "aux2", aux2, "aux3", aux3
!                             ! print*, "-aux3*[sine,-cos]", -aux3*[sine,-coss] 

!                             ! read(*,*)
!                             ! partlst(ptr%p%n - pa) = 1
                                
!                         else if (rs2 > 0 .and. rs1 == 0) then 
!                             gamma = theta(-pa +ptrn%p%n) + atan2(coss,sine)
!                             aux2 = -24*epsil*(1+A*sin(alpha*gamma+ph))*(1/r**2) * ((sigma*(1+B*sin(beta*gamma))/r)**6 *  &
!                                 (1-2* (sigma*(1+B*sin(beta*gamma))/r)**6 )) * &
!                                 [(x1(1)-x2(1)) - (rs1+rs2)*coss, (x1(2)-x2(2)) - (rs1+rs2)*sine]  
!                             aux3 = (1/r)* ( (4*epsil*(A*sin(alpha*gamma+ph) + 1 ))* (sigma/r)**6*(B*sin(beta*gamma)+1)**5 * &
!                                 ((sigma/r)**6 * 12*beta*B*cos(beta*gamma)*(B*sin(beta*gamma)+1)**6 - 6*beta*B*cos(beta*gamma)) + &
!                                 4*alpha*A*epsil*cos(alpha*gamma+ph)*(sigma/r)**6*(B*sin(beta*gamma)+1)**6 * &
!                                 ((sigma/r)**6 * (B*sin(beta*gamma)+1)**6 - 1)) 

                            
!                             Tor(-pa +ptrn%p%n) =  -aux3*(r+rs1+rs2) + Tor(-pa +ptrn%p%n)
!                             ptr%p%F = aux2 + ptr%p%F + aux3*[sine,-coss] 
!                             ptrn%p%F = -aux2 + ptrn%p%F 
!                             ! print*, "r2 > 0"
!                             ! print*, "ptr%p%F" , ptr%p%F
!                             ! print*, "ptrn%p%F" , ptrn%p%F, "particula"
!                             ! print*, "aux2", aux2, "aux3", aux3
!                             ! print*, "aux3*[sine,-cos]", aux3*[sine,-coss] 
!                             ! read(*,*)
!                             partlst(ptrn%p%n - pa) = 1
!                         else
!                             ! Não vou estudar aqui a interação entre duas partículas, elas vão se repelir como átomos
!                             aux2 = -(1/r**2)*(sigma*(1+B)/r)**6*(1-2*(sigma*(1+B)/r)**6)*24*epsil*(1+A)* & 
!                             [(x1(1)-x2(1)) - (rs1+rs2)*coss, (x1(2)-x2(2)) - (rs1+rs2)*sine]  
!                             ! print*, "L 395 r", r, "id",id
!                             ! implementar 24*ep0*(1/(r-rs)**2) * ((sig0/(r-rs))**6 *(1-2* (sig0/(r-rs))**6 )) * (r-rs) + (Np/(2*deltad))*((np.log((deltad*2 + sig0)/(r-rs) ) -1) + (deltad*2+sig0))*heaviside(r-(rs+sig0) , 1)*(1-heaviside(r-(deltad*2+rs+sig0) , 1))
!                             ! talvez o fricterm seja interessante em algum momento, vou deixar aqui
!                             fric_term = (fric_term1+fric_term2)/2
!                             if (fric_term > 0) then
!                                 fR = comp_fric([-(x1(1)-x2(1))+ (rs1+rs2)*coss,-(x1(2)-x2(2))+(rs1+rs2)*sine], &
!                                 [-(v1(1)-v2(1)),-(v1(2)-v2(2))],fric_term)
!                             end if

!                             ptr%p%F = aux2 + ptr%p%F +fR
!                             ptrn%p%F = -aux2 + ptrn%p%F - fR
!                         end if
!                     end if
!                     node => list_next(node) ! próxima partícula da célula
!                 end do
!                 ! print*, "L 871"
!                 !Células ao redor  !i é linha, j é coluna
!                 if (i /= mesh(2)+2) then ! se não for a última linha
!                     if (j == mesh(1)+2) then ! se for a última coluna
!                         node => list_next(malha(i+1,j)%list) !interagirá com a próxima linha apenas
!                         do while (associated(node))
!                             ptrn = transfer(list_get(node), ptrn) !outra particula selecionada
!                             sigma_b = propriedade(ptrn%p%grupo)%sigma
!                             epsil_b = propriedade(ptrn%p%grupo)%epsilon 
!                             rs2 = propriedade(ptrn%p%grupo)%rs 
!                             fric_term2 = propriedade(ptrn%p%grupo)%fric_term
!                             ! Lorenz-Betherlot rule for mixing epsilon sigma 
!                             if (sigma_a > sigma_b) then
!                                 rcut = r_cut*sigma_a + (rs1 + rs2)
!                                 sigma = 0.5*(sigma_a + sigma_b)
!                                 epsil = sqrt(epsil_a*epsil_b)
!                             else if (sigma_a < sigma_b) then 
!                                 rcut = r_cut*sigma_b + (rs1 + rs2)
!                                 sigma = 0.5*(sigma_a + sigma_b)
!                                 epsil = sqrt(epsil_a*epsil_b)
!                             else 
!                                 rcut = r_cut*sigma_a + (rs1 + rs2)
!                                 sigma = sigma_a
!                                 epsil = epsil_a
!                             end if 

!                             x2 = ptrn%p%x
!                             v2 = ptrn%p%v
!                             r = sqrt((x1(1)-x2(1))**2 + (x1(2)-x2(2))**2)
!                             coss = (x1(1)-x2(1))/r 
!                             sine = (x1(2)-x2(2))/r 
!                             r = r - rs1 - rs2 !raio
!                             ! print*, "L 666 r", r, "id",id
!                             if (r <= rcut) then
!                                 if (rs1 == 0 .and. rs2 == 0) then
!                                     aux2 = -(1/r**2)*(sigma/r)**6*(1-2*(sigma/r)**6)*24*epsil* & 
!                                     [(x1(1)-x2(1)) - (rs1+rs2)*coss, (x1(2)-x2(2)) - (rs1+rs2)*sine]  
!                                     ! print*, "L 395 r", r, "id",id
    
!                                     fric_term = (fric_term1+fric_term2)/2
!                                     if (fric_term > 0) then
!                                         fR = comp_fric([-(x1(1)-x2(1))+ (rs1+rs2)*coss,-(x1(2)-x2(2))+(rs1+rs2)*sine], &
!                                         [-(v1(1)-v2(1)),-(v1(2)-v2(2))],fric_term)
!                                     end if
!                                     ptr%p%F = aux2 + ptr%p%F +fR
!                                     ptrn%p%F = -aux2 + ptrn%p%F - fR
!                                 else if (rs1 > 0 .and. rs2 == 0) then
    
!                                     gamma = theta(-pa +ptr%p%n) + atan2(coss,sine)
!                                     aux2 = -24*epsil*(1+A*sin(alpha*gamma+ph))*(1/r**2) *  ( (sigma*(1+B*sin(beta*gamma))/r)**6 * &
!                                         (1-2* (sigma*(1+B*sin(beta*gamma))/r)**6 )) * & 
!                                         [(x1(1)-x2(1)) - (rs1+rs2)*coss, (x1(2)-x2(2)) - (rs1+rs2)*sine]  
!                                     aux3 = (1/r)* ( (4*epsil*(A*sin(alpha*gamma+ph) + 1))*(sigma/r)**6*(B*sin(beta*gamma)+1)**5 * &
!                                         ((sigma/r)**6 * 12*beta*B*cos(beta*gamma)* &
!                                         (B*sin(beta*gamma)+1)**6 - 6*beta*B*cos(beta*gamma)) + &
!                                         4*alpha*A*epsil*cos(alpha*gamma+ph)*(sigma/r)**6*(B*sin(beta*gamma)+1)**6 * &
!                                         ((sigma/r)**6 * (B*sin(beta*gamma)+1)**6 - 1)) 
    
!                                     Tor(-pa +ptr%p%n) =  -aux3*(r+rs1+rs2) + Tor(-pa +ptr%p%n)
!                                     ptr%p%F = aux2 + ptr%p%F 
!                                     ptrn%p%F = -aux2 + ptrn%p%F - aux3*[sine,-coss] 
!                                     ! partlst(ptr%p%n - pa) = 1
!                                 else if (rs2 > 0 .and. rs1 == 0) then 
!                                     gamma = theta(-pa +ptrn%p%n) + atan2(coss,sine)
!                                     aux2 = -24*epsil*(1+A*sin(alpha*gamma+ph))*(1/r**2) *  ( (sigma*(1+B*sin(beta*gamma))/r)**6 &
!                                     * (1-2* (sigma*(1+B*sin(beta*gamma))/r)**6 )) * & 
!                                     [(x1(1)-x2(1)) - (rs1+rs2)*coss, (x1(2)-x2(2)) - (rs1+rs2)*sine]  
!                                     aux3 = (1/r)* ( (4*epsil*(A*sin(alpha*gamma+ph) + 1 ))* (sigma/r)**6*(B*sin(beta*gamma)+1)**5 * &
!                                         ((sigma/r)**6 * 12*beta*B*cos(beta*gamma)*(B*sin(beta*gamma)+1)**6 - 6*beta*B*cos(beta*gamma)) + &
!                                         4*alpha*A*epsil*cos(alpha*gamma+ph)*(sigma/r)**6*(B*sin(beta*gamma)+1)**6 * &
!                                         ((sigma/r)**6 * (B*sin(beta*gamma)+1)**6 - 1)) 
    
!                                     Tor(-pa +ptrn%p%n) =  -aux3*(r+rs1+rs2) + Tor(-pa +ptrn%p%n)
!                                     ptr%p%F = aux2 + ptr%p%F + aux3*[sine,-coss] 
!                                     ptrn%p%F = -aux2 + ptrn%p%F 
!                                     partlst(ptrn%p%n - pa) = 1
!                                 else
!                                     ! Não vou estudar aqui a interação entre duas partículas, elas vão se repelir como átomos
!                                     aux2 = -(1/r**2)*(sigma*(1+B)/r)**6*(1-2*(sigma*(1+B)/r)**6)*24*epsil*(1+A)* & 
!                                     [(x1(1)-x2(1)) - (rs1+rs2)*coss, (x1(2)-x2(2)) - (rs1+rs2)*sine]  
!                                     ! print*, "L 395 r", r, "id",id
    
!                                     ! talvez o fricterm seja interessante em algum momento, vou deixar aqui
!                                     fric_term = (fric_term1+fric_term2)/2
!                                     if (fric_term > 0) then
!                                         fR = comp_fric([-(x1(1)-x2(1))+ (rs1+rs2)*coss,-(x1(2)-x2(2))+(rs1+rs2)*sine], &
!                                         [-(v1(1)-v2(1)),-(v1(2)-v2(2))],fric_term)
!                                     end if
    
!                                     ptr%p%F = aux2 + ptr%p%F +fR
!                                     ptrn%p%F = -aux2 + ptrn%p%F - fR
                                    
!                                 end if
!                             end if
!                             node => list_next(node) ! próxima partícula da célula
!                         end do
                        
!                         if (j /= 1) then !se não for a primeira coluna 
!                             node => list_next(malha(i+1,j-1)%list) !interagirá com a próxima linha apenas
!                             do while (associated(node))
!                                 ptrn = transfer(list_get(node), ptrn) !outra particula selecionada
!                                 sigma_b = propriedade(ptrn%p%grupo)%sigma
!                                 epsil_b = propriedade(ptrn%p%grupo)%epsilon 
!                                 rs2 = propriedade(ptrn%p%grupo)%rs 
!                                 fric_term2 = propriedade(ptrn%p%grupo)%fric_term
!                                 ! Lorenz-Betherlot rule for mixing epsilon sigma 
!                                 if (sigma_a > sigma_b) then
!                                     rcut = r_cut*sigma_a + (rs1 + rs2)
!                                     sigma = 0.5*(sigma_a + sigma_b)
!                                     epsil = sqrt(epsil_a*epsil_b)
!                                 else if (sigma_a < sigma_b) then 
!                                     rcut = r_cut*sigma_b + (rs1 + rs2)
!                                     sigma = 0.5*(sigma_a + sigma_b)
!                                     epsil = sqrt(epsil_a*epsil_b)
!                                 else 
!                                     rcut = r_cut*sigma_a + (rs1 + rs2)
!                                     sigma = sigma_a
!                                     epsil = epsil_a
!                                 end if 

!                                 x2 = ptrn%p%x
!                                 v2 = ptrn%p%v
!                                 r = sqrt((x1(1)-x2(1))**2 + (x1(2)-x2(2))**2)
!                                 coss = (x1(1)-x2(1))/r 
!                                 sine = (x1(2)-x2(2))/r 
!                                 r = r - rs1 - rs2 !raio
!                                 ! print*, "L 710 r", r, "id",id
!                                 if (r <= rcut) then
!                                     if (rs1 == 0 .and. rs2 == 0) then
!                                         aux2 = -(1/r**2)*(sigma/r)**6*(1-2*(sigma/r)**6)*24*epsil* & 
!                                         [(x1(1)-x2(1)) - (rs1+rs2)*coss, (x1(2)-x2(2)) - (rs1+rs2)*sine]  
!                                         ! print*, "L 395 r", r, "id",id
        
!                                         fric_term = (fric_term1+fric_term2)/2
!                                         if (fric_term > 0) then
!                                             fR = comp_fric([-(x1(1)-x2(1))+ (rs1+rs2)*coss,-(x1(2)-x2(2))+(rs1+rs2)*sine], &
!                                             [-(v1(1)-v2(1)),-(v1(2)-v2(2))],fric_term)
!                                         end if
!                                         ptr%p%F = aux2 + ptr%p%F +fR
!                                         ptrn%p%F = -aux2 + ptrn%p%F - fR
!                                 else if (rs1 > 0 .and. rs2 == 0) then
        
!                                         gamma = theta(-pa +ptr%p%n) + atan2(coss,sine)
!                                         aux2 = -24*epsil*(1+A*sin(alpha*gamma+ph))*(1/r**2) * ((sigma*(1+B*sin(beta*gamma))/r)**6 &
!                                             * (1-2* (sigma*(1+B*sin(beta*gamma))/r)**6 )) * &
!                                             [(x1(1)-x2(1)) - (rs1+rs2)*coss, (x1(2)-x2(2)) - (rs1+rs2)*sine]  
!                                         aux3 = (1/r)* ( (4*epsil*(A*sin(alpha*gamma+ph) + 1 ))* & 
!                                             (sigma/r)**6*(B*sin(beta*gamma)+1)**5 * &
!                                             ((sigma/r)**6 * 12*beta*B*cos(beta*gamma)* &
!                                             (B*sin(beta*gamma)+1)**6 - 6*beta*B*cos(beta*gamma)) + &
!                                             4*alpha*A*epsil*cos(alpha*gamma+ph)*(sigma/r)**6*(B*sin(beta*gamma)+1)**6 * &
!                                             ((sigma/r)**6 * (B*sin(beta*gamma)+1)**6 - 1)) 
        
!                                         Tor(-pa +ptr%p%n) =  -aux3*(r+rs1+rs2) + Tor(-pa +ptr%p%n)
!                                         ptr%p%F = aux2 + ptr%p%F 
!                                         ptrn%p%F = -aux2 + ptrn%p%F - aux3*[sine,-coss] 
!                                         ! partlst(ptr%p%n - pa) = 1
!                                     else if (rs2 > 0 .and. rs1 == 0) then 
!                                         gamma = theta(-pa +ptrn%p%n) + atan2(coss,sine)
!                                         aux2 = -24*epsil*(1+A*sin(alpha*gamma+ph))*(1/r**2) *  & 
!                                             ( (sigma*(1+B*sin(beta*gamma))/r)**6 * (1-2* (sigma*(1+B*sin(beta*gamma))/r)**6 )) * &
!                                             [(x1(1)-x2(1)) - (rs1+rs2)*coss, (x1(2)-x2(2)) - (rs1+rs2)*sine]  
!                                         aux3 = (1/r)* ( (4*epsil*(A*sin(alpha*gamma+ph) + 1 ))* & 
!                                             (sigma/r)**6*(B*sin(beta*gamma)+1)**5 * &
!                                             ((sigma/r)**6 * 12*beta*B*cos(beta*gamma)*(B*sin(beta*gamma)+1)**6 & 
!                                             - 6*beta*B*cos(beta*gamma)) + &
!                                             4*alpha*A*epsil*cos(alpha*gamma+ph)*(sigma/r)**6*(B*sin(beta*gamma)+1)**6 * &
!                                             ((sigma/r)**6 * (B*sin(beta*gamma)+1)**6 - 1)) 
        
!                                         Tor(-pa +ptrn%p%n) =  -aux3*(r+rs1+rs2) + Tor(-pa +ptrn%p%n)
!                                         ptr%p%F = aux2 + ptr%p%F + aux3*[sine,-coss] 
!                                         ptrn%p%F = -aux2 + ptrn%p%F 
!                                         partlst(ptrn%p%n - pa) = 1
!                                     else
!                                         ! Não vou estudar aqui a interação entre duas partículas, elas vão se repelir como átomos
!                                         aux2 = -(1/r**2)*(sigma*(1+B)/r)**6*(1-2*(sigma*(1+B)/r)**6)*24*epsil*(1+A)* & 
!                                         [(x1(1)-x2(1)) - (rs1+rs2)*coss, (x1(2)-x2(2)) - (rs1+rs2)*sine]  
!                                         ! print*, "L 395 r", r, "id",id
        
!                                         ! talvez o fricterm seja interessante em algum momento, vou deixar aqui
!                                         fric_term = (fric_term1+fric_term2)/2
!                                         if (fric_term > 0) then
!                                             fR = comp_fric([-(x1(1)-x2(1))+ (rs1+rs2)*coss,-(x1(2)-x2(2))+(rs1+rs2)*sine], &
!                                             [-(v1(1)-v2(1)),-(v1(2)-v2(2))],fric_term)
!                                         end if
        
!                                         ptr%p%F = aux2 + ptr%p%F +fR
!                                         ptrn%p%F = -aux2 + ptrn%p%F - fR
!                                     end if
!                                 end if
!                                 node => list_next(node) ! próxima partícula da célula
!                             end do                            
!                         end if
!                     else
!                          !interagirá com a próxima linha e coluna, e na diagonal
!                         node => list_next(malha(i,j+1)%list) 
                        
!                         do while (associated(node))
!                             ptrn = transfer(list_get(node), ptrn) !outra particula selecionada
!                             sigma_b = propriedade(ptrn%p%grupo)%sigma
!                             epsil_b = propriedade(ptrn%p%grupo)%epsilon 
!                             rs2 = propriedade(ptrn%p%grupo)%rs
!                             fric_term2 = propriedade(ptrn%p%grupo)%fric_term 
!                             ! Lorenz-Betherlot rule for mixing epsilon sigma 
!                             if (sigma_a > sigma_b) then
!                                 rcut = r_cut*sigma_a + (rs1 + rs2)
!                                 sigma = 0.5*(sigma_a + sigma_b)
!                                 epsil = sqrt(epsil_a*epsil_b)
!                             else if (sigma_a < sigma_b) then 
!                                 rcut = r_cut*sigma_b + (rs1 + rs2)
!                                 sigma = 0.5*(sigma_a + sigma_b)
!                                 epsil = sqrt(epsil_a*epsil_b)
!                             else 
!                                 rcut = r_cut*sigma_a + (rs1 + rs2)
!                                 sigma = sigma_a
!                                 epsil = epsil_a
!                             end if 

!                             x2 = ptrn%p%x
!                             v2 = ptrn%p%v
                            
!                             r = sqrt((x1(1)-x2(1))**2 + (x1(2)-x2(2))**2)
!                             coss = (x1(1)-x2(1))/r 
!                             sine = (x1(2)-x2(2))/r 
!                             r = r - rs1 - rs2 !raio
!                             ! print*, "L 757 r", r, "id",id
!                             if (r <= rcut) then
!                                 if (rs1 == 0 .and. rs2 == 0) then
!                                     aux2 = -(1/r**2)*(sigma/r)**6*(1-2*(sigma/r)**6)*24*epsil* & 
!                                     [(x1(1)-x2(1)) - (rs1+rs2)*coss, (x1(2)-x2(2)) - (rs1+rs2)*sine]  
!                                     ! print*, "L 395 r", r, "id",id
    
!                                     fric_term = (fric_term1+fric_term2)/2
!                                     if (fric_term > 0) then
!                                         fR = comp_fric([-(x1(1)-x2(1))+ (rs1+rs2)*coss,-(x1(2)-x2(2))+(rs1+rs2)*sine], &
!                                         [-(v1(1)-v2(1)),-(v1(2)-v2(2))],fric_term)
!                                     end if
!                                     ptr%p%F = aux2 + ptr%p%F +fR
!                                     ptrn%p%F = -aux2 + ptrn%p%F - fR
!                             else if (rs1 > 0 .and. rs2 == 0) then
    
!                                     gamma = theta(-pa +ptr%p%n) + atan2(coss,sine)
!                                     aux2 = -24*epsil*(1+A*sin(alpha*gamma+ph))*(1/r**2) * &
!                                         ( (sigma*(1+B*sin(beta*gamma))/r)**6 * (1-2* (sigma*(1+B*sin(beta*gamma))/r)**6 )) &
!                                         * [(x1(1)-x2(1)) - (rs1+rs2)*coss, (x1(2)-x2(2)) - (rs1+rs2)*sine]  
!                                     aux3 = (1/r)* ( (4*epsil*(A*sin(alpha*gamma+ph) + 1 ))* (sigma/r)**6*(B*sin(beta*gamma)+1)**5 * &
!                                         ((sigma/r)**6 * 12*beta*B*cos(beta*gamma)*(B*sin(beta*gamma)+1)**6 - 6*beta*B*cos(beta*gamma)) + &
!                                         4*alpha*A*epsil*cos(alpha*gamma+ph)*(sigma/r)**6*(B*sin(beta*gamma)+1)**6 * &
!                                         ((sigma/r)**6 * (B*sin(beta*gamma)+1)**6 - 1)) 
    
!                                     Tor(-pa +ptr%p%n) =  -aux3*(r+rs1+rs2) + Tor(-pa +ptr%p%n)
!                                     ptr%p%F = aux2 + ptr%p%F 
!                                     ptrn%p%F = -aux2 + ptrn%p%F - aux3*[sine,-coss] 
!                                     ! partlst(ptr%p%n - pa) = 1
!                                 else if (rs2 > 0 .and. rs1 == 0) then 
!                                     gamma = theta(-pa +ptrn%p%n) + atan2(coss,sine)
!                                     aux2 = -24*epsil*(1+A*sin(alpha*gamma+ph))*(1/r**2) * ( (sigma*(1+B*sin(beta*gamma))/r)**6 * & 
!                                             (1-2* (sigma*(1+B*sin(beta*gamma))/r)**6 )) * & 
!                                             [(x1(1)-x2(1)) - (rs1+rs2)*coss, (x1(2)-x2(2)) - (rs1+rs2)*sine]  
!                                     aux3 = (1/r)*( (4*epsil*(A*sin(alpha*gamma+ph) + 1 ))* (sigma/r)**6*(B*sin(beta*gamma)+1)**5 * &
!                                         ((sigma/r)**6 * 12*beta*B*cos(beta*gamma)*(B*sin(beta*gamma)+1)**6 - 6*beta*B*cos(beta*gamma)) + &
!                                         4*alpha*A*epsil*cos(alpha*gamma+ph)*(sigma/r)**6*(B*sin(beta*gamma)+1)**6 * &
!                                         ((sigma/r)**6 * (B*sin(beta*gamma)+1)**6 - 1)) 
    
!                                     Tor(-pa +ptrn%p%n) =  -aux3*(r+rs1+rs2) + Tor(-pa +ptrn%p%n)
!                                     ptr%p%F = aux2 + ptr%p%F + aux3*[sine,-coss] 
!                                     ptrn%p%F = -aux2 + ptrn%p%F 
!                                     partlst(ptrn%p%n - pa) = 1
!                                 else
!                                     ! Não vou estudar aqui a interação entre duas partículas, elas vão se repelir como átomos
!                                     aux2 = -(1/r**2)*(sigma*(1+B)/r)**6*(1-2*(sigma*(1+B)/r)**6)*24*epsil*(1+A)* & 
!                                     [(x1(1)-x2(1)) - (rs1+rs2)*coss, (x1(2)-x2(2)) - (rs1+rs2)*sine]  
!                                     ! print*, "L 395 r", r, "id",id
    
!                                     ! talvez o fricterm seja interessante em algum momento, vou deixar aqui
!                                     fric_term = (fric_term1+fric_term2)/2
!                                     if (fric_term > 0) then
!                                         fR = comp_fric([-(x1(1)-x2(1))+ (rs1+rs2)*coss,-(x1(2)-x2(2))+(rs1+rs2)*sine], &
!                                         [-(v1(1)-v2(1)),-(v1(2)-v2(2))],fric_term)
!                                     end if
    
!                                     ptr%p%F = aux2 + ptr%p%F +fR
!                                     ptrn%p%F = -aux2 + ptrn%p%F - fR
!                                 end if
!                             end if
!                             node => list_next(node) ! próxima partícula da célula                                
!                         end do

!                         node => list_next(malha(i+1,j)%list) !interagirá com a próxima linha 
!                         do while (associated(node))
!                             ptrn = transfer(list_get(node), ptrn) !outra particula selecionada
!                             sigma_b = propriedade(ptrn%p%grupo)%sigma
!                             epsil_b = propriedade(ptrn%p%grupo)%epsilon 
!                             rs2 = propriedade(ptrn%p%grupo)%rs 
!                             fric_term2 = propriedade(ptrn%p%grupo)%fric_term
!                             ! Lorenz-Betherlot rule for mixing epsilon sigma 
!                             if (sigma_a > sigma_b) then
!                                 rcut = r_cut*sigma_a + (rs1 + rs2)
!                                 sigma = 0.5*(sigma_a + sigma_b)
!                                 epsil = sqrt(epsil_a*epsil_b)
!                             else if (sigma_a < sigma_b) then 
!                                 rcut = r_cut*sigma_b + (rs1 + rs2)
!                                 sigma = 0.5*(sigma_a + sigma_b)
!                                 epsil = sqrt(epsil_a*epsil_b)
!                             else 
!                                 rcut = r_cut*sigma_a + (rs1 + rs2)
!                                 sigma = sigma_a
!                                 epsil = epsil_a
!                             end if 
!                             x2 = ptrn%p%x
!                             v2 = ptrn%p%v
!                             r = sqrt((x1(1)-x2(1))**2 + (x1(2)-x2(2))**2)
!                             coss = (x1(1)-x2(1))/r 
!                             sine = (x1(2)-x2(2))/r 
!                             r = r - rs1 - rs2 !raio
!                             ! print*, "L 798 r", r, "id",id
!                             if (r <= rcut) then
!                                 if (rs1 == 0 .and. rs2 == 0) then
!                                     aux2 = -(1/r**2)*(sigma/r)**6*(1-2*(sigma/r)**6)*24*epsil* & 
!                                     [(x1(1)-x2(1)) - (rs1+rs2)*coss, (x1(2)-x2(2)) - (rs1+rs2)*sine]  
!                                     ! print*, "L 395 r", r, "id",id
    
!                                     fric_term = (fric_term1+fric_term2)/2
!                                     if (fric_term > 0) then
!                                         fR = comp_fric([-(x1(1)-x2(1))+ (rs1+rs2)*coss,-(x1(2)-x2(2))+(rs1+rs2)*sine], &
!                                         [-(v1(1)-v2(1)),-(v1(2)-v2(2))],fric_term)
!                                     end if
!                                     ptr%p%F = aux2 + ptr%p%F +fR
!                                     ptrn%p%F = -aux2 + ptrn%p%F - fR
                                    
!                             else if (rs1 > 0 .and. rs2 == 0) then
    
!                                     gamma = theta(-pa +ptr%p%n) + atan2(coss,sine)
!                                     aux2 = -24*epsil*(1+A*sin(alpha*gamma+ph))*(1/r**2) * & 
!                                         ((sigma*(1+B*sin(beta*gamma))/r)**6 * (1-2* (sigma*(1+B*sin(beta*gamma))/r)**6 )) * &
!                                         [(x1(1)-x2(1)) - (rs1+rs2)*coss, (x1(2)-x2(2)) - (rs1+rs2)*sine]  
!                                     aux3 = (1/r)* ( (4*epsil*(A*sin(alpha*gamma+ph) + 1 ))* (sigma/r)**6*(B*sin(beta*gamma)+1)**5 * &
!                                         ((sigma/r)**6 * 12*beta*B*cos(beta*gamma)*(B*sin(beta*gamma)+1)**6 - 6*beta*B*cos(beta*gamma)) + &
!                                         4*alpha*A*epsil*cos(alpha*gamma+ph)*(sigma/r)**6*(B*sin(beta*gamma)+1)**6 * &
!                                         ((sigma/r)**6 * (B*sin(beta*gamma)+1)**6 - 1)) 
    
!                                     Tor(-pa +ptr%p%n) =  -aux3*(r+rs1+rs2) + Tor(-pa +ptr%p%n)
!                                     ptr%p%F = aux2 + ptr%p%F 
!                                     ptrn%p%F = -aux2 + ptrn%p%F - aux3*[sine,-coss] 
!                                     ! partlst(ptr%p%n - pa) = 1
!                                 else if (rs2 > 0 .and. rs1 == 0) then 
!                                     gamma = theta(-pa +ptrn%p%n) + atan2(coss,sine)
!                                     aux2 = -24*epsil*(1+A*sin(alpha*gamma+ph))*(1/r**2) *  ( (sigma*(1+B*sin(beta*gamma))/r)**6 * & 
!                                         (1-2* (sigma*(1+B*sin(beta*gamma))/r)**6 )) * &
!                                         [(x1(1)-x2(1)) - (rs1+rs2)*coss, (x1(2)-x2(2)) - (rs1+rs2)*sine]  
!                                     aux3 = (1/r)* ( (4*epsil*(A*sin(alpha*gamma+ph) + 1 ))* (sigma/r)**6*(B*sin(beta*gamma)+1)**5 * &
!                                         ((sigma/r)**6 * 12*beta*B*cos(beta*gamma)*(B*sin(beta*gamma)+1)**6 - 6*beta*B*cos(beta*gamma)) + &
!                                         4*alpha*A*epsil*cos(alpha*gamma+ph)*(sigma/r)**6*(B*sin(beta*gamma)+1)**6 * &
!                                         ((sigma/r)**6 * (B*sin(beta*gamma)+1)**6 - 1)) 
    
!                                     Tor(-pa +ptrn%p%n) =  -aux3*(r+rs1+rs2) + Tor(-pa +ptrn%p%n)
!                                     ptr%p%F = aux2 + ptr%p%F + aux3*[sine,-coss] 
!                                     ptrn%p%F = -aux2 + ptrn%p%F 
!                                     partlst(ptrn%p%n - pa) = 1
!                                 else
!                                     ! Não vou estudar aqui a interação entre duas partículas, elas vão se repelir como átomos
!                                     aux2 = -(1/r**2)*(sigma*(1+B)/r)**6*(1-2*(sigma*(1+B)/r)**6)*24*epsil*(1+A)* & 
!                                     [(x1(1)-x2(1)) - (rs1+rs2)*coss, (x1(2)-x2(2)) - (rs1+rs2)*sine]  
!                                     ! print*, "L 395 r", r, "id",id
    
!                                     ! talvez o fricterm seja interessante em algum momento, vou deixar aqui
!                                     fric_term = (fric_term1+fric_term2)/2
!                                     if (fric_term > 0) then
!                                         fR = comp_fric([-(x1(1)-x2(1))+ (rs1+rs2)*coss,-(x1(2)-x2(2))+(rs1+rs2)*sine], &
!                                         [-(v1(1)-v2(1)),-(v1(2)-v2(2))],fric_term)
!                                     end if
    
!                                     ptr%p%F = aux2 + ptr%p%F +fR
!                                     ptrn%p%F = -aux2 + ptrn%p%F - fR
!                                 end if
!                             end if
!                             node => list_next(node) ! próxima partícula da célula                                
!                         end do
                        
!                         node => list_next(malha(i+1,j+1)%list) !interagirá com a próxima linha e coluna
!                         do while (associated(node))
!                             ptrn = transfer(list_get(node), ptrn) !outra particula selecionada
!                             sigma_b = propriedade(ptrn%p%grupo)%sigma
!                             epsil_b = propriedade(ptrn%p%grupo)%epsilon 
!                             rs2 = propriedade(ptrn%p%grupo)%rs 
!                             fric_term2 = propriedade(ptrn%p%grupo)%fric_term
!                             ! Lorenz-Betherlot rule for mixing epsilon sigma 
!                             if (sigma_a > sigma_b) then
!                                 rcut = r_cut*sigma_a + (rs1 + rs2)
!                                 sigma = 0.5*(sigma_a + sigma_b)
!                                 epsil = sqrt(epsil_a*epsil_b)
!                             else if (sigma_a < sigma_b) then 
!                                 rcut = r_cut*sigma_b + (rs1 + rs2)
!                                 sigma = 0.5*(sigma_a + sigma_b)
!                                 epsil = sqrt(epsil_a*epsil_b)
!                             else 
!                                 rcut = r_cut*sigma_a + (rs1 + rs2)
!                                 sigma = sigma_a
!                                 epsil = epsil_a
!                             end if 

!                             x2 = ptrn%p%x
!                             v2 = ptrn%p%v
!                             r = sqrt((x1(1)-x2(1))**2 + (x1(2)-x2(2))**2)
!                             coss = (x1(1)-x2(1))/r 
!                             sine = (x1(2)-x2(2))/r 
!                             r = r - rs1 - rs2 !raio
!                             ! print*, "L 841 r", r, "id",id
!                             if (r <= rcut) then
!                                 if (rs1 == 0 .and. rs2 == 0) then
!                                     aux2 = -(1/r**2)*(sigma/r)**6*(1-2*(sigma/r)**6)*24*epsil* & 
!                                     [(x1(1)-x2(1)) - (rs1+rs2)*coss, (x1(2)-x2(2)) - (rs1+rs2)*sine]  
!                                     ! print*, "L 395 r", r, "id",id
    
!                                     fric_term = (fric_term1+fric_term2)/2
!                                     if (fric_term > 0) then
!                                         fR = comp_fric([-(x1(1)-x2(1))+ (rs1+rs2)*coss,-(x1(2)-x2(2))+(rs1+rs2)*sine], &
!                                         [-(v1(1)-v2(1)),-(v1(2)-v2(2))],fric_term)
!                                     end if
!                                     ptr%p%F = aux2 + ptr%p%F +fR
!                                     ptrn%p%F = -aux2 + ptrn%p%F - fR
!                                 else if (rs1 > 0 .and. rs2 == 0) then
    
!                                     gamma = theta(-pa +ptr%p%n) + atan2(coss,sine)
!                                     aux2 = -24*epsil*(1+A*sin(alpha*gamma+ph))*(1/r**2) *  ( (sigma*(1+B*sin(beta*gamma))/r)**6 * &
!                                         (1-2* (sigma*(1+B*sin(beta*gamma))/r)**6 )) * &
!                                         [(x1(1)-x2(1)) - (rs1+rs2)*coss, (x1(2)-x2(2)) - (rs1+rs2)*sine]  
!                                     aux3 = (1/r)* ( (4*epsil*(A*sin(alpha*gamma+ph) + 1 ))* (sigma/r)**6*(B*sin(beta*gamma)+1)**5 * &
!                                         ((sigma/r)**6 * 12*beta*B*cos(beta*gamma)*(B*sin(beta*gamma)+1)**6 - 6*beta*B*cos(beta*gamma)) + &
!                                         4*alpha*A*epsil*cos(alpha*gamma+ph)*(sigma/r)**6*(B*sin(beta*gamma)+1)**6 * &
!                                         ((sigma/r)**6 * (B*sin(beta*gamma)+1)**6 - 1)) 
    
!                                     Tor(-pa +ptr%p%n) =  -aux3*(r+rs1+rs2) + Tor(-pa +ptr%p%n)
!                                     ptr%p%F = aux2 + ptr%p%F 
!                                     ptrn%p%F = -aux2 + ptrn%p%F - aux3*[sine,-coss] 
!                                     ! partlst(ptr%p%n - pa) = 1
!                                 else if (rs2 > 0 .and. rs1 == 0) then 
!                                     gamma = theta(-pa +ptrn%p%n) + atan2(coss,sine)
!                                     aux2 = -24*epsil*(1+A*sin(alpha*gamma+ph))*(1/r**2) * &
!                                         ( (sigma*(1+B*sin(beta*gamma))/r)**6 * (1-2* (sigma*(1+B*sin(beta*gamma))/r)**6 )) * &
!                                         [(x1(1)-x2(1)) - (rs1+rs2)*coss, (x1(2)-x2(2)) - (rs1+rs2)*sine]  
!                                     aux3 = (1/r)* ( (4*epsil*(A*sin(alpha*gamma+ph) + 1 ))* (sigma/r)**6*(B*sin(beta*gamma)+1)**5 * &
!                                         ((sigma/r)**6 * 12*beta*B*cos(beta*gamma)*(B*sin(beta*gamma)+1)**6 - 6*beta*B*cos(beta*gamma)) + &
!                                         4*alpha*A*epsil*cos(alpha*gamma+ph)*(sigma/r)**6*(B*sin(beta*gamma)+1)**6 * &
!                                         ((sigma/r)**6 * (B*sin(beta*gamma)+1)**6 - 1)) 
    
!                                     Tor(-pa +ptrn%p%n) =  -aux3*(r+rs1+rs2) + Tor(-pa +ptrn%p%n)
!                                     ptr%p%F = aux2 + ptr%p%F + aux3*[sine,-coss] 
!                                     ptrn%p%F = -aux2 + ptrn%p%F 
!                                     partlst(ptrn%p%n - pa) = 1
!                                 else
!                                     ! Não vou estudar aqui a interação entre duas partículas, elas vão se repelir como átomos
!                                     aux2 = -(1/r**2)*(sigma*(1+B)/r)**6*(1-2*(sigma*(1+B)/r)**6)*24*epsil*(1+A)* & 
!                                     [(x1(1)-x2(1)) - (rs1+rs2)*coss, (x1(2)-x2(2)) - (rs1+rs2)*sine]  
!                                     ! print*, "L 395 r", r, "id",id
    
!                                     ! talvez o fricterm seja interessante em algum momento, vou deixar aqui
!                                     fric_term = (fric_term1+fric_term2)/2
!                                     if (fric_term > 0) then
!                                         fR = comp_fric([-(x1(1)-x2(1))+ (rs1+rs2)*coss,-(x1(2)-x2(2))+(rs1+rs2)*sine], &
!                                         [-(v1(1)-v2(1)),-(v1(2)-v2(2))],fric_term)
!                                     end if
    
!                                     ptr%p%F = aux2 + ptr%p%F +fR
!                                     ptrn%p%F = -aux2 + ptrn%p%F - fR
!                                 end if
!                             end if
!                             node => list_next(node) ! próxima partícula da célula                                
!                         end do
                        
!                         if (j /= 1) then !se não for a primeira coluna 
!                             node => list_next(malha(i+1,j-1)%list) !interagirá com a próxima linha apenas
!                             do while (associated(node))
!                                 ptrn = transfer(list_get(node), ptrn) !outra particula selecionada
!                                 sigma_b = propriedade(ptrn%p%grupo)%sigma
!                                 epsil_b = propriedade(ptrn%p%grupo)%epsilon 
!                                 rs2 = propriedade(ptrn%p%grupo)%rs 
!                                 fric_term2 = propriedade(ptrn%p%grupo)%fric_term
!                                 ! Lorenz-Betherlot rule for mixing epsilon sigma 
!                                 if (sigma_a > sigma_b) then
!                                     rcut = r_cut*sigma_a + (rs1 + rs2)
!                                     sigma = 0.5*(sigma_a + sigma_b)
!                                     epsil = sqrt(epsil_a*epsil_b)
!                                 else if (sigma_a < sigma_b) then 
!                                     rcut = r_cut*sigma_b + (rs1 + rs2)
!                                     sigma = 0.5*(sigma_a + sigma_b)
!                                     epsil = sqrt(epsil_a*epsil_b)
!                                 else 
!                                     rcut = r_cut*sigma_a + (rs1 + rs2)
!                                     sigma = sigma_a
!                                     epsil = epsil_a
!                                 end if 

!                                 x2 = ptrn%p%x
!                                 v2 = ptrn%p%v
!                                 r = sqrt((x1(1)-x2(1))**2 + (x1(2)-x2(2))**2)
!                                 coss = (x1(1)-x2(1))/r 
!                                 sine = (x1(2)-x2(2))/r 
!                                 r = r - rs1 - rs2 !raio
!                                 ! print*, "L 885 r", r, "id",id
!                                 if (r <= rcut) then
!                                     if (rs1 == 0 .and. rs2 == 0) then
!                                         aux2 = -(1/r**2)*(sigma/r)**6*(1-2*(sigma/r)**6)*24*epsil* & 
!                                         [(x1(1)-x2(1)) - (rs1+rs2)*coss, (x1(2)-x2(2)) - (rs1+rs2)*sine]  
!                                         ! print*, "L 395 r", r, "id",id
        
!                                         fric_term = (fric_term1+fric_term2)/2
!                                         if (fric_term > 0) then
!                                             fR = comp_fric([-(x1(1)-x2(1))+ (rs1+rs2)*coss,-(x1(2)-x2(2))+(rs1+rs2)*sine], &
!                                             [-(v1(1)-v2(1)),-(v1(2)-v2(2))],fric_term)
!                                         end if
!                                         ptr%p%F = aux2 + ptr%p%F +fR
!                                         ptrn%p%F = -aux2 + ptrn%p%F - fR
!                                 else if (rs1 > 0 .and. rs2 == 0) then
        
!                                         gamma = theta(-pa +ptr%p%n) + atan2(coss,sine)
!                                         aux2 = -24*epsil*(1+A*sin(alpha*gamma+ph))*(1/r**2) *  ( (sigma*(1+B*sin(beta*gamma))/r)**6 &
!                                             * (1-2* (sigma*(1+B*sin(beta*gamma))/r)**6 )) * &
!                                             [(x1(1)-x2(1)) - (rs1+rs2)*coss, (x1(2)-x2(2)) - (rs1+rs2)*sine]  
!                                         aux3 = (1/r)* ( (4*epsil*(A*sin(alpha*gamma+ph) + 1 ))* (sigma/r)**6*(B*sin(beta*gamma)+1)**5 * &
!                                             ((sigma/r)**6 * 12*beta*B*cos(beta*gamma)*(B*sin(beta*gamma)+1)**6 - 6*beta*B*cos(beta*gamma)) + &
!                                             4*alpha*A*epsil*cos(alpha*gamma+ph)*(sigma/r)**6*(B*sin(beta*gamma)+1)**6 * &
!                                             ((sigma/r)**6 * (B*sin(beta*gamma)+1)**6 - 1)) 
        
!                                         Tor(-pa +ptr%p%n) =  -aux3*(r+rs1+rs2) + Tor(-pa +ptr%p%n)
!                                         ptr%p%F = aux2 + ptr%p%F 
!                                         ptrn%p%F = -aux2 + ptrn%p%F - aux3*[sine,-coss] 
!                                         ! partlst(ptr%p%n - pa) = 1
!                                     else if (rs2 > 0 .and. rs1 == 0) then 
!                                         gamma = theta(-pa +ptrn%p%n) + atan2(coss,sine)
!                                         aux2 = -24*epsil*(1+A*sin(alpha*gamma+ph))*(1/r**2) * ((sigma*(1+B*sin(beta*gamma))/r)**6 * &
!                                         (1-2* (sigma*(1+B*sin(beta*gamma))/r)**6 )) * &
!                                         [(x1(1)-x2(1)) - (rs1+rs2)*coss, (x1(2)-x2(2)) - (rs1+rs2)*sine]  
!                                         aux3 = (1/r)* ( (4*epsil*(A*sin(alpha*gamma+ph) + 1 ))* (sigma/r)**6*(B*sin(beta*gamma)+1)**5 * &
!                                             ((sigma/r)**6 * 12*beta*B*cos(beta*gamma)*(B*sin(beta*gamma)+1)**6 - 6*beta*B*cos(beta*gamma)) + &
!                                             4*alpha*A*epsil*cos(alpha*gamma+ph)*(sigma/r)**6*(B*sin(beta*gamma)+1)**6 * &
!                                             ((sigma/r)**6 * (B*sin(beta*gamma)+1)**6 - 1)) 
        
!                                         Tor(-pa +ptrn%p%n) =  -aux3*(r+rs1+rs2) + Tor(-pa +ptrn%p%n)
!                                         ptr%p%F = aux2 + ptr%p%F + aux3*[sine,-coss] 
!                                         ptrn%p%F = -aux2 + ptrn%p%F 
!                                         partlst(ptrn%p%n - pa) = 1
!                                     else
!                                         ! Não vou estudar aqui a interação entre duas partículas, elas vão se repelir como átomos
!                                         aux2 = -(1/r**2)*(sigma*(1+B)/r)**6*(1-2*(sigma*(1+B)/r)**6)*24*epsil*(1+A)* & 
!                                         [(x1(1)-x2(1)) - (rs1+rs2)*coss, (x1(2)-x2(2)) - (rs1+rs2)*sine]  
!                                         ! print*, "L 395 r", r, "id",id
        
!                                         ! talvez o fricterm seja interessante em algum momento, vou deixar aqui
!                                         fric_term = (fric_term1+fric_term2)/2
!                                         if (fric_term > 0) then
!                                             fR = comp_fric([-(x1(1)-x2(1))+ (rs1+rs2)*coss,-(x1(2)-x2(2))+(rs1+rs2)*sine], &
!                                             [-(v1(1)-v2(1)),-(v1(2)-v2(2))],fric_term)
!                                         end if
        
!                                         ptr%p%F = aux2 + ptr%p%F +fR
!                                         ptrn%p%F = -aux2 + ptrn%p%F - fR
!                                     end if
!                                 end if
!                                 node => list_next(node) ! próxima partícula da célula
!                             end do                            
!                         end if
!                     end if
                    
!                 else ! se for a última lina, só interage com a celua ao lado 
!                     node => list_next(malha(i,j+1)%list) !interagirá com a próxima linha e coluna
!                     do while (associated(node))
!                         ptrn = transfer(list_get(node), ptrn) !outra particula selecionada
!                         sigma_b = propriedade(ptrn%p%grupo)%sigma
!                         epsil_b = propriedade(ptrn%p%grupo)%epsilon 
!                         rs2 = propriedade(ptrn%p%grupo)%rs 
!                         fric_term2 = propriedade(ptrn%p%grupo)%fric_term
!                         ! Lorenz-Betherlot rule for mixing epsilon sigma 
!                         if (sigma_a > sigma_b) then
!                             rcut = r_cut*sigma_a + (rs1 + rs2)
!                             sigma = 0.5*(sigma_a + sigma_b)
!                             epsil = sqrt(epsil_a*epsil_b)
!                         else if (sigma_a < sigma_b) then 
!                             rcut = r_cut*sigma_b + (rs1 + rs2)
!                             sigma = 0.5*(sigma_a + sigma_b)
!                             epsil = sqrt(epsil_a*epsil_b)
!                         else 
!                             rcut = r_cut*sigma_a + (rs1 + rs2)
!                             sigma = sigma_a
!                             epsil = epsil__a
!                         end if 

!                         x2 = ptrn%p%x
!                         v2 = ptrn%p%v
!                         r = sqrt((x1(1)-x2(1))**2 + (x1(2)-x2(2))**2)
!                         coss = (x1(1)-x2(1))/r 
!                         sine = (x1(2)-x2(2))/r 
!                         r = r - rs1 - rs2 !raio
!                         ! print*, "L 931 r", r, "id",id
!                         if (r <= rcut) then
!                             if (rs1 == 0 .and. rs2 == 0) then
!                                  aux2 = -(1/r**2)*(sigma/r)**6*(1-2*(sigma/r)**6)*24*epsil* & 
!                                  [(x1(1)-x2(1)) - (rs1+rs2)*coss, (x1(2)-x2(2)) - (rs1+rs2)*sine]  
!                                  ! print*, "L 395 r", r, "id",id
 
!                                  fric_term = (fric_term1+fric_term2)/2
!                                  if (fric_term > 0) then
!                                      fR = comp_fric([-(x1(1)-x2(1))+ (rs1+rs2)*coss,-(x1(2)-x2(2))+(rs1+rs2)*sine], &
!                                      [-(v1(1)-v2(1)),-(v1(2)-v2(2))],fric_term)
!                                  end if
!                                  ptr%p%F = aux2 + ptr%p%F +fR
!                                  ptrn%p%F = -aux2 + ptrn%p%F - fR
!                             else if (rs1 > 0 .and. rs2 == 0) then
!                                  gamma = theta(-pa +ptr%p%n) + atan2(coss,sine)
 
!                                  aux2 = -24*epsil*(1+A*sin(alpha*gamma+ph))*(1/r**2) * & 
!                                      ( (sigma*(1+B*sin(beta*gamma))/r)**6 * (1-2* (sigma*(1+B*sin(beta*gamma))/r)**6 )) * &
!                                      [(x1(1)-x2(1)) - (rs1+rs2)*coss, (x1(2)-x2(2)) - (rs1+rs2)*sine]  
!                                  aux3 = (1/r)* ( (4*epsil*(A*sin(alpha*gamma+ph) + 1 ))* (sigma/r)**6*(B*sin(beta*gamma)+1)**5 * &
!                                      ((sigma/r)**6 * 12*beta*B*cos(beta*gamma)*(B*sin(beta*gamma)+1)**6 - 6*beta*B*cos(beta*gamma)) + &
!                                      4*alpha*A*epsil*cos(alpha*gamma+ph)*(sigma/r)**6*(B*sin(beta*gamma)+1)**6 * &
!                                      ((sigma/r)**6 * (B*sin(beta*gamma)+1)**6 - 1)) 
 
!                                  ! print*, "r1 > 0 gamma", gamma
!                                  ! print*, "aux2, aux3",aux2, aux3
!                                  ! print*, "x1 =", x1
!                                  ! print*, "x2 =", x2
!                                  ! read(*,*)
!                                 Tor(-pa +ptr%p%n) =  -aux3*(r+rs1+rs2) + Tor(-pa +ptr%p%n)
!                                 ptr%p%F = aux2 + ptr%p%F 
!                                 ptrn%p%F = -aux2 + ptrn%p%F - aux3*[sine,-coss] 
!                                  ! partlst(ptr%p%n - pa) = 1
 
!                             else if (rs2 > 0 .and. rs1 == 0) then 
!                                 gamma = theta(-pa +ptrn%p%n) + atan2(coss,sine)
!                                 aux2 = -24*epsil*(1+A*sin(alpha*gamma+ph))*(1/r**2) * ((sigma*(1+B*sin(beta*gamma))/r)**6 *  &
!                                     (1-2* (sigma*(1+B*sin(beta*gamma))/r)**6 )) * &
!                                     [(x1(1)-x2(1)) - (rs1+rs2)*coss, (x1(2)-x2(2)) - (rs1+rs2)*sine]  
!                                 aux3 = (1/r)* ( (4*epsil*(A*sin(alpha*gamma+ph) + 1 ))* (sigma/r)**6*(B*sin(beta*gamma)+1)**5 * &
!                                     ((sigma/r)**6 * 12*beta*B*cos(beta*gamma)*(B*sin(beta*gamma)+1)**6 - 6*beta*B*cos(beta*gamma)) + &
!                                     4*alpha*A*epsil*cos(alpha*gamma+ph)*(sigma/r)**6*(B*sin(beta*gamma)+1)**6 * &
!                                     ((sigma/r)**6 * (B*sin(beta*gamma)+1)**6 - 1)) 

!                                 Tor(-pa +ptrn%p%n) =  -aux3*(r+rs1+rs2) + Tor(-pa +ptrn%p%n)
!                                 ptr%p%F = aux2 + ptr%p%F + aux3*[sine,-coss] 
!                                 ptrn%p%F = -aux2 + ptrn%p%F 
!                                 partlst(ptrn%p%n - pa) = 1
!                             else
!                                 ! Não vou estudar aqui a interação entre duas partículas, elas vão se repelir como átomos
!                                 aux2 = -(1/r**2)*(sigma*(1+B)/r)**6*(1-2*(sigma*(1+B)/r)**6)*24*epsil*(1+A)* & 
!                                 [(x1(1)-x2(1)) - (rs1+rs2)*coss, (x1(2)-x2(2)) - (rs1+rs2)*sine]  
!                                 ! print*, "L 395 r", r, "id",id

!                                 ! talvez o fricterm seja interessante em algum momento, vou deixar aqui
!                                 fric_term = (fric_term1+fric_term2)/2
!                                 if (fric_term > 0) then
!                                     fR = comp_fric([-(x1(1)-x2(1))+ (rs1+rs2)*coss,-(x1(2)-x2(2))+(rs1+rs2)*sine], &
!                                     [-(v1(1)-v2(1)),-(v1(2)-v2(2))],fric_term)
!                                 end if

!                                 ptr%p%F = aux2 + ptr%p%F +fR
!                                 ptrn%p%F = -aux2 + ptrn%p%F - fR
!                             end if
!                         end if
!                         node => list_next(node) ! próxima partícula da célula                                
!                     end do
!                 end if
!                 node => next_node
!             end do
!         end do
!     end do
!      !print*, 'Linha 143'
!     ! if (id == 0) read(*,*)
!     ! call MPI_barrier(MPI_COMM_WORLD, ierr) 
! end subroutine comp_FT    