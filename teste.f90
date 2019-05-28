program teste
    use mod1
   
    call system('mkdir temp')
    open(10,file='temp/test.txt',status="replace")
    write(10,*) 'lalalala'
    print*,'a'
    
    
    
    
end program teste
