# Start of the makefile
# Defining variables
objects = m_config.o linkedlist.o mod1.o saida.o lennard.o data.o mod0.o matprint.o randnormal.o 
f90comp = mpiifort
switch = -O3 -assume protect_parens,minus0 -prec-div -prec-sqrt 
# Makefile
execname: $(objects)
	$(f90comp) -o lennard $(switch) $(objects)
randnormal.mod: randnormal.o randnormal.f90
	$(f90comp) -c $(switch) randnormal.f90
mod1.mod: mod1.o mod1.f90
	$(f90comp) -c $(switch) mod1.f90
matprint.mod: matprint.o matprint.f90
	$(f90comp) -c $(switch) matprint.f90
data.mod: data.o data.f90 
	$(f90comp) -c $(switch) data.f90
linkedlist.mod: linkedlist.o linkedlist.f90
	$(f90comp) -c $(switch) linkedlist.f90
mod0.mod: mod0.o mod0.f90
	$(f90comp) -c $(switch) mod0.f90
saida.mod: saida.o saida.f90
	$(f90comp) -c $(switch) saida.f90
m_config.mod: m_config.o m_config.f90
	$(f90comp) -c $(switch) m_config.f90
randnormal.o: randnormal.f90
	$(f90comp) -c $(switch) randnormal.f90
mod1.o: mod1.f90
	$(f90comp) -c $(switch) mod1.f90
matprint.o: matprint.f90 
	$(f90comp) -c $(switch) matprint.f90
data.o: data.f90 
	$(f90comp) -c $(switch) data.f90
linkedlist.o: linkedlist.f90
	$(f90comp) -c $(switch) linkedlist.f90
mod0.o: mod0.f90
	$(f90comp) -c $(switch) mod0.f90
saida.o: saida.f90
	$(f90comp) -c $(switch) saida.f90
m_config.o: saida.f90
	$(f90comp) -c $(switch) m_config.f90
lennard.o: mod1.mod linkedlist.mod mod0.mod saida.mod data.mod matprint.mod m_config.mod randnormal.mod lennard.f90
	$(f90comp) -c $(switch) lennard.f90
# Cleaning everything
clean:
	rm $(objects)
#	rm -f *.o *.mod *.MOD
# End of the makefile
