# Start of the makefile
# Defining variables
objects = m_config.o linkedlist.o lennard.o data.o matprint.o randnormal.o mod1.o
f90comp = mpiifort
switch = -O3
# Makefile
execname: $(objects)
	$(f90comp) -o lennard.out $(switch) $(objects)
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
m_config.o: m_config.f90
	$(f90comp) -c $(switch) m_config.f90
lennard.o: mod1.mod linkedlist.mod data.mod matprint.mod m_config.mod randnormal.mod lennard.f90
	$(f90comp) -c $(switch) lennard.f90
# Cleaning everything
clean:
	rm $(objects)
	rm -f *.o *.mod *.MOD
# End of the makefile
