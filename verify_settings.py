
import configparser

change = False
aux = True
config = configparser.ConfigParser()
config.optionxform=str

while (aux):
    try:
        config.read('settings.ini')
        aux = False
    except configparser.DuplicateSectionError as a:
        print(a)
        a = input("\nChange the settings file then enter to continue.\n")
        pass


N = int(config['global']['N'].split()[0])
ntype = int(config['global']['Ntype'].split()[0])
quant = []

i = ntype 
aux = True
while aux: 
    try:
        quant.append(int(config['par_'+str(i)]['quantidade'].split()[0])) 
        print("Number of particle types greater then expected (ntype = {} but there is par_0 to par_{}).\n".format(ntype,i))
        nty = input("Enter the right number of particle types.\n")
        i = int(nty) 
        config.set('global','Ntype',nty)
    except KeyError:
        aux = False 
        pass 


j = 0

aux = True
while aux:
    try:
        for i in range(ntype):
            quant.append(int(config['par_'+str(i)]['quantidade'].split()[0]))
            j+= 1
        aux = False 
    except:
        print("Number of particle types (ntype) smaller then expected. It will correct to proceed with ntype = {}".format(j))
        config.set('global','Ntype',str(j))
        ntype = j 
        j = 0
        change = True
        

if sum(quant) != N: 
    print("The particle groups have a different number of particles then specified {}. \nIt will be corrected to {}".format(N,sum(quant)))
    config.set('global','N',str(sum(quant)))
    change = True

if change:
    print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n")
    print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n")
    a = input("Changes will be made to the settings file. Enter to proceed or Ctrl+C to abort.")
    a = input()
    with open('settings.ini', 'w') as configfile:
        config.write(configfile)

try:
    from evtk.hl import pointsToVTK 
except ModuleNotFoundError:
    print("Evtk module not found! You will not be able to convert the .CSV files to .VTU. Try:\n")
    print("conda install -c e3sm evtk     OR\npip install pyevtk")

with open('settings.ini','r') as file:
    data = file.readlines()

for line in range(len(data)):
    if data[line][0] != "[" and data[line][0] != "#" and data[line][0] != " ":
        data[line] = '\t' + data[line]

with open('settings.ini','w') as file:
    file.writelines(data)



    

