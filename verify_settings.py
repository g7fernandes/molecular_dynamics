
import configparser
import os

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
        os.system('xdg-open settings.ini')
        input("Fix the config file, then enter to continue...")
        config.read('settings.ini')
        # nty = input("Enter the right number of particle types.\n")
        # i = int(nty) 
        # config.set('global','Ntype',nty)
    except KeyError:
        aux = False 
        pass 


j = 0

aux = True
while aux:
    try:
        for i in range(ntype):
            quant.append(int(config['par_'+str(i)]['quantidade'].split()[0]))
        aux = False 
    except:
        print("Number of particle types (ntype) smaller then expected. It maybe correct to proceed with ntype = {}".format(j))
        os.system('xdg-open settings.ini')
        input("Fix the config file, then enter to continue...")
        config.read('settings.ini')

aux = True
while aux:
    for i in range(ntype):
        x_file = config['par_'+str(i)]['x'].split()[0]
        x_file = x_file.replace("'","")
        x_file = x_file.replace('"','')
        if os.path.isfile(x_file):
            aux = False
        else:
            print("Position file {} not found.".format(x_file))
            os.system('xdg-open settings.ini')
            input("Fix the config file, then enter to continue...")
            config.read('settings.ini')
        v_file = config['par_'+str(i)]['v_file'].split()[0]
        v_file = v_file.replace("'","")
        v_file = v_file.replace('"','')
        if v_file[0] == '%' or os.path.isfile(x_file):
            aux = False
        else:
            print("Velocity file {} not found.".format(v_file))
            os.system('xdg-open settings.ini')
            input("Fix the config file, then enter to continue...")
            config.read('settings.ini')
        

if sum(quant) != N: 
    print("The particle groups have a different number of particles then specified {}. \nIt should maybe be corrected to {}".format(N,sum(quant)))
    os.system('xdg-open settings.ini')
    input("Fix the config file, then enter to continue...")
    config.read('settings.ini')


try:
    from evtk.hl import pointsToVTK 
except ModuleNotFoundError:
    print("Evtk module not found! You will not be able to convert the .CSV files to .VTU. Try:\n")
    print("conda install -c e3sm evtk     OR\npip install pyevtk")





    

