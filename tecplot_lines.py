import numpy as np
import tecplot as tp
import pandas as pd
from tecplot.constant import PlotType, Color,Units, LinePattern, AxisTitleMode
import sys
# if '-c' in sys.argv:
#     tp.session.connect()

tp.session.connect()

# IMPORTA LINHAS DO CSV 
csv_input = pd.read_csv('Distributions.csv')
hl = list(csv_input) #headers list

## NOME DAS ZONAS ##
zonas = \
[['A','neon_3sep_stat'], \
['B','carlos_2sep_stat'], \
['C','neon4sep_stat'], \
['D','kubuntu_3sep_stat'], \
['E','kubuntu_5sep_stat'], \
['F','suse_6sep_stat'], \
['G','rafa_3sep_stat'], \
['H','carlos_4sep_stat'], \
['I','rafa_4sep_stat'], \
['J','kubuntu_4sep_stat'], \
['K','rafa_6sep_stat'], \
['L','suse_3sep_stat'], \
['M','rafa_5sep_stat'], \
['N','carlos_3sep_stat'], \
['O','kubuntu_7_sep_stat'], \
['P','suse_7sep_STAT'], \
['Q','rafa_7sep_stat']]


for i in range(len(hl)):
    print("{} | {} | {}\n".format(i,[zonas[j][0] for j in range(len(zonas)) if zonas[j][1] == hl[i]][0],hl[i]))
linhas = input("Entre os numeros dos dados que quer plotar separados por espaço:\n").split()
linhas = [int(x) for x in linhas]

frame = tp.active_frame()
frame.activate()
frame.height = 8
frame.width = 11
dataset = frame.create_dataset("Data", ['x', 'y'])
for dado in linhas:
    col = [zonas[j][0] for j in range(len(zonas)) if zonas[j][1] == hl[dado]][0]
    x = csv_input[hl[dado]].to_numpy()
    zone = dataset.add_ordered_zone(col, len(x))
    zone.values('x')[:] = np.linspace(0,1,len(x))   
    zone.values('y')[:] = x

# Set plot type to XYLine
plot = frame.plot(PlotType.XYLine)
plot.activate()

############# EIXOS #############################
x_axis = plot.axes.y_axis(0)
x_axis.title.title_mode = AxisTitleMode.UseText
x_axis.title.text = 'x'
x_axis.fit_range_to_nice()

y_axis = plot.axes.y_axis(0)
y_axis.title.title_mode = AxisTitleMode.UseText
y_axis.title.text = 'Fraction of Particles'
y_axis.fit_range_to_nice()

for ax in [x_axis, y_axis]:
    ax.title.font.typeface = 'Times'
    ax.title.font.bold = False
    ax.title.font.italic = False
    ax.title.font.size_units = Units.Frame
    ax.title.font.size = 7

################ LINHAS #####################

# Line patters list 
# DashDot = 2 DashDotDot = 5
# Dashed = 1 Dotted = 3
# LongDash = 4 Solid = 0

# Show all linemaps and make the lines a bit thicker
cont = 1
for lmap in plot.linemaps():
    lmap.show = True
    lmap.line.line_thickness = 0.6
    lmap.line.line_pattern = LinePattern(cont)
    lmap.line.color = Color.Black
    cont += 1

################# LEGENDA #####################
legend = plot.legend
legend.font.typeface = 'Times'
legend.show = True


################ IMAGEM SALVA ######################
tp.export.save_png('add_ordered_zones.png', 600, supersample=3)

# # clear plot
# plot.delete_linemaps()