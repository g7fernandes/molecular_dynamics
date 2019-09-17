import numpy as np
import tecplot as tp
import pandas as pd
from tecplot.constant import PlotType, Color
import sys
if '-c' in sys.argv:
    tp.session.connect()

# IMPORTA LINHAS DO CSV 
csv_input = pd.read_csv('Distributions.csv')
hl = list(csv_input) #headers list

frame = tp.active_frame()
frame.activate()
frame.height = 8
frame.width = 11
dataset = frame.create_dataset("Data", ['x', 'y'])
for col in hl:
    x = csv_input[col].to_numpy()
    zone = dataset.add_ordered_zone(col, len(x))
    zone.values('x')[:] = np.linspace(0,1,len(x))   
    zone.values('y')[:] = x


# Set plot type to XYLine
plot = frame.plot(PlotType.XYLine)
plot.activate()

# Show all linemaps and make the lines a bit thicker
for lmap in plot.linemaps():
    lmap.show = True
    lmap.line.line_thickness = 0.6

plot.legend.show = True

tp.export.save_png('add_ordered_zones.png', 600, supersample=3)