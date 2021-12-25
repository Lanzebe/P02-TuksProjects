import pandas as pd 
import numpy as np

a= pd.read_excel(r'E:\OneDrive\ENGINEERING\01-ENGINEERING TOOLS\03-VIBRATION\MVR group project\RoadProfile.xlsx', sheet_name='Sheet2')
x = np.array(a).flatten()
print('x:',x)

b= pd.read_excel(r'E:\OneDrive\ENGINEERING\01-ENGINEERING TOOLS\03-VIBRATION\MVR group project\RoadProfile.xlsx', sheet_name='Sheet3')
y = np.array(b).flatten()
print('y:',y)

xy_array=np.array([x,
                  y])
print('xy_array:',xy_array)


np.save('X-Y random road', xy_array)