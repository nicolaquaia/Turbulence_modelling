# modify planeChannel tutorials to apply to cylinder case

import pandas as pd
import matplotlib.pyplot as plt



def open_boundary(file_path):
    
    # Open and process the file
    with open(file_path, 'r') as file:
        lines = file.readlines()[1:-1]  # Skip the first line

    # Process each line to extract values
    data = []
    for line in lines:
        line = line.strip()[1:-1]  # Remove '(' and ')'
        values = list(map(float, line.split()))  # Convert to floats
        data.append(values)
    
    n = len(data[0])
    print(n)

    if n == 1:
        name_list = ['l']
    if n == 3:
        name_list = ['x', 'y', 'z']

    if n == 6:
        name_list = ['xx', 'yy', 'zz', 'xy', 'xz', 'yz']

    # Create a DataFrame from the processed data
    df = pd.DataFrame(data, columns=name_list)

    return df


points = open_boundary('LES_smaTurbBoundary/constant/boundaryData/left/points')
U      = open_boundary('LES_smaTurbBoundary/constant/boundaryData/left/0/U')
L      = open_boundary('LES_smaTurbBoundary/constant/boundaryData/left/0/L')
R      = open_boundary('LES_smaTurbBoundary/constant/boundaryData/left/0/R')

plt.figure()
plt.plot(points['y'], U['x'], label='U')
plt.plot(points['y'], L['l'], label='L')
plt.plot(points['y'], R['yy'], label='Ryy')
plt.plot(points['y'], R['xx'], label='Rxx')
plt.legend()
plt.show()

# actual mesh on y direction
# region 7: top left 
nQuarter = 14
grading = 1
# region 8: central left
nyOuter = 28
grading = 2


