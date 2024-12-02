# modify planeChannel tutorials to apply to cylinder case

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np


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

    if n == 1:
        name_list = ['l']
    if n == 3:
        name_list = ['x', 'y', 'z']

    if n == 6:
        name_list = ['xx', 'yy', 'zz', 'xy', 'xz', 'yz']

    # Create a DataFrame from the processed data
    df = pd.DataFrame(data, columns=name_list)

    return df

def computeCentroidGrid(nQuarter = 14, nyOuter = 28, rOuter45 = np.sqrt(0.5),  ymax = 2.5 ):

    y7 = np.linspace(0, rOuter45, nQuarter+1)
    y7 = y7 + rOuter45/nQuarter/2

    y8 = np.logspace(np.log10(y7[-1]), np.log10(ymax-0.05), num=nyOuter)
    points_new = np.concatenate((y7, y8[1:]))

    flipped_points = -np.flip(points_new)
    points_final = np.concatenate((flipped_points, points_new))
    return points_final

def interpolateCentralFlow(F, points, points_mesh, scale, velocity = False):
    points = np.array(points['y'])
    n = F.shape[1]
    name_list = F.columns.tolist()

    start = np.argmax(points > 0.5)
    stop = np.argmin(points < 1.51)

    points = (points[start:stop] - 1)*5
    
    data = {}
    for name in name_list:
        array = np.array(F[name].iloc[start:stop])
        interpolated_array = np.interp(points_mesh, points, array) * scale
        data[name] = interpolated_array
    
    if velocity == True:
        data = {
            'x': np.ones(len(points_mesh)),
            'y': np.zeros(len(points_mesh)),
            'z': np.zeros(len(points_mesh))
        }
    
    return data



def save_output(data, filepath, length=False):
    # Check if the input is a dictionary
    if not isinstance(data, dict):
        raise ValueError("Input data must be a dictionary.")

    # Transpose the data so that rows correspond to points
    rows = zip(*data.values())  # Transpose the dictionary values

    # Save data to file
    with open(filepath, "w") as f:
        f.write("(\n")  # Open parenthesis for the block
        for row in rows:
            # Format each row (triplet of x, y, z)
            formatted_row = "({})".format(" ".join(f"{val:.4e}" for val in row))
            if length:  # Handle the `length` parameter condition
                formatted_row = "{}".format(" ".join(f"{val:.4e}" for val in row))
            f.write(f"{formatted_row}\n")
        f.write(")\n")  # Close parenthesis for the block


def plot_original(points, U, L, R):

    fig1, axs1 = plt.subplots(1, 3, num=2, figsize=(10,3.5), clear=True)
    axs1[0].plot(points['y'], U['x'])
    axs1[0].set_ylabel("horizontal velocity U")
    axs1[0].axvline(x=0.5, color='k', linestyle='--')
    axs1[0].axvline(x=1, color='k', linestyle='--')
    axs1[0].set_xlabel("y")
    axs1[0].grid(True)

    axs1[1].plot(points['y'], L['l'])
    axs1[1].axvline(x=0.5, color='k', linestyle='--')
    axs1[1].axvline(x=1, color='k', linestyle='--')
    axs1[1].set_ylabel("length scale L")
    axs1[1].set_xlabel("y")
    axs1[1].grid(True)

    axs1[2].plot(points['y'], R['xx'], label='xx')
    axs1[2].plot(points['y'], R['yy'], label='yy')
    axs1[2].plot(points['y'], R['zz'], label='zz')
    axs1[2].plot(points['y'], R['xy'], label='xy')
    axs1[2].plot(points['y'], R['xz'], label='xz')
    axs1[2].plot(points['y'], R['yz'], label='yz')
    axs1[2].axvline(x=0.5, color='k', linestyle='--')
    axs1[2].axvline(x=1, color='k', linestyle='--')
    axs1[2].set_ylabel("Reynolds Stress Tensor R")
    axs1[2].set_xlabel("y")
    axs1[2].legend()
    axs1[2].grid(True)

    plt.tight_layout()
    #plt.savefig('given_points.pdf', format='pdf')
    plt.show()


def main():
    import pandas as pd
    import matplotlib.pyplot as plt
    import numpy as np

    scale = 1
    
    points = open_boundary('points_orig')
    U      = open_boundary('0/U_orig')
    L      = open_boundary('0/L_orig')
    R      = open_boundary('0/R_orig')

    #plot_original(points, U, L, R)

    points_mesh = computeCentroidGrid()

    U_new = interpolateCentralFlow(U, points, points_mesh, scale, velocity=True)
    L_new = interpolateCentralFlow(L, points, points_mesh, scale)
    R_new = interpolateCentralFlow(R, points, points_mesh, scale)
    points_out = {
        'x':np.zeros(len(points_mesh)),
        'y':points_mesh,
        'z':np.zeros(len(points_mesh))
    }

    save_output(points_out, "points")
    save_output(U_new, "0/U")
    save_output(L_new, "0/L", length = True)
    save_output(R_new, "0/R")

    #plot_original(points_out, U_new, L_new, R_new)


main()

