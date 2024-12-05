# compute turbulence intensity
'''
defined as U_rms / U
with U_rms = np.sqrt((u_prime_x**2 + u_prime_y**2 + u_prime_z**2)/3)
'''
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

plt.rcParams.update({'font.size': 14})  # Change 14 to your desired font size



def open_boundary(file_path):

    with open(file_path, 'r') as file:
        lines = file.readlines()
    data_lines = lines[:-1]

    last = data_lines[-1]
    last = last.strip().split()
    nprobe = int((len(last)-1)/3)

    data = []
    for line in data_lines:
        parts = line.strip().split()
        if parts[0] == '#':
            continue  # Skip comment lines

        # Extract time
        time = float(parts[0])
        entry = [time]

        for probe_idx in range(nprobe):
            x = list(map(float, parts[probe_idx*3+1].strip('()').split()))
            y = list(map(float, parts[probe_idx*3+2].strip('()').split()))
            z = list(map(float, parts[probe_idx*3+3].strip('()').split()))

            entry.extend(x + y + z)
        data.append(entry)
    
    name_list = ['time']
    for probe_idx in range(nprobe):
        name_list.extend([f'x{probe_idx}', f'y{probe_idx}', f'z{probe_idx}'])

    # Create a DataFrame from the processed data
    df = pd.DataFrame(data, columns=name_list)

    return df


    

def compute_intensity(U, start):

    filtered_U = U[U['time'] > start]
    time = np.array(filtered_U['time'])

    result = {
        'time':time
    }

    nprobe = int((U.shape[1]-1)/3)
    for i in range(nprobe):
        xname = f'x{i}'
        yname = f'y{i}'
        zname = f'z{i}'
        iname = f'probe{i}'
        U_x = np.array(filtered_U[xname])
        U_y = np.array(filtered_U[yname])
        U_z = np.array(filtered_U[zname])
        u_prime_x = U_x - np.mean(U_x)
        u_prime_y = U_y - np.mean(U_y)
        u_prime_z = U_z - np.mean(U_z)

        U_rms = np.sqrt((u_prime_x**2 + u_prime_y**2 + u_prime_z**2)/3)
        intensity = U_rms / U_x

        result[iname] = intensity*100

    df = pd.DataFrame(result)
    return df


def main():
    start_time = 5

    # U = open_boundary('U_sma_inter')
    # intensity = compute_intensity(U, start_time)
    # print(f'intensity central: {np.mean(intensity['probe0']):.3f}%')
    # print(f'intensity top    : {np.mean(intensity['probe1']):.3f}%')
    # print(f'intensity bottom : {np.mean(intensity['probe2']):.3f}%')

    #data_name = ['U_keq_inter', 'U_keq_value', 'U_sma_inter', 'U_sma_value', 'U_keq_final', 'u_sma_final']
    data_name = ['U_keq_inter', 'U_sma_inter', 'U_keq_value', 'U_smaRef']

    fig, axs = plt.subplots(3, len(data_name), num=2, figsize=(12,6), clear=True)
    for idx, name in enumerate(data_name):
        U = open_boundary(name)
        intensity = compute_intensity(U, start_time)

        mean0 = np.mean(intensity['probe0'])
        mean1 = np.mean(intensity['probe1'])
        mean2 = np.mean(intensity['probe2'])

        axs[1,idx].plot(intensity['time'], intensity['probe0'], color='b')
        axs[1,idx].axhline(y=mean0, color='k', linestyle='--')
        axs[1,idx].text(
            x=0.5 * (intensity['time'].min() + intensity['time'].max()),  # Midpoint of the x-range
            y=mean0 + 0.1,  # Slightly above the line
            s=f'{mean0:.3f}%',
            color='k',
            ha='center',  # Horizontal alignment
            bbox=dict(facecolor='white', boxstyle='round,pad=0.5')  # White background box
        )

        axs[0,idx].plot(intensity['time'], intensity['probe1'], color='g')
        axs[0,idx].axhline(y=mean1, color='k', linestyle='--')
        axs[0,idx].text(
            x=0.5 * (intensity['time'].min() + intensity['time'].max()),  # Midpoint of the x-range
            y=mean1 + 0.1,  # Slightly above the line
            s=f'{mean1:.3f}%',
            color='k',
            ha='center',  # Horizontal alignment
            bbox=dict(facecolor='white', boxstyle='round,pad=0.5')  # White background box
        )
        axs[2,idx].plot(intensity['time'], intensity['probe2'], color='r')
        axs[2,idx].axhline(y=mean2, color='k', linestyle='--')
        axs[2,idx].text(
            x=0.5 * (intensity['time'].min() + intensity['time'].max()),  # Midpoint of the x-range
            y=mean2+ 0.1,  # Slightly above the line
            s=f'{mean2:.3f}%',
            color='k',
            ha='center',  # Horizontal alignment
            bbox=dict(facecolor='white', boxstyle='round,pad=0.5')  # White background box
        )

        #axs[idx].plot(intensity['time'], np.ones(len(intensity['time']))*0, color='k')
        #axs[idx].plot(intensity['time'], np.ones(len(intensity['time']))*(+offset), color='k')
        #axs[idx].plot(intensity['time'], np.ones(len(intensity['time']))*(-offset), color='k')

        axs[0,idx].grid()
        axs[1,idx].grid()
        axs[2,idx].grid()

        axs[2,idx].set_xlabel("time")
        if idx == 0:
            axs[1,idx].set_ylabel("turbulence intensity [%]")
        axs[0,idx].set_title(name)

    plt.tight_layout()
    #plt.savefig('turbulent_intensity.pdf', format='pdf')
    plt.show()


if __name__ == "__main__":
    main()