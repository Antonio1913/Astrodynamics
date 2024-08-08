import matplotlib.pyplot as plt
import numpy as np
import matplotlib.animation as animation
from mpl_toolkits.mplot3d import Axes3D


def orbitplot(body_data, names, AU = False):

    #Defining the Conversion between AU and Km
    km2AU = 1.496*10**8

    # Define desired dimensions in pixels
    width_in_pixels = 1200
    height_in_pixels = 600
    desired_dpi = 100  # Increase DPI for better resolution

    # Convert pixels to inches
    width_in_inches = width_in_pixels / desired_dpi
    height_in_inches = height_in_pixels / desired_dpi

    # Setting the Global Properties of the Figure
    # Setting Surrounding Background of Plot to Black
    plt.rcParams['axes.facecolor'] = 'black'

    # Setting the Type of Font and Font Size
    plt.rcParams['font.family'] = 'serif'
    plt.rcParams['font.serif'] = 'Times New Roman'
    plt.rcParams['font.size'] = 12

    # Adjusted the Characteristics of the Figure
    fig = plt.figure(facecolor='black', figsize=(width_in_inches, height_in_inches), dpi=desired_dpi)
    ax = fig.add_subplot(111, projection='3d')
    ax.set_facecolor('black')

    #Initializing Limits
    max_xlim = [-float('inf'), float('inf')]
    max_ylim = [-float('inf'), float('inf')]
    max_zlim = [-float('inf'), float('inf')]

    # Plotting the Inputted Data with Label and Marking Initial Position
    if isinstance(body_data, list):

        if AU:
            for i, bodies in enumerate(body_data):
                ax.plot(bodies[:, 0] / km2AU, bodies[:, 1] / km2AU, bodies[:, 2] / km2AU, label=names[i])
                ax.plot([bodies[0, 0] / km2AU], [bodies[0, 1] / km2AU], [bodies[0, 2] / km2AU], 'o')

                # Finding Limit of One Set of Data
                set_xlim, set_ylim, set_zlim = equalaxis(bodies)

                # Updating Limit To Account For All Groups of Data
                max_xlim = [max(max_xlim[0], set_xlim[0]), min(max_xlim[1], set_xlim[1])]
                max_ylim = [max(max_ylim[0], set_ylim[0]), min(max_ylim[1], set_ylim[1])]
                max_zlim = [max(max_zlim[0], set_zlim[0]), min(max_zlim[1], set_zlim[1])]


        else:
            for i, bodies in enumerate(body_data):
                ax.plot(bodies[:, 0], bodies[:, 1], bodies[:, 2], label=names[i])
                ax.plot([bodies[0, 0]], [bodies[0, 1]], [bodies[0, 2]], 'o')

                # Finding Limit of One Set of Data
                set_xlim, set_ylim, set_zlim = equalaxis(bodies)

                # Updating Limit To Account For All Groups of Data
                max_xlim = [max(max_xlim[0], set_xlim[0]), min(max_xlim[1], set_xlim[1])]
                max_ylim = [max(max_ylim[0], set_ylim[0]), min(max_ylim[1], set_ylim[1])]
                max_zlim = [max(max_zlim[0], set_zlim[0]), min(max_zlim[1], set_zlim[1])]

    else:
        if AU:
            ax.plot(body_data[:, 0] / km2AU, body_data[:, 1] / km2AU, body_data[:, 2] / km2AU, label=names)
            ax.plot([body_data[0, 0] / km2AU], [body_data[0, 1] / km2AU], [body_data[0, 2] / km2AU], 'o')

            # Finding Limit of One Set of Data
            set_xlim, set_ylim, set_zlim = equalaxis(body_data)

            # Updating Limit To Account For All Groups of Data
            max_xlim = [max(max_xlim[0], set_xlim[0]), min(max_xlim[1], set_xlim[1])]
            max_ylim = [max(max_ylim[0], set_ylim[0]), min(max_ylim[1], set_ylim[1])]
            max_zlim = [max(max_zlim[0], set_zlim[0]), min(max_zlim[1], set_zlim[1])]


        else:
            ax.plot(body_data[:, 0], body_data[:, 1], body_data[:, 2], label=names)
            ax.plot([body_data[0, 0]], [body_data[0, 1]], [body_data[0, 2]], 'o')

            # Finding Limit of One Set of Data
            set_xlim, set_ylim, set_zlim = equalaxis(body_data)

            # Updating Limit To Account For All Groups of Data
            max_xlim = [max(max_xlim[0], set_xlim[0]), min(max_xlim[1], set_xlim[1])]
            max_ylim = [max(max_ylim[0], set_ylim[0]), min(max_ylim[1], set_ylim[1])]
            max_zlim = [max(max_zlim[0], set_zlim[0]), min(max_zlim[1], set_zlim[1])]

    # Setting Axis Limits
    ax.set_xlim(max_xlim)
    ax.set_ylim(max_ylim)
    ax.set_zlim(max_zlim)

    # Adjusts the Characteristics of the Plot
    if AU:
        ax.set_xlabel('X (AU)', color='white')
        ax.set_ylabel('Y (AU)', color='white')
        ax.set_zlabel('Z (AU)', color='white')

    else:
        ax.set_xlabel('X (Km)', color='white')
        ax.set_ylabel('Y (Km)', color='white')
        ax.set_zlabel('Z (Km)', color='white')

    ax.tick_params(axis='x', colors='white')
    ax.tick_params(axis='y', colors='white')
    ax.tick_params(axis='z', colors='white')

    # Adjusted the Characteristics of the Legend
    ax.legend(facecolor='white', bbox_to_anchor=(1.05, 0.5))

    # Displays the Plot
    plt.show(block=True)

def equalaxis(body_data):
    # Set equal scaling
    max_range = np.array([body_data[:, 0].max() - body_data[:, 0].min(),
                          body_data[:, 1].max() - body_data[:, 1].min(),
                          body_data[:, 2].max() - body_data[:, 2].min()]).max() / 2.0

    mid_x = (body_data[:, 0].max() + body_data[:, 0].min()) * 0.5
    mid_y = (body_data[:, 1].max() + body_data[:, 1].min()) * 0.5
    mid_z = (body_data[:, 2].max() + body_data[:, 2].min()) * 0.5

    set_xlim = (mid_x - max_range, mid_x + max_range)
    set_ylim = (mid_y - max_range, mid_y + max_range)
    set_zlim = (mid_z - max_range, mid_z + max_range)

    return set_xlim, set_ylim, set_zlim






