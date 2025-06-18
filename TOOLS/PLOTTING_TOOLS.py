
import matplotlib.pyplot as plt
import numpy as np


# body_data - should be a list

def orbitplot(body_data, names, AU=False):

    # Defining the Conversion between AU and Km
    km2AU = 1.496*10**8

    # Define desired dimensions in pixels
    width_in_pixels = 1200
    height_in_pixels = 600
    desired_dpi = 100  # Increase DPI for better resolution

    # Convert pixels to inches
    width_in_inches = width_in_pixels / desired_dpi
    height_in_inches = height_in_pixels / desired_dpi

    # Setting the Global Properties of the Figure
    # Setting the Type of Font and Font Size
    plt.rcParams['font.family'] = 'serif'
    plt.rcParams['font.serif'] = 'Times New Roman'
    plt.rcParams['font.size'] = 12

    # SETTING AXES CONFIGURATIONS
    plt.rcParams['axes3d.xaxis.panecolor'] = "black"
    plt.rcParams['axes3d.yaxis.panecolor'] = "black"
    plt.rcParams['axes3d.zaxis.panecolor'] = "black"
    plt.rcParams['axes.facecolor'] = "black"
    plt.rcParams['axes.labelcolor'] = "white"

    # SETTING LEGEND CONFIGURATIONS
    plt.rcParams['legend.labelcolor'] = "white"
    plt.rcParams['legend.facecolor'] = "black"

    # SETTING FIGURE CONFIGURATIONS
    plt.rcParams['figure.facecolor'] = "black"

    # CREATING THE FIGURE
    fig = plt.figure(figsize=(width_in_inches, height_in_inches), dpi=desired_dpi)
    ax = fig.add_subplot(111, projection='3d')

    # Initializing Limits
    max_xlim = [float('inf'), -float('inf')]
    max_ylim = [float('inf'), -float('inf')]
    max_zlim = [float('inf'), -float('inf')]

    # Plotting the Inputted Data with Label and Marking Initial Position
    if len(names) > 1:

        for i, bodies in enumerate(body_data):
            if AU:
                bodies /= km2AU
                ax.plot(bodies[:, 0], bodies[:, 1], bodies[:, 2], label=names[i])
                ax.plot(bodies[0, 0], bodies[0, 1], bodies[0, 2], 'o')
            else:
                ax.plot(bodies[:, 0], bodies[:, 1], bodies[:, 2], label=names)
                ax.plot(bodies[0, 0], bodies[0, 1], bodies[0, 2], 'o')

            # Finding Limit of One Set of Data
            set_xlim, set_ylim, set_zlim = equalaxis(bodies)

            # Updating Limit To Account For All Groups of Data
            max_xlim = [min(max_xlim[0], set_xlim[0]), max(max_xlim[1], set_xlim[1])]
            max_ylim = [min(max_ylim[0], set_ylim[0]), max(max_ylim[1], set_ylim[1])]
            max_zlim = [min(max_zlim[0], set_zlim[0]), max(max_zlim[1], set_zlim[1])]


    else:
        if AU:
            ax.plot(body_data[:, 0] / km2AU, body_data[:, 1] / km2AU, body_data[:, 2] / km2AU, label=names)
            ax.plot([body_data[0, 0] / km2AU], [body_data[0, 1] / km2AU], [body_data[0, 2] / km2AU], 'o')
        else:
            ax.plot(body_data[:, 0], body_data[:, 1], body_data[:, 2] / km2AU, label=names)
            ax.plot(body_data[0, 0], body_data[0, 1], body_data[0, 2], 'o')

        # Finding Limit of One Set of Data
        set_xlim, set_ylim, set_zlim = equalaxis(body_data)

        # Updating Limit To Account For All Groups of Data
        max_xlim = (max(max_xlim[0], set_xlim[0]), min(max_xlim[1], set_xlim[1]))
        max_ylim = (max(max_ylim[0], set_ylim[0]), min(max_ylim[1], set_ylim[1]))
        max_zlim = (max(max_zlim[0], set_zlim[0]), min(max_zlim[1], set_zlim[1]))

    # Setting Axis Limits
    ax.set_xlim(max_xlim)
    ax.set_ylim(max_ylim)
    ax.set_zlim(max_zlim)

    # Adjusts the Characteristics of the Plot
    if AU:
        ax.set_xlabel('X (AU)')
        ax.set_ylabel('Y (AU)')
        ax.set_zlabel('Z (AU)')

    else:
        ax.set_xlabel('X (Km)')
        ax.set_ylabel('Y (Km)')
        ax.set_zlabel('Z (Km)')

    # SETTING TICK MARK CONFIGURATIONS
    ax.tick_params(axis='x', colors='white')
    ax.tick_params(axis='y', colors='white')
    ax.tick_params(axis='z', colors='white')

    # SETTING LEGEND POSITION
    ax.legend(bbox_to_anchor=(1.65, 0.35))

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
