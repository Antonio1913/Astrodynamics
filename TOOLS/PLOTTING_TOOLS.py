import matplotlib.pyplot as plt
import matplotlib.animation as animation
from mpl_toolkits.mplot3d import Axes3D


def orbitplot(body_data, names, AU = False, animate = False):

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

    # Plotting the Inputted Data with Label and Marking Initial Position
    if AU:
        for i, bodies in enumerate(body_data):
            ax.plot(bodies[:, 0] / km2AU, bodies[:, 1] / km2AU, bodies[:, 2] / km2AU, label=names[i])
            ax.plot([bodies[0, 0] / km2AU], [bodies[0, 1] / km2AU], [bodies[0, 2] / km2AU], 'o')
    else:
        for i, bodies in enumerate(body_data):
            ax.plot(bodies[:, 0], bodies[:, 1], bodies[:, 2], label=names[i])
            ax.plot([bodies[0, 0]], [bodies[0, 1]], [bodies[0, 2]], 'o')

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








