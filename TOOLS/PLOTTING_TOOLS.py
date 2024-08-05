import matplotlib as plt


def orbitplot(body_data):

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    for bodies in body_data:
        ax.plot(bodies[:, 0], bodies[:, 1], bodies[:, 2], label=names(bodies)
        ax.plot([bodies[0, 0]], [bodies[0, 1]], [bodies[0, 2]], 'o')

    return plt.show()