import matplotlib.pyplot as plt
import numpy as np
from matplotlib import animation
from matplotlib import colors


TOTAL_TIME = 100

from simulator.global_variables import TOTAL_TIME

figure, axis = plt.subplots(figsize=(5, 8))
colormap = colors.ListedColormap(['blue', 'green', 'red', 'darkred'])


def animate(i):
    norm = results[i]
    axis.imshow(norm, cmap=colormap, vmin=0, vmax=3)
    axis.set_title("Plot2 - Time step: {} = {} min".format(i, elapsed_time_collection[i]),
                   fontsize=20)  # TODO: elapsed time on the simulation, check it


axis.set_xticks([])
axis.set_yticks([])

plt.xlabel("67 m")
plt.ylabel("31 m")

animation = animation.FuncAnimation(figure, animate, frames=np.arange(0, TOTAL_TIME), interval=5)
animation.save('fire_simulation.gif', dpi=80, writer='imagemagick')

plt.close()
print("finished")
