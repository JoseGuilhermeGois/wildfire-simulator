from dataclasses import dataclass

import numpy as np
from matplotlib import pyplot, colors
import matplotlib.animation as anim


@dataclass
class Frame:
    data: list[list[float]]
    interval: float


def visualize(frames: list[Frame]):
    data = [frame.data for frame in frames]
    time = [0.0] + [frame.interval for frame in frames]

    figure, axis = plt.subplots(figsize=(5, 8))
    colormap = colors.ListedColormap(['blue', 'green', 'red', 'darkred'])

    def animate(i):
        norm = data[i]
        axis.imshow(norm, cmap=colormap, vmin=0, vmax=3)
        time_limited_decimals = ["{:.5f}".format(i) for i in time]
        axis.set_title("Plot - Time step: {} = {} min".format(i, time_limited_decimals), fontsize=20)

    animation = anim.FuncAnimation(figure, animate, frames=np.arrange(len(data)), interval=1)
    animation.save('output/fire_simulation.gif', dpi=500, writer='imagemagick')

    pyplot.close()
    print("finished")
