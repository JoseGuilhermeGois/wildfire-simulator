from dataclasses import dataclass

from matplotlib import colors
import matplotlib.pyplot as plt
import matplotlib.animation as anim


@dataclass
class Frame:
    data: list[list[int]]
    interval: float


def visualize(frames: list[Frame]):
    time_counter = 0
    time = []
    data = [frame.data for frame in frames]

    time = [sum(frame.interval for frame in frames[:i]) for i, _ in enumerate(frames)]

    fig, ax = plt.subplots(figsize=(5, 8))
    colormap = colors.ListedColormap(['blue', 'green', 'red', 'darkred'])

    def animate(i):
        norm = data[i]
        ax.imshow(norm, cmap=colormap, vmin=0, vmax=3)
        time_limited_decimals = ["{:.5f}".format(i) for i in time]
        ax.set_title("Plot - Time step: {} = {} min".format(i, time_limited_decimals[i]), fontsize=20)

    ax.set_xticks([])
    ax.set_yticks([])

    plt.xlabel("100 m")
    plt.ylabel("50 m")

    print("Doing the animation! Zzzz")

    animation = anim.FuncAnimation(fig, animate, frames=len(data), interval=5)

    animation.save('output/fire_simulation.gif', dpi=80, writer='pillow')

    print("finished")
