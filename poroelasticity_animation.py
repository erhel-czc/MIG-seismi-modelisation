import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from matplotlib.widgets import Slider
from result import Result

# --- Paramètres du cercle ---
R = 1e3             # rayon



theta = np.linspace(0, 2 * np.pi, 100)

# --- Préparation de la figure ---
fig, axes = plt.subplots(2, 2, figsize=(10, 8))
plt.subplots_adjust(bottom=0.25)
sliprate, = axes[0,0].plot([], [], '-k')
sigma, = axes[0,1].plot([], [], '-k')
pressure, = axes[1,0].plot([], [], '-k')
position, = axes[1,1].plot([], [], 'bo', markersize=8)
plt.subplots_adjust(bottom=0.25)

axes[0,0].set_title("Slip rate evolution")
axes[0,1].set_title("Sigma evolution")
axes[1,0].set_title("Pressure evolution")
axes[1,1].set_title("Position par rapport à la source")

axes[1,1].set_xlim(-R-1, R+1)
axes[1,1].set_ylim(-R-1, R+1)
axes[1,1].set_aspect("equal")

slider_ax = fig.add_axes([0.2, 0.1, 0.6, 0.05])
angle_slider = Slider(
    slider_ax,
    "Angle (rad)",
    0,
    99,
    valinit=0,
    valstep=1,
)
axes[1,1].plot(R*np.cos(theta), R*np.sin(theta), 'gray', alpha=0.5)

def update(index):
    path = f'Results/Poroelasticity/100/{theta[index]}'
    data = Result.load_results(path)
    T =data.T * data.pd.dc / data.pd.V_p
    sliprate.set_data(T, data.Phi)
    sigma.set_data(T, data.Sigma_n)
    pressure.set_data(T, data.P)
    x = R * np.cos(theta[index])
    y = R * np.sin(theta[index])
    position.set_data([x], [y])
    for ax in axes.flatten():
        ax.relim()
        ax.autoscale_view()
    fig.canvas.draw_idle()


angle_slider.on_changed(update)
plt.show()


