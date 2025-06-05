import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from numba import njit

L = 2e4
h_estuary = 10
h_river = 3
w_river = 200
w_estuary = 2000
Q_river = 800.0 

tide_amplitude = 3.0
tide_period = 12 * 3600 + 25 * 60


n_manning = 0.01

nx = 500
dx = L / nx
x = np.linspace(0, L, nx)

g = 9.81
cfl = 0.2
dt = cfl * dx / (np.sqrt(g * h_estuary + tide_amplitude))
tmax = 2 * tide_period
t = 0.0


w = w_river + (w_estuary - w_river) * np.exp(-4 * x / L)
zb = -h_river + (h_river - h_estuary) * (1 - x / L)**2
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8))

# ax1.plot(x/1000, w, 'b-', linewidth=2)
# ax1.set_ylabel("Largeur (m)")
# ax1.set_title("Évolution de la largeur de la Garonne (estuaire → Podensac)")
# ax1.grid(True)

# ax2.plot(x/1000, zb, 'r-', linewidth=2)
# ax2.set_ylabel("Profondeur (m)")
# ax2.set_xlabel("Distance (km)")
# ax2.set_title("Évolution de la profondeur")
# ax2.grid(True)

# plt.tight_layout()
# plt.show()



fig, ax = plt.subplots(figsize=(12, 5))
ax.set_xlim(0, L / 1000)
ax.set_ylim(np.min(zb)-1, 2*tide_amplitude)
ax.set_xlabel("Distance (km)")
ax.set_ylabel("Elevation (m)")
ax.grid(True, alpha=0.3)

line_water, = ax.plot(x / 1000, zb + h, 'b-', lw=2, label="Surface")
line_bed, = ax.plot(x / 1000, zb, 'k-', lw=1.5, label="Bed")
ax.legend()
ax.set_title("t = 0 h")

def update(frame):
    global h, u, h_new, u_new, t

    
    t += dt

    line_water.set_ydata(h)
    ax.set_title(f"Garonne Tidal Bore - t = {t/3600:.2f} h")
    return line_water

anim = FuncAnimation(fig, update, frames=int(tmax / (dt*5) ), interval=50, blit=True)
plt.show()
