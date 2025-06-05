import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from numba import njit

Lx = 1e5
h_sea = 10
h_riv = 3
h_static = h_riv
tide_amplitude = 3.0
tide_period = 12 * 3600 + 25 * 60
Q_river = 1000.0
river_width = 200.0
n_manning = 0.01

nx = 200
dx = Lx / nx
x = np.linspace(0, Lx, nx)

g = 9.81
cfl = 0.9
dt = cfl * dx / (np.sqrt(g * h_sea + tide_amplitude) + 10)
tmax = 2 * tide_period
t = 0.0

@njit
def bathymetry(x):
    zb = (h_sea - h_riv) * (x / Lx)
    zb += 0.5 * np.exp(-((x - 0.7 * Lx) / (0.1 * Lx)) ** 2)
    return -zb

zb = bathymetry(x)

h = h_static * np.ones(nx)
q = np.zeros(nx)


fig, ax = plt.subplots(figsize=(12, 5))
ax.set_xlim(0, Lx / 1000)
ax.set_ylim(np.min(zb)-1, h_static + 2*tide_amplitude)
ax.set_xlabel("Distance (km)")
ax.set_ylabel("Elevation (m)")
ax.grid(True, alpha=0.3)

line_water, = ax.plot(x / 1000, zb + h, 'b-', lw=2, label="Surface")
line_bed, = ax.plot(x / 1000, zb, 'k-', lw=1.5, label="Bed")
ax.legend()
title = ax.set_title("t = 0 h")

@njit
def step_lax_friedrichs(h, q, zb, t, dx, dt):
    nx = len(h)
    h_new = np.copy(h)
    q_new = np.copy(q)

    F1 = q
    F2 = q**2 / h + 0.5 * g * h**2

    for i in range(1, nx - 1):
        h_new[i] = 0.5 * (h[i+1] + h[i-1]) - dt/(2*dx) * (F1[i+1] - F1[i-1])
        q_new[i] = 0.5 * (q[i+1] + q[i-1]) - dt/(2*dx) * (F2[i+1] - F2[i-1])
        dzb = (zb[i+1] - zb[i-1]) / (2*dx)
        q_new[i] -= dt * g * h[i] * dzb

        h_f = max(h[i], 0.1)
        q_new[i] -= dt * g * n_manning**2 * q[i] * abs(q[i]) / h_f**(10/3)

    eta_tide = tide_amplitude * np.sin(2 * np.pi * t / tide_period)
    h_new[-1] = h_static + eta_tide
    q_new[-1] = q_new[-2]
    q_new[0] = q_new[1]
    h_new[0] = h_new[1]
    return h_new, q_new

def update(frame):
    global h, q, t
    for _ in range(5):
        h, q = step_lax_friedrichs(h, q, zb, t, dx, dt)
        t += dt

    line_water.set_ydata(h)
    title.set_text(f"Garonne Tidal Bore - t = {t/3600:.2f} h")
    return line_water, title

anim = FuncAnimation(fig, update, frames=int(tmax / dt / 5), interval=50, blit=True)
plt.show()
