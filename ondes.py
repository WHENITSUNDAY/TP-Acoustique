import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, FFMpegWriter
from itertools import chain

# Paramètres physiques
Lx, Ly = 0.10, 0.08
piezo_diameter = 0.025
alu_thickness = 0.0083
alu_width = 0.067

c_w = 1500
c_a = 5000

rho_w = 1000
rho_a = 2700

f_source = 2e6

# Paramètres numériques
t = 0
nx, ny = 100, 80
dx = Lx / nx
dy = Ly / ny
dt = 0.9 * dx / (c_a * np.sqrt(2))  # CFL

piezo_width = int(ny * (piezo_diameter / Ly))
piezo_left = ny // 2 - piezo_width // 2
piezo_right = ny // 2 - piezo_width // 2

# Paramètres de la PML
#sponge_width = nx // 10
#damping_strength = 0.2

pml_width = 15 
sigma_max = 1e6 

sigma = np.zeros((nx, ny))

for j in range(ny):
    if j < piezo_left or j >= piezo_left + piezo_width:
        for i in range(pml_width):
            sigma[i, j] = sigma_max * ((pml_width - i) / pml_width)**3
            sigma[nx-1-i, j] = sigma_max * ((pml_width - i) / pml_width)**3

for j in range(pml_width):
    sigma[:, j] = sigma_max * ((pml_width - j) / pml_width)**3
    sigma[:, ny-1-j] = sigma_max * ((pml_width - j) / pml_width)**3

# Initialisation des tableaux
p = np.zeros((nx, ny))
p_prev = np.zeros((nx, ny))

print(sigma)
celerity_matrix = np.full((nx, ny), c_w)  # Par défaut, on est uniquement dans l'eau

has_aluminium = False  # True si on veut ajouter la plaque d'alu, False sinon

if has_aluminium:
    alu_distance = 0.05
    alu_x0 = int(alu_distance / Lx * nx)
    alu_y0 = int((Ly - alu_width) / 2 / Ly * ny)
    alu_x1 = alu_x0 + int(alu_thickness / Lx * nx)
    alu_y1 = alu_y0 + int(alu_width / Ly * ny)

    celerity_matrix[alu_x0:alu_x1, alu_y0:alu_y1] = c_a


fig, ax = plt.subplots(figsize=(10, 5))
im = ax.imshow(p.T, cmap='viridis', extent=[0, nx * dx, 0, ny * dy], vmin=-2, vmax=2)
plt.colorbar(im, label='Pression acoustique (Pa)')
ax.set_xlabel('x (m)')
ax.set_ylabel('y (m)')
ax.set_title('Propagation d\'ondes acoustiques avec PML')

if has_aluminium:
    alu_patch = plt.Rectangle((alu_x0 * dx, alu_y0 * dy), alu_thickness, alu_width, fill=False, color='black', alpha=1, label='Plaque d\'Aluminium')
    ax.add_patch(alu_patch)

emitter = plt.Line2D([0, 0], [piezo_left * dy, (piezo_left + piezo_width) * dy], color='yellow', linewidth=6, label='Émetteur')
receiver = plt.Line2D([nx * dx, nx * dx], [piezo_left * dy, (piezo_left + piezo_width) * dy], color='green', linewidth=6, label='Récepteur')
ax.add_line(emitter)
ax.add_line(receiver)

ax.legend()

def source(t):
    if t < 0.5 / f_source:
        return 20 * np.sin(2 * np.pi * f_source * t)
    else:
        return 0
def update(frame):
    global p, p_prev, t

    p_next = np.zeros((nx, ny))
    if frame == 0:
        for j in range(piezo_left, piezo_left + piezo_width):
            p_next[0, j] = source(t)

    for i in range(1, nx - 1):
        for j in range(1, ny - 1):
            laplacian = (p[i + 1, j] - 2 * p[i, j] + p[i - 1, j]) / dx**2 + (p[i, j + 1] - 2 * p[i, j] + p[i, j - 1]) / dy**2
            celerity = celerity_matrix[i, j]
            p_next[i, j] = (2 * p[i, j] - p_prev[i, j] + (celerity**2 * dt**2) * laplacian)

    # Conditions aux limites
    for i in range(nx):
        for j in range(ny):
            if sigma[i,j] > 0:
                p_next[i, j] *= np.exp(-sigma[i,j] * dt)

    p_prev, p = p, p_next
    t += dt

    im.set_array(p.T)
    im.set_clim(vmin=-0.5, vmax=0.5)
    return im

ani = FuncAnimation(fig, update, frames=200, interval=40, blit=False)

plt.show()
#writer = FFMpegWriter(fps=20, metadata=dict(artist='Me'), bitrate=1800)
#ani.save('acoustic_wave_propagation.mp4', writer=writer)
#plt.close()
