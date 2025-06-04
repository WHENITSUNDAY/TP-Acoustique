import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, FFMpegWriter
import time
from numba import njit, prange
"""Ce code a pour but de simuler la propagation d'ondes acoustiques dans un milieu bidimensionnel entre deux dispositifs piézoélectriques, l'un émetteur et l'autre récepteur.
On considère un premier cas où le milieu est seulement composé d'où, puis un second cas où on ajoute une plaque d'aluminium entre les deux dispositifs piézoélectriques.
On peut donc visualiser la propagation complexe de l'onde acoustique dans l'eau, puis dans l'eau et l'aluminium (ce que l'on ne peut pas faire pendant le TP ! Car on a seulement accès au signal reçu)
En parallèle, on modélise notamment le signal électrique reçu par le récepteur piézoélectrique, qui est proportionnel à la pression acoustique moyenne dans la zone du récepteur.
Cela permet donc de comparer (visuellement) avec les résultats expérimentaux obtenus pendant le TP au Fablab.
Le but est donc de montrer aux élèves que la simulation numérique est un outil puissant pour comprendre et visualiser des phénomènes physiques complexes, 
parfois difficilement accessibles par l'expérimentation directe, et ce avec un simple code Python d'une centaine de lignes.

La propagation des ondes acoustiques est modélisée par les équations d'Alembert. On considère les hypothèses suivantes :
-Milieu homogène, isotrope, linéaire et élastique.
-Propagation bidimensionnelle de l'onde acoustique
-Milieu non-dispersif : l'eau est considérée comme un fluide parfait (pas de viscosité), l'aluminium comme un solide parfaitement élastique uniquement soumis à des ondes longitudinales.
Il faut bien prendre en compte que la propagation de l'onde dans l'aluminium n'est pas fidèlement modélisée par les équations d'Alembert, mais cela simplifie grandement la simulation (plus accessible pour les élèves).

Pour le côté numérique, on utilise la méthode des différences finies pour résoudre les équations d'Alembert, avec un schéma explicite le plus simple possible (Ordre 2 en espace et temps).
On utilise des conditions aux limites absorbantes de type PML pour éviter les réflexions aux bords du domaine (qu'on peut théoriquement considérer comme infini en prenant en compte la dimension de l'aquarium relativement au domaine)
L'utilisation d'un filtre passe-bas est primordiale pour lisser les artefacts de dispersion numérique, ce qui clarifie la visualisation de l'onde."""

start = time.time()

# Paramètres physiques
has_aluminium = False # True si on veut ajouter la plaque d'alu, False sinon
alu_distance = 0.03 # Distance entre l'émetteur et la plaque d'aluminium (m)
Lx, Ly = 0.08, 0.08 # Dimensions du domaine (m)
piezo_diameter = 0.025 # Diamètre des dispositifs piézoélectriques (m)
alu_thickness = 0.0083 # Épaisseur de la plaque d'aluminium (m)
alu_width = 0.067 # Largeur de la plaque d'aluminium (m)
c_w = 1482 # Vitesse du son dans l'eau (m/s)
c_a = 6320 # Vitesse du son (longitudinal) dans l'aluminium (m/s) (on néglige le cisaillement)

f_source = 2e6 # Fréquence de la source (Hz), 2MHz -> Ultrasons.
alpha = 0.010
rho_w = 1000 
rho_a = 2700

# Paramètres numériques
t = 0
tmax = 1.8e-4
nx, ny = 128, 128
dx = Lx / nx
dy = Ly / ny
dt = 0.95 * dx / (c_a * np.sqrt(2))  # CFL


piezo_width = int(ny * (piezo_diameter / Ly))
piezo_bottom = ny // 2 - piezo_width // 2

if has_aluminium:
    ampli = 2.5
else:
    ampli = 1.5

# Paramètres de la PML
pml_width = nx // 5
sigma_max = 5e5 

sigma = np.zeros((nx, ny), dtype=np.float32)

i = np.arange(pml_width)
j_top = np.arange(0, piezo_bottom)
j_bot = np.arange(piezo_bottom + piezo_width, ny)

if piezo_bottom > 0:
    sigma[np.ix_(i, j_top)] = sigma_max * ((pml_width - i)[:, None] / pml_width) * ((piezo_bottom - j_top)[None, :] / piezo_bottom)
    sigma[np.ix_(nx-1-i, j_top)] = sigma_max * ((pml_width - i)[:, None] / pml_width) * ((piezo_bottom - j_top)[None, :] / piezo_bottom)

if ny - (piezo_bottom + piezo_width) > 0:
    sigma[np.ix_(i, j_bot)] = sigma_max * ((pml_width - i)[:, None] / pml_width) * ((j_bot[None, :] - (piezo_bottom + piezo_width - 1)) / (ny - (piezo_bottom + piezo_width)))
    sigma[np.ix_(nx-1-i, j_bot)] = sigma_max * ((pml_width - i)[:, None] / pml_width) * ((j_bot[None, :] - (piezo_bottom + piezo_width - 1)) / (ny - (piezo_bottom + piezo_width)))

j = np.arange(pml_width)
sigma[:, j] += sigma_max * ((pml_width - j) / pml_width) ** 3
sigma[:, ny-1-j] += sigma_max * ((pml_width - j) / pml_width) ** 3

"""plt.figure(figsize=(6, 5))
plt.imshow(sigma.T, origin='lower', cmap='viridis', extent=[0, nx*dx, 0, ny*dy])
plt.colorbar(label='Intensité de sigma')
plt.title('Carte de la matrice sigma (PML)')
plt.xlabel('x (m)')
plt.ylabel('y (m)')
plt.tight_layout()
plt.show()"""

# Initialisation des tableaux
p = np.zeros((nx, ny), dtype=np.float32)
p_prev = np.zeros((nx, ny), dtype=np.float32)
p_next = np.zeros((nx, ny), dtype=np.float32)
celerity_matrix = np.full((nx, ny), c_w**2, dtype=np.float32)  # Par défaut, on est uniquement dans l'eau...

inv_dx2 = 1.0 / dx**2
inv_dy2 = 1.0 / dy**2
dt2 = dt**2

#Si on ajoute l'aluminium, on modifie la matrice de célérité aux points où l'aluminium est présent.
if has_aluminium:
    alu_x0 = int(alu_distance / Lx * nx)
    alu_y0 = int((Ly - alu_width) / 2 / Ly * ny)
    alu_x1 = alu_x0 + int(alu_thickness / Lx * nx)
    alu_y1 = alu_y0 + int(alu_width / Ly * ny)
    celerity_matrix[alu_x0:alu_x1, alu_y0:alu_y1] = c_a**2

# Configuration du plot (partie pas intéressante pour les élèves)
fig = plt.figure(figsize=(8, 8))
gs = fig.add_gridspec(2, 1, height_ratios=[5, 1], width_ratios=[1], hspace=0.2, top=0.92, bottom=0.08, left=0.15, right=0.98)

ax1 = fig.add_subplot(gs[0])
im = ax1.imshow(p.T, cmap='magma', extent=[0, nx * dx , 0, ny * dy], vmin=-1, vmax=1)
plt.colorbar(im, ax=ax1, label='Pression acoustique (Pa)')
ax1.set_xlabel('x (m)')
ax1.set_ylabel('y (m)')
ax1.set_title('Propagation d\'ondes acoustiques avec PML')


if has_aluminium:
    alu_patch = plt.Rectangle((alu_x0 * dx, alu_y0 * dy), alu_thickness, alu_width, fill=False, color='black', label='Plaque d\'Aluminium')
    ax1.add_patch(alu_patch)
emitter = plt.Line2D([0, 0], [piezo_bottom * dy, (piezo_bottom + piezo_width) * dy], color="#FFDD6D", linewidth=6, label='Émetteur')
receiver = plt.Line2D([nx * dx, nx * dx], [piezo_bottom * dy, (piezo_bottom + piezo_width) * dy], color="#3E177E", linewidth=6, label='Récepteur')

ax1.add_line(emitter)
ax1.add_line(receiver)
ax1.legend(loc='upper right')

ax2 = fig.add_subplot(gs[1])
ax2.set_xlim(0, tmax*1e6)
ax2.set_ylim(-1, 1)
ax2.set_xlabel('Temps (µs)')
ax2.set_ylabel('Signal reçu (normalisé)')
ax2.grid(True)


pos1 = ax1.get_position()
pos2 = ax2.get_position()
ax2.set_position([pos1.x0, pos2.y0, pos1.width, pos2.height])

time_points = []
signal_points = []
line, = ax2.plot([], [], color='purple', lw=2, markersize=4)

print("Init:", time.time() - start)

""" Fonctions principales du code"""
def source(t):
    pulse_duration = 3 / f_source
    if t < pulse_duration:
        window = 0.5 * (1 - np.cos(2 * np.pi * t / pulse_duration))
        return 700 * window * np.sin(2 * np.pi * f_source * t)
    return 0

@njit(fastmath=True, cache=True)
def step_numba(p, p_prev, p_next, celerity_matrix, sigma, alpha, dt, dt2, inv_dx2, inv_dy2, piezo_bottom, piezo_width, source_val):
    nx, ny = p.shape
    # Source
    p[0, piezo_bottom:piezo_bottom + piezo_width] = source_val

    # Différences finies et PML
    for i in range(1, nx-1):
        for j in range(1, ny-1):
            laplacian = (
                (p[i+1, j] + p[i-1, j] - 2 * p[i, j]) * inv_dx2 +
                (p[i, j+1] + p[i, j-1] - 2 * p[i, j]) * inv_dy2
            )
            p_next[i, j] = (
                (2 * p[i, j] - p_prev[i, j] +
                 (celerity_matrix[i, j] * dt2) * laplacian) *
                np.exp(-sigma[i, j] * dt)
            )
    # Passe-bas 
    for i in range(1, nx-1):
        for j in range(1, ny-1):
            p_next[i, j] += alpha * (
                p_next[i+1, j] + p_next[i-1, j] +
                p_next[i, j+1] + p_next[i, j-1] -
                4 * p_next[i, j]
            )
    received = np.mean(p_next[-2, piezo_bottom:piezo_bottom + piezo_width])
    return p_next, received

def update(frame):
    global p, p_prev, p_next, t, time_points, signal_points

    source_val = source(t)
    p_next.fill(0)
    p_next, received = step_numba(
        p, p_prev, p_next, celerity_matrix, sigma, alpha, dt, dt2, inv_dx2, inv_dy2,
        piezo_bottom, piezo_width, source_val
    )

    np.copyto(p_prev, p)
    np.copyto(p, p_next)

    time_points.append(t*1e6)
    signal_points.append(received * ampli)

    if frame % frame_step == 0:
        im.set_array(p.T)
        ax1.set_title(f"Propagation d'ondes acoustiques avec PML\nTemps : {t*1e6:.2f} µs")
        line.set_data(np.array(time_points), np.array(signal_points))

    t += dt

    return im, line

nframes = int(tmax / dt)
fps = 100
frame_step = 1

anim = FuncAnimation(fig, update, frames=nframes, interval=1000/fps, blit=False)
writer = FFMpegWriter(fps=fps, bitrate=1800)

if has_aluminium:
    anim.save('onde_acoustique_eau_alu.mp4', writer=writer, dpi = 200)
else:
    anim.save('onde_acoustique_eau.mp4', writer=writer, dpi = 200)

print("Avec Numba:", time.time() - start)
plt.close()