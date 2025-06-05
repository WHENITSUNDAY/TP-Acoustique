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
has_aluminium = True # True si on veut ajouter la plaque d'alu, False sinon

Lx, Ly = 0.08, 0.08 # Dimensions du domaine (m)
piezo_diameter = 0.025 # Diamètre des dispositifs piézoélectriques (m)
alu_thickness = 0.0083 # Épaisseur de la plaque d'aluminium (m)
alu_width = 0.067 # Largeur de la plaque d'aluminium (m)
alu_distance = 0.03 # Distance entre l'émetteur et la plaque d'aluminium (m)
c_w = 1500 # Vitesse du son dans l'eau (m/s) (ou dans un fluide en général)
c_a = 6000 # Vitesse du son (longitudinal) dans l'aluminium (m/s) (on néglige le cisaillement) (ou dans un solide en général)


f_source = 2e6 # Fréquence de la source (Hz), 2MHz -> Ultrasons.
t = 0
tmax = 3e-4 # Durée de la simulation (s)

c_max = c_w if not has_aluminium else c_a

nx, ny = int(128*(Lx/Ly)), 128
dx = Lx / nx
dy = Ly / ny
dt = 0.8 * dx / (c_max * np.sqrt(2))  # CFL
alpha = 0.008 * dt / dx**2 
print (alpha)
#Optimisation : on précalcule les coefficients pour éviter de les recalculer à chaque itération
inv_dx2 = 1.0 / dx**2
inv_dy2 = 1.0 / dy**2
dt2 = dt**2

piezo_width = int(ny * (piezo_diameter / Ly))
piezo_bottom = ny // 2 - piezo_width // 2

nframes = int(tmax / dt)

all_p = np.zeros((nframes, nx, ny), dtype=np.float32)
all_signal = np.zeros(nframes, dtype=np.float32)
all_time = np.zeros(nframes, dtype=np.float32)
# Initialisation des tableaux
p = np.zeros((nx, ny), dtype=np.float32)
p_prev = np.zeros((nx, ny), dtype=np.float32)
p_next = np.zeros((nx, ny), dtype=np.float32)
celerity_matrix = np.full((nx, ny), c_w**2, dtype=np.float32)  # Par défaut, on est uniquement dans l'eau...

#Si on ajoute l'aluminium, on modifie la matrice de célérité aux points où l'aluminium est présent.
if has_aluminium:
    alu_x0 = int(alu_distance / Lx * nx)
    alu_y0 = int((Ly - alu_width) / 2 / Ly * ny)
    alu_x1 = alu_x0 + int(alu_thickness / Lx * nx)
    alu_y1 = alu_y0 + int(alu_width / Ly * ny)
    celerity_matrix[alu_x0:alu_x1, alu_y0:alu_y1] = c_a**2

has_low_pass_filter = True # True si on veut ajouter un filtre passe-bas, False sinon
has_pml = True # True si on veut ajouter des conditions aux limites absorbantes de type PML, False sinon

sigma = np.zeros((nx, ny), dtype=np.float32)

# Paramètres de la PML (à ne pas toucher)
if has_pml:
    pml_width = nx // 5
    sigma_max = 0.2 * c_w / dx

    

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

assert alu_distance + alu_thickness < Lx, "La plaque d'aluminium doit être entièrement dans le domaine !"
assert piezo_diameter < Ly, "Le diamètre des dispositifs piézoélectriques doit être inférieur à la hauteur du domaine !"

@njit(fastmath=True, cache=True)
def compute_wave_propagation(nframes, nx, ny, p, p_prev, p_next, celerity_matrix, sigma, alpha, dt, dt2, inv_dx2, inv_dy2, piezo_bottom, piezo_width, f_source):
    all_p = np.zeros((nframes, nx, ny), dtype=np.float32)
    all_signal = np.zeros(nframes, dtype=np.float32)
    t = 0.0
    for frame in range(nframes):
        pulse_duration = 3 / f_source
        if t < pulse_duration:
            window = 0.5 * (1 - np.cos(2 * np.pi * t / pulse_duration))
            source_val = 700 * window * np.sin(2 * np.pi * f_source * t)
        else:
            source_val = 0.0

        p_next.fill(0)

        # Source
        p[0, piezo_bottom:piezo_bottom + piezo_width] = source_val

        # Différences finies et PML avec Laplacien d'ordre 4 à l'intérieur, ordre 2 près des bords
        for i in range(1, nx-1):
            for j in range(1, ny-1):
                if i < 2 or i > nx-3 or j < 2 or j > ny-3:
                    laplacian = (
                    (p[i+1, j] + p[i-1, j] - 2 * p[i, j]) * inv_dx2 +
                    (p[i, j+1] + p[i, j-1] - 2 * p[i, j]) * inv_dy2
                    )
                else:
                    laplacian = (
                    (-p[i+2, j] + 16*p[i+1, j] - 30*p[i, j] + 16*p[i-1, j] - p[i-2, j]) * (inv_dx2/12) +
                    (-p[i, j+2] + 16*p[i, j+1] - 30*p[i, j] + 16*p[i, j-1] - p[i, j-2]) * (inv_dy2)/12)

                if has_pml:
                    p_next[i, j] = (
                        (2 * p[i, j] - p_prev[i, j] +
                         (celerity_matrix[i, j] * dt2) * laplacian) *
                        np.exp(-sigma[i, j] * dt)
                    )
                else:
                    p_next[i, j] = (
                        (2 * p[i, j] - p_prev[i, j] +
                        (celerity_matrix[i, j] * dt2) * laplacian)
                    )
        if has_low_pass_filter:
            for i in range(1, nx - 1):
                for j in range(1, ny - 1):
                    laplacian = (
                        p_next[i+1, j] + p_next[i-1, j] +
                        p_next[i, j+1] + p_next[i, j-1] -
                        4 * p_next[i, j]
                    )
                    p_next[i, j] += alpha * laplacian

        received = np.mean(p_next[-2, piezo_bottom:piezo_bottom + piezo_width])
        all_p[frame, :, :] = p_next
        all_signal[frame] = received

        tmp = p_prev
        p_prev = p
        p = p_next
        p_next = tmp
        t += dt

    return all_p, all_signal

print("Init:", time.time() - start)

all_p, all_signal = compute_wave_propagation(
    nframes, nx, ny, p, p_prev, p_next, celerity_matrix, sigma, alpha, dt, dt2, inv_dx2, inv_dy2,
    piezo_bottom, piezo_width, f_source
)


print("Simulation complète:", time.time() - start)
fig = plt.figure(figsize=(8, 8))
gs = fig.add_gridspec(2, 1, height_ratios=[5, 1], width_ratios=[1], hspace=0.2, top=0.92, bottom=0.08, left=0.15, right=0.98)

ax1 = fig.add_subplot(gs[0])
im = ax1.imshow(all_p[0].T, cmap='magma', extent=[0, nx*dx, 0, ny*dy], vmin=-1, vmax=1, origin='lower', aspect='auto')
plt.colorbar(im, ax=ax1, label='Pression acoustique (Pa)')
ax1.set_xlabel('x (m)')
ax1.set_ylabel('y (m)')
ax1.set_title('Propagation d\'ondes acoustiques avec PML')

if has_aluminium:
    alu_patch = plt.Rectangle((alu_x0*dx, alu_y0*dy), alu_thickness, alu_width, fill=False, color='black', label='Plaque d\'Aluminium')
    ax1.add_patch(alu_patch)
    
emitter = plt.Line2D([0, 0], [piezo_bottom*dy, (piezo_bottom + piezo_width)*dy], color="#FFDD6D", linewidth=6, label='Émetteur')
receiver = plt.Line2D([nx*dx, nx*dx], [piezo_bottom*dy, (piezo_bottom + piezo_width)*dy], color="#3E177E", linewidth=6, label='Récepteur')
ax1.add_line(emitter)
ax1.add_line(receiver)
ax1.legend(loc='upper right')

ax2 = fig.add_subplot(gs[1])
ax2.set_xlim(0, tmax*1e6)
all_signal = all_signal / np.max(np.abs(all_signal))  # Normalisation du signal
ax2.set_ylim(1.2*np.min(all_signal), 1.2*np.max(all_signal))
ax2.set_xlabel('Temps (µs)')
ax2.set_ylabel('Signal reçu (normalisé)')
ax2.grid(True)

time_points = np.linspace(0, tmax, nframes) * 1e6
line, = ax2.plot([], [], color='purple', lw=2)

def init():
    line.set_data([], [])
    return im, line

def update(frame):
    data = all_p[frame].T
    vmax = np.max(np.abs(data))
    if vmax == 0:
        vmax = 1e-12
    im.set_data(data / vmax)
    im.set_clim(-1, 1) 
    line.set_data(time_points[:frame+1], all_signal[:frame+1])
    return im, line

frames_step = 10
frames_to_show = np.arange(0, nframes, frames_step)
anim = FuncAnimation(fig, update, frames=frames_to_show, init_func=init, interval=50, blit=False, repeat=False)

plt.show()
# writer = FFMpegWriter(fps=int(len(frames_to_show)/20), bitrate=1800)
# if has_aluminium:
#     anim.save('onde_acoustique_eau_alu_parallel.mp4', writer=writer, dpi=100)
# else:
#     anim.save('onde_acoustique_eau_parallel.mp4', writer=writer, dpi=100)

# print("Ecriture vidéo:", time.time() - start)
