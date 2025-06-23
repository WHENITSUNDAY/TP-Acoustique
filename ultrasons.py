import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import time
from numba import njit
from matplotlib.colors import ListedColormap
import tkinter as tk
from shape import *
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
show_energy = True #True si on veut montrer le profil énergetique, False sinon
has_aluminium = True # True si on veut ajouter la plaque d'alu, False sinon
# Paramètres physiques

aluminium_shape = 1 # 1 pour une plaque rectangulaire, 2 pour un cercle, 3 pour hexagone, 4 pour triangle équilatéral, 5 fente classique, 6 fentes de Young
Lx, Ly = 0.08, 0.08 # Dimensions du domaine (m)
piezo_diameter = 0.025 # Diamètre des dispositifs piézoélectriques (m)
f_source = 2e6 # Fréquence de la source (Hz), 2MHz -> Ultrasons.
tmax = 2e-4 # Durée de la simulation (s)
c_w = 1500 # Vitesse du son dans l'eau (m/s) (ou dans un fluide en général)
c_a = 6000 # Vitesse du son (longitudinal) dans l'aluminium (m/s) (ou dans un solide en général)

rho_w = 1000
rho_a = 2700
c_max = c_w if not has_aluminium else c_a

nx, ny = int(128*(Lx/Ly)), 128
dx = Lx / nx
dy = Ly / ny
dt = 0.8 * dx / (c_max * np.sqrt(2))  # CFL
alpha = 0.008 * dt / dx**2 

inv_dx2 = 1.0 / dx**2
inv_dy2 = 1.0 / dy**2
dt2 = dt**2

# Matrice de matériau, 1 si alu, 0 sinon
material_matrix = np.zeros((nx, ny), dtype=np.float32)

patch = None
patches = []

if has_aluminium:
    if aluminium_shape == 1:  # Rectangle
        alu_distance = 0.03  # Distance entre l'émetteur et la plaque d'aluminium (m)
        alu_thickness = 0.008  # Épaisseur de la plaque d'aluminium (m)
        alu_width = 0.067  # Largeur de la plaque d'aluminium (m)
        i = 0
        material_matrix, patch = create_rectangle(material_matrix, nx, ny, Lx, Ly, alu_distance, alu_width, alu_thickness, i)
    elif aluminium_shape == 2:  # Cercle
        alu_radius = 0.02  # Rayon du cercle
        alu_center = (Lx/2, Ly/2)  # Centre du cercle d'aluminium (m)
        material_matrix, patch = create_circle(material_matrix, nx, ny, Lx, Ly, alu_center, alu_radius)

    elif aluminium_shape == 3:  # Hexagone
        alu_radius = 0.02 
        alu_center = (Lx/2, Ly/2)
        material_matrix, patch = create_hexagon(material_matrix, nx, ny, Lx, Ly, alu_center, alu_radius)

    elif aluminium_shape == 4:  # Triangle équilatéral
        tri_height = 0.04
        tri_base = 0.04
        alu_center = (Lx/2, Ly/2)
        i = 0
        material_matrix, patch = create_triangle(material_matrix, nx, ny, Lx, Ly, alu_center, tri_base, tri_height, i)

    elif aluminium_shape == 5:  # Plaque avec fente centrale
        slit_width = 0.001
        alu_distance = 0.03
        alu_thickness = 0.004
        alu_width = 0.067
        material_matrix, patches = create_slit(material_matrix, nx, ny, Lx, Ly, alu_distance, alu_width, alu_thickness, slit_width)

    elif aluminium_shape == 6:  # Fentes de Young
        slit_width = 0.003
        slit_sep = 0.012
        alu_distance = 0.03
        alu_thickness = 0.004
        alu_width = 0.067
        material_matrix, patches = create_double_slit(material_matrix, nx, ny, Lx, Ly, alu_distance, alu_width, alu_thickness, slit_width, slit_sep)

celerity_matrix = np.where(material_matrix == 1, c_a, c_w).astype(np.float32)
rho_matrix = np.where(material_matrix == 1, rho_a, rho_w).astype(np.float32)        

piezo_width = int(ny * (piezo_diameter / Ly))
piezo_bottom = ny // 2 - piezo_width // 2

nframes = int(tmax / dt)

# Initialisation des tableaux
all_p = np.zeros((nframes, nx, ny), dtype=np.float32)
all_signal = np.zeros(nframes, dtype=np.float32)
all_time = np.zeros(nframes, dtype=np.float32)
p = np.zeros((nx, ny), dtype=np.float32)
p_prev = np.zeros((nx, ny), dtype=np.float32)
p_next = np.zeros((nx, ny), dtype=np.float32)

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
    j_top = ny - 1 - j

    sigma_bot = sigma_max * ((pml_width - j) / pml_width)**3
    sigma_top = sigma_max * ((pml_width - j) / pml_width)**3

    sigma[:, j] = np.maximum(sigma[:, j], sigma_bot[None, :])
    sigma[:, j_top] = np.maximum(sigma[:, j_top], sigma_top[None, :])

if has_aluminium:
    if aluminium_shape == 1:
        assert alu_distance + alu_thickness < Lx, "La plaque d'aluminium doit être entièrement dans le domaine !"
    elif aluminium_shape == 2 or aluminium_shape == 3 :
        assert alu_radius < Lx / 2, "Le rayon de la plaque d'aluminium doit être inférieur à la moitié de la largeur du domaine !"

assert piezo_diameter < Ly, "Le diamètre des dispositifs piézoélectriques doit être inférieur à la hauteur du domaine !"

cmap = ListedColormap(['#1f77b4', '#b0b0b0'])

root = tk.Tk()
screen_height = root.winfo_screenheight()
root.destroy()
dpi = plt.rcParams['figure.dpi']

fig_h = (screen_height * 0.8) / dpi


plt.figure(figsize=(1.2*fig_h, fig_h))
plt.imshow(material_matrix.T, cmap=cmap, extent=[0, Lx, 0, Ly], origin='lower', aspect='auto', vmin=0, vmax=1)
plt.title('Matrice de matériau (0: Eau, 1: Aluminium)')
plt.xlabel('x (m)')
plt.ylabel('y (m)')
cbar = plt.colorbar(ticks=[0, 1], label='Matériau')
cbar.ax.set_yticklabels(['Eau', 'Aluminium'])
plt.tight_layout()
plt.show()

@njit(fastmath=True, cache=True)
def compute_wave_propagation(nframes, nx, ny, p, p_prev, p_next, celerity_matrix, rho_matrix, sigma, alpha, dt, dt2, inv_dx2, inv_dy2, piezo_bottom, piezo_width, f_source):
    all_p = np.zeros((nframes, nx, ny), dtype=np.float32)
    all_signal = np.zeros(nframes, dtype=np.float32)
    all_energy = np.zeros((nframes, nx, ny), dtype=np.float32)
    energy_sum = np.zeros((nx, ny), dtype=np.float32)
    all_time = np.zeros(nframes, dtype=np.float32)
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

                c = celerity_matrix[i, j]
                rho = rho_matrix[i, j]

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
                         (c**2 * dt2) * laplacian) *
                        np.exp(-sigma[i, j] * dt)
                    )
                else:
                    p_next[i, j] = (
                        (2 * p[i, j] - p_prev[i, j] +
                        (c**2 * dt2) * laplacian)
                    )

                #calcul de l'energie cumulative
                dpdt = (p[i, j] - p_prev[i, j]) / dt
                v = dpdt / (rho * c)
                e_instant = 0.5 * rho * v**2 + 0.5 * p[i, j]**2 / (rho * c**2)
                energy_sum[i, j] += e_instant
                all_energy[frame, i, j] = np.log2(energy_sum[i, j] / (frame + 1) + 1) 

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
        all_time[frame] = t
        tmp = p_prev
        p_prev = p
        p = p_next
        p_next = tmp
        t += dt

    return all_p, all_signal, all_energy, all_time

print("Init:", time.time() - start)

all_p, all_signal, all_energy, all_time = compute_wave_propagation(
    nframes, nx, ny, p, p_prev, p_next, celerity_matrix, rho_matrix, sigma, alpha, dt, dt2, inv_dx2, inv_dy2,
    piezo_bottom, piezo_width, f_source
)

if np.max(np.abs(all_signal)) != 0:
    all_signal = all_signal / np.max(np.abs(all_signal))  # Normalisation du signal
print("Simulation complète:", time.time() - start)


#Figure principale : pression acoustique et signal reçu
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=((0.9*fig_h, fig_h)), gridspec_kw={'height_ratios': [5, 1]})
plt.subplots_adjust(hspace=0.35, top=0.95, bottom=0.08, left=0.08, right=0.98)

im = ax1.imshow(all_p[0].T, cmap='magma', extent=[0, nx*dx, 0, ny*dy], vmin=-1, vmax=1, origin='lower', aspect='auto')
plt.colorbar(im, ax=ax1, label='Pression acoustique (Pa)')
ax1.set_xlabel('x (m)')
ax1.set_ylabel('y (m)')
ax1.set_title('Propagation d\'ondes acoustiques avec PML')

if has_aluminium:
    if patch:
        ax1.add_patch(patch)
    if patches:
        for patch in patches:
            ax1.add_patch(patch)


emitter = plt.Line2D([0, 0], [piezo_bottom*dy, (piezo_bottom + piezo_width)*dy], color="#FFDD6D", linewidth=6, label='Émetteur')
receiver = plt.Line2D([nx*dx, nx*dx], [piezo_bottom*dy, (piezo_bottom + piezo_width)*dy], color="#3E177E", linewidth=6, label='Récepteur')
ax1.add_line(emitter)
ax1.add_line(receiver)
ax1.legend(loc='upper right')

# Plot signal reçu
ax2.set_xlim(0, tmax*1e6)
ax2.set_ylim(1.2*np.min(all_signal), 1.2*np.max(all_signal))
ax2.set_xlabel('Temps (µs)')
ax2.set_ylabel('Signal reçu (normalisé)')
ax2.grid(True)
time_points = np.linspace(0, tmax, nframes) * 1e6
line, = ax2.plot([], [], color='purple', lw=2)

paused = False
def on_press(event):
    global paused
    if event.key == ' ':
        paused = not paused
        if paused:
            anim.event_source.stop()
        else:
            anim.event_source.start()

def update(frame):
    p_data = all_p[frame].T
    vmax_p = np.max(np.abs(p_data))
    if vmax_p == 0:
        vmax_p = 1e-12
    im.set_data(p_data / vmax_p)
    im.set_clim(-1, 1)
    line.set_data(time_points[:frame+1], all_signal[:frame+1])
    ax1.set_title(f'Propagation d\'ondes acoustiques avec PML\nTemps : {all_time[frame]*1e6:.2f} µs')
    return im, line

interval = 30
step = max(1, nframes // (20 * 1000 // interval))

fig.canvas.mpl_connect('key_press_event', on_press)
anim = FuncAnimation(fig, update, frames=range(0, nframes, step), interval=interval, blit=False, repeat=True)

plt.show()
plt.close(fig)

if show_energy:
    fig_energy, ax_energy = plt.subplots(figsize=(1.2*fig_h, fig_h))
    im_energy = ax_energy.imshow(all_energy[0].T, cmap='inferno', extent=[0, nx*dx, 0, ny*dy], origin='lower', aspect='auto', vmax=1)
    plt.colorbar(im_energy, ax=ax_energy, label='Énergie (J)')
    ax_energy.set_xlabel('x (m)')
    ax_energy.set_ylabel('y (m)')
    ax_energy.set_title('Distribution de l\'énergie')

    paused = False

    def on_press_energy(event):
        global paused
        if event.key == ' ':
            paused = not paused
            if paused:
                anim_energy.event_source.stop()
            else:
                anim_energy.event_source.start()

    def update_energy(frame):
        energy_data = all_energy[frame].T
        vmax_e = np.max(np.abs(energy_data))
        if vmax_e == 0:
            vmax_e = 1e-12
        im_energy.set_data(energy_data / vmax_e)
        im_energy.set_clim(0, 1)
        ax_energy.set_title(f'Distribution de l\'énergie cumulée\nTemps : {all_time[frame]*1e6:.2f} µs')
        return im_energy,

    fig_energy.canvas.mpl_connect('key_press_event', on_press_energy)
    anim_energy = FuncAnimation(fig_energy, update_energy, frames=range(0, nframes, step), interval=interval, blit=False, repeat=True)
    plt.show()
    plt.close(fig_energy)