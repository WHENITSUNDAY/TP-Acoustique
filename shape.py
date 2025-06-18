import numpy as np
import matplotlib.pyplot as plt
from matplotlib.path import Path
from matplotlib.patches import Polygon, RegularPolygon

def create_rectangle(material_matrix, nx, ny, Lx, Ly, distance, width, thickness, material=1):

    dx, dy = Lx/nx, Ly/ny

    x0 = int(distance / Lx * nx)
    y0 = int((Ly - width) / 2 / Ly * ny)
    x1 = x0 + int(thickness / Lx * nx)
    y1 = y0 + int(width / Ly * ny)
    material_matrix[x0:x1, y0:y1] = material

    patch = plt.Rectangle((x0*dx, y0*dy), thickness, width, fill=False, color='black', label='Plaque d\'Aluminium')

    return material_matrix, patch

def create_circle(material_matrix, nx, ny, Lx, Ly, center, radius, material=1):
    cx = center[0] * nx / Lx
    cy = center[1] * ny / Ly
    r = radius * nx / Lx
    
    for i in range(nx):
        for j in range(ny):
            if np.sqrt((i - cx)**2 + (j - cy)**2) <= r:
                material_matrix[i, j] = material

    patch = plt.Circle((center[0], center[1]), radius, fill=False, color='black', label='Boule d\'Aluminium')

    return material_matrix, patch

def create_hexagon(material_matrix, nx, ny, Lx, Ly, center, radius, material=1):
    cx = center[0] * nx / Lx
    cy = center[1] * ny / Ly
    r = radius * nx / Lx
    
    for i in range(nx):
        for j in range(ny):
            x = i - cx
            y = j - cy
            x_abs = abs(x)
            y_abs = abs(y)
            if (y_abs <= np.sqrt(3) * r / 2 and 
                x_abs <= r and 
                np.sqrt(3) * x_abs + y_abs <= np.sqrt(3) * r):
                material_matrix[i, j] = material

    patch = RegularPolygon((center[0], center[1]), numVertices=6, radius=radius, orientation=np.pi/6, fill=False, color='black', label='Hexagone d\'Aluminium')

    return material_matrix, patch

def create_triangle(material_matrix, nx, ny, Lx, Ly, center, base, height, material=1):

    v1 = (center[0] - base/2, center[1] - height/2)
    v2 = (center[0] + base/2, center[1] - height/2)
    v3 = (center[0], center[1] + height/2)
    
    pts = [
        (int(v1[0] / Lx * nx), int(v1[1] / Ly * ny)),
        (int(v2[0] / Lx * nx), int(v2[1] / Ly * ny)),
        (int(v3[0] / Lx * nx), int(v3[1] / Ly * ny)),
    ]
    
    X, Y = np.meshgrid(np.arange(nx), np.arange(ny), indexing='ij')
    points = np.vstack((X.flatten(), Y.flatten())).T
    path = Path(pts)
    mask = path.contains_points(points).reshape((nx, ny))
    material_matrix[mask] = material

    patch = Polygon([v1, v2, v3], closed=True, fill=False, color='black', label="Triangle d'Aluminium")
    return material_matrix, patch

def create_slit(material_matrix, nx, ny, Lx, Ly, distance, width, thickness, slit_width, plate_material=1, slit_material=0):
    patches = []
    dx, dy = Lx/nx, Ly/ny
    material_matrix, patch_1 = create_rectangle(material_matrix, nx, ny, Lx, Ly, distance, width, thickness, plate_material)

    slit_y0 = int((Ly / 2 - slit_width / 2) / Ly * ny)
    slit_y1 = int((Ly / 2 + slit_width / 2) / Ly * ny)
    x0 = int(distance / Lx * nx)
    y0 = int((Ly - width) / 2 / Ly * ny)
    x1 = x0 + int(thickness / Lx * nx)
    material_matrix[x0:x1, slit_y0:slit_y1] = slit_material

    alu_patch = plt.Rectangle((x0*dx, y0*dy), thickness, width, fill=False, color='black', label='Plaque fendue')
    slit_patch = plt.Rectangle((x0*dx, Ly/2 - slit_width/2), thickness, slit_width, fill=False, color='black', linestyle='dashed')

    patches = [alu_patch, slit_patch]

    return material_matrix, patches

def create_double_slit(material_matrix, nx, ny, Lx, Ly, distance, width, thickness, slit_width, slit_sep, plate_material=1, slit_material=0):
    
    dx, dy = Lx/nx, Ly/ny

    material_matrix, patch_1 = create_rectangle(material_matrix, nx, ny, Lx, Ly, distance, width, thickness, plate_material)
    center1 = Ly / 2 - slit_sep / 2
    center2 = Ly / 2 + slit_sep / 2
    
    x0 = int(distance / Lx * nx)
    x1 = x0 + int(thickness / Lx * nx)
    y0 = int((Ly - width) / 2 / Ly * ny)
    for center in [center1, center2]:
        slit_y0 = int((center - slit_width / 2) / Ly * ny)
        slit_y1 = int((center + slit_width / 2) / Ly * ny)
        material_matrix[x0:x1, slit_y0:slit_y1] = slit_material
    alu_patch = plt.Rectangle((x0*dx, y0*dy), thickness, width, fill=False, color='black', label='Plaque (fentes de Young))')
    slit1_patch = plt.Rectangle((x0*dx, center1 - slit_width/2), thickness, slit_width, fill=False, color='black', linestyle='dashed')
    slit2_patch = plt.Rectangle((x0*dx, center2 - slit_width/2), thickness, slit_width, fill=False, color='black', linestyle='dashed')
    
    patches = [alu_patch, slit1_patch, slit2_patch]

    return material_matrix, patches