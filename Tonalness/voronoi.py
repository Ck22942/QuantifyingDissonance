import numpy as np
from scipy.spatial import Voronoi
from shapely.geometry import Polygon
import matplotlib.pyplot as plt
import math
from math import gcd
from functools import reduce

OCTAVE = 600

LIMIT = 50000
def gcd3(A, B, C):
    return reduce(gcd, [A, B, C])

def make_farey_coordinates(n):
    coordinates = set()  
    for i in range(n, -1, -1):
        for j in range(n, i, -1):
            for k in range(n, j, -1):
                if j != 0 and k != 0:  
                    r1 = i / float(j)
                    r2 = j / float(k)
                    if (r1 <= 0.95) and (r2 <= 0.95) and (1 >= (r2 * r1) >= 0.5) and (gcd3(i,j,k) == 1) and (i * j * k < LIMIT):
                        coordinates.add((r1, r2))  
    return np.array(list(coordinates))  



def calculate_voronoi_areas(vor):
    areas = []
    polygons = []
    for region_idx in vor.point_region:
        region = vor.regions[region_idx]
        if -1 in region or len(region) == 0:  
            areas.append(None)
            polygons.append(None)
        else:
            polygon = Polygon([vor.vertices[i] for i in region])
            areas.append(polygon.area)
            polygons.append(polygon)
    return areas, polygons

points = make_farey_coordinates(37)  
vor = Voronoi(points)
print(len(points))




areas, polygons = calculate_voronoi_areas(vor)


finite_areas = [area for area in areas if area is not None]
if finite_areas:
    min_area, max_area = min(finite_areas), max(finite_areas)
else:
    min_area, max_area = 0, 1  



fig, ax = plt.subplots(figsize=(8, 8))
for polygon, area in zip(polygons, areas):
    if polygon is not None:  
        normalized_area = (area - min_area) / (max_area - min_area) 
        color = plt.cm.viridis(normalized_area*20)  
        x, y = polygon.exterior.xy
        ax.fill(x, y, color=color, edgecolor='black', alpha=0.6)



ax.set_xlim(points[:, 0].min()/1.25, points[:, 0].max()*1.25)
ax.set_ylim(points[:, 1].min()/1.25, points[:, 1].max()*1.25)

points2 = np.array([[8/9, 9/10], [7/8, 8/13], [3/4, 4/6] , [2/3,3/4]])

# Plot the points
ax.scatter(points2[:, 0], points2[:, 1], color='blue', label='Points')


def custom_xscale(x, y):
    x1 = -np.log2(x)
    y1 = -np.log2(y) / 2
    return x1 + y1

def custom_yscale(y):
    y1 = -np.log2(y)
    return y1 * -np.sqrt(3) / 2

# Inverse transformations (needed for FuncScale)
def inverse_custom_xscale(x_scaled, y):
    # Reconstruct x from x_scaled and y using the reverse transformation logic
    return 2 ** (-x_scaled - (-np.log2(y) / 2))

def inverse_custom_yscale(y_scaled):
    return 2 ** (-y_scaled / (-np.sqrt(3) / 2))


# Add grid with more minor ticks
ax.grid(True, which="both", linestyle='--', linewidth=0.5)

# Add title
plt.title("Plot with Both Axes as Negative Logarithmic Scales")
plt.show()

