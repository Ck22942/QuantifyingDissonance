import numpy as np
from scipy.spatial import Voronoi
from shapely.geometry import Polygon, Point
import matplotlib.pyplot as plt

def make_farey_coordinates(n):
    coordinates = set()  # Use a set to avoid duplicates
    for i in range(n, -1, -1):
        for j in range(n, i, -1):
            for k in range(n, j, -1):
                if j != 0 and k != 0:  # Avoid division by zero
                    r1 = i / float(j)
                    r2 = j / float(k)
                    if r2 * r1 >= 0.5:
                        coordinates.add((r1, r2))  # Add to the set
    return np.array(list(coordinates))  # Convert to a NumPy array for compatibility

# Function to calculate Voronoi cell areas
def calculate_voronoi_areas(vor):
    areas = []
    polygons = []
    for region_idx in vor.point_region:
        region = vor.regions[region_idx]
        if -1 in region or len(region) == 0:  # Ignore infinite regions
            areas.append(None)
            polygons.append(None)
        else:
            polygon = Polygon([vor.vertices[i] for i in region])
            areas.append(polygon.area)
            polygons.append(polygon)
    return areas, polygons

# Generate Farey coordinates and compute Voronoi diagram
points = make_farey_coordinates(37)  # Generate Farey coordinates for n=10
vor = Voronoi(points)

def gaussian_2d(x, y, x_i, y_i, sigma):
    coefficient = 1 / (2 * np.pi * sigma**2)
    exponent = -((x - x_i)**2 + (y - y_i)**2) / (2 * sigma**2)
    return coefficient * np.exp(exponent)

def integrate_gaussian_over_voronoi_cell(polygon, x_i, y_i, sigma, num_samples=1000):
  
    minx, miny, maxx, maxy = polygon.bounds
    
    total_value = 0
    for _ in range(num_samples):

        x, y = np.random.uniform(minx, maxx), np.random.uniform(miny, maxy)
        point = Point(x, y)
        if polygon.contains(point):
            total_value += gaussian_2d(x, y, x_i, y_i, sigma)
    
    polygon_area = polygon.area
    bounding_box_area = (maxx - minx) * (maxy - miny)
    return total_value * (polygon_area / bounding_box_area)



sigma = 0.2  
probabilities = []

for region_idx, point in zip(vor.point_region, points):
    region = vor.regions[region_idx]
    if -1 in region or len(region) == 0:  
        probabilities.append(None)
    else:
        polygon = Polygon([vor.vertices[i] for i in region])
        x_i, y_i = point 
        prob = integrate_gaussian_over_voronoi_cell(polygon, x_i, y_i, sigma)
        probabilities.append(prob)

for i, prob in enumerate(probabilities):
    print(f"Voronoi cell {i}: Probability = {prob if prob is not None else 'Infinite'}")


fig, ax = plt.subplots(figsize=(8, 8))
for region_idx, prob in zip(vor.point_region, probabilities):
    region = vor.regions[region_idx]
    if -1 in region or len(region) == 0:
        continue  # Skip infinite regions
    polygon = Polygon([vor.vertices[i] for i in region])
    if prob is not None:
        color = plt.cm.viridis(prob / max(filter(None, probabilities)))  # Normalize for color mapping
        x, y = polygon.exterior.xy
        ax.fill(x, y, color=color, edgecolor='black', alpha=0.7)

plt.scatter(points[:, 0], points[:, 1], c='black', s=10)
plt.title("Voronoi Diagram with Gaussian Cell Probabilities")
plt.axis('off')
plt.show()