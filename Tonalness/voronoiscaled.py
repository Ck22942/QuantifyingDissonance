import pickle
import numpy as np
from scipy.spatial import Voronoi
from shapely.geometry import Polygon , Point
import matplotlib.pyplot as plt
import math
from math import gcd
from functools import reduce

#######
gmem = dict()
OCTAVE = 900
LIMIT = 50000
seeds = 100000
s_percent = 0.01
s = np.log2(1+s_percent)*OCTAVE*np.sqrt(3/2)


half_g_width = round(s*4)
hgw_ratio_down = 1/(2**(half_g_width/OCTAVE))
hgw_ratio_up = 2 / hgw_ratio_down


##sorting
sqrt3_2 = np.sqrt(3)/2
_L = 100000
_R= _L + ((1325*1325)/2)
local_data = _R+ ((1325*1325)/2)
color_idx = _L -  16




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

    return (total_value / num_samples) * (polygon_area / bounding_box_area)



def gcd3(A, B, C):
    return reduce(gcd, [A, B, C])

def make_farey_coordinates(n):
    coordinates = set()  
    C_total_loops=0
    idx=0
    A=1
    A_loop=0
    counter=0
    while(A_loop<270):
        A_loop+=1
        for i in range(n, -1, -1):
            B_ = max(math.floor(A*hgw_ratio_down),1)
            B_loop = 0
            while(B_loop<270 and (B_/A)<hgw_ratio_up):
                B_loop+=1
                C=max(math.floor(A*hgw_ratio_down),1)
                C_loop = 0
                while(C_loop<270 and (C/A)<hgw_ratio_up):
                    C_loop+=1
                    GCD = gcd3(A,B_,C)
                    if (C / A > hgw_ratio_down) and (C / B_ > hgw_ratio_down) and (B_ / A > hgw_ratio_down) and (C / B_ < hgw_ratio_up)  and (A * B_ * C < LIMIT):
                        x1 = np.log2(B_/A)*OCTAVE
                        y1 = np.log2(C/B_)*OCTAVE
                        x=round(x1+(y1/2))
                        y=round(y1*sqrt3_2)
                        x+= half_g_width
                        y+= half_g_width
                        gmem[seeds+x+(y*1325)]=max(1/((A*B_*C)^(1/3)),gmem[seeds+x+(y*1325)])
                        counter+=1
                        


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

def transform_x(x, y):
    x1 = -np.log2(x)
    y1 = -np.log2(y) / 2
    return x1 + y1

def transform_y(y):
    y1 = -np.log2(y)
    return y1 * np.sqrt(3) / 2


points = make_farey_coordinates(50)
vor = Voronoi(points)
areas , polygons = calculate_voronoi_areas(vor)

finite_areas = [area for area in areas if area is not None]
min_area, max_area = min(finite_areas), max(finite_areas)

transformed_vertices = np.array([
    [transform_x(v[0], v[1]), transform_y(v[1])] if v[0] > 0 and v[1] > 0 else [np.nan, np.nan]
    for v in vor.vertices
])



fig, ax = plt.subplots(figsize=(8, 8))
idx = 0
for region_idx, area in zip(vor.point_region, areas):
    region = vor.regions[region_idx]
    if -1 in region or len(region) == 0:  
        continue
    polygon = [transformed_vertices[i] for i in region if i != -1]
    polygon = np.array(polygon)

    if (area is not None):
        
        color = plt.cm.rainbow_r((log_normalized[idx] - 1) / 1.75 )  
        
       
        ax.fill(polygon[:, 0], polygon[:, 1], color=color, edgecolor='None', alpha=0.7)
    idx += 1





ax.set_xlim(0, 1)  
ax.set_ylim(0, 1)




plt.title("Transformed Voronoi Diagram Colored by Original Area")
plt.axis('off')
plt.show()

sorted_list = sorted(log_normalized)


plt.plot(sorted_list, marker='o', linestyle='-', color='b', label='Sorted Data')
plt.xlabel('Index')
plt.ylabel('Value')
plt.title('Plot of Sorted List')
plt.legend()
plt.grid()
plt.show()





# n = len(polygons)
# print("-------" + str(n) + " --------- ")
# entropies = []

# for i in range(n):
#     HE_j = 0
#     x_i , y_i = points[i]
#     for j  in range(n):
#         if areas[j] is not None:  
#             x_0 , y_0 = points[j]
#             if (x_0 - x_i)**2 + (y_0 - y_i)**2 < 9 *(sigma**2): 
#                 p = integrate_gaussian_over_voronoi_cell(polygons[j] , x_i , y_i , sigma, 1000)
#                 if p != 0:
#                     HE_j -= p * np.log(p)

#     entropies.append(HE_j)
#     print(i)  

# saveEntropies = [float(x) for x in entropies]

# print(saveEntropies)

# import pickle
# with open("newdate.pkl", "wb") as f:
#     pickle.dump(saveEntropies, f)

# print("   ")