import pickle
import numpy as np
from scipy.spatial import Voronoi
from shapely.geometry import Polygon , Point
import matplotlib.pyplot as plt
import math
from math import gcd
from functools import reduce
import pygame

#############################################################
##START STUFF
#############################################################


BLACK = (0, 0, 0)
WHITE = (255, 255, 255)
RED = (255, 0, 0)
GREEN = (0, 255, 0)
BLUE = (0, 0, 255)


gmem = np.zeros(2*6_850_000, dtype=np.float32)

OCTAVE = 900
LIMIT = 50000
s_percent = 0.01
s = np.log2(1+s_percent)*OCTAVE*np.sqrt(3/2)
how_much_octaves = 1.5


#gmem
seeds = 100000
data = seeds + (1500*1500)
data2 = data + (1500*1500)
outline_data = data + (1500*1500)


##right size of kernel based on sd
half_g_width = round(s*4)
hgw_ratio_down = 1/(2**(half_g_width/OCTAVE))
hgw_ratio_up = 2 / hgw_ratio_down
hgw_ratio_up *= how_much_octaves

##sorting
_L = 100000
_R= _L + ((1325*1325)/2)
local_data = np.zeros(1325 * 1325, dtype=float)
color_idx = np.zeros(16, dtype=float)

##rest
counter = 0
done = 0
big_sum = 0
sqrt3_2 = np.sqrt(3)/2
sum = 0
count_pass2 = 0

#########################################################
##GAUSS STUFF
#########################################################

gauss_kernel = {}
def gauss(distance):
    return math.exp((distance / s) ** 2 * -0.5) / (s * math.sqrt(2 * math.pi))

def make_gauss_kernel():
    d = -half_g_width
    for _ in range(half_g_width * 2 + 1):
        gauss_kernel[d] = gauss(d)
        d += 1


def apply_gauss(w, h, src, dest):
    i = -1
    for _ in range(h):
        left_edge = i
        right_edge = i + w
        for _ in range(w):
            i += 1
            weight = gmem[src + i]
            if weight > 0:
                start_i = max(left_edge, i - half_g_width)
                stop_i = min(right_edge, i + half_g_width)
                g = start_i
                for _ in range(stop_i - start_i + 1):
                    gmem[dest + g] += weight * gauss_kernel[g - i]
                    g += 1

################################################################
##SEED STUFF
################################################################

def gcd3(A, B, C):
    return reduce(gcd, [A, B, C])


def do_stuff():
    coordinates = set()  
    C_total_loops=0
    idx=0
    A=1
    A_loop=0
    global counter
    counter = 0
    while(A_loop<270):
        A_loop+=1
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
                    gmem[seeds + x + (y * 1325)] = max( 1 / ((A * B_ * C) ** (1 / 3)), gmem[seeds + x + (y * 1325)])
                    counter+=1  
                C+=1
                C_total_loops+1
            B_+=1
        A+=1


############################################################
##SOME STUFF!!!!!
############################################################



def finish(ary):
    z = -1
    big_sum = 0
    big_count = 0
    sdev_sum = 0

    # Loop through 1325x1325
    for _ in range(1325):
        for _ in range(1325):
            z += 1
            big_sum += gmem[ary + z]
            if gmem[ary + z] > 0:
                big_count += 1

    # Calculate areas and mean
    mean1 = big_sum / big_count
    circle_area = (half_g_width ** 2) * math.pi
    full_area = ((half_g_width * 2 + OCTAVE) ** 2) * (sqrt3_2 / 2)
    avg_local_sum = (circle_area / full_area) * counter * mean1

    sum_pass2 = 0
    count_pass2 = 0
    z = -1
    sort_ary_idx = -1

    # Normalize and calculate entropy
    for _ in range(1325):
        for _ in range(1325):
            z += 1
            point = gmem[ary + z]
            if point:
                point /= avg_local_sum
                gmem[ary + z] = -point * math.log(point)
                sum_pass2 += gmem[ary + z]
                count_pass2 += 1
                local_data[count_pass2] = gmem[ary + z]
            else:
                gmem[ary + z] = 0

    mean2 = sum_pass2 / count_pass2

    # Calculate standard deviation
    z = -1
    for _ in range(1325):
        for _ in range(1325):
            z += 1
            sdev_sum += (gmem[ary + z] - mean2) ** 2

    plot_sdev = math.sqrt(sdev_sum / count_pass2)



def transpose_matrix(src, dest, rows, cols):
    c_i = 1
    d_i = -1
    for _ in range(cols):
        c_i += 1
        r_i = -1
        for _ in range(rows):
            r_i += 1
            d_i += 1
            gmem[dest + d_i] = gmem[src + c_i + (cols * r_i)]


def wipe_gmem(wipe_start):
    gm_idx = 0
    for _ in range(1325):
        for _ in range(1325):
            gmem[wipe_start + gm_idx] = 0
            gm_idx += 1

#################################################
### SORTINGG 
#################################################

def partition(arr, low, high):
    i = low - 1
    for j in range(low, high):
        if arr[j] <= arr[high]:
            i += 1
            arr[i], arr[j] = arr[j], arr[i]
    arr[i + 1], arr[high] = arr[high], arr[i + 1]
    return i + 1

def quicksort(arr, low, high):
    global color_idx  # Declare color_idx as global
    stack = []
    stack.append(low)
    stack.append(high)
    while stack:
        high = stack.pop()
        low = stack.pop()
        pivot_index = partition(arr, low, high)
        if pivot_index - 1 > low:
            stack.append(low)
            stack.append(pivot_index - 1)
        if pivot_index + 1 < high:
            stack.append(pivot_index + 1)
            stack.append(high)

    chunk = len(local_data) // 16
    color_idx[:] = [local_data[i * chunk] for i in range(16)]


########################################
## TOGETHER?!?!??!?
#########################################


def recompute():
    make_gauss_kernel()

    do_stuff()
    apply_gauss(1325,1325,seeds,data)
    transpose_matrix(data,data2,1325,1325)
    wipe_gmem(data)
    apply_gauss(1325,1325,data2,data)
    wipe_gmem(data2)
    transpose_matrix(data,data2,1325,1325)
    print("hello")
    finish(data2)
    print("got past this bit")
    last_s = s_percent

    quicksort(local_data, 0, count_pass2)
    print("past here too")
    done= 0
    last_s = s_percent


# def defer_recompute():
#     defer_recompute -= 1
#     if defer_recompute == 0:
#         recompute()

#####################################################
####GUI STUFF
#####################################################
recompute()
width, height = 1500, 1500
pygame.init()
screen = pygame.display.set_mode((width, height))


def cubehelix3_16(in_val):
    if in_val < color_idx[1]:
        r, gg, b = 0, 0, 0
    elif in_val < color_idx[2]:
        r, gg, b = 0, 39, 12
    elif in_val < color_idx[3]:
        r, gg, b = 0, 68, 60
    elif in_val < color_idx[4]:
        r, gg, b = 0, 80, 131
    elif in_val < color_idx[5]:
        r, gg, b = 3, 75, 202
    elif in_val < color_idx[6]:
        r, gg, b = 72, 60, 252
    elif in_val < color_idx[7]:
        r, gg, b = 156, 43, 255
    elif in_val < color_idx[8]:
        r, gg, b = 235, 36, 244
    elif in_val < color_idx[9]:
        r, gg, b = 255, 45, 194
    elif in_val < color_idx[10]:
        r, gg, b = 255, 73, 134
    elif in_val < color_idx[11]:
        r, gg, b = 255, 115, 86
    elif in_val < color_idx[12]:
        r, gg, b = 255, 164, 67
    elif in_val < color_idx[13]:
        r, gg, b = 235, 209, 85
    elif in_val < color_idx[14]:
        r, gg, b = 211, 241, 135
    elif in_val < color_idx[15]:
        r, gg, b = 215, 255, 200
    else:
        r, gg, b = 255, 255, 255

    return r , gg , b 




def draw_from_mem():
    gfx_y = -1
    gfx_r = gfx_g = gfx_b = 0


    for _ in range(1325):
        gfx_y += 1
        gfx_x = -1
        for _ in range(1325):
            gfx_x += 1
            v = gmem[data2 + gfx_x + (gfx_y * 1325)]
            v = max(0, v)
            rr , ggg , bb = cubehelix3_16(v)
            screen.set_at((gfx_x, gfx_y), (rr, ggg, bb))

    ##outline()
    global done
    done = 1
running = True
while running:
    for event in pygame.event.get():
        if event.type == pygame.QUIT: 
            running = False
    if not running :
        break
    
    screen.fill(WHITE)
    draw_from_mem()

    pygame.image.save(screen, "output_image.png")
    pygame.display.flip()


pygame.quit()
