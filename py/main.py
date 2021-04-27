import numpy as np
import time
import math
from progress.bar import ChargingBar    

large_motion = False
b_nim_error = False
N = 10
DIM = 3
SIZE = 3
L_FREE_MOTION = np.power(SIZE * SIZE * SIZE, 1/3.0) / (2*N)
dt = 0.001
SPEED = 1
eps = 1
sigma = 1
rcut = 2.5 * sigma
rmin = 0.0001
fulltime = dt

r = np.array([np.zeros([DIM]) for _ in np.arange(N)])
dr = np.array([np.zeros([DIM]) for _ in np.arange(N)])
v = np.array([np.zeros([DIM]) for _ in np.arange(N)])
f = np.array([np.zeros([DIM]) for _ in np.arange(N)])
m = np.array([np.zeros([DIM]) for _ in np.arange(N)])

def clear_file():
     with open('coords.xmol', 'w', encoding='utf-8') as file:
        pass

def new_frame(file_handle):
    file_handle.write(f'{N}\n')

def print_to_file():
    global v
    global r
    global fulltime
    with open('coords.xmol', 'a', encoding='utf-8') as file:
        new_frame(file)
        file.write(f'***** time = * {fulltime:.4} *****\n')
        for i in np.arange(N):
            file.write(f'Ar {r[i][0]:.4f} {r[i][1]:.4f} {r[i][2]:.4f} {v[i][0]:.4f} {v[i][1]:.4f} {v[i][2]:.4f}\n')

def init_system(zero_v=False):
    global m
    global list_coords
    if not zero_v:
        for i in np.arange(N):
            v[i] =  np.array([((SPEED - (-SPEED)) * np.random.random() + (-SPEED)) for _ in np.arange(DIM)])
    k = np.ceil(np.power(N, 1.0 / 3))
    dh = SIZE / k
    m = np.ones(N)
    counter = 0
    for x in np.arange(k):
        for y in np.arange(k):
            for z in np.arange(k):
                if counter < N:
                    r[counter] = np.array([(x + 1.0 / 2) * dh, (y + 1.0 / 2) * dh, (z + 1.0 / 2) * dh])
                    dr[counter] = np.array([ v[counter][0] * 2 * dt, v[counter][1] * 2 * dt, v[counter][2] * 2 * dt])
                    if lenght(dr[counter]) > L_FREE_MOTION:
                        print('init error')
                    counter += 1
    print_to_file()

def force_LD(r):
    if r > rcut:
        return 0
    if r < rmin:
        return force_LD(rmin)
    x = sigma / r
    return -48 * eps * (np.power(x, 13, dtype=np.float64) - 0.5 * np.power(x, 7, dtype=np.float64))

def verle_r(r, dr, f, m, dt): 
    return r + (dr + (f / (2 * m)) * np.square(dt))

def verle_v(dr, dt):
    return dr / (2 * dt)

def nim_fix(coord, size):
    if coord >= size / 2.0:
        coord = size - coord
    elif coord <= -size / 2.0:
        coord = coord + size
    return coord

def nim_vanilla(r1, r2, size):
    crds = [(-(r1[i]-r2[i])) for i in range(len(r1))]
    dist = [nim_fix(crd, size) for crd in crds]
    return dist

# nearist image method
def nim(r1, r2, size):
    global b_nim_error
    crds = np.array([(-(r1[i]-r2[i])) for i in np.arange(np.size(r1))])
    dist = np.array([nim_fix(crd, size) for crd in crds])
    if lenght(dist) > lenght(np.array([size, size, size])):
        if not b_nim_error:
            b_nim_error = True
            print('nim error')
    return dist

def lenght_vanilla(r):
    return math.sqrt(sum([r[i]**2 for i in range(len(r))]))

def lenght(r):
    return np.sqrt(np.sum([np.square(r[i]) for i in np.arange(np.size(r))]))

def normalaize(vec):
    return vec / lenght(vec)

def calc_forces():
    for i in np.arange(N):
        f[i] = np.zeros(DIM)
        for j in np.arange(N):
            if i != j:
                rij = nim(r[i], r[j], SIZE)
                _rij = lenght(rij)
                ff = force_LD(_rij)
                _dr = rij.copy()
                _dr = normalaize(_dr)
                f[i] += _dr * ff

def correct_coord(coord, left_boundary, right_boundary):
    l = right_boundary - left_boundary
    d = 0
    if coord >= right_boundary:
        d = coord - left_boundary
        coord = coord - l * np.floor(d / l)
    elif coord < left_boundary:
        d = left_boundary - coord
        coord = right_boundary - l * (d / l - np.floor(d / l))
    return coord

#periodic boundary condition
def pbc(r, size):
    return np.array([correct_coord(r[i], 0, size) for i in np.arange(r.size)]) 

def integrate():
    global large_motion
    for i in np.arange(N):
        r_tmp = r[i].copy() # хранит t-dt
        r[i] = verle_r(r[i], dr[i], f[i], m[i], dt) # вычсиляем координаты по алгоритму верле для t+dt
        dr[i] = r[i] - r_tmp # вычисляем разность координат между t+dt и t
        if lenght(dr[i]) > L_FREE_MOTION : # если слишком большая разноть, то значть что-то не так
            if not large_motion:
                large_motion = True
                print(' - too large motion detected')
        v[i] = verle_v(dr[i], dt) # вычисялем скорость по алгоритму верле
        r[i] = pbc(r[i], SIZE)

def main():
    global fulltime
    steps = 1000
    clear_file()
    st = time.time()
    init_system(zero_v=False)
    with ChargingBar('Steps', max=steps, suffix='%(percent)d%%') as bar:
        for i in np.arange(1, steps+1):
            calc_forces()
            integrate()
            print_to_file()
            fulltime += dt
            bar.next()
    ed = time.time()
    print(f'Done! Time: {(ed - st):.3f}')

if __name__ == '__main__':
    #main()

    a = [0.1, 0.1, 0.1]
    b = [9.9, 9.8, 9.7]
    dist = nim_vanilla(a, b, 10)
    lenght_vanilla(dist)

