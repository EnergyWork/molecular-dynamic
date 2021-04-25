import numpy as np
from progress.bar import ChargingBar

large_motion = False
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
time = dt

r = np.array([np.zeros([DIM]) for _ in np.arange(N)])
dr = np.array([np.zeros([DIM]) for _ in np.arange(N)])
v = np.array([np.zeros([DIM]) for _ in np.arange(N)])
f = np.array([np.zeros([DIM]) for _ in np.arange(N)])
m = np.array([np.zeros([DIM]) for _ in np.arange(N)])

list_coords = np.array([])

def clear_file():
     with open('coords.xmol', 'w', encoding='utf-8') as file:
        pass

def new_frame(file_handle):
    file_handle.write(f'{N}\n')

def print_to_file(ti):
    global v
    global r
    global time
    with open('coords.xmol', 'a', encoding='utf-8') as file:
        new_frame(file)
        file.write(f'***** time = * {time:.4} *****\n')
        for i in np.arange(N):
            file.write(f'Ar {r[i][0]:.4f} {r[i][1]:.4f} {r[i][2]:.4f} {v[i][0]:.4f} {v[i][1]:.4f} {v[i][2]:.4f}\n')

def init_system():
    global m
    global list_coords
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
    list_coords = np.append(list_coords, r)
    print_to_file(0)

def force_LD(r):
    if r > rcut:
        return 0
    if r < rmin:
        return force_LD(rmin)
    x = sigma / r
    #return -48 *  eps / sigma * (np.power(x, 13, dtype=np.float64) - 0.5 * np.power(x, 7, dtype=np.float64))
    return -4 * eps / sigma * (12 * np.power(x, 13, dtype=np.float64) - 6 * np.power(x, 7, dtype=np.float64))

def potential_LD(r):
    if r > rcut:
        return 0
    if r < rmin:
        return potential_LD(rmin)
    x = sigma / r
    return 4 * eps * (np.power(x, 12, dtype = np.float64) - np.power(x, 6, dtype = np.float64))

def verle_r(r, dr, f, m, dt): 
    return r + (dr + (f / (2 * m)) * np.square(dt))
    #return 2*r-dr+(f/m)*np.square(dt)

def verle_v(r, dr, dt):
    return dr / (2 * dt)
    #return (r-dr)/(2*dt)

def nim_fix(coord):
    if coord >= SIZE / 2.0:
        coord = SIZE - coord
    elif coord <= -SIZE / 2.0:
        coord = coord + SIZE
    return coord

# nearist image method
def nim(r1, r2, size):
    x = -(r1[0] - r2[0])
    y = -(r1[1] - r2[1])
    z = -(r1[2] - r2[2])

    dist = np.zeros(DIM)
    dist[0] = nim_fix(x)
    dist[1] = nim_fix(y)
    dist[2] = nim_fix(z)

    return dist

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
    l = np.abs(right_boundary - left_boundary)
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
        r_tmp = r[i].copy()
        r[i] = verle_r(r[i], dr[i], f[i], m[i], dt)
        dr[i] = r[i] - r_tmp
        if lenght(dr[i]) > L_FREE_MOTION :
            if not large_motion:
                large_motion = True
                print(' - too large motion detected')
        v[i] = verle_v(r[i], dr[i], dt)
        r[i] = pbc(r[i], SIZE)

def main():
    global list_coords
    global time
    steps = 1000
    clear_file()
    init_system()
    with ChargingBar('Steps', max=steps, suffix='%(percent)d%%') as bar:
        for i in np.arange(1, steps+1):
            calc_forces()
            integrate()
            list_coords = np.append(list_coords, r)
            print_to_file(i)
            time += dt
            bar.next()
    print('Done!')

if __name__ == '__main__':
    main()
