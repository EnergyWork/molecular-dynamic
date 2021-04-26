import numpy as np
import time
from progress.bar import ChargingBar

class MDSystem:
    def __init__(self, n_atoms, cube_size, dim, speed):
        self.N = n_atoms
        self.SIZE = cube_size
        self.DIM = dim
        self.SPEED = speed
        self.__init_vars()
    
    def __init_vars(self):
        self.large_motion = False
        self.b_nim_error = False
        self.L_FREE_MOTION = np.power(self.SIZE * self.SIZE * self.SIZE, 1/3.0) / (2*self.N)
        self.dt = 0.001
        self.eps = 1
        self.sigma = 1
        self.rcut = 2.5 * self.sigma
        self.rmin = 0.0001
        self.fulltime = self.dt

    def __init_arrays(self):
        self.r = np.array([np.zeros([self.DIM]) for _ in np.arange(self.N)])
        self.dr = np.array([np.zeros([self.DIM]) for _ in np.arange(self.N)])
        self.v = np.array([np.zeros([self.DIM]) for _ in np.arange(self.N)])
        self.f = np.array([np.zeros([self.DIM]) for _ in np.arange(self.N)])
        self.m = np.array([np.zeros([self.DIM]) for _ in np.arange(self.N)])

    def __clear_file(self):
        with open('coords.xmol', 'w', encoding='utf-8') as file:
            pass

    def __new_frame(self, file_handle):
        file_handle.write(f'{self.N}\n')

    def print_to_file(self):    
        with open('coords.xmol', 'a', encoding='utf-8') as file:
            self.__new_frame(file)
            file.write(f'***** time = * {self.fulltime:.4} *****\n')
            for i in np.arange(self.N):
                file.write(f'Ar {self.r[i][0]:.4f} {self.r[i][1]:.4f} {self.r[i][2]:.4f} {self.v[i][0]:.4f} {self.v[i][1]:.4f} {self.v[i][2]:.4f}\n')

    def init_system(self, zero_v=False):
        self.__clear_file()
        self.__init_arrays()
        if not zero_v:
            for i in np.arange(self.N):
                self.v[i] = np.array([((self.SPEED - (-self.SPEED)) * np.random.random() + (-self.SPEED)) for _ in np.arange(self.DIM)])
        k = np.ceil(np.power(self.N, 1.0 / 3))
        dh = self.SIZE / k
        self.m = np.ones(self.N)
        counter = 0
        for x in np.arange(k):
            for y in np.arange(k):
                for z in np.arange(k):
                    if counter < self.N:
                        self.r[counter] = np.array([(x + 1.0 / 2) * dh, (y + 1.0 / 2) * dh, (z + 1.0 / 2) * dh])
                        self.dr[counter] = np.array([self.v[counter][0] * 2 * self.dt, self.v[counter][1] * 2 * self.dt, self.v[counter][2] * 2 * self.dt])
                        if self.__lenght(self.dr[counter]) > self.L_FREE_MOTION:
                            print(' - init error')
                        counter += 1
        self.print_to_file()

    def __force_LD(self, r):
        if r > self.rcut:
            return 0
        if r < self.rmin:
            return self.__force_LD(self.rmin)
        x = self.sigma / r
        return -48 * self.eps * (np.power(x, 13, dtype=np.float64) - 0.5 * np.power(x, 7, dtype=np.float64))

    def __verle_r(self, r, dr, f, m, dt): 
        return r + (dr + (f / (2 * m)) * np.square(dt))

    def __verle_v(self, dr, dt):
        return dr / (2 * dt)

    def __nim_fix(self, coord, size):
        if coord >= size / 2.0:
            coord = size - coord
        elif coord <= -size / 2.0:
            coord = coord + size
        return coord

    # nearist image method
    def __nim(self, r1, r2, size):
        crds = np.array([(-(r1[i]-r2[i])) for i in np.arange(np.size(r1))])
        dist = np.array([self.__nim_fix(crd, size) for crd in crds])
        if self.__lenght(dist) > self.__lenght(np.array([size, size, size])):
            if not self.b_nim_error:
                self.b_nim_error = True
                print(' - nim error')
        return dist

    def __lenght(self, r):
        return np.sqrt(np.sum([np.square(r[i]) for i in np.arange(np.size(r))]))

    def __normalaize(self, vec):
        return vec / self.__lenght(vec)

    def calc_forces(self):
        for i in np.arange(self.N):
            self.f[i] = np.zeros(self.DIM)
            for j in np.arange(self.N):
                if i != j:
                    rij = self.__nim(self.r[i], self.r[j], self.SIZE)
                    _rij = self.__lenght(rij)
                    ff = self.__force_LD(_rij)
                    _dr = rij.copy()
                    _dr = self.__normalaize(_dr)
                    self.f[i] += _dr * ff

    def __correct_coord(self, coord, left_boundary, right_boundary):
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
    def __pbc(self, r, size):
        return np.array([self.__correct_coord(r[i], 0, size) for i in np.arange(r.size)]) 

    def integrate(self):
        for i in np.arange(self.N):
            r_tmp = self.r[i].copy() # хранит t-dt
            self.r[i] = self.__verle_r(self.r[i], self.dr[i], self.f[i], self.m[i], self.dt) # вычсиляем координаты по алгоритму верле для t+dt
            self.dr[i] = self.r[i] - r_tmp # вычисляем разность координат между t+dt и t
            if self.__lenght(self.dr[i]) > self.L_FREE_MOTION : # если слишком большая разноть, то значть что-то не так
                if not self.large_motion:
                    self.large_motion = True
                    print(' - too large motion detected')
            self.v[i] = self.__verle_v(self.dr[i], self.dt) # вычисялем скорость по алгоритму верле
            self.r[i] = self.__pbc(self.r[i], self.SIZE)
            self.fulltime += self.dt
    

def main():
    steps = 1000
    st = time.time()
    system = MDSystem(n_atoms=10, cube_size=3, dim=3, speed=1)
    system.init_system(zero_v=False)
    with ChargingBar('Steps', max=steps, suffix='%(percent)d%%') as bar:
        for i in np.arange(1, steps+1):
            system.calc_forces()
            system.integrate()
            system.print_to_file()
            bar.next()
    ed = time.time()
    print(f'Done! Time: {(ed - st):.3f}')

if __name__ == '__main__':
    main()