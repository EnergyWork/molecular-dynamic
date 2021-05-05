import numpy as np
import time
from progress.bar import ChargingBar

class MDSystem:
    def __init__(self, n_atoms, cube_size, dim, speed):
        self.__N = n_atoms
        self.__SIZE = cube_size
        self.__DIM = dim
        self.__SPEED = speed
        self.__init_vars()
    
    def __init_vars(self):
        self.__large_motion = False
        self.__b_nim_error = False
        self.__L_FREE_MOTION = np.power(self.__SIZE * self.__SIZE * self.__SIZE, 1/3.0) / (2*self.__N)
        self.__dt = 0.001
        self.__eps = 1
        self.__sigma = 1
        self.__rcut = 2.5 * self.__sigma
        self.__rmin = 0.0001
        self.__fulltime = self.__dt

    def __init_arrays(self):
        self.__r = np.array([np.zeros([self.__DIM]) for _ in np.arange(self.__N)])
        self.__dr = np.array([np.zeros([self.__DIM]) for _ in np.arange(self.__N)])
        self.__v = np.array([np.zeros([self.__DIM]) for _ in np.arange(self.__N)])
        self.__f = np.array([np.zeros([self.__DIM]) for _ in np.arange(self.__N)])
        self.__m = np.array([np.zeros([self.__DIM]) for _ in np.arange(self.__N)])

    def __clear_file(self):
        with open('.\\py\\coords.xmol', 'w', encoding='utf-8') as file:
            pass
        with open('.\\py\\debug.txt', 'w', encoding='utf-8') as fdbg:
            pass

    def __new_frame(self, file_handle):
        file_handle.write(f'{self.__N}\n')

    def print_to_file(self):    
        with open('.\\py\\coords.xmol', 'a', encoding='utf-8') as file:
            self.__new_frame(file)
            file.write(f'***** time = * {self.__fulltime:.4} *****\n')
            for i in np.arange(self.__N):
                file.write(f'Ar {self.__r[i][0]:.4f} {self.__r[i][1]:.4f} {self.__r[i][2]:.4f} {self.__v[i][0]:.4f} {self.__v[i][1]:.4f} {self.__v[i][2]:.4f}\n')

    def init_system(self, zero_v=False):
        self.__clear_file()
        self.__init_arrays()
        if not zero_v:
            for i in np.arange(self.__N):
                self.__v[i] = np.array([((self.__SPEED - (-self.__SPEED)) * np.random.random() + (-self.__SPEED)) for _ in np.arange(self.__DIM)])
        k = np.ceil(np.power(self.__N, 1.0 / 3))
        dh = self.__SIZE / k
        self.__m = np.ones(self.__N)
        counter = 0
        for x in np.arange(k):
            for y in np.arange(k):
                for z in np.arange(k):
                    if counter < self.__N:
                        self.__r[counter] = np.array([(x + 1.0 / 2) * dh, (y + 1.0 / 2) * dh, (z + 1.0 / 2) * dh])
                        self.__dr[counter] = np.array([self.__v[counter][0] * 2 * self.__dt, self.__v[counter][1] * 2 * self.__dt, self.__v[counter][2] * 2 * self.__dt])
                        with open('.\\py\\debug.txt', 'a', encoding='utf-8') as fdbg:
                            strdbg = f"{ self.__dr[counter][0]:.9f}  { self.__dr[counter][1]:.9f}  { self.__dr[counter][2]:.9f}\n"
                            fdbg.write(strdbg)
                        if self.__lenght(self.__dr[counter]) > self.__L_FREE_MOTION:
                            print(' - init error')
                        counter += 1
                    else:
                        break
        self.print_to_file()

    def __force_LD(self, r):
        if r > self.__rcut:
            return 0
        if r < self.__rmin:
            return self.__force_LD(self.__rmin)
        x = self.__sigma / r
        return -48 * self.__eps / self.__sigma * (np.power(x, 13, dtype=np.float64) - 0.5 * np.power(x, 7, dtype=np.float64))

    def __verle_r(self, r, dr, f, m, dt):
        return r + (dr + (f / (2 * m)) * np.square(dt))

    def __verle_v(self, dr, dt):
        return dr / (2 * dt)

    def __nim_fix(self, coord, size):
        if coord >= size / 2.0:
            coord = coord - size
        elif coord <= -size / 2.0:
            coord = coord + size
        return coord

    # nearist image method
    def __nim(self, r1, r2, size):
        crds = np.array([(-(r1[i]-r2[i])) for i in np.arange(np.size(r1))]) # -(r1[i]-r2[i])
        dist = np.array([self.__nim_fix(crd, size) for crd in crds])
        if self.__lenght(dist) > self.__lenght(np.array([size, size, size])):
            if not self.__b_nim_error:
                self.__b_nim_error = True
                print(' - nim error')
        return dist

    def __lenght(self, r):
        return np.sqrt(np.sum([np.square(r[i]) for i in np.arange(np.size(r))]))

    def __normalaize(self, vec):
        return vec / self.__lenght(vec)

    def calc_forces(self):
        for i in np.arange(self.__N):
            self.__f[i] = np.zeros(self.__DIM)
            for j in np.arange(self.__N):
                if i != j:
                    rij = self.__nim(self.__r[i], self.__r[j], self.__SIZE)
                    _rij = self.__lenght(rij)
                    ff = self.__force_LD(_rij)
                    # with open('.\\py\\debug.txt', 'a', encoding='utf-8') as fdbg:
                    #     strdbg = f"{ff:.9f}\n"
                    #     fdbg.write(strdbg)
                    _dr = rij.copy()
                    _dr = self.__normalaize(_dr)
                    self.__f[i] += _dr * ff
            # with open('.\\py\\debug.txt', 'a', encoding='utf-8') as fdbg:
            #     strdbg = f"{ self.__f[i][0]:.9f}  { self.__f[i][1]:.9f}  { self.__f[i][2]:.9f}\n"
            #     fdbg.write(strdbg)

    def __correct_coord(self, coord, left_boundary, right_boundary):
        l = right_boundary - left_boundary
        d = 0
        if coord >= right_boundary:
            d = np.floor(coord / l)
            coord = coord - d * l
        elif coord < left_boundary:
            d = np.floor(np.abs(coord) / l)
            coord = right_boundary - (np.abs(coord) - d * l)
        return coord

    #periodic boundary condition
    def __pbc(self, r, size):
        return np.array([self.__correct_coord(r[i], 0, size) for i in np.arange(self.__DIM)]) 

    def integrate(self):
        for i in np.arange(self.__N):
            r_tmp = self.__r[i].copy() # хранит t-dt
            self.__r[i] = self.__verle_r(self.__r[i], self.__dr[i], self.__f[i], self.__m[i], self.__dt) # вычсиляем координаты по алгоритму верле для t+dt
            self.__dr[i] = self.__r[i] - r_tmp # вычисляем разность координат между t+dt и t
            # with open('.\\py\\debug.txt', 'a', encoding='utf-8') as fdbg:
            #     strdbg = f"{ self.__r[i][0]:.9f}  { self.__r[i][1]:.9f}  { self.__r[i][2]:.9f}\n"
            #     fdbg.write(strdbg)
            if self.__lenght(self.__dr[i]) > self.__L_FREE_MOTION : # если слишком большая разноть, то значть что-то не так
                if not self.__large_motion:
                    self.__large_motion = True
                    print(' - too large motion detected')
            self.__v[i] = self.__verle_v(self.__dr[i], self.__dt) # вычисялем скорость по алгоритму верле
            self.__r[i] = self.__pbc(self.__r[i], self.__SIZE)
            self.__fulltime += self.__dt
    

def main():
    #np.random.seed(42)
    steps = 1000
    st = time.time()
    system = MDSystem(n_atoms=10, cube_size=3, dim=3, speed=2)
    system.init_system(zero_v=True)
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