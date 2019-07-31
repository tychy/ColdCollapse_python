import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import math
matplotlib.rcParams['figure.figsize'] = [8.0, 8.0]

G = 1
M = 1
R = 1
eps = 0.03125
def spherical():
    x = np.random.randn()
    y = np.random.randn()
    z = np.random.randn()
    if (x**2 + y**2 + z**2)<= 1:
        return np.array([x, y, z])
    else:
        return spherical()


def calc_W(N, xarray):
    W = 0
    m_particle = M / N
    for i in range(N-1):
        for j in range(i+1,N):
            dist = np.linalg.norm(xarray[i] - xarray[j])
            W += -1 * m_particle**2 * math.pow(dist**2 + eps**2, -1/2)
    return W

def make_spherical_df(N, rv):
    # 一様乱数で位置の配列を作成
    xarray = []
    for i in range(N):
        xarray.append(spherical())
    
    # 初期速度分布を計算するための分散を計算
    W = calc_W(N, xarray)
    std = math.sqrt(2*rv*abs(W)/(3*M))
    
    varray = []
    for i in range(N):
        varray.append(np.random.normal(0, std, (3)))
    return xarray, varray


def calc_K(N, varray):
    K = 0
    for i in range(N):
        K += (M / N) * (np.linalg.norm(varray[i])**2) / 2
    return K


def calc_force(N, xarray, varray):
    m_p = M / N # m_particle
    acarray = [np.array([0.0, 0.0, 0.0])] * N
    for i in range(N-1):
        for j in range(i+1, N):
            xi = xarray[i]
            xj = xarray[j]
            r_vec = xj - xi
            dist = np.linalg.norm(r_vec)
            r3inv = math.pow(dist**2 + eps**2, - 3 / 2)
            acarray[i] = acarray[i] + m_p * r3inv * r_vec
            acarray[j] = acarray[j] - m_p * r3inv * r_vec
    return acarray


def leap_frog(N, dt, x0, v0, a0):
    v_half = [v0[i] + a0[i] * dt / 2 for i in range(N)]
    x1 = [x0[i] + v_half[i] * dt for i in range(N)]
    a1 = calc_force(N, x0, v0)
    v1 = [v_half[i] + a1[i] * dt / 2 for i in range(N)]
    return x1, v1, a1

def main():
    N = 1024 #　粒子数
    dt = 0.03125  #　刻み幅
    t_end = 10  # 終了時刻
    rv = 0.5  # ビリアル係数
    
    xarray, varray = make_spherical_df(N, rv)
    acarray = [np.array([0.0, 0.0, 0.0])] * N
    KW = []
    loss = []
    KWt = []
    E = calc_K(N, varray) + calc_W(N ,xarray)
    
    idx = 0 
    print("steps", t_end/dt)
    for t in np.linspace(0.0, t_end, t_end/dt):
        xarray, varray, acarray = leap_frog(N, dt, xarray, varray, acarray)
        
        KW.append(calc_K(N, varray)/abs(calc_W(N ,xarray)))
        loss.append((calc_K(N, varray)+calc_W(N ,xarray))/E)
        KWt.append(t)
        
        
        if idx % 30 == 0:
            print(t)
            xn = np.array(xarray)
            plt.scatter(xn[:,0], xn[:, 1])
            plt.title("Main")
            plt.xlabel("x axis")
            plt.ylabel("y axis")
            plt.xlim(-2, 2)
            plt.ylim(-2, 2)
            plt.show()
            
            # aa
            plt.xlabel("t")
            plt.ylabel("virial ratio")
            plt.plot(KWt, KW)
            plt.show()
            
            plt.xlabel("t")
            plt.ylabel("energy loss")
            plt.plot(KWt, loss)
            plt.show()
            
            plt.title("r hist")
            plt.hist(np.linalg.norm(np.array(xarray), axis=1))
            plt.show()
            plt.title("sum(v) hist")
            plt.hist(np.sum(np.array(varray), axis=1))
            plt.show()
            plt.title("sum(a) hist")
            plt.hist(np.sum(np.array(acarray),axis=1))
            plt.show()
        idx += 1


if __name__ == '__main__':
    main()