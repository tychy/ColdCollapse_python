import numpy as np
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import math
matplotlib.rcParams['figure.figsize'] = [15.0, 5.0]
from utils import *
from const import *

def main():
    N = 200 #　粒子数
    t_end = 5  # 終了時刻
    step = 500 * t_end
    rv = 0.1  # ビリアル係数
    
    xarray, varray = make_spherical_df(N, rv)
    acarray = [np.array([0.0, 0.0, 0.0])] * N
    KW = []
    loss = []
    KWt = []
    E = calc_K(N, varray) + abs(calc_W(N ,xarray))
    
    idx = 0 
    for i in range(step):
        dt = t_end / step
        t = dt * i
        xarray, varray, acarray = leap_frog(N, dt, xarray, varray, acarray)
        # print(varray)
        KW.append(calc_K(N, varray)/abs(calc_W(N ,xarray)))
        loss.append((calc_K(N, varray)+abs(calc_W(N ,xarray)))/E)
        KWt.append(t)
        if idx % 30 == 0:
            print(t)
            xn = np.array(xarray)
            plt.subplot(1, 3, 1)
            plt.scatter(xn[:,0], xn[:, 1])
            plt.title("Main")
            plt.xlabel("x axis")
            plt.ylabel("y axis")
            plt.xlim(-10, 10)
            plt.ylim(-10, 10)
            plt.plot()
            
            # a
            plt.subplot(1, 3, 2)
            plt.xlabel("t")
            plt.ylabel("virial ratio")
            plt.plot(KWt, KW)
            plt.plot()
            plt.subplot(1, 3, 3)
            plt.xlabel("t")
            plt.ylabel("energy loss")
            plt.plot(KWt, loss)
            plt.plot()
            plt.savefig("result/{:.2f}time.png".format(t))
            plt.close()
        idx += 1
    return


if __name__ == "__main__":
    main()
