import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import math
matplotlib.rcParams['figure.figsize'] = [16.0, 8.0]
from utils import *
from const import *


def test_spherical(n_sample):
    x_ls, y_ls, z_ls = [], [], []
    for i in range(n_sample):
        x, y, z = spherical()
        x_ls.append(x)
        y_ls.append(y)
        z_ls.append(z)
    fig = plt.figure(figsize=(16, 8))
    ax1 = fig.add_subplot(121)
    ax2 = fig.add_subplot(122)
    ax1.title.set_text('xy Plot')
    ax2.title.set_text('xz Plot')

    ax1.set_xlim([-1,1])
    ax1.set_ylim([-1,1])
    ax1.scatter(x_ls, y_ls)
    ax1.plot()
    
    ax2.set_xlim([-1,1])
    ax2.set_ylim([-1,1])
    ax2.scatter(x_ls, z_ls)
    ax2.plot()
    plt.show()


def test_v(n_sample):
    a = make_spherical_df(n_sample, 0.5)
    x_ls = []
    y_ls = []
    for i in range(n_sample):
        x_ls.append(a[1][i][0])
        y_ls.append(a[1][i][1])
    plt.scatter(x_ls, y_ls)
    plt.show()
    return


def test_rv(n_sample):
    a =  make_spherical_df(n_sample, 0.4)
    w = calc_W(n_sample, a[0])
    k = calc_K(n_sample, a[1])
    print(k / abs(w))
    return