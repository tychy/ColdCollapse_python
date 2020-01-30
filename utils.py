import numpy as np
import math
from math import pow, sqrt
from const import *


"""
fail
def rad_spherical(number_of_particles):
    radius = np.random.uniform(0.0,1.0, (number_of_particles)) 
    theta = np.random.uniform(0.,1.,(number_of_particles)) * np.pi
    phi = np.random.uniform(0.,1.,(number_of_particles)) * 2 * np.pi
    x = radius * np.sin(theta) * np.cos(phi)
    y = radius * np.sin(theta) * np.sin(phi)
    z = radius * np.cos(theta)
    return (x,y,z)
"""
def spherical():
    x = np.random.randn()
    y = np.random.randn()
    z = np.random.randn()
    if (np.linalg.norm((x, y, z)))<= 1:
        return np.array([x, y, z])
    else:
        return spherical()


def calc_W(N, xarray):
    W = 0
    m_particle = M / N
    for i in range(N-1):
        for j in range(i+1,N):
            dist = np.linalg.norm(xarray[i] - xarray[j])
            W += -1 * pow(m_particle, 2) * pow(pow(dist + eps, 2), -1/2)
    return W


def make_spherical_df(N, rv):
    # 一様乱数で位置の配列を作成
    xarray = []
    for i in range(N):
        xarray.append(spherical())
    
    # 初期速度分布を計算するための分散を計算
    W = calc_W(N, xarray)
    std = math.sqrt(2*rv*abs(W)/(3*M))
    #print(std)
    varray = []
    for i in range(N):
        varray.append(np.random.normal(0, std, (3)))
    print(varray)
    return xarray, varray


def calc_K(N, varray):
    K = 0
    for i in range(N):
        K += (M / N) * pow(np.linalg.norm(varray[i]), 2) / 2
    return K


def calc_force(N, xarray, varray):
    m_p = M / N # m_particle
    acarray = [np.array([0.0, 0.0, 0.0])] * N
    for i in range(N-1):
        for j in range(i+1, N):
            xi = xarray[i]
            xj = xarray[j]
            r_vec = xj - xi
            # print(r_vec)
            dist = np.linalg.norm(r_vec)
            r3inv = pow(pow(dist, 2) + pow(eps, 2), - 3 / 2)
            acarray[i] = acarray[i] + m_p * r3inv * r_vec
            acarray[j] = acarray[j] - m_p * r3inv * r_vec
    return acarray


def leap_frog(N, dt, x0, v0, a0):
    v_half = [v0[i] + a0[i] * dt / 2 for i in range(N)]
    # print(v_half)
    x1 = [x0[i] + v_half[i] * dt for i in range(N)]
    a1 = calc_force(N, x0, v0)
    v1 = [v_half[i] + a1[i] * dt / 2 for i in range(N)]
    return x1, v1, a1