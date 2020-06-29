#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 15 14:13:57 2020

@author: rdc
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator


def loss_huber(x, k):
    res = (np.abs(x) <= k) * x ** 2 * 0.5
    res += (np.abs(x) > k) * (k * np.abs(x) - k**2/2)
    return(res)


def score_huber(x, k):
    res = (np.abs(x) <= k) * x
    res += (np.abs(x) > k) * (k * np.sign(x))
    return(res)


def score_huber_deriv(x, k):
    res = x*0
    res = (np.abs(x) <= k) * 1
    return(res)


k = 0.7
x = np.arange(-4, 4, 0.0001)


fig, ax = plt.subplots(ncols=3, figsize=(8, 4), dpi=150)
ax[0].plot(x, loss_huber(x, k))
ax[0].yaxis.set_major_locator(MaxNLocator(integer=True))
ax[0].set(xlabel=r'$u$', ylabel=r"$\rho_H(u,k={0})$".format(k))


ax[1].plot(x, score_huber(x, k))
ax[1].set_yticks([-k, 0, k])
ax[1].set(xlabel=r'$u$', ylabel=r"$\psi_H(u,k={0})$".format(k))


ax[2].plot(x, score_huber_deriv(x, k))
ax[2].yaxis.set_major_locator(MaxNLocator(integer=True))
ax[2].set(xlabel=r'$u$', ylabel=r"$\psi^\prime_H(u,k={0})$".format(k))

plt.subplots_adjust(wspace=1)
plt.savefig(
        fname='huberized.png',
        dpi=150,
        transparent=True
 )
plt.show()
