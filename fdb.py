import matplotlib.pyplot as plt
import math
import numpy as np

#------------------------------------------------------------------------
# start condition
f0        = 1.143 * 10**6
p0        = 0.0
r0        = 0.0
Sr0       = 0.0
ra0       = 0.0

# etha calculation
mp        = 938.256*10**6 
E         = 640 * 10**6
gamma     = 1 + E/mp
alpha     = 0.0513
etha      = abs(alpha - 1 / gamma**2)

# syncrotron frequency
f_s       = 701.458
omega_s   = 2 * math.pi * f_s

# deflection caused by magnetic field non-uniformity 
rB        = 0.0001

# damping parameters
omega_ref = 6 * 10**3
thau_I    = 0.25 * 10**(-4)
g_i       = 0.05
g_a       = 3.0

# run parameters
dt       = 10**(-6)
Tcycle   = 2*10**(-5)
Ncycle   = 1000
Niternal = math.floor(Tcycle/dt)

#------------------------------------------------------------------------

def oscillator(f, df, p, r, step):
    r_next = r + p * step
    p_next = p - omega_s**2 * (r-rB - df / f / etha) * step
    return r_next, p_next

def compensation(Sr, ra, r, r_diff, step):
    Sr_next = Sr + r * omega_ref * step
    ra_next = ra - ra / thau_I * step + r_diff
    return Sr_next, ra_next

def stimulus(f, Sr, ra):
    df = - f * (g_i * Sr + g_a * ra)
    return df

def run(Ncycle, Niternal, f0, p0, r0, Sr0, ra0):
    Sr           = Sr0
    ra           = ra0
    r            = r0
    p            = p0
    f            = f0
    df           = 0

    t_data       = [0 ]
    r_data       = [r ]
    p_data       = [p ]
    Sr_data      = [Sr]
    ra_data      = [ra]
    r_diff_data  = [0 ]
    f_data       = [f0]

    r_cycle_prev = 0
    r_cycle_diff = 0

    t = 0
        
    for i in range(Ncycle):
        
        for j in range(Niternal):

            t = t + dt
            r, p = oscillator(f, df, p, r, dt)
            
            t_data.append(t)
            r_data.append(r)
            p_data.append(p)
            Sr_data.append(Sr)
            ra_data.append(ra)
            r_diff_data.append(r_cycle_diff)
            f_data.append(f)

        df = stimulus(f, Sr, ra)
        f  = f + df

        r_cycle_diff = r - r_cycle_prev
        r_cycle_prev = r
        Sr, ra = compensation(Sr, ra, r, r_cycle_diff, Tcycle)

    return t_data, r_data, p_data, r_diff_data, Sr_data, ra_data, f_data

#------------------------------------------------------------------------
t, r, p, r_diff, Sr, ra, f = run(Ncycle, Niternal, f0, p0, r0, Sr0, ra0)

fig, axs = plt.subplots(6, 1)

axs[0].plot(t, r)
axs[0].plot(t[::Niternal], r[::Niternal], 'b.')
axs[0].set_xlabel("t")
axs[0].set_ylabel("r")
axs[0].grid()

axs[1].plot(t, p)
axs[1].set_xlabel("t")
axs[1].set_ylabel("p")
axs[1].grid()

axs[2].plot(t, r_diff)
axs[2].set_xlabel("t")
axs[2].set_ylabel("r_diff")
axs[2].grid()

axs[3].plot(t, ra)
axs[3].set_xlabel("t")
axs[3].set_ylabel("ra")
axs[3].grid()

axs[4].plot(t, Sr)
axs[4].set_xlabel("t")
axs[4].set_ylabel("Sr")
axs[4].grid()

axs[5].plot(t, f)
axs[5].set_xlabel("t")
axs[5].set_ylabel("f")
axs[5].grid()

#------------------------------------------------------------------------
# reference, no damping

g_i       = 0.0
g_a       = 0.0

t, r, p, r_diff, Sr, ra, f = run(Ncycle, Niternal, f0, p0, r0, Sr0, ra0)

axs[0].plot(t, r)
axs[0].set_xlabel("t")
axs[0].set_ylabel("r")

#------------------------------------------------------------------------

plt.show()
