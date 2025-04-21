import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import h as h_SI, c as c_SI, G as G_SI, physical_constants

# — convert constants to cgs —
h    = h_SI * 1e7               # Planck’s constant [erg·s]
c    = c_SI * 1e2               # speed of light [cm/s]
G    = G_SI * 1e3               # gravitational constant [cm³/(g·s²)]
m_n  = physical_constants['neutron mass'][0] * 1e3   # neutron mass [g]

# — define our natural length, mass, and energy‐density scales —
r_tilde   = (1/np.pi) * (h/(m_n*c))**1.5 * (c/np.sqrt(m_n*G))
M_tilde   = (1/np.pi) * (h/(m_n*c))**1.5 * (c**3/np.sqrt(m_n*G**3))
eps_scale = np.pi**2 * m_n**4 * c**5 / h**3

# — assign scales for use in integration →
r_scale = r_tilde    # [cm]
m_scale = M_tilde    # [g]
Msun    = 1.989e33   # convert back to solar masses

# — unit conversions for EOS →
MeV_to_erg = 1.602176634e-6
fm3_to_cm3 = 1e39
conv_P     = MeV_to_erg * fm3_to_cm3

# — polynomial coefficients for each EOS fit —
eos_coefs = {
    'AV14+UVII': [2.6511,  76.744,  -183.611, 459.906,  -122.832],
    'UV14+UVII': [7.57891, -1.23275, 227.384, -146.596, 324.823, -120.355],
    'UV14+TNI':  [6.33041, -28.1793, 288.397,  -65.2281]
}

def build_eos_table(coefs, n_max=2.5, N=1000):
    """
    Precompute P(n), ε(n), and dP/dn on a grid of densities.
    Returns arrays (n, P_dimless, eps_dimless, dPdn_dimless).
    """
    n     = np.linspace(0, n_max, N)
    P_t   = np.zeros(N)
    eps_t = np.zeros(N)
    dpdn  = np.zeros(N)

    for i, ni in enumerate(n):
        # Evaluate E(n) polynomial and its derivatives
        E_n     = sum(c * ni**j for j, c in enumerate(coefs))
        dE_dn   = sum(j * c * ni**(j-1) for j, c in enumerate(coefs) if j>0)
        d2E_dn2 = sum(j*(j-1) * c * ni**(j-2)
                      for j, c in enumerate(coefs) if j>1)

        # Physical pressure, energy density, and slope
        P_phys   = conv_P * ni**2 * dE_dn
        eps_phys = (ni * fm3_to_cm3) * (E_n + 939.0) * MeV_to_erg
        dpdn_phys= conv_P * (2*ni*dE_dn + ni**2 * d2E_dn2)

        # Convert to dimensionless units
        P_t[i]   = P_phys   / eps_scale
        eps_t[i] = eps_phys / eps_scale
        dpdn[i]  = dpdn_phys/ eps_scale

    return n, P_t, eps_t, dpdn

def integrate_tov(nc, n_grid, P_t, eps_t, dpdn_t, h=1e-4):
    """
    Simple RK4 integration of the TOV eqns given central density nc.
    Returns (M in M☉, R in km).
    """
    # starting radius and mass
    R = 1e-6
    eps0 = np.interp(nc, n_grid, eps_t)
    u = (4/3)*np.pi * eps0 * R**3
    t = nc

    def derivs(R, u, t):
        # interpolate values at current t
        P   = np.interp(t, n_grid, P_t)
        eps = np.interp(t, n_grid, eps_t)
        dp  = np.interp(t, n_grid, dpdn_t)
        # stop if out of bounds or dPdn≤0
        if t<=0 or t>=n_grid[-1] or dp<=0:
            return None
        du  = 4*np.pi * R**2 * eps
        dPdr= - (eps+P)*(u+4*np.pi*R**3*P)/(R*(1-2*u/R))
        dt  = dPdr / dp
        return du, dt

    # RK4 loop
    while True:
        if np.interp(t, n_grid, P_t) < 0:
            break
        k1 = derivs(R, u, t)
        k2 = derivs(R+0.5*h, u+0.5*h*k1[0], t+0.5*h*k1[1]) if k1 else None
        k3 = derivs(R+0.5*h, u+0.5*h*k2[0], t+0.5*h*k2[1]) if k2 else None
        k4 = derivs(R+  h,   u+  h*k3[0], t+  h*k3[1])   if k3 else None
        if None in (k1, k2, k3, k4):
            break
        u += (h/6)*(k1[0]+2*k2[0]+2*k3[0]+k4[0])
        t += (h/6)*(k1[1]+2*k2[1]+2*k3[1]+k4[1])
        R += h

    # convert to M☉ and km
    M    = u * m_scale / Msun
    R_km = R * r_scale / 1e5
    return M, R_km

# — main driver: sweep over central densities and plot —
nc_vals = np.linspace(0.01, 2.3, 100)
results = {}
for name, coefs in eos_coefs.items():
    n, P_t, eps_t, dpdn_t = build_eos_table(coefs)
    Ms, Rs = [], []
    for nc in nc_vals:
        M, R_km = integrate_tov(nc, n, P_t, eps_t, dpdn_t)
        Ms.append(M)
        Rs.append(R_km)
    results[name] = (np.array(Ms), np.array(Rs))

# — observational constraints for two pulsars —
obs = {
    'J0030+0451': {'M':1.44, 'M_err':[0.14,0.15], 'R':13.02, 'R_err':[1.06,1.24]},
    'J0740+6620': {'M':2.08, 'M_err':[0.07,0.07], 'R':13.7,  'R_err':[1.5, 2.6 ]},
}

# 1) Mass vs. density
plt.figure()
for name, (Ms, _) in results.items():
    plt.plot(nc_vals, Ms, label=name)
plt.xlabel('n₀ (fm⁻³)')
plt.ylabel('M (M☉)')
plt.title('Mass vs Central Density')
plt.ylim(0,3)
plt.grid(); plt.legend()

# 2) Radius vs. density
plt.figure()
for name, (_, Rs) in results.items():
    plt.plot(nc_vals, Rs, label=name)
plt.xlabel('n₀ (fm⁻³)')
plt.ylabel('R (km)')
plt.title('Radius vs Central Density')
plt.ylim(0,15)
plt.grid(); plt.legend()

# 3) Mass–Radius with error bars
plt.figure()
for name, (Ms, Rs) in results.items():
    plt.plot(Rs, Ms, label=name)
for key, d in obs.items():
    plt.errorbar([d['R']], [d['M']],
                 xerr=np.array(d['R_err']).reshape(2,1),
                 yerr=np.array(d['M_err']).reshape(2,1),
                 fmt='o', capsize=5, label=key)
plt.xlabel('R (km)')
plt.ylabel('M (M☉)')
plt.title('Mass–Radius Relation')
plt.ylim(0,3)
plt.grid(); plt.legend()

plt.tight_layout()
plt.show()
