# Heat Transfert polar convection 

# Content

* [Objectives](#objectives)
* [About this repository](#about-this-repository)
* [Code description](#code-description)
* [Résultats](#résultats)

# Objectives
expliquer le but de ce travail = traduction python et voir comment se comporte différente intrusion de chaleur mais aussi que se passe t'il si l on change la tolérence

# About this repository
In this repository you can find:
- the main Python code [Python_porous convection.py](link)
- the Matlab code base [Matlab_porous convection.m](link)
- the code that allows the differences between the two main codes to be perceived [Python_Matlab_Comparator.py](link)

# Code description
* [Darcy velocity](#link)
* [Compressible Mass conservation](#link)
* [Fluid pressure time update](#link)
* [Temperature diffusion](#link)
* [Temperature advection](#link)

## Darcy velocity
description

$q_x = -\frac{k}{\eta_f}\frac{\partial P_f}{\partial x}$

$q_y = -\frac{k}{\eta_f} \left( \frac{\partial P_f}{\partial y} + \rho g \right)$

```md
qx[1:-1, 1:-1] = -k_etaf*(np.diff(Pf[:, 1:-1], axis=0)/dx)
qy[1:-1, 1:-1] = -k_etaf*(np.diff(Pf[1:-1, :], axis=1)/dy 
+ (rho_g[1:-1, :-1] + rho_g[1:-1, 1:])/2)
```

⤴️ [_back to content_](#content)

## Compressible Mass conservation
description

$\frac{\partial P_f}{\partial t} = - \left( \frac{\partial q_x}{\partial x} + \frac{\partial q_y}{\partial y} \right) \frac{1}{\beta_f}$

```md
dPfdt   = -(np.diff(qx, axis=0)/dx + np.diff(qy, axis=1)/dy)/betaf
```

⤴️ [_back to content_](#content)

##  Fluid pressure time update
description

$\Delta P_f = \frac{\Delta P_f}{\Delta t} \Delta t
P_f^{new} = P_f^{old} + \Delta P_f$

$\frac{\partial P_f}{\partial t} \approx \frac{\Delta P_f}{\Delta t} = \frac{P_f^{new} - P_f^{old}}{\Delta t}$

```md
deltaPf = dt*dPfdt
Pf      = Pf + deltaPf 
max_delta_Pf.append(np.max(np.abs(deltaPf)))
if max_delta_Pf[len(max_delta_Pf)-1]/np.max(np.abs(Pf)) < tol:
    break
```

⤴️ [_back to content_](#content)

## Temperature diffusion
description

$\frac{\partial T}{\partial t} = \frac{\partial}{\partial x} \left( \frac{\lambda}{\rho c_p} \frac{\partial T}{\partial x} \right) + \frac{\partial}{\partial y} \left( \frac{\lambda}{\rho c_p} \frac{\partial T}{\partial y} \right)$

$T_f^{new} = T_f^{old} + \frac{\Delta T}{\Delta t} \Delta t$

$\frac{\partial T}{\partial t} \approx \frac{\Delta T}{\Delta t} = \frac{T^{new} - T^{old}}{\Delta t}$

```md
dTdt = np.diff(lam_rhoCp * np.diff(T[:, 1:-1], axis=0) / dx, axis=0) / dx \
+ np.diff(lam_rhoCp * np.diff(T[1:-1, :], axis=1) / dy, axis=1) / dy
T[1:-1, 1:-1] += dTdt * dt
```

⤴️ [_back to content_](#content)

## Temperature advection
description

$\frac{\partial T}{\partial t} = - q_x \frac{\partial T}{\partial x} - q_y \frac{\partial T}{\partial y}$

$T_f^{new} = T_f^{old} + \frac{\Delta T}{\Delta t} \Delta t$

```md
T[:, 1:] -= dt * np.maximum(0, qy[:,1:-1]) * np.diff(T, axis=1) / dy
T[:,:-1] -= dt * np.minimum(0, qy[:,1:-1]) * np.diff(T, axis=1) / dy
T[1:,:]  -= dt * np.maximum(0, qx[1:-1,:]) * np.diff(T, axis=0) / dx
T[:-1,:] -= dt * np.minimum(0, qx[1:-1,:]) * np.diff(T, axis=0) / dx
```

⤴️ [_back to content_](#content)
# Résults
* [Difference between Python and Matlab code](Difference-between-Python-and-Matlab-code)
* [Different heat intrusion situation](link doc)

## Difference between Python and Matlab code

Here are some screenshots at different time steps of the Python (left) and Matlab (right) visualizations

image de la comparaison des codes

## Different heat intrusion situation

* variante A
vidéo
* variante B
vidéo
* variante C
vidéo
* variante D
vidéo
* variante E
vidéo

⤴️ [_back to code description_](#code-description)
