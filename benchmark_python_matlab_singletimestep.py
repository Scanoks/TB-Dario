import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import os

# chose which file
it      = 0     # time step

# import
abspath = os.path.abspath(__file__) # https://stackoverflow.com/questions/1432924/python-change-the-scripts-working-directory-to-the-scripts-own-directory
dname = os.path.dirname(abspath)
os.chdir(dname)

qx_python = np.genfromtxt("output_P/output_it_"+ str(it) +"_qx.csv"  , delimiter=',')
qy_python = np.genfromtxt("output_P/output_it_"+ str(it) +"_qy.csv"  , delimiter=',')
Pf_python = np.genfromtxt("output_P/output_it_"+ str(it) +"_Pf.csv"  , delimiter=',')
T_python  = np.genfromtxt("output_P/output_it_"+ str(it) +"_T.csv"   , delimiter=',')

qx_matlab = np.genfromtxt("output_M/output_it_"+ str(it) +"_qx_M.csv", delimiter=',')
qy_matlab = np.genfromtxt("output_M/output_it_"+ str(it) +"_qy_M.csv", delimiter=',')
Pf_matlab = np.genfromtxt("output_M/output_it_"+ str(it) +"_Pf_M.csv", delimiter=',')
T_matlab  = np.genfromtxt("output_M/output_it_"+ str(it) +"_T_M.csv" , delimiter=',')

# ...needed for plotting...
Ly = 1e3                        # Length in y, m
deltaT = 100                    # Temperature difference, K
Lx_Ly       = 2                 # pour définir après que Lx est 2* plus grand que Ly
Lx = Lx_Ly * Ly                 # Length in x, m
nx = 150                         # Number of grid points in x 
ny = int(nx * Ly / Lx)          # Number of grid points in y  ## fix = int sur python
dx = Lx/(nx-1)                  # défini la longueur de dx
dy = Ly/(ny-1)                  # défini la longueur de dy
x, y = np.meshgrid(np.arange(-Lx/2, Lx/2 + dx, dx), np.arange(-Ly/2, Ly/2 + dy, dy), indexing='ij')
xqx, yqx = np.meshgrid(np.arange(-Lx/2 - 0.5*dx, Lx/2 + 0.5*dx + dx, dx), np.arange(-Ly/2, Ly/2 + dy, dy), indexing='ij')
xqy, yqy = np.meshgrid(np.arange(-Lx/2, Lx/2 + dx, dx), np.arange(-Ly/2- 0.5*dy, Ly/2 + 0.5*dy + dy, dy), indexing='ij')

lam_rhoCp = 1e-6                # Heat diffusivity, m^2/s
k_etaf = 1e-12                  # Permeability(m^2)/fluid viscosity (Pa*s) = m^2/(Pa*s)
time_char = Ly**2/lam_rhoCp     # m^2/(m^2/s) = s
Pr_char = lam_rhoCp/k_etaf      # (m^2/s)/(m^2/(Pa*s)) = Pa
q_char = Ly/time_char           # m/s

# Non-dimensionalization
# We need 4 independent characteristic scales, we chose:
# -> characteristic time scale,        tc = time_char
# -> characteristic length scale,      Lc = Ly
# -> characteristic pressure scale,    Pc = Pr_char
# -> characteristic temperature scale, Tc = deltaT
# A useful derived characteristic scale:
# -> characteristic (Darcy) velocity scale,    qc = Lc/tc = q_char


fig = plt.figure()

# plot qx
ax = fig.add_subplot(4,3,1)
pc = ax.pcolor(xqx/Ly,yqx/Ly,qx_python/q_char,cmap='jet')
ax.set_aspect('equal','box')
plt.locator_params(axis='both', nbins=3)
divider = make_axes_locatable(ax) # https://stackoverflow.com/questions/18195758/set-matplotlib-colorbar-size-to-match-graph
cax = divider.append_axes("right", size="5%", pad=0.05)
cb = cb = plt.colorbar(pc,cax=cax,format='%.0e') # https://stackoverflow.com/questions/22012096/how-to-set-number-of-ticks-in-plt-colorbar, 2nd solution
cb.ax.locator_params(nbins=3)
ax.set_title('Python')
ax.set_ylabel(r'$\dfrac{y}{L_y}$')
textstr = r'$\dfrac{q_x}{q_c}$' # https://matplotlib.org/stable/gallery/text_labels_and_annotations/placing_text_boxes.html
#textstr = r'$\dfrac{q_x}{q_c} = q_x \cdot \dfrac{\rho c_p}{\lambda} \cdot L_y$' # https://matplotlib.org/stable/gallery/text_labels_and_annotations/placing_text_boxes.html
ax.text(0.025, 0.95, textstr, transform=ax.transAxes, fontsize=14, verticalalignment='top', 
    bbox=dict(boxstyle='square', facecolor='white', alpha=0.5))

ax = fig.add_subplot(4,3,2)
pc = ax.pcolor(xqx/Ly,yqx/Ly,qx_matlab/q_char,cmap='jet')
ax.set_aspect('equal','box')
plt.locator_params(axis='both', nbins=3)
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.05)
cb = plt.colorbar(pc,cax=cax,format='%.0e')
cb.ax.locator_params(nbins=3)
ax.set_title('Matlab')

ax = fig.add_subplot(4,3,3)
pc = ax.pcolor(xqx/Ly,yqx/Ly,(qx_python-qx_matlab)/q_char,cmap='jet')
ax.set_aspect('equal','box')
plt.locator_params(axis='both', nbins=3)
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.05)
cb = plt.colorbar(pc,cax=cax,format='%.0e')
cb.ax.locator_params(nbins=3)
ax.set_title('Python-Matlab')

# plot qy
ax = fig.add_subplot(4,3,4)
pc = ax.pcolor(xqy/Ly,yqy/Ly,qy_python/q_char,cmap='jet')
ax.set_aspect('equal','box')
plt.locator_params(axis='both', nbins=3)
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.05)
cb = plt.colorbar(pc,cax=cax,format='%.0e')
cb.ax.locator_params(nbins=3)
ax.set_ylabel(r'$\dfrac{y}{L_y}$')
textstr = r'$\dfrac{q_y}{q_c}$'
#textstr = r'$\dfrac{q_y}{q_c} = q_y \cdot \dfrac{\rho c_p}{\lambda} \cdot L_y$'
ax.text(0.025, 0.95, textstr, transform=ax.transAxes, fontsize=14, verticalalignment='top', 
    bbox=dict(boxstyle='square', facecolor='white', alpha=0.5))

ax = fig.add_subplot(4,3,5)
pc = ax.pcolor(xqy/Ly,yqy/Ly,qy_matlab/q_char,cmap='jet')
ax.set_aspect('equal','box')
plt.locator_params(axis='both', nbins=3)
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.05)
cb = plt.colorbar(pc,cax=cax,format='%.0e')
cb.ax.locator_params(nbins=3)

ax = fig.add_subplot(4,3,6)
pc = ax.pcolor(xqy/Ly,yqy/Ly,(qy_python-qy_matlab)/q_char,cmap='jet')
ax.set_aspect('equal','box')
plt.locator_params(axis='both', nbins=3)
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.05)
cb = plt.colorbar(pc,cax=cax,format='%.0e')
cb.ax.locator_params(nbins=3)

# plot Pf
ax = fig.add_subplot(4,3,7)
pc = ax.pcolor(x/Ly,y/Ly,Pf_python/Pr_char,cmap='jet')
ax.set_aspect('equal','box')
plt.locator_params(axis='both', nbins=3)
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.05)
cb = plt.colorbar(pc,cax=cax,format='%.0e')
cb.ax.locator_params(nbins=3)
ax.set_ylabel(r'$\dfrac{y}{L_y}$')
textstr = r'$\dfrac{P_f}{P_c}$'
#textstr = r'$\dfrac{P_f}{P_c}  = P_f \cdot \dfrac{\rho c_p}{\lambda} \cdot \dfrac{k}{\eta_f}$'
ax.text(0.025, 0.95, textstr, transform=ax.transAxes, fontsize=14, verticalalignment='top', 
    bbox=dict(boxstyle='square', facecolor='white', alpha=0.5))

ax = fig.add_subplot(4,3,8)
pc = ax.pcolor(x/Ly,y/Ly,Pf_matlab/Pr_char,cmap='jet')
ax.set_aspect('equal','box')
plt.locator_params(axis='both', nbins=3)
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.05)
cb = plt.colorbar(pc,cax=cax,format='%.0e')
cb.ax.locator_params(nbins=3)

ax = fig.add_subplot(4,3,9)
pc = ax.pcolor(x/Ly,y/Ly,(Pf_python-Pf_matlab)/Pr_char,cmap='jet')
ax.set_aspect('equal','box')
plt.locator_params(axis='both', nbins=3)
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.05)
cb = plt.colorbar(pc,cax=cax,format='%.0e')
cb.ax.locator_params(nbins=3)

# plot T
ax = fig.add_subplot(4,3,10)
pc = ax.pcolor(x/Ly,y/Ly,T_python/deltaT,cmap='jet')
ax.set_aspect('equal','box')
plt.locator_params(axis='both', nbins=3)
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.05)
cb = plt.colorbar(pc,cax=cax,format='%.0e')
cb.ax.locator_params(nbins=3)
ax.set_xlabel(r'$\dfrac{x}{L_y}$')
ax.set_ylabel(r'$\dfrac{y}{L_y}$')
textstr = r'$\dfrac{T}{T_c}$'
#textstr = r'$\dfrac{T}{T_c}  = \dfrac{T}{\Delta T}$'
ax.text(0.025, 0.95, textstr, transform=ax.transAxes, fontsize=14, verticalalignment='top', 
    bbox=dict(boxstyle='square', facecolor='white', alpha=0.5))

ax = fig.add_subplot(4,3,11)
pc = ax.pcolor(x/Ly,y/Ly,T_matlab/deltaT,cmap='jet')
ax.set_aspect('equal','box')
plt.locator_params(axis='both', nbins=3)
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.05)
cb = plt.colorbar(pc,cax=cax,format='%.0e')
cb.ax.locator_params(nbins=3)
ax.set_xlabel(r'$\dfrac{x}{L_y}$')

ax = fig.add_subplot(4,3,12)
pc = ax.pcolor(x/Ly,y/Ly,(T_python-T_matlab)/deltaT,cmap='jet')
ax.set_aspect('equal','box')
plt.locator_params(axis='both', nbins=3)
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.05)
cb = plt.colorbar(pc,cax=cax,format='%.0e')
cb.ax.locator_params(nbins=3)
ax.set_xlabel(r'$\dfrac{x}{L_y}$')

#fig.tight_layout()
fig.suptitle("Time step, it = "+ str(it))
plt.draw()
plt.show()
