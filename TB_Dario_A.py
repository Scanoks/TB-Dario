import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import os

# Physics ¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦
# Dimension-independent variables ==========================================================================
Ly          = 1e3                                   # Length in y, m
lam_rhoCp   = 1e-6                                  # Heat diffusivity, m^2/s
k_etaf      = 1e-12                                 # Permeability(m^2)/fluid viscosity (Pa*s) = m^2/(Pa*s)
k_etafmin   = 1e-12                                 # Perméabilité si T < 360C°
k_etafmax   = 1e-18                                 # Perméabilité si T > 500C°
deltaT      = 100                                   # Temperature difference, K

# Useful scales ============================================================================================
time_char   = Ly**2/lam_rhoCp                       # m^2/(m^2/s) = s ## ** = ^ sur matlab
Pr_char     = lam_rhoCp/k_etaf                      # (m^2/s)/(m^2/(Pa*s)) = Pa

# nondim ===================================================================================================
Lx_Ly       = 2                                     # pour définir après que Lx est 2* plus grand que Ly
w_Ly        = 0.1                                   # permet de gérer la taille de T
Ra          = 1e6                                   # Ra = alphrhofg*deltaT*k_etaf*Ly/ lam_rhoCp
tt_nondim   = 1e-2                                  # permet de définir mon tt

# Dimension-dependent variables =============================================================================
Lx      = Lx_Ly * Ly                                # Length in x, m
alphrhofg = Ra * lam_rhoCp / k_etaf / deltaT / Ly   # Pa*s/m^2/K*m^2/s/m = Pa/m/K  % thermal expansion*density*gravity, Pa/m/K
tt      = tt_nondim * time_char                     # Total time, s
w       = w_Ly * Ly                                 # Défini mon w pour gérer la taille de mon T 

# Numerical variables =======================================================================================
dt      = 2*3e-9 * time_char                   # Time step, s
betaf   = 1e-4 / Pr_char                            # Fluid compressibility, 1/Pa
nx      = 110                                       # Number of grid points in x 
ny      = int(nx * Ly / Lx)                         # Number of grid points in y  ## fix = int sur python
nout    = 10                                       # Plot every nout time step
st      = 3                                         # Plot every st velocity, nombre de flèche
niter   = 1000000                                         # défini la taille de ma boucle inter
tol     = 1e-6                                      # tolérence pour break ma boucle inter
output  = 1                                         # write output files every nout time step

# Prétraitement ==============================================================================================
nt      = 3e6                                       # pas de temps total, int et // permet d'arrondir à l'entier le plus proche
dx      = Lx/(nx-1)                                 # défini la longueur de dx
dy      = Ly/(ny-1)                                 # défini la longueur de dy
x, y    = np.meshgrid(np.arange(-Lx/2, Lx/2 + dx, dx), np.arange(-Ly/2, Ly/2 + dy, dy), indexing='ij')
    # np.meshgrid va créer une grille pour x et y et np.arrange créer le tableau avec des valeurs allant de -Lx/2 à Lx/2 avec un pas de temps de taille dx idem y, indexing='ij' permet d'avoir un système de coordonnées cartésiennes car sinon par defaut on a 'xy' qui correspond à un systeme de vecteur

# Initialisation des T ========================================================================================
# Conditions initiales
time    = 0                                         # s
phi     = 90/180*np.pi                              # angle de rotation
x0      = 0                                         # défini position en x
y0      = -np.max(y)/2                              # défini position en y
a       = Lx/20                                     # défini la déformation de la longueur 
b       = a/1                                       # défini la déformation de la largeur
T       = -deltaT/2+(y-Ly/2)/Ly*-deltaT
T      += deltaT * np.exp(-((x-x0)*np.cos(phi) + (y-y0)*np.sin(phi))**2/a**2 - (-(y-y0)*np.cos(phi) + (x-x0)*np.sin(phi))**2/b**2) # gaussien en forme d'une ellipse tournée // formule dans matlab #T =  (deltaT * np.exp((-x**2/wn_x - (y/wn_y + np.max(y/wn_y)/2)**2)/w**2))      # np.exp pour exponentielle ajouter T+= ........ pour ajouter des forme a ma T

T[:, 0]     = deltaT / 2    # heating at the base by deltaT/2
T[:, -1]    = -deltaT / 2

Pf      = 0 * T                                     # défini Pf comme une matric de 0 de taille T()

# q normal ========================================================================================================
qx      = np.zeros((nx+1, ny))                      # création d'un matrice de 0 de la taille de nos nx et ny, pour qx on demarre à la deuxiemme ligne
qy      = np.zeros((nx, ny+1))                      # idem mais ici on enlève la derniere colonne   

# Action ¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦
arr      = np.arange(0, nt+1, dtype
                     =int)             # défini la range de ma boucle # cet écriture permet de corriger une erreur dans la prise en compte de range dans python contrairement à matlab
plt.ion()                                           # sets interactive mode on (needed for plot update)
fig, axs = plt.subplots(3, 2 )

if output==1:
    abspath = os.path.abspath(__file__) # https://stackoverflow.com/questions/1432924/python-change-the-scripts-working-directory-to-the-scripts-own-directory
    dname   = os.path.dirname(abspath)
    os.chdir(dname)
    if not os.path.exists('output_P'):
        os.makedirs('output_P')
for it in arr:
    
    # Hydro ¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦
    rho_g           = -alphrhofg*T      # flottabilité
    max_delta_Pf    = []
    
    for inter in range(1,niter):
        
        # q normal =================================================================================================
        qx[1:-1, 1:-1] = -k_etaf*(np.diff(Pf[:, 1:-1], axis=0)/dx)  # dans python on commence avec l'indice 0 contrairement à matlab il faut donc adapter ça en indiquant 1:-1 pour enlever la premiere et derniere valeur
        qy[1:-1, 1:-1] = -k_etaf*(np.diff(Pf[1:-1, :], axis=1)/dy 
        + (rho_g[1:-1, :-1] + rho_g[1:-1, 1:])/2)
        # qx dépend de la pression notamment la différence et aussi influencé par la permeabilité pour qy il y a aussi l influence de la gravité

        #Pf normal =================================================================================================
        dPfdt   = -(np.diff(qx, axis=0)/dx + np.diff(qy, axis=1)/dy)/betaf                      # calcul de la pression en fonction du mouvement des fluide et du temps
        deltaPf = dt*dPfdt
        Pf      = Pf + deltaPf 
        max_delta_Pf.append(np.max(np.abs(deltaPf)))

        if max_delta_Pf[len(max_delta_Pf)-1]/np.max(np.abs(Pf)) < tol:
            break

    # Thermo ¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦
    # diffusion normal =============================================================================================
    dTdt = np.diff(lam_rhoCp * np.diff(T[:, 1:-1], axis=0) / dx, axis=0) / dx \
    + np.diff(lam_rhoCp * np.diff(T[1:-1, :], axis=1) / dy, axis=1) / dy
    T[1:-1, 1:-1] += dTdt * dt

    # advection normal   dTdt = - Vx*dTdx - Vy*dTdy ================================================================
    T[:, 1:] -= dt * np.maximum(0, qy[:,1:-1]) * np.diff(T, axis=1) / dy   # up
    T[:,:-1] -= dt * np.minimum(0, qy[:,1:-1]) * np.diff(T, axis=1) / dy    # down
    T[1:,:]  -= dt * np.maximum(0, qx[1:-1,:]) * np.diff(T, axis=0) / dx    # right
    T[:-1,:] -= dt * np.minimum(0, qx[1:-1,:]) * np.diff(T, axis=0) / dx    # left

    # boundary conditions normal ===================================================================================
    T[:, 0]     = deltaT / 2    # heating at the base by deltaT/2
    T[:, -1]    = -deltaT / 2   # cooling at the top by -deltaT/2
    T[-1, :]     = T[-2, :]
    T[0, :]     = T[1, :]
    time        += dt
    
    # Visualisation ¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦
    if np.mod(it, nout) == 0:   # postprocessing
        print(it)
        plt.clf()

    # quiver =======================================================================================================
        x_quiver    = x[1:-1:st, 1:-1:st]/Ly
        y_quiver    = y[1:-1:st, 1:-1:st]/Ly
        qx_quiver   = qx[2:-1:st, 1:-1:st]
        qy_quiver   = qy[1:-1:st, 2:-1:st]

        norm = np.sqrt(qx_quiver**2 + qy_quiver**2)
        qx_quivers = qx_quiver/norm
        qy_quivers = qy_quiver/norm

    # plot T normal ================================================================================================
        plt.clf()
        ax = fig.add_subplot(1,1,1)
        im=ax.pcolormesh(x/Ly, y/Ly, T/deltaT,cmap='jet',shading='gouraud')
        ax.quiver(x_quiver, y_quiver, qx_quivers, qy_quivers)
        ax.set_title('it = ' + str(int(it)) + '   temps = ' + str(int(dt/3600*it/24)) + ' jours ')
        ax.set_aspect('equal', 'box')
        plt.xlabel('X')
        plt.ylabel('Y')
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.05)

        plt.colorbar(im, cax=cax)
        
        plt.pause(0.1)
        plt.draw()
        plt.show()
        
        it_long = str(it*10).zfill(7)
        plt.savefig('./A_img/A_img_'+it_long+'.png')

# Visualisaiton finale =============================================================================================            

plt.ioff()
print('=== please close the figure window to continue ===')
plt.show() # keeps last frame blocked, so it doesn't disappear - only works when 'ioff()'
print('\n=== finished :) ===\n')

