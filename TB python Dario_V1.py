import numpy as np
import matplotlib.pyplot as plt
# Physics
# Dimension-independent variables
Ly = 1e3            # Length in y, m
lam_rhoCp = 1e-6    # Heat diffusivity, m^2/s
k_etaf = 1e-12      # Permeability(m^2)/fluid viscosity (Pa*s) = m^2/(Pa*s)
deltaT = 100        # Temperature difference, K
# Useful scales
time_char = Ly**2/lam_rhoCp # m^2/(m^2/s) = s ## ** = ^ sur matlab
Pr_char = lam_rhoCp/k_etaf  # (m^2/s)/(m^2/(Pa*s)) = Pa
# nondim
Lx_Ly       = 2
w_Ly        = 0.1
Ra          = 1e4           # Ra = alphrhofg*deltaT*k_etaf*Ly/ lam_rhoCp
tt_nondim   = 1e-2
# Dimension-dependent variables
Lx = Lx_Ly * Ly             # Length in x, m
alphrhofg = Ra * lam_rhoCp / k_etaf / deltaT / Ly   # Pa*s/m^2/K*m^2/s/m = Pa/m/K  % thermal expansion*density*gravity, Pa/m/K
tt = tt_nondim * time_char  # Total time, s
w = w_Ly * Ly
# Numerical variables
dt = 3e-9 * time_char       # Time step, s
betaf = 1e-4 / Pr_char      # Fluid compressibility, 1/Pa
nx = 150                    # Number of grid points in x 
ny = int(nx * Ly / Lx)      # Number of grid points in y  ## fix = int sur python
nout = 100                  # Plot every nout time step
st = 5               # Plot every st velocity
niter = 1000
#print('niter = ',niter)
tol = 1e-6
# Prétraitement
nt = int(tt // dt)         # int et // permet d'arrondir à l'entier le plus proche
#nt = np.fix(tt/dt)
#nt = 1e6
print(nt)
dx = Lx/(nx-1)              # défini la longueur de dx
dy = Ly/(ny-1)
#, y = np.meshgrid(np.arange(-Lx/2, Lx/2, dx), np.arange(-Ly/2, Ly/2, dy), indexing='ij')
x, y = np.meshgrid(np.arange(-Lx/2, Lx/2 + dx, dx), np.arange(-Ly/2, Ly/2 + dy, dy), indexing='ij')
                           # np.meshgrid va créer une grille pour x et y et np.arrange créer le tableau avec des valeurs allant de -Lx/2 à Lx/2 avec un pas de temps de taille dx idem y
# Conditions initiales
time = 0                    # s
T = (deltaT * np.exp((-x**2 - (y + np.max(y)/2)**2)/w**2)) # np.exp pour exponentielle
#T = deltaT * (np.random.rand(nx, ny) - 0.5) ## crée des valeur aléatoir entre -0.5 et 0.5 pour la matrice de T
Pf = 0 * T
qx = np.zeros((nx+1, ny))   # création d'un matrice de 0 de la taille de nos nx et ny, pour qx on demarre à la deuxiemme ligne
qy = np.zeros((nx, ny+1))   # idem mais ici on enlève la derniere colonne
# Action
#tic = time.time() # formnule équivalente à tic dans matlab
arr = np.arange(0, nt+1)
#print(arr)
#arr = arr.astype(int)
#print(range(arr))
# arr permet de corriger une erreur dans la prise en compte de range

plt.ion() # sets interactive mode on (needed for plot update)
fig = plt.figure()
ax = fig.add_subplot(111)
#data = (T/deltaT)
#print('data shape = ',data.shape)
#im = ax.imshow(data, origin='lower')
plt.draw()
#plt.colorbar()
#t.plot(y)
#plt.pause(0.001)
plt.show()

for it in arr:
    rho_g = -alphrhofg*T  # flottabilité
    # Hydro
    max_delta_Pf = []
    for inter in range(1,niter):
        qx[1:-1, 1:-1] = -k_etaf*(np.diff(Pf[:, 1:-1], axis=0)/dx) # dans python on commence avec l'indice 0 contrairement à matlab il faut donc adapter ça en indiquant 1:-1 pour enlever la premiere et derniere valeur
        qy[1:-1, 1:-1] = -k_etaf*(np.diff(Pf[1:-1, :], axis=1)/dy 
            + (rho_g[1:-1, :-1] + rho_g[1:-1, 1:])/2)
        # betaf*dPfdt = 0 = -( diff(qx,1,0)/dx + diff(qy,1,1)/dy )
        dPfdt = -(np.diff(qx, axis=0)/dx + np.diff(qy, axis=1)/dy)/betaf # axis=1 correspond à la deuxième ligne et axis=0 la première
        deltaPf = dt*dPfdt
        Pf = Pf + deltaPf
        max_delta_Pf.append(np.max(np.abs(deltaPf)))
        #if max_delta_Pf[-1]/np.max(np.abs(Pf)) < tol:               # arret lorsque qu'on dépasse la tolérence
        #    break
        #print(max_delta_Pf[len(max_delta_Pf)-1]/np.max(np.abs(Pf)))
        if max_delta_Pf[len(max_delta_Pf)-1]/np.max(np.abs(Pf)) < tol:
            #print('got out, broke')
            break
    #print('got out, it = ',it,' err = ',max_delta_Pf[len(max_delta_Pf)-1]/np.max(np.abs(Pf)))

    # Thermo
    # diffusion
    dTdt = np.diff(lam_rhoCp * np.diff(T[:, 1:-1], axis=0) / dx, axis=0) / dx \
    + np.diff(lam_rhoCp * np.diff(T[1:-1, :], axis=1) / dy, axis=1) / dy
    T[1:-1, 1:-1] += dTdt * dt

    # advection    dTdt = - Vx*dTdx - Vy*dTdy
    T[:, 1:] -= dt * np.maximum(0, qy[:, 1:-1]) * np.diff(T, axis=1) / dy  # up
    T[:,:-1] -= dt * np.minimum(0, qy[:,1:-1]) * np.diff(T, axis=1) / dy  # down
    T[1:,:] -= dt * np.maximum(0, qx[1:-1,:]) * np.diff(T, axis=0) / dx # right
    T[:-1,:] -= dt * np.minimum(0, qx[1:-1,:]) * np.diff(T, axis=0) / dx # left

    # boundary conditions
    T[:, 0] = deltaT / 2  # heating at the base by deltaT/2
    T[:, -1] = -deltaT / 2  # cooling at the top by -deltaT/2
    T[0, :] = T[1, :]
    #time += dt
    if np.mod(it, nout) == 0: # postprocessing
        print(it,'it')
        plt.clf()
        #data = np.transpose(T/deltaT)
        #im.set_data(data)
        plt.pcolor(x/Ly, y/Ly, T/deltaT)
        x_quiver = x[1:-1:st, 1:-1:st]/Ly
        y_quiver = y[1:-1:st, 1:-1:st]/Ly
        qx_quiver = qx[2:-1:st, 1:-1:st]
        qy_quiver = qy[1:-1:st, 2:-1:st]
        """print(np.shape(y/Ly), 'y/Ly')
        print(np.shape(x/Ly), 'x/Ly')
        print(np.shape(qx_quiver), 'qx_quiver')
        print(np.shape(T/deltaT), 'T/deltaT')
        print(np.shape(qx_quiver))
        print(np.shape(qx_quiver))
        print(np.shape(x_quiver))
        print(np.shape(y_quiver))"""
        plt.quiver(x_quiver, y_quiver, qx_quiver, qy_quiver, color='k')
      # plt.quiver(x_quiver, y_quiver)
      # plt.quiver(x[1:-1:st, 1:-1:st]/Ly, y[1:-1:st, 1:-1:st]/Ly,#qy[1:-1:st, 2::st], qx[2::st, 1:-1:st])
      # plt.quiver(x[1:-1:st, 1:-1:st]/Ly, y[1:-1:st, 1:-1:st]/Ly)
      # plt.quiver(y[1:-1:st, 1:-1:st]/Lx, x[1:-1:st, 1:-1:st]/Lx)
      # Tracer le champ de vecteurs avec les axes inversés
      # q = ax.quiver(y[1:-1:st, 1:-1:st]/Lx, x[1:-1:st, 1:-1:st]/Lx)
        
        plt.title(it)
        plt.xlabel('X')
        plt.ylabel('Y')
        plt.title(f'Temperature at time: {it}')
        plt.draw()
        plt.colorbar()
        plt.pause(0.001)
        plt.show()

        
"""         plt.figure(figsize=(10, 6))
        plt.pcolormesh(x/Ly, y/Ly, np.transpose(T/deltaT), shading='gouraud')
        plt.colorbar()  
        plt.axis('image')
        plt.quiver(x[1:-1:st, 1:-1:st]/Ly, y[1:-1:st, 1:-1:st]/Ly, qx[2::st, 1:-1:st], qy[1:-1:st, 2::st], color='black')
        plt.xlabel('X')
        plt.ylabel('Y')
        plt.draw()
        plt.pause(0.001)
        plt.show() """

plt.ioff()
print('=== please close the figure window to continue ===')
plt.show() # keeps last frame blocked, so it doesn't disappear - only works when 'ioff()'
print('\n=== finished :) ===\n')

    #if it % nout == 0: # postprocessing
    #    break
    # subplot(122),
    #    plt.pcolor(x/Ly,y/Ly,),T/deltaT
    #    plt.show()
    #plt.colorbar(), plt.title('Temperature at time: {:.2f}'.format(it*dt/time_char))
    #plt.axis('image')
    #plt.quiver(x[::st,::st]/Ly,y[::st,::st]/Ly,qx[1::st,::st],qy[::st,1::st],color='k')
    #plt.xlabel('X'), plt.ylabel('Y')
    #plt.draw()
    # pour plot ddans python on ajoute l'extension matplotlib.pyplot
    # save('step' + str(int(it/nout)))

#plt.show()
# formule pour faire le toc
#end_time = time.time()

#elapsed_time = end_time - start_time

#print(f"Elapsed time: {elapsed_time:.2f} seconds")
