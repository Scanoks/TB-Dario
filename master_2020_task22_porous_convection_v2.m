clear variables,close all,clf,colormap(jet),delete step*.mat
%physics
% dim independent
Ly          = 1e3;      % Length in y, m
lam_rhoCp   = 1e-6;     % heat diffusivity, m^2/s  
k_etaf      = 1e-12;    % permeability(m^2)/ fluid viscosity (Pa*s) = m^2/(Pa*s) 
deltaT      = 100;      % temperature difference, K 
% usefull scales
time_char   = Ly^2/lam_rhoCp;   % m^2/(m^2/s) = s
Pr_char     = lam_rhoCp/k_etaf; % (m^2/s)/ ( m^2/(Pa*s)) = Pa
% nondim
Lx_Ly       = 2;
w_Ly        = 0.1;
Ra          = 1e6;   % Ra = alphrhofg*deltaT*k_etaf*Ly/ lam_rhoCp
tt_nondim   = 1e-2;
% dim dependent
Lx          = Lx_Ly*Ly; % Length in x, m
alphrhofg   = Ra*lam_rhoCp/k_etaf/deltaT/Ly; %   Pa*s/m^2/K*m^2/s/m = Pa/m/K  % thermal expansion*density*gravity, Pa/m/K
tt          = tt_nondim*time_char;      % total time, s
w           = w_Ly*Ly;   
%numerics
dt          = 3e-9*time_char;% time step, s
betaf       = 1e-4/Pr_char;  % fluid compressibility, 1/Pa
nx          = 150;           % number of grid points in x 
ny          = fix(nx*Ly/Lx); % number of grid points in y 
nout        = 100;         	 % plot every nout time step
st          = 5;             % plot every st velocity
niter       = 1e7;
tol         = 1e-6;
output      = 1;             % write output files every nout time step
%preprocessing
nt          = fix(tt/dt);
dx          = Lx/(nx-1);
dy          = Ly/(ny-1);
[x,y]       = ndgrid(-Lx/2:dx:Lx/2,-Ly/2:dy:Ly/2);
%initial conditions
time        = 0;        % s
T           = deltaT*exp( (-x.^2 -(y+ max(y(:))/2).^2 )/w^2 );
%T           = deltaT*(rand(nx,ny)-0.5);
Pf          = 0*T;
qx          = zeros(nx+1, ny);
qy          = zeros(nx, ny+1);
%action
if output==1
    if ~exist('output_M','dir'); mkdir('output_M'); end
end
tic
for it = 0:nt
    rho_g    = -alphrhofg*T;  % buoyancy
    % Hydro 
    max_delta_Pf = [];
    for iter = 1:niter
        qx(2:end-1,2:end-1) = - k_etaf*(diff(Pf(:,2:end-1),1,1)/dx );
        qy(2:end-1,2:end-1) = - k_etaf*(diff(Pf(2:end-1,:),1,2)/dy ...
            + (rho_g(2:end-1,1:end-1) + rho_g(2:end-1,2:end))/2 );
        % betaf*dPfdt = 0 = -( diff(qx,1,1)/dx + diff(qy,1,2)/dy )
        dPfdt               = -( diff(qx,1,1)/dx + diff(qy,1,2)/dy )/betaf;
        deltaPf             = dt*dPfdt;
        Pf                  = Pf + deltaPf;
        max_delta_Pf(iter)  = max(abs(deltaPf(:)));
        if max_delta_Pf(iter)/max(abs(Pf(:))) < tol,break,end
    end
    % Thermo
    % diffusion
    dTdt = diff( lam_rhoCp*diff(T(:,2:end-1),1,1)/dx ,1,1)/dx ...
        +  diff( lam_rhoCp*diff(T(2:end-1,:),1,2)/dy ,1,2)/dy;
    T(2:end-1,2:end-1)= T(2:end-1,2:end-1) + dTdt*dt;
    % advection    dTdt = - Vx*dTdx - Vy*dTdy
    T(:,2:end  ) = T(:,2:end  ) - dt*max(0,qy(:,2:end-1)).*diff(T,1,2)/dy; % up
    T(:,1:end-1) = T(:,1:end-1) - dt*min(0,qy(:,2:end-1)).*diff(T,1,2)/dy; % down        
    T(2:end  ,:) = T(2:end  ,:) - dt*max(0,qx(2:end-1,:)).*diff(T,1,1)/dx; % right
    T(1:end-1,:) = T(1:end-1,:) - dt*min(0,qx(2:end-1,:)).*diff(T,1,1)/dx; % left    
    % boundary conditions
    T(:,  1) = deltaT/2; % heating at the base by deltaT/2
    T(:,end) =-deltaT/2; % cooling at the top by -deltaT/2    
    T(1,:)   = T(2,:);  
    T(end,:) = T(end-1,:);
    time         = time + dt;
    if mod(it,nout)==0%postprocessing
      %  subplot(122),
        pcolor(x/Ly,y/Ly,T/deltaT),shading flat
        colorbar,title(['Temperature at time: ',num2str(it*dt/time_char,2)]),axis image
        hold on,quiver(x(1:st:end,1:st:end)/Ly,y(1:st:end,1:st:end)/Ly ...
            ,qx(2:st:end,1:st:end),qy(1:st:end,2:st:end),'k'),hold off         
        %caxis([-1 1]/10)
       % subplot(121),loglog(max_delta_Pf)
       xlabel('X'), ylabel('Y')
        drawnow
%         save(['step',int2str(it/nout)])
        if output==1% save variables
            writematrix(qx,['output_M/output_it_',num2str(it),'_qx_M.csv']) 
            writematrix(qy,['output_M/output_it_',num2str(it),'_qy_M.csv']) 
            writematrix(Pf,['output_M/output_it_',num2str(it),'_Pf_M.csv']) 
            writematrix(T,['output_M/output_it_',num2str(it),'_T_M.csv']) 
        end
    end
end
cpus = toc