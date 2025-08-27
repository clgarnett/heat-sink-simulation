tic

xmin = 0; xmax = 2; ymin = 0; ymax = 3; tmax = 1.5;             % Heat sink dimensions.
S_heat = 25; S_x_min = 0.5; S_x_max = 1.5; S_y_min = 1.75; S_y_max = 2.25;  % Source heat and dimensions.
temperature_min = 0; temperature_max = 3;                         % Color scale range.
Nx = 40; Ny = 60; Nt = 15;         % Number of points in the horizontal, vertical and time discretizations. 

Dx = xmax/(Nx - 1); Dy = ymax/(Ny - 1);

% Our matrices Ax and Ay
Ax = zeros(Nx, Nx);
Ax(1:Nx-1, 2:Nx) = Ax(1:Nx-1, 2:Nx) + diag(ones(Nx-1, 1));
Ax(1:Nx, 1:Nx) = Ax(1:Nx, 1:Nx) - 2*diag(ones(Nx, 1));
Ax(2:Nx, 1:Nx-1) = Ax(2:Nx, 1:Nx-1) + diag(ones(Nx-1, 1));
Ax = Ax / (Dx)^2;

Ay = zeros(Ny, Ny);
Ay(1:Ny-1, 2:Ny) = Ay(1:Ny-1, 2:Ny) +diag(ones(Ny-1, 1));
Ay(1:Ny, 1:Ny) = Ay(1:Ny, 1:Ny) - 2*diag(ones(Ny, 1));
Ay(2:Ny, 1:Ny-1) = Ay(2:Ny, 1:Ny-1) +diag(ones(Ny-1, 1));
Ay = Ay / (Dy)^2;

dltx = (xmax-xmin) / (Nx-1); dlty = (ymax-ymin) / (Ny-1);
[x, y] = ndgrid(xmin + (0:Nx-1)*dltx, ymin + (0:Ny-1)*dlty);

% Source term
S = S_heat * (x > S_x_min) .* (x < S_x_max) .* (y > S_y_min) .* (y < S_y_max);  % Source dimensions
figure; pcolor(x,y,S); xlabel('x'); ylabel('y'); title('Source')
colorbar; shading interp; colormap hot; axis equal; axis tight

U = zeros(Ny,Nx); % Initial solution
dltt = tmax/Nt;   % Time step
f = figure(2); hold on

for i = 1:Nt
    t = i*dltt;

    % Find the matrix U at the next time point:     (A*dt - I)U + U(B*dt) = -dt*S - U_{-1}
    U_next = lyap(dltt*Ay-eye(Ny), Ax*dltt, U + dltt*rot90(S));
    U = U_next
    

    % Plot solution
    figure(2)
    pcolor(x,y,rot90(U, -1)); shading interp; colormap hot; colorbar;
    xlabel('x'); ylabel('y'); title(sprintf('t = %g', t))
    clim([temperature_min temperature_max]) % Fix the colors
    axis equal; axis tight
    M(i) = getframe(f); % Store the figure to make movie
end

toc

figure(2)
movie(f, M, 1) 