clear, figure(1), clf,colormap(jet)
load('mnt2_sample25.mat'); % mnt2
B = mnt2 - min(mnt2(:)); B = B/max(B(:));
%physics
Lx     = 2;    % length of the computational domain in x
Ly     = 2;    % length of the computational domain in y
g      = 9.8;  % gravity
mu_f   = 1;    % viscosity
Ke     = 10e-3;% coefficient
%numerics
nx     = 181;  % number of grid points in x 
ny     = 184;  % number of grid points in y 
nt     = 5000; % number of time stepd
%preprocessing
dx     = Lx/(nx-1); % grid spacing in x 
dy     = Ly/(ny-1); % grid spacing in y 
x      = 0:dx:Lx;   % coordinates of grid points
y      = 0:dy:Ly;   % coordinates of grid points
dt     = 4e-5;      % time step
%initial
[x2d y2d] = ndgrid(x,y);
H         =  0.5*exp( -((x2d-0.5)-0.5*Lx).^2/0.1 - ((y2d-0.5)-0.5*Ly).^2/0.1)...
           + 0.5*exp( -(x2d-0.5*Lx).^2/0.1 - (y2d-0.5*Ly).^2/0.1);
surf(x2d,y2d,B),shading flat,light,view(123,32),title('B(x,y)')
surf(x2d,y2d,H),shading flat,light,view(123,32),title('H')
surf(x2d,y2d,H+B),shading flat,light,view(123,32),title('S')
%action
for it=1:nt
Hx     = (H(1:end-1,:)+H(2:end,:))/2;
Hy     = (H(:,1:end-1)+H(:,2:end))/2;
qx     =  -g/(3*mu_f)* Hx.^3.*diff(H + B,1,1)/dx;
qy     =  -g/(3*mu_f)* Hy.^3.*diff(H + B,1,2)/dy;
qxa    = (qx(:,1:end-1)+qx(:,2:end))/2;
qya    = (qy(1:end-1,:)+qy(2:end,:))/2;
Ha     = (Hx(:,1:end-1)+Hx(:,2:end))/2;
Vel    = sqrt(qxa.^2 + qya.^2)./Ha; 
dBdt   = Ke*Vel;
B(2:end,2:end) = B(2:end,2:end) + dBdt*dt;  
dHdt   =  -diff(qx(:,2:end-1),1,1)/dx-diff(qy(2:end-1,:),1,2)/dy;
H(2:end-1,2:end-1)= H(2:end-1,2:end-1)+dt*dHdt;
    if mod(it,100)==0
    surf(x ,y ,(B+H)');
    view(123,32),axis image, colorbar
    shading interp,light,lighting phong
    xlabel('distance in x direction');
    ylabel('distance in y direction');
    zlabel('distance in z direction');
    title(it),drawnow
    end
end