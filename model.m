clear, figure(1), clf,colormap(jet)
load('mnt2_sample25.mat'); % mnt2
B = mnt2 - min(mnt2(:)); B = B/max(B(:));
%load('snow_mountain','snow_mountain');
%colormap(snow_mountain)
%physics
Lx     = 2;  % length of the computational domain in x
Ly     = 2;  % length of the computational domain in y
g      = 9.8;
mu_f   = 1;
Ke     = 10e-3;
%numerics
nx     = 181;  % number of grid points in x 
ny     = 184;  % number of grid points in y 
nt     = 5000; % number of time stepd
%preprocessing
dx     = Lx/(nx-1);     % grid spacing in x 
dy     = Ly/(ny-1);     % grid spacing in y 
x      = 0:dx:Lx; % coordinates of grid points
y      = 0:dy:Ly; % coordinates of grid points
dt     = 4e-5;          % time step
%initial
[x2d y2d] = ndgrid(x,y);
%Hg       =  0.5*exp( -(x2d-0.5*Lx).^2/0.1 - (y2d-0.5*Ly).^2/0.1);
Hg       =  0.5*exp( -((x2d-0.5)-0.5*Lx).^2/0.1 - ((y2d-0.5)-0.5*Ly).^2/0.1)...
           + 0.5*exp( -(x2d-0.5*Lx).^2/0.1 - (y2d-0.5*Ly).^2/0.1);
surf(x2d,y2d,B),shading flat,light,view(123,32),title('B(x,y)')
surf(x2d,y2d,Hg),shading flat,light,view(123,32),title('Hg')
surf(x2d,y2d,Hg+B),shading flat,light,view(123,32),title('S')
B0 = B;
%action
for it=1:nt
Hgx     = (Hg(1:end-1,:)+Hg(2:end,:))/2;
Hgy     = (Hg(:,1:end-1)+Hg(:,2:end))/2;
qgx     =  -g/(3*mu_f)* Hgx.^3.*diff(Hg + B,1,1)/dx;
qgy     =  -g/(3*mu_f)* Hgy.^3.*diff(Hg + B,1,2)/dy;
qxa     = (qgx(:,1:end-1)+qgx(:,2:end))/2;
qya     = (qgy(1:end-1,:)+qgy(2:end,:))/2;
Hga     = (Hgx(:,1:end-1)+Hgx(:,2:end))/2;
Vel     = sqrt(qxa.^2 + qya.^2)./Hga; 
dBdt    = Ke*Vel;
B(2:end,2:end) = B(2:end,2:end) + dBdt*dt;  
contourf(B-B0);
dHgdt   =  -diff(qgx(:,2:end-1),1,1)/dx-diff(qgy(2:end-1,:),1,2)/dy;
Hg(2:end-1,2:end-1)= Hg(2:end-1,2:end-1)+dt*dHgdt;
    if mod(it,100)==0
    surf(x ,y ,(B+Hg)');
    view(123,32),axis image
    shading interp,light,lighting phong
    xlabel('distance in x direction');
    ylabel('distance in y direction');
    zlabel('distance in z direction');
    title(it),drawnow
    end
end