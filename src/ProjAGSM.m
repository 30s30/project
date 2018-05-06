%% Bryan Soto UH:1284525
%% Final Project SC- Gauss-seidel method
clc; 
clear all; 

% Set the variables 

nx = 10; % nodes at the x axis; all dirichlet boundaries
ny = nx; % nodes at y axis; bottom is neumann boundary so, ghost nodes 

tol = 1e-9; %set the tolerance for gauss approx. 
err = 1000 ; % initialize error to check for tolerance every iteration 
itr = 0; %initialize gaussian iteration 

% define parameters; 
ax= -pi;ay= -pi;bx= pi;by= pi;
lambda= pi; 
Lx = 2*pi; % (-pi,pi) 
Ly = 2*pi; % (-pi,pi) 

dx = linspace(ax,bx,nx);
dy = linspace(ay,by,ny)'; 

hx = dx(2)-dx(1); hy=dy(2)-dy(1);

[x2,y2] = meshgrid(dx,dy);

y2 = flipud(y2); % flip the grid to make it look like my grid

% create the grid for interior points and the boundary conditions 
u = zeros(ny,nx);
u_z=u; 



%% Dirichlet Boundary Conditions 

% left boundary (D)
l_b = y2.*(y2-ay).^2;
% right boundary (D)
r_b = ((y2-ay).^2).*cos(pi*y2/ay);
% upper boundary (D) 
up_b = by*(by-ay)^2+((x2-ax)/(bx-ax))*((((by-ay)^2)*cos(pi*by/ay))-(by*(by-ay)^2));


% add boundaries to matrix 
u_z(2:ny-1,1) = l_b(2:ny-1,1); 
u_z(2:ny-1,nx) =r_b(2:ny-1,nx); 
u_z(1,2:nx-1) = up_b(1,2:nx-1); 
% fix corner boundaries where they are an addition of each other
u_z(1,1) = u_z(1,2)+u_z(2,1);  %upper left corner
u_z(ny,1)= u_z(ny-1,1) + u_z(ny,2); %lower left corner 
u_z(ny,nx) = u_z(ny-1,nx) + u_z(ny,nx-1); %lower right
u_z(1,nx) = u_z(2,nx)+u_z(1,nx-1); %upper right 

% create f matrix of inner values 

f = cos((pi/2)*(2*((x2-ax)/(bx-ax))+1)).*sin(pi*(y2-ay)/(by-ay));

% Fix bondary conditions 

dudy = 0;   %Neumann boundary value
u_0=u_z; 
%% Gauss-Seidel 

while err>tol

    
 for j=2:ny-1 %for loop for the y variable
        for i=2:nx-1 %for loop for the x variable
     
         
         if i == nx-1   % Neumann condition  
         u_z(ny,j) = (1/(4-(lambda*hx^2)))*(2*u_z(ny-1,2)+u_z(ny,j-1)+u_z(ny,j+1)-(hx^2)*f(ny,j)-2*hy*dudy);
         
         end
         
         u_z(i,j)=(1/(4-(lambda*hx^2)))*(u_0(i-1,j)+u_z(i+1,j)+u_0(i,j-1)+u_z(i,j+1)-f(i,j)*(hx^2)); %numerical algorithm for the inner nodes at k+1 using the outside nodes plus the source
             
        end
 end
 
err = max(max(abs(u_z(2:nx-1,2:ny-1)-u_0(2:nx-1,2:ny-1))));
itr =itr+1;
  
     
X(:,itr) = max(err);

end 











