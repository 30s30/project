%% Scientific Computing Final Project
%% Helmholtz Equation-AcH
%% Bryan Soto Salinas|UH ID: 12781


clc; 
clear all; 

% Set the variables 

nx = 5*2^1; % nodes at the x axis; all dirichlet boundaries
ny = nx; % nodes at y axis; bottom is neumann boundary so, ghost nodes 

tol = 1e-9; %set the tolerance for gauss approx. 
err = 1; % initialize error to check for tolerance every iteration 

% Define parameters; 
ax= -pi;ay= -pi;bx= pi;by= pi;
 
Lx = bx-ax; % (-pi,pi) 
Ly = by-ay; % (-pi,pi) 

dx = linspace(ax,bx,nx);
dy = linspace(ay,by,ny)'; 

% since same x-y grid length and assume equal spacing dx and # of nodes. 

h = Lx/(nx+1);
hx=h; 
hy=h;

[x2,y2] = meshgrid(dx,dy); % create 2 x

y2 = flipud(y2); % flip the grid to make it look like my grid

% create the grid for interior points and the boundary conditions 
u = zeros(ny,nx);

%% Dirichlet Boundary Conditions 


%Vectorizing the boundary values 

% left boundary (D)
l_b = y2.*(y2-ay).^2;
% right boundary (D)
r_b = ((y2-ay).^2).*cos(pi*y2/ay);
% upper boundary (D) 
up_b = by*(by-ay)^2+((x2-ax)/(bx-ax))*((((by-ay)^2)*cos(pi*by/ay))-(by*(by-ay)^2));


% add boundaries to main matrix U
u(2:ny-1,1) = l_b(2:ny-1,1); 
u(2:ny-1,nx) =r_b(2:ny-1,nx); 
u(1,2:nx-1) = up_b(1,2:nx-1); 
% fix corner boundaries where they are an addition of each other
u(1,1) = u(1,2)+u(2,1);  %upper left corner
u(ny,1)= u(ny-1,1) + u(ny,2); %lower left corner 
u(ny,nx) = u(ny-1,nx) + u(ny,nx-1); %lower right
u(1,nx) = u(2,nx)+u(1,nx-1); %upper right 

% create f matrix of inner values and add value for lambda depending on the

%%  Gauss-Seidel Method
t = input('Type number of case you want to solve: #1) [F=F(x,y);lambda=pi] OR  #2) [F=0;lambda=0] OR  #3) [F=F(x,y);lambda=0] \n#  '); 

if t ==1   
        f = cos((pi/2)*(2*((x2-ax)/(bx-ax))+1)).*sin(pi*(y2-ay)/(by-ay));
        lambda=pi;
elseif t == 2
       f = zeros(ny,nx);
       lambda=0;
elseif t == 3
        f = cos((pi/2)*(2*((x2-ax)/(bx-ax))+1)).*sin(pi*(y2-ay)/(by-ay));
        lambda=0;   
end
          
%The while loop that evaluates the only Neumann BC, the discretization and also deals with the corners 
%by computing the average of the two boxes that are adjacent to corner
c=1; itr = 0; % initialize iteration count
dudy=0; % Neumann condition boundary value
tic   %initialize timer
while  err > tol
u_p = u;    % Define u_k+1 value as u_p to be the future value
    for j = 2:nx-1
        for i = 2:ny-1   % Indexing is applied for better performance 
            if c == 1
            c= c +1 ;
             u(1,1)=(u(1,2)+u(2,1))/2;              %Fix the corner values by taking average
             u(1,nx)= (u(1,nx-1)+u(2,nx))/2;
             u(ny,1)= (u(ny-1,1)+u(ny,2))/2;
             u(ny,nx)= (u(ny,nx-1)+u(ny-1,nx))/2;
             
            end;
            
            u(i,j)=(1/(4-(lambda*hx^2)))*(u(i-1,j)+u_p(i+1,j)+u(i,j-1)+u_p(i,j+1)-f(i,j)*(hx^2)); %numerical algorithm for the inner nodes at k+1 using the outside nodes plus the source
        end
        
        
            u(ny,j) = (1/(4-(lambda*hx^2)))*(2*u(ny-1,j)+u(ny,j-1)+u(ny,j+1)-(hx^2)*f(ny,j)-2*hy*dudy); %Dicretized Neumann condition implemented to matrix
        
    end
   
  
   u(1,1)= (u(1,2)+u(2,1))/2;           %Fix corner values 
   u(1,nx)= (u(1,nx-1)+u(2,nx))/2;
   u(ny,1)= (u(ny-1,1)+u(ny,2))/2;
   u(ny,nx)= (u(ny,nx-1)+u(ny-1,nx))/2;
   E = u - u_p;
   err = sqrt(sum(sum(E).^2)); %L_2 form error
   itr = itr+1;
   Eeval(:,itr) = [itr,err];
   
end 
toc
u_avg = mean(mean(u(2:nx-1,2:ny-1)))
disp(['The number of time steps taken to converge using GS Method is: ',num2str(itr)]); 
Eeval= Eeval'; 

figure 
plot(Eeval(:,1),Eeval(:,2),'-r'), axis tight, xlabel(' # of iterations'),ylabel('Error'), title(['# iterations vs L2 Form error for a tolerance=',num2str(tol),', for N=',num2str(nx)])

figure 
surf(x2,y2,u),title(['3D plot for Gauss-Seidel Method solution',',  tolerance=',num2str(tol),', for N=',num2str(nx)])

figure 
contour(u),title(['Contour plot for SOR Method solution',',  tolerance=',num2str(tol),', for N=',num2str(nx)])











