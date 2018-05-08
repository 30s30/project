clear all; clc; close all;
addpath('./')

%                                                 f_a(by)+....
ax= -pi;%                        (-pi,pi)     ____________________ (pi,pi)  
ay= -pi;   %                                 |                    |                          
bx= pi;%                                 f_a |                    |  g_a
by= pi;%                                     |                    |     
%                                            |____________________|
lambda = pi; %                   (-pi,-pi)           du/dy =0       (pi,-pi)  

NUM=5; 

for iter=1:NUM 
N=5*2^iter;	
M=5*2^iter;
hx(iter)=(bx-ax)/(N+1);  
hy(iter)=(by-ay)/(M+1);


[X{iter},Y{iter}]= meshgrid(ax:hx(iter):bx,ay:hy(iter):by);  % create grid to do the boundary conditions

B{iter} = BVP(X{iter},Y{iter},ax,ay,bx,by,N,M);  % call the function to solve boundary conditions
f{iter} = rhs(X{iter},Y{iter},ax,ay,bx,by);


% SOLVING RIGHT-HAND SIDE

%% Neumann conditions 
% In-between points
f{iter}(1,3:M)= (f{iter}(1,(3:M))*hx(iter)^2)+(2*hy(iter)*B{iter}(1,3:M));   % TOP (my bottom)

% Corners 
f{iter}(1,2) = (f{iter}(1,2)*hx(iter)^2)+(2*hy(iter)*B{iter}(1,2))-B{iter}(1,1); %TOPLEFT
f{iter}(1,M+1) = (f{iter}(1,M+1)*hx(iter)^2)+(2*hy(iter)*B{iter}(1,M+1))- B{iter}(1,M+2); %TOPRIGHT

%% Dirichlet conditions 
f{iter}(2:N,2)= (f{iter}(2:N,2)*hx(iter)^2)-B{iter}(2:N,1);  %LEFT 
f{iter}(N+1,3:M) = (f{iter}(N+1,3:M)*hx(iter)^2)-B{iter}(N+2,3:M); %BOTTOM 
f{iter}(2:N,M+1) = (f{iter}(2:N,M+1)*hx(iter)^2)-B{iter}(2:N,M+2); %RIGHT 

% Corners 

f{iter}(N+1,2) = (f{iter}(N+1,2)*hx(iter)^2)- B{iter}(N+2,2)- B{iter}(N+1,1); %BOTTOMLEFT
f{iter}(N+1,M+1) = (f{iter}(N+1,M+1)*hx(iter)^2) -B{iter}(N+2,M+1)-B{iter}(N+1,M+2); %BOTTOMRIGHT

% Setting right hand side as vector 

f{iter} = [f{iter}(1:N+1,2:M+1)]';   % only the N points of matter 
f{iter} = f{iter}(:);   

% Build the pentadiagonal matrix 

A{iter} = pentad(N,M,hx(iter)); 

% Find solution 

u{iter} = A{iter}\f{iter}; 


end 

