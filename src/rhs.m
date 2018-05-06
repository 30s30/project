function f = rhs(X,Y,ax,ay,bx,by) 

f = cos(pi/2*(2*((X-ax)/(bx-ax))+1)).*sin(pi*(Y-ay)/(by-ay));


end 