function B = BVP(X,Y,ax,ay,bx,by,N,M) 

B(:,1) = Y(:,1).*(Y(:,1)-ay).^2;                        %left 
B(:,M+2)= (Y(:,M+2)-ay).^2 .*cos(pi*Y(:,M+2)/ay);       %right
B(1,:) = zeros(N+2,1);                                % top (my bottom)
B(N+2,:) = (by*(by-ay)^2)+((X(N+2,:)-ax)/(bx-ax))*(((by-ay)^2*cos(pi*by/ay))-(by*(by-ay)^2));  % bottom (my top)

end 

