function A = pentad(N,M,hx)
% Build matrix B
r2 = (-(4-pi)*hx^2)*ones(N+1,1);
r = ones(N,1);
B = diag(r2,0) + diag(r,1) + diag(r,-1);
% Sparse matrix B
B = sparse(B);
% Build sparse identity matrix
I = speye(N);
% Build tridiagonal block matrix A
A = kron(B,I) + kron(I,B);
% Modify upper right diagonal 

A(1:N-1,M+1:2*M)= A(1:N-1,M+1:2*M)*2; 
A(N,2*M) = A(N,2*M)*2;
end 



