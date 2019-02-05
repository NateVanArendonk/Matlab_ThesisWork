function M=dup(x,n);
% M=dup(x,n)
% Duplicates a column vector n times to return a matrix

%x=(1:5)'
M=x(:,ones(n,1));
