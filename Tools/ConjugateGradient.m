function xminres = ConjugateGradient(A,b,tol,x0)
%
% Conjugate gradient algorithm: basic implementation
%
% A is a symmetric positive definite matrix
% b is a vector
% tol is a tolerance for the residual vector
% x0 is an initialization - zero by default
% 
% the algorithm terminates after a number of iterations
% and keeps track of the minimum residual vector found

% maximum number of iterations
maxiter = 3000;

% default tolerance (if not provided)
if nargin<3
	tol=1e-12;
end

% size of the matrix
n = size(A,1);

% initialization: zero by default	
if nargin<4
	x = zeros(n,1);
else
	x = x0;
end

% print info when residual decreases enough
kprint  = 1;

% store vector with minimal residual
minres  = 1e16;
xminres = x;
itminres = 1;

Converged = 0;

% initial gradient and descent direction
g = A*x-b;
d = -g;

% Conjugate gradient loop
for i=1:maxiter
	% one matrix vector product
	vv    = A*d; 
	% classical formulas
	gamma = -dot(d,g)/dot(d,vv);
	x     = x+gamma*d;
	g     = g+gamma*vv;
	beta  = dot(g',vv)/dot(d',vv);
	d     = -g+beta*d; 
	
	ng = norm(g);
	% compare with minimal residual
	if ng<minres
		minres   = ng;
		xminres  = x;
		itminres = i;
	end
	% print if residual decreases enough
	if minres<10^(-kprint)
		fprintf("Iter %4d | norm %.3e  | min %.3e \n",i,norm(g),minres);
		kprint = kprint+1;
	end
	
	% stop if tolerance is reached
	if ng <tol
		Converged = 1;
		break
	end
end

if(Converged==1)
	fprintf("==========Desired tolerance reached\n");
else
	fprintf("==========Maximum number of iterations reached\n");
end
fprintf("Minimal residual %.3e found at iteration %d\n",minres,itminres);
% output vector of minimal residual
res = xminres;
