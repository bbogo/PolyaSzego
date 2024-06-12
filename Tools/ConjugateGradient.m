function xminres = ConjugateGradient(A,b,tol,x0)

maxiter = 3000;

if nargin<3
	tol=1e-10;
end

n = size(A,1);


if nargin<4
x = zeros(n,1);
else
x = x0;
end

kprint  = 1;

% store vector with minimal residual
minres  = 1e16;
xminres = x;

g = A*x-b;
d = -g;

for i=1:maxiter
	vv    = A*d; % mat vect mult
	%dot(d,vv)
	%gamma = -Dot_(d',g)/Dot_(d',vv);
	gamma = -dot(d,g)/dot(d,vv);
	x     = x+gamma*d;
	g     = g+gamma*vv;%A*x-b;
	%beta  = Dot_(g',vv)/Dot_(d',vv);
	beta  = dot(g,vv)/dot(d,vv);
	
	d     = -g+beta*d; 
	
	
	ng = norm(g);
	if ng<minres
		minres = ng;
		xminres = x;
	end
	
	if minres<10^(-kprint)
		fprintf("Iter %4d | norm %.3e  | min %.3e \n",i,norm(g),minres);
		kprint = kprint+1;
	end
	
	
	if ng <tol
		break
	end
end

res = xminres;
