function result = verifylss_sparse(A,b,opts)

% Sparse version of verifylss
% for positive definite A
%
% the opts structure should contain a lower bound
% for the eigenvalue of A

minEig = opts.minEig; % estimate for smallest eigenvalue of A

disp("Verifylss Sparse");

xt   = mldivide(A,b.mid);           % floating point system

% compute residual in interval arithmetic
res  = A*intval(xt)-b; 

nres = norm(res);

ub   = 1/minEig*nres; 

result = midrad(xt,ub.sup*ones(size(xt)));
