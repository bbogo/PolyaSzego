function PolyaHessInterval(n,m)
%
%
% PolyaHessInterval(n,m)
%
% Inputs: n - number of vertices of the polygon
%         m - number of division points on a ray
%             for mesh construction
% 
% 1. Solves and validates bounds for first eigenvalue
%    and eigenvector using INTLAB. 
%    In particular, a constant sign is proved for the
%    first eigenvector. 
%
% 2. For the second eigenvalue an enclosure based on the
%    residual is computed. The residual is computed
%    using interval arithmetic.
%
% All matrices involved, including the Rigidity and
% Mass matrices are computed analytically,
% the numerical computations use intervals
%
% References: 
% [B. Bogosel, D. Bucur, Polygonal Faber-Krahn inequality
% Local minimality via Validated Computing, 2024]
%
% [B. Bogosel, D. Bucur, On the Polygonal Faber-Krahn inequality,
% Journal de l'Ecole Polytechnique, 2024]


h = 1/intval(m);       % mesh size
Pi = intval('pi');     % interval enclosure for pi
Theta = 2*Pi/n;        % central angle
ArTri = 0.5*h^2*sin(Theta);   % area of small triangle

% Mesh generation for a slice
% indices of nodes inside the polygon are also given
[pts,tri,Inside,ends] = MeshSlice(n,m);
nv = size(pts,1);



% Computing the first eigenvalue
% disp("Assembly rigidity and mass")
opts.arM = ArTri;
opts.n   = n;

% exact assembly using interval matrices
[K,M]       = dir_assemKM_polya(pts,tri,opts);

% Select only nodes for interior points
% to solve the Dirichlet problem
K1 = K(Inside,Inside);
M1 = M(Inside,Inside);

% compute the first eigenvalue on the slice 
% floating point computation
disp("Floating point Generalized eigenvalue problem");
tic
[v,d] = eigs(K1.mid,M1.mid,1,'sm','Tolerance',1e-14);
toc

fprintf("     First eigenvalue on the slice = %.6f\n",d(1));


fprintf("Residual Method\n");
% distance between eigenvalues (to be quantified)

% compute the residual with Intlab
residual1 = K1*v(:,1)-d(1)*M1*v(:,1);

%invMr     = verifylss(M1,residual1);

%resnorm1  = sqrt(abs(dot(residual1,invMr)))

resnorm1_sq  = abs(dot(residual1,residual1))
resnorm1  = (intval(2)/(min(n,6)*ArTri))*resnorm1_sq



% for normalization error
v1       = v(:,1);
vecnorm1 = (dot(v1.',M1*v1));

vecnorm1.rad
vecnorm1.mid


lb1inf = d(1)-sqrt(resnorm1/vecnorm1);
lb1sup = d(1)+sqrt(resnorm1/vecnorm1);
LB1 = infsup(lb1inf.inf,lb1sup.sup);

fprintf("     First eigenvalue enclosure = [%.12f,%.12f]\n",LB1.inf,LB1.sup);

fprintf("     Interval radius - first eigenvalue  %.2e\n",rad(LB1));
v1 = v(:,1); % store first eigenvector


% Compute an enclosure for the second eigenvalue
% The full mesh needs to be used since there is no
% symmetry anymore. 
%
% Theoretical results imply that the second eigenvalue
% is double. Comparing with the fourth eigenvalue
% on the unit disk assures that the enclosed 
% eigenvalue is indeed the second one!

disp(" ");

disp("Second eigenvalue enclosure using residual");
% Compute an enclosure for the second eigenvalue
% The full mesh needs to be used since there is no
% symmetry anymore. 
%
% Theoretical results imply that the second eigenvalue
% is double. Comparing with the fourth eigenvalue
% on the unit disk assures that the enclosed 
% eigenvalue is indeed the second one!


% construct full symmetric mesh for the regular polygon
res = PolyaMesh(n,m);


pts    = res.pts;
tri    = res.tri;
Inside = res.Inside;
Index0 = res.Index0;


% Assembly of full mass and rigidity
% explicit matrices using intervals
[K,M] = dir_assemKM_polya(pts,tri,opts);

hy = h*sin(pi/n);

% index  of diagonal elements to modify in K,M 
ImodUp = find(and(pts(:,2)>1e-6,pts(:,2)<=(hy+1e-6)));



for i=1:length(ImodUp)
	jj = ImodUp(i);
	K(jj,jj) = 2*tan(Theta/2)+1/sin(Theta)+1/tan(Theta/2);
	M(jj,jj) = ArTri*5/6;
end

Pos   = find(pts(:,2)>1e-6);
InPos = intersect(Inside,Pos);

K1 = K(InPos,InPos);
M1 = M(InPos,InPos);

size(K)
size(K1)

% Solve the floating point generalized eigenvalue problem
[v,d] = eigs(K1.mid,M1.mid,2,'sm','Tolerance',1e-12,...
             'IsSymmetricDefinite',true);



d = diag(d)

% Compute the residual vector for the second
% eigenvalue
residual2 = K1*v(:,1)-d(1)*M1*v(:,1);


% Compute norms with respect to the mass matrix
% allowing to find a guaranteed enclosure
% See Proposition 4.4 in the paper


	% Residual norm estimate
	% Remark 4.5
disp("     Direct estimation - eigenvalues of mass matrix");
resnorm2 = abs(dot(residual2,residual2));
resnorm2  = (intval(2)/(min(n,6)*ArTri))*resnorm2;

vecnorm2 = (dot(v(:,1),M1*v(:,1)));

lb2inf = d(1)-sqrt(resnorm2/vecnorm2);
lb2sup = d(1)+sqrt(resnorm2/vecnorm2);

% Construct the enclosure for the second eigenvalue
LB2 = infsup(lb2inf.inf,lb2sup.sup);
%LB2 = LB2/ArTri;
fprintf("     Second eigenvalue enclosure = [%.12f,%.12f]\n",LB2.inf,LB2.sup);

fprintf("     Interval radius - second eigenvalue  %.2e\n",rad(LB2));


% ======================= Eigenvector estimate ===================

a = LB2.inf-LB1.sup; % difference between first two eigenvalues
a = intval(a);
Qub = resnorm1/a^2;

ub1 = Qub+(sqrt(vecnorm1)-1)^2;
ub2 = Qub+((intval(1)-vecnorm1+Qub)/(1+sqrt(vecnorm1-Qub)))^2;

ev_radM = max(ub1,ub2);
mmm = ev_radM.sup


%abs(vecnorm1-1)
ev_rad  = 2/(min(n,6)*ArTri)*ev_radM+abs(vecnorm1-1)
ev_rad  = sqrt(ev_rad)
% compute an enclosure for the first eigenvalue
% and the first eigenvector using INTLAB

X1 = midrad(v1,ev_rad.sup*ones(size(v1)));




fprintf("     First eigenvalue enclosure = [%.12f,%.12f]\n",LB1.inf,LB1.sup);

fprintf("     Interval radius - eigenvalue  %.2e\n",rad(LB1));
fprintf("     Interval radius - eigenvector %.2e\n",max(rad(X1)));

% Check that the first eigenvector has a constant
% sign: this guarantees that the eigenvalue enclosed
% is indeed the first one
disp("Check if first eigenvector has constant sign");
fprintf("     Min lower bound=%f\n",min(X1.inf));
fprintf("     Max lower bound=%f\n",max(X1.inf));

fprintf("     Min upper bound=%f\n",min(X1.sup));
fprintf("     Max upper bound=%f\n",max(X1.sup));


% ================================================================




     
save_string = ['IntervalPolya_',num2str(n),'_',num2str(m),'.mat'];

difflam = LB2.inf-LB1.sup

% Save the result
save(['IntervalPolya_',num2str(n),'_',num2str(m),'.mat'],...
      'LB1','LB2','X1');

fprintf(['Results saved in file: ',save_string,'\n']);

