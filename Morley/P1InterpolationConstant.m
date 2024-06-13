function C4 = P1InterpolationConstant(a,b,m)
% Computing an upper bound for the
% interpolation constant on the triangle with vertices
% (0,0), (1,0), (a,b)
%
% Input: third vertex coordinates (a,b)
% Mesh parameter m: divide each side into m segments
%
% m rather small should be used, e.g. m=10
%
% Output: certified upper bound for the interpolation constant
%
% Requires a working installation of Intlab
%
% based on the paper: K. Kobayashi, On the interpolation
% constants over triangular elements
% 

%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Area of the triangle
ArT = b/2;

% Precomputing P1 rigidity matrix on the triangle (scale invariant)
	KT = intval(zeros(3,3));
	KT(1,1) = b/2+(a-1)^2/b/2;
	KT(2,2) = b/2+a^2/b/2;
	KT(3,3) = 1/b/2;

	KT(1,2) = -b/2-a*(a-1)/b/2;
	KT(2,1) = KT(1,2);
	KT(1,3) = (a-1)/b/2;
	KT(3,1) = KT(1,3);
	KT(2,3) = -a/b/2;
	KT(3,2) = KT(2,3);


% store data regarding the triangle
data.KT  = KT;
data.ArT = ArT/m^2;
data.l12 = 1/m;
data.l13 = sqrt(a^2+b^2)/m;
data.l23 = sqrt((a-1)^2+b^2)/m;
data.m   = m;
data.a   = a;
data.b   = b;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% switch to floating point if a and b are not intervals
% the mesh is computed to have triangle indices
% but all information is explicit using intervals
try
	[pts,tri,Inside,ends] = MeshTriangle(a.mid,b.mid,m);
catch
	[pts,tri,Inside,ends] = MeshTriangle(a,b,m);
end

% Assemble matrices for the H^2 semi-norm and the gradient norm
% for the Morley element
[Kxx,Kx] = AssemblyMorleyTriInterval(pts,tri,data);

% Dirichlet condition: only indices of the corners
Dind = [1,m+1,size(pts,1)];

ndof = size(Kxx,1);

% Find DoF indices for non-Dirichlet conditions
inside = setdiff(1:ndof,Dind);
% Find DoF indices for points different from corners (plotting)
inside0 = setdiff(1:size(pts,1),Dind);

% Submatrices excluding Dirichlet nodes (corners)
Kxx0 = Kxx(inside,inside);
Kx0  = Kx(inside,inside);

% Make matrices symmetric
K1 = 0.5*(Kxx0+Kxx0');
K2 = 0.5*(Kx0 +Kx0'); 

% Floating point estimation for smallest eigenvalue
[v,d] = eigs(K1.mid,K2.mid,1,'sm','Tolerance',1e-12);

% Lower bound tentative - first eigenvalue of floating point problem
% minus something small

LB = intval(d(1)-1e-6);

% Define a difference matrix
AA    = Kxx0-LB*Kx0;

% Check if this matrix is symmetric positive definite
rr = isspd(AA,1,1,1);

% if positivity checks out then a lower bound is validated
if rr
	LB;
	C4 = sqrt(1/LB);
	% Kobayashi upper bound formula
	C4 = sqrt(m^2/(m^2-1))*C4;
	fprintf("The desired bound was found\n");
else
	fprintf("The desired bound was not found\n");
end


cc = zeros(size(pts,1),1);
cc(inside0) = v(1:(length(inside0)),1);

clf
patch('Faces',tri,'Vertices',pts,'FaceVertexCData',cc,'FaceColor','Interp','Edgecolor','k');

colorbar

axis equal

