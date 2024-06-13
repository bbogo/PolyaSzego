function [Kxx,Kx] = AssemblyMorley(pts,tri,data)

% Assembly for Morley finite element
% H2 semi norm and gradient norm
%
% Input: mesh parameters
% pts     points:    npt  x 2
% tri     triangles: ntri x 3

KT  = data.KT;  % explicit interval triangle rigidity matrix
arK = data.ArT; % area of small triangle in mesh (same for all, interval)
% edge lengths
l12 = data.l12;
l13 = data.l13;
l23 = data.l23;
m   = data.m;
a   = data.a;
b   = data.b;


np = size(pts,1); % number of points
nt = size(tri,1); % number of triangles

% all edges
AllEdges = [tri(:,[2,3]);
            tri(:,[1,3]);
            tri(:,[1,2]);];
            
% Sort edges
[AllEdges,IS] = sort(AllEdges,2);
% Put a sign related to original orientation
Sign = -2*IS(:,1)+3;
Sign = reshape(Sign,[],3); % orientation of edges
            
% Find unique edges
[Edges,IA,IC] = unique(AllEdges,'rows');

% number of edges
nbe = size(Edges,1);

% number of degrees of freedom
ndof = np+nbe;

% Degrees of freedom: [points values, edge values]

% Initialization
% Rigidity matrix for second derivatives
Kxx    = intval(sparse(ndof,ndof));
% Rigidity matrix for first derivatives
Kx     = intval(sparse(ndof,ndof));

% iterate on triangles and add to matrices

R = [0 -1; 1 0]; % 90 degrees rotation matrix

% initialize matrix for rotating coordinate systems aligned to edges
% same for all triangles, orientation flips after m*(m+1)/2

M2init = intval(zeros(6,6));
M2init(5:6,5:6) = eye(2);

cc = (a-1)/sqrt((a-1)^2+b^2);
ss = b/sqrt((a-1)^2+b^2);

M2init(1:2,1:2) = [cc -ss;
                   ss cc];
                  
cc = -a/sqrt(a^2+b^2);
ss = -b/sqrt(a^2+b^2);

M2init(3:4,3:4) = [cc -ss;
                   ss cc];

% Assembly loop
for j = 1:nt
	% vertex indices
	act    = tri(j,:); % actual triangle

	subpts = pts(act,:); % actual points
	pt1    = pts(tri(j,1),:);
	pt2    = pts(tri(j,2),:);
	pt3    = pts(tri(j,3),:);


	% edge indices: 23-13-12
	ic = IC([j,j+nt,j+2*nt]);

	% corresponding edges
	Edge = Edges(ic,:);

	% edge sign orientations
	SS = Sign(j,:);
	
	% find gradients as affine functions
	M1 = kron([-1 1 1; 1 -1 1; 1 1 -1],eye(2));
	
	% Orientation flip for triangle in second half
	if j<=(m*(m+1)/2)
		M2 = M2init;
	else
		M2 = -M2init;
	end
	
	% third matrix
	M3 = [0     -1/l23  1/l23  0      0     0;
		  0      0      0      SS(1)  0     0;
		  1/l13  0     -1/l13  0      0     0;
		  0      0      0      0      SS(2) 0;
		 -1/l12  1/l12  0      0      0     0;
		  0      0      0      0      0     SS(3)];
		  
	% permutation matrix (bring px onto first three,
	%                     py onto last three)

	P  = [1 0 0 0 0 0;
		  0 0 1 0 0 0;
		  0 0 0 0 1 0;
		  0 1 0 0 0 0;
		  0 0 0 1 0 0;
		  0 0 0 0 0 1];
	
	% intermediary matrix	  
	Interm = P*M1*M2*M3;

	% build P1 mass and rigidity matrices for T

	Mt = arK/12*[2 1 1;
	             1 2 1;
	             1 1 2];

	% rigidity matrix - precomputed
	Kt = KT;
		          
	M2 = kron(eye(2),Mt);
	K2 = kron(eye(2),Kt);

	% Compute contribution to current DoFs
	fullK = Interm'*K2*Interm;
	fullM = Interm'*M2*Interm;
	
	% actual dofs [vert; edges]
	actdof = [act(:); np+ic(:)];
	
	% add contribution to the assembled matrices
	Kxx(actdof,actdof) = Kxx(actdof,actdof)+fullK;
	Kx(actdof,actdof)  = Kx(actdof,actdof) +fullM;
end

subplot(1,2,1)
spy(Kx)

subplot(1,2,2)
spy(Kxx)





