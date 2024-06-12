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

pl = 0;   % plotting

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
[v,d] = eigs(K1.mid,M1.mid,1,'sm','Tolerance',1e-12);
toc

fprintf("     First eigenvalue on the slice = %.6f\n",d(1));


% compute an enclosure for the first eigenvalue
% and the first eigenvector using INTLAB

disp("Validated Generalized eigenvalue problems");
disp("     First eigenvalue: on slice, exploiting symmetry");
tic
[LB1,X1] = verifyeig_BB(K1,d(1),v(:,1),M1);
toc
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

u = zeros(nv,2);
u = intval(zeros(nv,2));
u(Inside)=X1;

if(pl==1)
	clf
	patch('Faces',tri,'Vertices',pts,'FaceVertexCData',u(:,1).mid,'FaceColor','Interp','Edgecolor','none');
	axis equal
	colorbar 
end


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

% construct full symmetric mesh for the regular polygon
res = PolyaMesh(n,m);


pts    = res.pts;
tri    = res.tri;
Inside = res.Inside;
Index0 = res.Index0;


% Assembly of full mass and rigidity
% explicit matrices using intervals
[K,M] = dir_assemKM_polya(pts,tri,opts);

K1 = K(Inside,Inside);
M1 = M(Inside,Inside);

size(K1)

% Solve the floating point generalized eigenvalue problem
[v,d] = eigs(K1.mid,M1.mid,2,'sm','Tolerance',1e-12,...
             'IsSymmetricDefinite',true);
d

centers = 1/3*(pts(tri(:,1),:)+pts(tri(:,2),:)+pts(tri(:,3),:));

ItriPos = find(centers(:,2)>=-1e-6);
TriPos  = tri(ItriPos,:);

hy = h*sin(pi/n);
Ip = find(pts(:,2)>=(-hy-1e-6));
Imod = find(and(pts(:,2)<-1e-6,pts(:,2)>=(-hy-1e-6)));


clf
hold on

plot(pts(:,1),pts(:,2),'.b','MarkerSize',10);

pts(Imod,2) = 0;
%plot(pts(Imod,1),pts(Imod,2),'.m','MarkerSize',10);

Posp = find(pts(:,2)>=(1e-6));
InPos = intersect(Posp,Inside);

plot(pts(InPos,1),pts(InPos,2),'.k','MarkerSize',15);

patch('Faces',TriPos,'Vertices',pts,'FaceColor','none','Edgecolor','k');




plot(pts(Ip,1),pts(Ip,2),'.','color',[0,0.5,0],'MarkerSize',10);
%plot(pts(Imod,1),pts(Imod,2),'.m','MarkerSize',10);

ImodUp = find(and(pts(:,2)>1e-6,pts(:,2)<=(hy+1e-6)));

plot(pts(ImodUp,1),pts(ImodUp,2),'.r','MarkerSize',25);

axis equal
axis off
axis tight

K0 = K(Inside,Inside);
M0 = M(Inside,Inside);

K1 = K0(Ip,Ip);
M1 = M0(Ip,Ip);

[K,M] = dir_assemKM(pts,TriPos);

Kl = K(m+3,:);
Kl(abs(Kl)>0)

pause


ImodUp = find(and(pts(:,2)>1e-6,pts(:,2)<=(hy+1e-6)));

Kd = diag(K);
Kd(ImodUp)


plot(pts(ImodUp,1),pts(ImodUp,2),'.g');


%lk = K(ImodUp(2),:)
%lk(abs(lk)>1e-5)

%lk = K(InPos(1),:)
%lk(abs(lk)>1e-5)



K1 = K(InPos,InPos);
M1 = M(InPos,InPos);

K11 = K1(abs(K1)>0);
M11 = M1(abs(M1)>0);

uniquetol(full(K11(:)))
uniquetol(full(M11(:)))/ArTri.mid

%full(M1(1:10,1:10))/ArTri.mid

size(K1)

% Solve the floating point generalized eigenvalue problem
[v,d] = eigs(K1,M1,2,'sm','Tolerance',1e-12,...
             'IsSymmetricDefinite',true);
hold off
pause


d = diag(d)

% Compute the residual vector for the second
% eigenvalue
residual = K1*v(:,2)-d(2)*M1*v(:,2);


% Compute norms with respect to the mass matrix
% allowing to find a guaranteed enclosure
% See Proposition 4.4 in the paper

meth = 2;

if meth==1
	% First method: using Intlab
	disp("Verifylss for residual");
	tic
	ressol  = verifylss(M1,residual);
	toc
	resnorm = sqrt(abs(dot(residual,ressol)))
else
	% Second method: residual norm estimate
	% Remark 4.5
	disp("     Direct estimation - eigenvalues of mass matrix");
	resnorm2 = norm(residual);
	resnorm  = sqrt(intval(2)/(min(n,6)*ArTri))*resnorm2;
end
vecnorm = sqrt(dot(v(:,2),M1*v(:,2)));

lb2inf = d(2)-resnorm/vecnorm;
lb2sup = d(2)+resnorm/vecnorm;

% Construct the enclosure for the second eigenvalue
LB2 = infsup(lb2inf.inf,lb2sup.sup);
fprintf("     Second eigenvalue enclosure = [%.12f,%.12f]\n",LB2.inf,LB2.sup);

fprintf("     Interval radius - second eigenvalue  %.2e\n",rad(LB2));

     
save_string = ['IntervalPolya_',num2str(n),'_',num2str(m),'.mat'];



% Save the result
%save(['IntervalPolya_',num2str(n),'_',num2str(m),'.mat'],...
%      'LB1','LB2','X1');

%fprintf(['Results saved in file: ',save_string,'\n']);

