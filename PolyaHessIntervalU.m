function res = PolyaHessIntervalU(n,m,hh)

% PolyaHessInterval(n,m,hh)
%
% Inputs: n - number of vertices of the polygon
%         m - number of division points on a ray 
%             for mesh construction
%         hh - an alternate mesh size to check if 
%              computations can validate the local minimality
%
% Computation and validation of positivity for
% the eigenvalues of the Hessian matrix of the first
% Dirichlet-Laplacian eigenvalue on an n-gon
%
% 
% References: 
% [1] [B. Bogosel, D. Bucur, Polygonal Faber-Krahn inequality
% Local minimality via Validated Computing, 2024]
%
% [2] [B. Bogosel, D. Bucur, On the Polygonal Faber-Krahn inequality,
% Journal de l'Ecole Polytechnique, 2024]


% plotting parameter: only used for testing
pl = 0;

Pi = intval('pi');          % interval enclosure for pi
theta = 2*Pi/n;             % central angle
Area  = n*0.5*sin(theta);   % area of the regular n-gon
h = 1/intval(m);            % mesh parameter
ArTri = 0.5*h^2*sin(theta); % area of a small triangle in the mesh

opts.arM = ArTri;
opts.n   = n;

gamfact  = 10;            % factor for scaling the normalization condition


% Mesh generation for a slice of the regular n-gon
[pts0,tri0,Inside0,ends0] = MeshSlice(n,m);

% loading eigenvalue and eigenfunction information
% if the information is not available then it is computed
str = ['./Results/IntervalPolya_',num2str(n),'_',num2str(m),'.mat'];

try
	load(str);
catch
	if m<=400
		PolyaHessInterval(n,m);
	else
		PolyaHessInterval_Sparse(n,m);
	end
	load(str);
end

% The first eigenvector X1 is normalized on T_+
% need to divide by sqrt(n) to have normalization
% on the whole polygon 

sqrtn = intval(sqrt(intval(n)));
X1 = X1/sqrtn;

% set the first eigenfunction on the inside nodes

% set the first and second eigenvalues as interval enclosures
lam1 = LB1;
lam2 = LB2;

% create mapping from slice to the whole mesh
% to transporte the first eigenfunction on a slice to the whole mesh


% Construct the full mesh of the regular polygon
res = PolyaMesh(n,m);

pts    = res.pts;
tri    = res.tri;
Inside = res.Inside;
% index of the corresponding point in the first slice
% to transfer information from the first slice to all
% the mesh

Index0 = res.Index0;
nv = size(pts,1);
nt = size(tri,1);

% Finite Element Assembly
% this routine computes analytical values
% and evaluates them in interval matrices
disp("Assembly Mass and rigidity for the whole mesh using intervals")
tic
[K_int,M_int] = dir_assemKM_polya(pts,tri,opts);
toc

% ======================= 
% region detecting using centroids for mesh triangles
reg = zeros(size(tri,1),1);
centers = 1/3*(pts(tri(:,1),:)+pts(tri(:,2),:)+pts(tri(:,3),:));

if(pl==1)
   plot(centers(:,1),centers(:,2),'.b');
end

ang = atan2(centers(:,2),centers(:,1));
ang(ang<0) = ang(ang<0)+2*pi;
reg = floor(ang/(2*pi/n));

% plotting for checking region detection
if(pl==1)
	patch('Faces',tri,'Vertices',pts,'FaceVertexCData',reg(:),'FaceColor','flat','Edgecolor','none');
	axis equal
	axis tight
end

% ========================
% solving equations for material derivatives associated to 
% the movement of the point a_0 = (1,0)
% U^1,U^2 in the references

% Find triangles near the ray [0,0]->[1,0]
% Triangle indices for the union Tplus U Tminus
Tplus = find(or(reg==0,reg==(n-1)));

% Triangle indices for the region T_+
Tpp   = find(reg==0);
% Triangle indices for the region T_-
Tmm   = find(reg==(n-1));

% See Proposition 4.6 in [1] for more details
% regarding the assembly
disp("Assembly For Saddle point systems using interval matrices")
tic
Kxx_int = dir_assemKMxx_polya(pts,tri(Tplus,:),opts);
Kyy_int = dir_assemKMyy_polya(pts,tri(Tplus,:),opts);
Kxy_p = dir_assemKMxy_polya_plus(pts,tri(Tpp,:),opts);
Kxy_m = dir_assemKMxy_polya_minus(pts,tri(Tmm,:),opts);

Kxy_int = Kxy_p-Kxy_m;
Kyx_int = Kxy_p+Kxy_m;
toc

% Generate eigenfunction on the whole mesh using intervals
% u0_int on the first slice
u0_int = intval(zeros(size(pts0,1),1));
u0_int(Inside0) = X1;

% u_int on the whole mesh
u_int = intval(zeros(size(pts,1),1));
u_int = u0_int(Index0);

% number of points on segment [0,0]-->[1,0]
nseg = m+1;

if(pl==1)
	% plot points on the first ray if needed
	plot(pts(1:nseg,1),pts(1:nseg,2),'.');
	axis equal
end

% Store values of u on the segment
useg = u_int(1:nseg);

% form two vectors to apply the trapeze 
% quadratule rule on the first ray
usm = useg(1:end-1);
usp = useg(2:end);

% compute upper bound for L^2 norm of dx(u) on S_0
% see [2] Remark 5.5
% an additional term involving |grad(u)-grad(u_h)| is added later on
intuseg = sqrt(lam1*h*0.5*sum(usm.^2+usp.^2))

fprintf("Estimate for L^2 norm of dx(u) on S_0 = %f\n",intuseg.mid);

% build Right hand sides for U^1, U^2
% Equations (49), (50) in [1]
rhs1_int = 2*Kxx_int*u_int-1/tan(theta)*(Kxy_int)*u_int-2*LB1/n*M_int*u_int;
rhs2_int = -2/tan(theta)*Kyy_int*u_int+Kyx_int*u_int;

% build constraint vector
con0 = M_int*u_int;
gam = gamfact*1/norm(con0);

% re-normalize constraint vector
con = gam*con0;


% consider submatrices corresponding to the inside nodes
% Dirichlet boundary conditions
K0_int = K_int(Inside,Inside);
M0_int = M_int(Inside,Inside);

% Matrix for the linear system for U^1, U^2
mat = K0_int-LB1*M0_int; 

% Matrix including the constraint vectors
Mat = [mat, con(Inside);
       con(Inside).', 0];
       
% Right hand sides including the constraints
R1 = [rhs1_int(Inside);0];
R2 = [rhs2_int(Inside);0];

% radii for the right hand sides
fprintf("	Radius R1 = %.3e\n",max(R1.rad));
fprintf("	Radius R2 = %.3e\n",max(R2.rad));

% Check orthogonality on the first eigenfunction
dotR1 = dot(R1,[u_int(Inside);0]);
dotR2 = dot(R2,[u_int(Inside);0]);

fprintf("Orthogonality RHS on u_1: %f %f\n",dotR1.mid,dotR2.mid);

disp("Solve the Linear systems: Conjugte Gradient with zero initialization");
disp("Floating point residual estimate");

U1 = zeros(nv,1);
U2 = zeros(nv,1);

% tolerance for the conjugate gradient
% algorithm will try to reach this and return
% the vector with smallest residual
tolCG = 1e-12; 
CG = 1;    % Use the conjugate gradient implementation or something else

if (CG==1)
	% Apply the conjugate gradient method for the two systems
	tic
	U1(Inside) = ConjugateGradient(mat.mid,rhs1_int(Inside).mid,tolCG);
	toc
	% test orthogonality on the first eigenfunction
	dotU1=dot(U1,con);
	fprintf("Orthogohality of U1 on u1: %f %f\n",dotU1.inf,dotU1.sup);

	tic
	U2(Inside) = ConjugateGradient(mat.mid,rhs2_int(Inside).mid,tolCG);
	toc
	% test orthogonality on the first eigenfunction
	dotU2=dot(U2,con);
	fprintf("Orthogohality of U2 on u1: %f %f\n",dotU2.inf,dotU2.sup);
else
    %setup = struct('type','ilutp','droptol',1e-6);
	%[L,U] = ilu(Mat.mid,setup);
	U1t  = bicg(Mat.mid,[rhs1_int(Inside).mid;0],1e-15,2000);
	U2t  = bicg(Mat.mid,[rhs2_int(Inside).mid;0],1e-15,2000);
	U1(Inside) = U1t(1:end-1);
	U2(Inside) = U2t(1:end-1);
end

% Compute the residual for the "big system" including the constraint
% using interval arithmetics
Residual1 = Mat*[U1(Inside);0]-[rhs1_int(Inside);0];
Residual2 = Mat*[U2(Inside);0]-[rhs2_int(Inside);0];

% Eventual plotting for U_1,U_2
if(pl==1)
	figure(1)
	subplot(1,2,1)
	patch('Faces',tri,'Vertices',pts,'FaceVertexCData',U1(:),'FaceColor','Interp','Edgecolor','none');
	axis equal
	colorbar 
	axis tight
	subplot(1,2,2)
	patch('Faces',tri,'Vertices',pts,'FaceVertexCData',U2(:),'FaceColor','Interp','Edgecolor','none');
	axis equal
	colorbar 
	axis tight
end

% Error estimation for saddle point systems
% Lemma 4.7 in [1], computation of quantities involved
% is detailed below 
% Error estimation using norms of residuals:

normR1 = norm(Residual1);
normR2 = norm(Residual2);

% lower bound for the difference lam2-lam1
lbDiff = LB2-LB1;

bb     = con(Inside);
normbb = norm(bb);
normA  = max(mat(:));

w      = (LB2-LB1)/gam^2;  

% upper bound for ||A(w)^{-1}||_2 from the estimate
ubInvAw = max(1/lbDiff,1/(gam^2*w));

% upper bound for |A(w)||b^Tb|^(-1)
ubAw   = max(mat(:))/norm(bb)^2+w;
sqrt5    = sqrt(intval('5'));

% upper bound for error U
ConstU   = intval(2/(sqrt5-1))*max(ubInvAw,ubAw);

% Finding error starting from the residual
ErrorU1  = ConstU*normR1;
ErrorU2  = ConstU*normR2;

fprintf("  Error estimate U1 - discrete: %.2e\n",ErrorU1.sup);
fprintf("  Error estimate U2 - discrete: %.2e\n",ErrorU2.sup);


% Set interval solution starting from the floating point
% solution and enclosing with the residual error
U1_int = midrad(U1,ErrorU1.sup);
U2_int = midrad(U2,ErrorU2.sup);

% plotting U1, U2 if needed
if(pl==1)
	% this part uses floating point computations
	% only used for plotting ==============================================
	ar = abs(sb_simpvol(pts,tri));
	% compute gradients as P0 functions
	[g1x,g1y,g1z,g2x,g2y,g2z,g3x,g3y,g3z]  = sb_pdetrg(pts,tri,ar);
	% ========================================

	it1                                    = tri(:,1);
	it2                                    = tri(:,2);
	it3                                    = tri(:,3);

	dxu = u_int(it1).*g1x+u_int(it2).*g2x+u_int(it3).*g3x;
	dyu = u_int(it1).*g1y+u_int(it2).*g2y+u_int(it3).*g3y;


	dxU1 = U1_int(it1).*g1x+U1_int(it2).*g2x+U1_int(it3).*g3x;
	dyU1 = U1_int(it1).*g1y+U1_int(it2).*g2y+U1_int(it3).*g3y;


	dxU2 = U2_int(it1).*g1x+U2_int(it2).*g2x+U2_int(it3).*g3x;
	dyU2 = U2_int(it1).*g1y+U2_int(it2).*g2y+U2_int(it3).*g3y;

	figure(3)
	patch('Faces',tri,'Vertices',pts,'FaceVertexCData',dxU1(:),'FaceColor','flat','Edgecolor','none');
	axis equal
	colorbar 

	figure(4)
	patch('Faces',tri,'Vertices',pts,'FaceVertexCData',dyU1(:),'FaceColor','flat','Edgecolor','none');
	axis equal
	colorbar 
	
	figure(5)
	cts1 = 1/3*(pts(tri(:,1),1)+pts(tri(:,2),1)+pts(tri(:,3),1));
	cts2 = 1/3*(pts(tri(:,1),2)+pts(tri(:,2),2)+pts(tri(:,3),2));
	quiver(cts1,cts2,dxU1,dyU1);
end

% In what follows, the a-priori estimates for u_1, U^1, U^2 are used
% to generate enclosures for the eigenvalues of the Hessian matrix

% ============ Hessian terms =================

% compute all integrals of the form 
% d(x,y)u * d(x,y) U

% Keep region indices in memory for each one of the slices
for i=1:n
	regions{i} = find(reg==(i-1));
end

% initialize vectors for storing integrals of the form
% integral(slice)(d{x,y}(u) d{x,y}(U^i))

T1xx = intval(zeros(n,1));   % int(slice_i)(dx(u) dx(U^1))
T1yy = intval(zeros(n,1));   % int(slice_i)(dy(u) dy(U^1))
T1xy = intval(zeros(n,1));   % int(slice_i)(dx(u) dy(U^1))
T1yx = intval(zeros(n,1));   % int(slice_i)(dy(u) dx(U^1))

T2xx = intval(zeros(n,1));   % int(slice_i)(dx(u) dx(U^2))
T2yy = intval(zeros(n,1));   % int(slice_i)(dy(u) dy(U^2))
T2xy = intval(zeros(n,1));   % int(slice_i)(dx(u) dy(U^2))
T2yx = intval(zeros(n,1));   % int(slice_i)(dy(u) dx(U^2))

T1xy_full = intval(zeros(n,1)); % int(slice_i)(dx(u) dy(U^1)+dy(u) dx(U^1))
T2xy_full = intval(zeros(n,1)); % int(slice_i)(dx(u) dy(U^2)+dy(u) dx(U^2))

MU1        = intval(zeros(n,1));% int(slice_i)(u_1*U^1)
MU2        = intval(zeros(n,1));% int(slice_i)(u_1*U^2)

% Assemble analytically partial rigidity matrices on slices
opts.reg = reg;
[Kxx,Kyy,Kxy,MU] = dir_assemKpartial_polya(pts,tri,opts);

% Compute the desired integrals using the partial
% rigidity matrices (everything with intervals)
arT = ArTri;
for i=1:n
	T1xx(i) = dot(Kxx{i}*u_int,U1_int); % int(slice_i)(dx(u) dx(U^1))
	T2xx(i) = dot(Kxx{i}*u_int,U2_int); % int(slice_i)(dx(u) dx(U^2))
	
	T1yy(i) = dot(Kyy{i}*u_int,U1_int); % int(slice_i)(dy(u) dy(U^1))
	T2yy(i) = dot(Kyy{i}*u_int,U2_int); % int(slice_i)(dy(u) dy(U^2))
	
	T1xy_full(i) = 	dot(Kxy{i}*u_int,U1_int); % int(slice_i)(dx(u) dy(U^1)+dy(u) dx(U^1))
	T2xy_full(i) = 	dot(Kxy{i}*u_int,U2_int); % int(slice_i)(dx(u) dy(U^2)+dy(u) dx(U^2))
	
	MU1(i)    = 	dot(MU{i}*u_int,U1_int); % int(slice_i)(u*U^1)
	MU2(i)    = 	dot(MU{i}*u_int,U2_int); % int(slice_i)(u*U^1)
end

fprintf("Interval size comparison: using mass vs rigidity matrix\n");
fprintf("   for evaluating quadratures involving u_1\n");

rig_comp = T1xx+T1yy;  % int(T_i) grad(U)*grad(u_1)
mas_comp = LB1*MU1;    % int(T_i) lambda_1*U*u_1
fprintf("Midpoint comparison= %.3e\n",max(rig_comp.mid-mas_comp.mid));
fprintf("Radius   comparison= %.3e\n",max(rig_comp.rad./mas_comp.rad));

% integrals on T0 of dxu and dyu
Ax   = dot(Kxx{1}*u_int,u_int);  % = int(T0)(dx(u)^2);
Ay   = dot(Kyy{1}*u_int,u_int);  % = int(T0)(dy(u)^2);

% add factor appearing in the formula for the eigenvalue of the Hessian
aaa = Ax*sin(theta)/(1-cos(theta))/2/n;
bbb = Ay*sin(theta)/(1-cos(theta))/2/n;

% initialize alpha, beta, gamma
% corresponding to A_k, B_k, C_k in Theorem 2.1 in [1] 
alphas = intval(zeros(n,1));
betas  = intval(zeros(n,1));
gammas = intval(zeros(n,1));

jj = 0:n-1;

ss = sin((2*jj+1)*theta).';
cc = cos((2*jj+1)*theta).';

% Compute alpha, beta, gamma following Theorem 2.1 in [1] 
% All elements are computed using intervals
for k=0:(n-1)
    vv = cos((jj+1)*k*theta)+cos(jj*k*theta);
    ww = (cos((jj+1)*k*theta)-cos(jj*k*theta))/sin(theta);
	alphas(k+1)= 2*n*(1-cos(k*theta))/sin(theta)*Ax...
	          -2*Area*(LB1*dot(vv,MU1)...
	          +dot(ww,-ss.*(T1xx-T1yy)+cc.*(T1xy_full)));
	          
	vv = cos(theta)*( cos((jj+1)*k*theta)-cos(jj*k*theta))/sin(theta);
    ww = (cos((jj+1)*k*theta)-cos(jj*k*theta))/sin(theta);
	betas(k+1)= 2*n*(1-cos(k*theta))/sin(theta)*Ay...
	          -2*Area*(LB1*dot(vv,MU2)...
	          +dot(ww,-cc.*(T2xx-T2yy)-ss.*(T2xy_full)));  

	vv = cos(theta)/sin(theta)*( sin((jj+1)*k*theta)-sin(jj*k*theta));
    ww = (sin((jj+1)*k*theta)-sin(jj*k*theta))/sin(theta);
	gammas(k+1)= -2*Area*(LB1*dot(vv,MU1)...
	             +dot(ww,-cc.*(T1xx-T1yy)-ss.*(T1xy_full)));                     
end

% Resulting alphas,betas,gammas and eigenvalues as intervals
% before computing error estimates

lambdas1 = 0.5*(alphas+betas-sqrt((alphas-betas).^2+4*gammas.^2));
lambdas2 = 0.5*(alphas+betas+sqrt((alphas-betas).^2+4*gammas.^2));

disp("A_k, B_k, C_k, midpoint and radius of each interval");
for k=1:n
	fprintf("k=%2d\n",k);
	fprintf("mid  A_k=%7.3f B_k=%7.3f C_k=%7.3f mu_1=%7.3f mu_2=%7.3f\n",...
	         alphas(k).mid,betas(k).mid,gammas(k).mid,lambdas1(k).mid,lambdas2(k).mid); 
	fprintf("rad  A_k=%6.1e B_k=%6.1e C_k=%6.1e mu_1=%6.1e mu_2=%6.1e\n",...
	         alphas(k).rad,betas(k).rad,gammas(k).rad,lambdas1(k).rad,lambdas2(k).rad);    
	disp(" ");      
end

% ================= A priori estimates ================

h    = 1/intval(m);      % corresponding h
if nargin>2
	% set alternate h if needed for testing purposes
    h   = intval(hh);
end

% H^1 extension constant from the n-gon to R^2
% Lemma 3.2 in [1]
if n<6
	C_ext = intval(4);
else
	C_ext = sqrt(4+24*cos(theta)^2);
end

% constants for estimations: Morley finite element computations
% see Appendix A in [1]
% see the code regarding Morley assembly for the computations

% load from file
try
	load ConstantsP1;
catch
	OptimalC = [0,0,0.8499,0.4912,0.3697,0.3200,0.3146,0.3107,0.3104,0.3127];
end


% Choose the constant corresponding to the isosceles triangles
% with top angle 2pi/n 
C1 = OptimalC(n);

% estimates from the Xuefeng Liu article
% estimate for lambda 1
difflam = lam1^3*C1^2/(1+C1^2*h^2*lam1^2)*h^2;

% increase interval for analytical value of lambda1 to include the 
% error estimate. When the continuous quantity is required use lam1, otherwise use the discrete LB1
lam1 = infsup(LB1.inf-difflam.sup,LB1.sup+difflam.sup);

% estimate for lambda2
difflam2 = lam2^3*C1^2/(1+C1^2*h^2*lam2^2)*h^2;

% increase interval for analytical value of lambda2 to include the 
% error estimate. When the continuous quantity is required use lam2, otherwise use the discrete LB2
lam2 = infsup(LB2.inf-difflam2.sup,LB2.sup+difflam2.sup);

% compare lambda 2 with the fourth eigenvalue on the ball = 26.31
% if smaller, then the second eigenvalue is enclosed

fprintf("       A priori bounds for lambda_2 [%.3f %.3f]\n",lam2.inf,lam2.sup);
fprintf("       Fourth eigenvalue unit disk  =26.31\n ");

% Interpolation estimates [Liu, Oishi]
% gradient of u and projection P(u) on the FE space
diffgradup = C1*h*lam1;
% L2 estimate between u and P(u)
diffL2up   = C1^2*h^2*lam1;

% first component estimate ||grad(u)-grad(u_h)|| - interpolation error
Cu1 = C1*lam1*h;

% In the following use (3.14) from Cances Dusson Maday Stamm Vohralik
% |grad(u_h)|^2-lam1 = ErrorGradL2^2-lam1*ErrorUL2^2
% difflam = ErrorGradL2^2-lam1*ErrorUL2^2<=ErrorGradL2^2+lam1*ErrorUL2^2
%
% Therefore if two estimates are known a third one follows
% Replace successively difflam with the estimate, if the latter is better
% 
% sometimes estimates are improved. If not, nothing changes.

% See Section 2.2 for notations
for i=1:5
	% estimate for p-bar 
	GradPbar = sqrt(LB2)/(LB2-LB1)*(difflam+LB1*diffL2up);
	L2Pbar   = 1/sqrt(LB2)*GradPbar;

	% bounds for alpha
	OneMinusAlpha = L2Pbar^2+diffL2up*(2+diffL2up);
	%LbAlpha       = 1-OneMinusAlpha;

	% second and third component
	Cu2 = OneMinusAlpha*sqrt(LB1); % estimate for second term
	Cu3 = GradPbar;           % estimate for third term
		  
	% Cu1 is largest, of order h. The other two are of order h^2     
		  
	% L2 estimate for gradient of u      
	Cgradu = Cu1+Cu2+Cu3;

	% L2 estimate for u
	Cu = lam1*C1^2*h^2+OneMinusAlpha+L2Pbar;

	% identity from Vohralik et al
	% replace one of them with the smaller one. 
	
	term1 = Cgradu^2+lam1*Cu^2;
	term2 = difflam;
	difflam = min(term1,term2);
	fprintf("Iteration %d: L2grad = %.3e | L2u = %.3e | difflam = %.3e\n",i,Cgradu.sup,Cu.sup,difflam.sup);
end

% estimate for H^1 norm
uH1 = sqrt(Cu^2+Cgradu^2);

% estimate for |f-f_h| in H^{-1}. 
% See Section 5 from [2], Proposition 5.6

diffh = 2*sqrt(2)/sin(theta)/sqrt(n)*Cgradu+1/sqrt(1+lam1)*(2/n*difflam+2*LB1/n*Cu);

% put all parameters in one place 
opts.C1      = C1;
opts.lam1    = lam1;
opts.lam2    = lam2;
opts.LB1     = LB1;
opts.LB2     = LB2;
opts.h       = h;
opts.difflam = difflam;
opts.Cgradu  = Cgradu;
opts.Cu      = Cu;

% error when computing int(T0)(dx(u)^2) and (dy(u)^2)
% Estimate for the gradient error + Cauchy-Schwarz + symmetry
errorab = Cgradu/n*sqrt(lam1+LB1);

% error bound for integral of \partial_r u on a ray
% See Remark 5.5 from [2]
if mod(n,2) ==0
	intuseg = intuseg+0.5*Cgradu; 
else
    k = (n-1)/2;
	intuseg = intuseg+sqrt((k+1)/(2*n))*Cgradu;
end

% bounds for US0 from Section 3 in [1] 
% see proof of Theorem 3.4
GradUS0 = lam2/(lam2-lam1)*intuseg;
US0     = 1/sqrt(lam2)*GradUS0;
c0      = lam1/4;  % L infinity bound on u1

% ineq:  X^2+c<=b+a*sqrt(X^2+c)
aa = C_ext^2/2/sqrt(2)*(sqrt(lam1^2+lam1));
cc = GradUS0^2;
bb = (lam1^2*US0^2+0.5*c0^2)+cc;
 
Q = 0.5*(aa+sqrt(aa^2+4*bb));
D2Using = sqrt(Q^2-cc);
opts.D2Using = D2Using;

% =========================================
% estimates needed for problems of the form
% a(U,v) = (f_reg,v)+(f_sing,v)
% U orthogonal to the first eigenfunction
% Theorem 3.6 in [1], Section 5 in [2]
% =========================================

% see [2] Section 5 for details regarding the
% practical estimates below. Straightforward most of the time

% data for U1 
% estimate coming from |(a tensor b)c|<=|a||b||c|
f1Hm1        = 2/sin(theta)*sqrt(2*lam1/n);
% finer alternative 
f1Hm1        = 1/sin(theta)*sqrt(2*lam1/n)+1/sin(theta)*sqrt(aaa+errorab);

f1reg        = 2*lam1/sin(theta)*sqrt(2/n)+2*lam1/n; % to check if the last term can be removed
% finer alternative
f1reg        = 2*sqrt(2*lam1)/sin(theta)*sqrt(aaa+errorab)+2*lam1/n;

% Sum of coefficients before the norm of U0
f1sing       = 4/tan(theta);
% ===============================

% ||f1_h-f1||
diffh1       = diffh;
% ==================================

[gradUV1,L2UV1,gradV1,V1L2,diffgradtildeV1,diffL2tildeV1,gradUh1,difftV1Uh,U1L2] = estimate(f1Hm1,f1reg,f1sing,diffh1,opts);

% data for U2
% estimate coming from |(a tensor b)c|<=|a||b||c|
f2Hm1        = 2/sin(theta)*sqrt(2*lam1/n);

% better alternative
f2Hm1        = 2*sqrt(2)/tan(theta)*sqrt(bbb+errorab)+sqrt(2*lam1/n);


f2reg        = 2*lam1/sin(theta)*sqrt(2/n);
% alternative
f2reg        = 2*sqrt(2)/sin(theta)*sqrt(bbb+errorab);

% Sum of coefficients before the norm of U0 
f2sing       = 2;
% ===============================

% |f2-f2h|
diffh2       = diffh;%2/sin(theta)*sqrt(2/n)*Cgradu+4*sqrt(lam1)/(n*sin(theta)*sqrt(1+lam1))*(Cgradu+sqrt(lam1)*Cu);
% ==================================

[gradUV2,L2UV2,gradV2,V2L2,diffgradtildeV2,diffL2tildeV2,gradUh2,difftV2Uh,U2L2] = estimate(f2Hm1,f2reg,f2sing,diffh2,opts);

% intval array for storing the results
res = intval(zeros(2*n,1));

% structures for storing the results
% FEM interval results
AsFEM = alphas;
BsFEM = betas;
CsFEM = gammas;
l1FEM = lambdas1;
l2FEM = lambdas2;

% Final results
AsFin = intval(zeros(n,1));
BsFin = intval(zeros(n,1));
CsFin = intval(zeros(n,1));
l1Fin = intval(zeros(n,1));
l2Fin = intval(zeros(n,1));

% A priori estimates
AsEstimate = intval(zeros(n,1));
BsEstimate = intval(zeros(n,1));
CsEstimate = intval(zeros(n,1));

% loop for adding a priori estimates to the finite element
% approximations for the Hessian eigenvalues 
% See [1] Theorem 2.1 for the formulas and 
% [2] Section 4 for proofs
for k=0:n-1
	a = alphas(k+1);
   	b = betas(k+1);
   	c = gammas(k+1);
   	
   	fprintf("a.mid = %f  |  rad = %f \n",a.mid,a.rad);
   	fprintf("b.mid = %f  |  rad = %f \n",b.mid,b.rad);
   	fprintf("c.mid = %f  |  rad = %f \n",c.mid,c.rad);
   
   	lambda2 = 0.5*(a+b+sqrt((a-b)^2+4*c^2));
   	lambda1 = 0.5*(a+b-sqrt((a-b)^2+4*c^2));
   	
   	res(2*k+1) = lambda1;
   	res(2*k+2) = lambda2;
	
	fprintf("Lambda1: %f\n",lambda1.mid);
    fprintf("Lambda2: %f\n",lambda2.mid);

	% estimate for A_k ===============================
	
	falphaHm1        = lam1*sqrt((1+cos(k*theta))/(1+lam1))+sqrt(lam1*(1-cos(k*theta)))/sin(theta);
	falphareg        = sqrt(1+cos(k*theta))*lam1+sqrt(2)/sin(theta)*sqrt(1-cos(k*theta))*lam1; % to check if the last term can be removed
	
	jj = 0:(n-1);
	coef = 2*(1-cos(k*theta))*cos(theta)/sin(theta)*sum(abs(cos(jj*k*theta)));
		
	falphasing       = coef;
	diffhalpha       = sqrt((1+cos(k*theta))/(1+lam1))*(difflam+lam1*Cu)+sqrt(1-cos(k*theta))/sin(theta)*Cgradu;
	
	[gradUValpha,L2UValpha,gradValpha,ValphaL2,diffgradtildeValpha,diffL2tildeValpha,gradUhalpha,difftValphaUh,UalphaL2] = estimate(falphaHm1,falphareg,falphasing,diffhalpha,opts);
	
	% See Section 3.1 in [1]
	% Remark: Vtilde is a projection of V on the orthogonal of u{1,h}
	% therefore the L^2 and gradient norm for Vtilde will be smaller than those for V
	
	fprintf("Estimate A_k   k=%d: ",k);
	% first  term
	erralpha = gradUV1.*gradUValpha;%+lam1*ValphaL2*L2UV1+lam1*U1L2*L2UValpha;
	% second term
	erralpha = erralpha+gradValpha*diffgradtildeV1+gradV1*diffgradtildeValpha+difflam*V1L2*ValphaL2+lam1*ValphaL2*diffL2tildeV1+lam1*V1L2*diffL2tildeValpha;
	% third  term
	erralpha = erralpha+gradV1*difftValphaUh+gradUhalpha*difftV1Uh;

	if(k==0)
		% for k=0 the discrete and exact value of A_0 is 0
		erralpha = intval(0);
	end
	
	fprintf(" %f \n",erralpha.sup);

	% estimate for B_k ===============================

	fbetaHm1        = lam1*cos(theta)/sin(theta)*sqrt((1-cos(k*theta))/(1+lam1))+sqrt(lam1*(1-cos(k*theta)))/sin(theta); 
	fbetareg        = cos(theta)/sin(theta)*sqrt(1-cos(k*theta))*lam1+sqrt(2)/sin(theta)*sqrt(1-cos(k*theta))*lam1; 
	
	jj = 0:(n-1);
	coef = 2*(1-cos(k*theta))*sum(abs(cos(jj*k*theta)));
		
	fbetasing       = coef;
	
	% U2 is odd w.r.t y, singular part for B_k is even w.r.t y
	% therefore the singular component is zero analytically and numerically
	% it can be ignored in the estimates! 
	fbetasing = 0; 
	diffhbeta       = cos(theta)/sin(theta)*sqrt((1-cos(k*theta))/(1+lam1))*(difflam+lam1*Cu)+sqrt((1-cos(k*theta)))/sin(theta)*Cgradu;
	
	[gradUVbeta,L2UVbeta,gradVbeta,VbetaL2,diffgradtildeVbeta,diffL2tildeVbeta,gradUhbeta,difftVbetaUh,UbetaL2] = estimate(fbetaHm1,fbetareg,fbetasing,diffhbeta,opts);
	
	fprintf("Estimate B_k   k=%d: ",k);

	% See Section 3.1 in [1]
	% first term
	errbeta = gradUV2.*gradUVbeta;
	
	% second term
	errbeta = errbeta+gradVbeta*diffgradtildeV2+gradV2*diffgradtildeVbeta+difflam*V2L2*VbetaL2+lam1*VbetaL2*diffL2tildeV2+lam1*V2L2*diffL2tildeVbeta;
	
	% third term
	errbeta = errbeta+gradV2*difftVbetaUh+gradUhbeta*difftV2Uh;
	
	fprintf(" %f \n",errbeta.sup);
	
	% estimate for C_k ===============================
	fgamma1Hm1        = lam1*cos(theta)/sin(theta)*sqrt((1-cos(k*theta))/(1+lam1))+sqrt(lam1*(1-cos(k*theta)))/sin(theta); 
	fgamma1reg        = cos(theta)/sin(theta)*sqrt(1-cos(k*theta))*lam1+sqrt(2)/sin(theta)*sqrt(1-cos(k*theta))*lam1; % to check if the last term can be removed
	
	jj = 0:(n-1);
	coef = 2*(1-cos(k*theta))*sum(abs(sin(jj*k*theta)));
		
	fgamma1sing       = coef;
	
	% U1 is even w.r.t y, singular part for f_C is odd w.r.t y
	% therefore the singular component is zero analytically and numerically
	% it can be ignored in the estimates! 
	fgamma1sing = 0;  
	diffhgamma1       = cos(theta)/sin(theta)*sqrt((1-cos(k*theta))/(1+lam1))*(difflam+lam1*Cu)+sqrt(1-cos(k*theta))/sin(theta)*Cgradu;
	
	[gradUVgamma1,L2UVgamma1,gradVgamma1,Vgamma1L2,diffgradtildeVgamma1,diffL2tildeVgamma1,gradUhgamma1,difftVgamma1Uh,Ugamma1L2] = estimate(fgamma1Hm1,fgamma1reg,fgamma1sing,diffhgamma1,opts);
	
	fprintf("Estimate C_k   k=%d: ",k);
    % first term
	errgamma1 = gradUVgamma1*gradUV1;
	
	% second term
	errgamma1 = errgamma1+gradVgamma1*diffgradtildeV1+gradV1*diffgradtildeVgamma1+difflam*V1L2*Vgamma1L2+lam1*Vgamma1L2*diffL2tildeV1+lam1*V1L2*diffL2tildeVgamma1;
	
	% third term
	errgamma1 = errgamma1+gradV1*difftVgamma1Uh+gradUhgamma1*difftV1Uh;

	fprintf(" %f \n",errgamma1.sup);
	
	% Combine all a priori errors regarding Hessian eigenvalues
	% See formulas in Theorem 2.1 in [1]
	erralpha  = 2*Area*erralpha+2*n*(1-cos(k*theta))/sin(theta)*errorab;
	errbeta   = 2*Area*errbeta +2*n*(1-cos(k*theta))/sin(theta)*errorab;
	errgamma1 = 2*Area*errgamma1;
	
	fprintf("Err alpha %f\n",erralpha.sup);
	fprintf("Err beta  %f\n",errbeta.sup);
	fprintf("Err gamma %f\n",errgamma1.sup);
	
	% for k=1, n-1 alpha=beta=gamma
	% take best estimate among the given ones
	if or(k==1,k==n-1)
		err = min([erralpha.sup,errbeta.sup,errgamma1.sup]);
		err = intval(err);
		erralpha = err;
		errbeta  = err;
		errgamma1= err;
	end
	
	rada = intval(a.rad)+intval(erralpha.sup);
	radb = intval(b.rad)+intval(errbeta.sup);
	radc = intval(c.rad)+intval(errgamma1.sup);
	
	ia = midrad(a.mid,rada.sup);
	ib = midrad(b.mid,radb.sup);
	ic = midrad(c.mid,radc.sup);
		
	% formulas for the Hessian eigenvalues
	% Theorem 2.1 in [1], Section 4 in [2]
	ilam1 = 0.5*(ia+ib-sqrt((ia-ib)^2+4*ic^2));
	ilam2 = 0.5*(ia+ib+sqrt((ia-ib)^2+4*ic^2));
	
	AsFin(k+1) = ia;
	BsFin(k+1) = ib;
	CsFin(k+1) = ic;
	l1Fin(k+1) = ilam1;
	l2Fin(k+1) = ilam2;
	
	AsEstimate(k+1) = erralpha.sup;
	BsEstimate(k+1) = errbeta.sup;
	CsEstimate(k+1) = errgamma1.sup;
	
	display(' ');
	fprintf("===============\n");
	fprintf("Interval alpha:   %f %f rad=%f\n",ia.inf,ia.sup,ia.rad);
	fprintf("Interval beta :   %f %f rad=%f\n",ib.inf,ib.sup,ib.rad);
	fprintf("Interval gamma:   %f %f rad=%f\n",ic.inf,ic.sup,ic.rad);
	fprintf("===============\n");
	fprintf("Interval Lambda1: %f %f\n",ilam1.inf,ilam1.sup);
    fprintf("Interval Lambda2: %f %f\n",ilam2.inf,ilam2.sup);
	fprintf("===============\n");
	
	res(2*k+1) = ilam1;
	res(2*k+2) = ilam2;
	disp(" ");
end

[~,I] = sort(res.mid);
res = res(I);

fprintf("Number of positive eigenvalues = %d\n",sum(res.inf>0));

if sum(res.inf>0)==(2*n-4)
	fprintf("Proof of local minimality succeeded!\n");
end

% if a different mesh size is considered, 
% estimate the required number of degrees of freedom
if nargin>3
   display("Estimated number of points");
   fprintf("new meshp %d\n",floor(1/hh)+1);
   nmeshp = floor(1/hh)+1;
end

fprintf("Degrees of Freedom (full mesh)%d\n", 1+n*m*(m+1)/2);

struc.AsFin = AsFin;
struc.BsFin = BsFin;
struc.CsFin = CsFin;
struc.l1Fin = l1Fin;
struc.l2Fin = l2Fin;
struc.AsEstimate = AsEstimate;
struc.BsEstimate = BsEstimate;
struc.CsEstimate = CsEstimate;

struc.AsFEM = AsFEM;
struc.BsFEM = BsFEM;
struc.CsFEM = CsFEM;
struc.l1FEM = l1FEM;
struc.l2FEM = l2FEM;

struc.rig_comp = rig_comp;  % int(T_i) grad(U)*grad(u_1)
struc.mas_comp = mas_comp;

save_string = ['./Results/CompPolya_',num2str(n),'_',num2str(m),'.mat'];

% Save the result for further processing
save(save_string,'struc');

function [gradUV,L2UV,gradV,VL2,diffgradtildeV,diffL2tildeV,gradUh,difftVUh,UL2] = estimate(fHm1,freg,fsing,diffh,o)

	% construct function to evaluate estimates from Theorem 3.6. 
	% this is applied two times for the two functions appearing
	% in the formulas 
	% A_k = ... + a(U^1,W^{A_k})
	% B_k = ... + a(U^2,W^{B_k})
	% C_k = ... + a(U^1,W^{C_k})
	
	% fsing is the coefficient for the singular part multiplying
	% the norm of U_{S_0} in Theorem 3.4
	C1 = o.C1;
	lam1      = o.lam1; % enclosure for continuous eigenvalues
	lam2      = o.lam2;
	LB1       = o.LB1;  % enclosure for discrete eigenvalues
	LB2       = o.LB2;
	h         = o.h;
	difflam   = o.difflam;
	Cgradu    = o.Cgradu;
	Cu        = o.Cu;
	D2Using   = o.D2Using;
	uH1       = sqrt(Cu^2+Cgradu^2);

	gradU = sqrt(lam2*(1+lam2))/(lam2-lam1)*fHm1; 
	UL2   = 1/sqrt(lam2)*gradU;

	% ||U-V|| in H1 regular and singular parts

	% Estimate for singular part: Lemma 3.5 in [1]
	gradUVsing = (C1*h)*sqrt(2)*fsing*D2Using;

	% ===========================================
	% Classical L^2 estimate
	gradUVreg  = C1*h*(lam1*UL2+freg);

	% L2 estimates
	L2UVsing   = C1*h*gradUVsing;
	L2UVreg    = C1*h*gradUVreg; 

	% estimates for ||U-V||
	gradUV = gradUVsing+gradUVreg;
	L2UV   = L2UVsing+L2UVreg;

	% estimates for V
	gradV      = sqrt(lam1)*UL2+sqrt(1+1/lam1)*fHm1;
	VL2        = 1/sqrt(lam1)*gradV;

	% estimates for tilde V
	diffgradtildeV = sqrt(LB1)*(L2UV+VL2*Cu);
	diffL2tildeV   = L2UV+VL2*Cu;   

	gradUh         = sqrt(LB2*(LB2+1))/(LB2-LB1)*(fHm1+diffh);
	difftVUh       = sqrt(LB2)/(LB2-LB1)*(difflam*UL2+LB1*L2UV+sqrt(1+LB2)*diffh);




