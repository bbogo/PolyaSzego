function [K,M,W] = dir_assemKM(p,t,Ibord,ww)

% penalisation pour Dirichlet
% K   = K + sparse(Ibord,Ibord,1e12*ones(1,length(Ibord)),size(K,1),size(K,2));
% comment trouver Ibord?????????

%
%  Assemblage P1 2D/3D
%  [K,M,W] = sb_assembKM(p,t,ww)
%
%  NB : ww tient compte d un poid attache aux integrales. Il faut que numel(ww) = size(p,1)
%
%

global solverindex


nt                                       = size(t,1);
np                                       = size(p,1);
dim                                      = size(p,2);
ar                                       = abs(sb_simpvol(p,t(:,1:(dim+1))));
arK                                      = ar;
arM                                      = ar;
W                                        = 0;
% ajout des poids
if(exist('ww','var'))
  if(numel(ww)==size(p,1))
    ww                                   = sum(ww(t(:,1:(dim+1))),2)/(dim+1);
    arM                                  = ar(:).*ww(:);
  elseif(numel(ww)~=size(t,1))
    error('Problem of dimension with the weights in sb_assembKM.m');
  end
  if(dim==3)
    error('l ajout des poids pour les integrales de masse et de rigidite n est que prevu en 2D pour l instant (cf sb_assembKM.m)')
  end
end


if(size(t,2)<=4) % 2D ou 2D en 3D
  if(size(p,2)>2) % cas beltrami
    disp('mesh beltrami???')
    pause
    [K,M]                                = sb_laplacian_beltrami(p,t);
    return
  end

  it1                                    = t(:,1);
  it2                                    = t(:,2);
  it3                                    = t(:,3);
  [g1x,g1y,g1z,g2x,g2y,g2z,g3x,g3y,g3z]  = sb_pdetrg(p,t,arK);
  c3                                     = (g1x.*g2x+g1y.*g2y+g1z.*g2z).*arK;
  c1                                     = (g2x.*g3x+g2y.*g3y+g2z.*g3z).*arK;
  c2                                     = (g3x.*g1x+g3y.*g1y+g3z.*g1z).*arK;
  %c3
  %c1
  %c2
  %pause
  K                                      = sparse(it1,it2,c3,np,np);
  K                                      = K+sparse(it2,it3,c1,np,np);
  K                                      = K+sparse(it3,it1,c2,np,np);
  K                                      = K+K.';
  K                                      = K+sparse(it1,it1,-c2-c3,np,np);
  K                                      = K+sparse(it2,it2,-c3-c1,np,np);
  K                                      = K+sparse(it3,it3,-c1-c2,np,np);

  if nargin>2
  K = K + sparse(Ibord,Ibord,1e12*ones(1,length(Ibord)),size(K,1),size(K,2));
  end
  aod                                    = arM/12; % Off diagonal element
  ad                                     = 2*aod; % Diagonal element
  M                                      = sparse(it1,it2,aod,np,np);
  M                                      = M+sparse(it2,it3,aod,np,np);
  M                                      = M+sparse(it3,it1,aod,np,np);
  M                                      = M+M.';
  M                                      = M+sparse(it1,it1,ad,np,np);
  M                                      = M+sparse(it2,it2,ad,np,np);
  M                                      = M+sparse(it3,it3,ad,np,np);
  return
end



%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% version vectorialisee de ce qui suit (utilise uniquement en 3D)

% pour la matrice de rigidite
vertices                                 = p(t(:,1:(dim+1))',:);
vertices                                 = [ones(size(vertices,1),1),vertices]';
Mt                                       = eye(dim);
GG                                       = zeros(nt*(dim+1),dim);

if(dim==3)
  if(isempty(solverindex))
    solverindex                          = SolverIndex(4,1e5);
  end
end

for k = 1:dim
  %B                                      = invmat4x4(vertices')*[0;Mt(:,k)]; 
  if(isempty(solverindex))
    GG(:,k)                              = SolverNxN(vertices,repmat([0;Mt(:,k)],nt,1));
  else
    GG(:,k)                              = SolverNxN(vertices,repmat([0;Mt(:,k)],nt,1),solverindex);
  end
  %disp(num2str(nt))
  %max(abs(GG(:,k)-B(:)))
end
[IGG,JGG,SGG]                            = find(GG);
JGG                                      = JGG + (ceil(IGG/(dim+1)) - 1)*dim ;
GGsp                                     = sparse(IGG,JGG,SGG)';
Gsp                                      = GGsp'*GGsp;
[IGG,JGG,SGG]                            = find(Gsp);
JGG                                      = mod(JGG - 1,dim+1) + 1;
GG                                       = (2/prod(1:dim))*full(sparse(IGG,JGG,SGG));

% Essai trop lent
%vertices                                 = p(t(:,1:(dim+1))',:);
%vertices                                 = [ones(size(vertices,1),1),vertices]';
%Mt                                       = [zeros(1,dim);eye(dim)];
%groups                                   = mat2cell(vertices,size(vertices,1),(dim+1)*ones(size(vertices,2)/(dim+1),1))';
%GG2                                      = cell2mat(cellfun(@(x) x\Mt,groups,'UniformOutput',false));
%groups                                   = mat2cell(GG2,(dim+1)*ones(size(GG2,1)/(dim+1),1),size(GG2,2));
%GGt                                      = (2/prod(1:dim))*cell2mat(cellfun(@(x) x*transpose(x),groups,'UniformOutput',false));

art                                      = repmat(ar(:),1,dim+1)';
GG                                       = bsxfun(@times,GG,art(:));
if(dim==3)
  GG                                     = 3*GG;
end




% pour la matrice de masse
if(dim==2)
  pa                                     = p(t(:,2),:) - p(t(:,1),:);                                     
  pb                                     = p(t(:,3),:) - p(t(:,1),:);  
  Mdets                                  = abs(pa(:,1).*pb(:,2) - pa(:,2).*pb(:,1));
  %Mdets                                  =          p(t(:,2),1).*p(t(:,3),2)-p(t(:,2),2).*p(t(:,3),1);
  %Mdets                                  = Mdets - (p(t(:,1),1).*p(t(:,3),2)-p(t(:,1),2).*p(t(:,3),1));
  %Mdets                                  = Mdets + (p(t(:,1),1).*p(t(:,2),2)-p(t(:,1),2).*p(t(:,2),1));
  Mdets                                  = Mdets/24;
  Mdets                                  = repmat(Mdets(:),1,dim+1)';
  GM                                     = bsxfun(@times,repmat([2 1 1;1 2 1;1 1 2],nt,1),Mdets(:));
elseif(dim==3)
  pa                                     = p(t(:,2),:) - p(t(:,1),:);                                     
  pb                                     = p(t(:,3),:) - p(t(:,1),:);  
  pc                                     = p(t(:,4),:) - p(t(:,1),:);  
  Mdets                                  =         pa(:,1).*(pb(:,2).*pc(:,3) - pb(:,3).*pc(:,2));
  Mdets                                  = Mdets - pa(:,2).*(pb(:,1).*pc(:,3) - pb(:,3).*pc(:,1));
  Mdets                                  = Mdets + pa(:,3).*(pb(:,1).*pc(:,2) - pb(:,2).*pc(:,1));
  Mdets                                  = abs(Mdets)/120;
  Mdets                                  = repmat(Mdets(:),1,dim+1)';
  GM                                     = bsxfun(@times,repmat([2 1 1 1;1 2 1 1;1 1 2 1;1 1 1 2],nt,1),Mdets(:));
end




if(dim==2)
  It                                     = repmat(t(:,1:(dim+1)),1,dim+1)';
  It                                     = It(:);
  Jt                                     = t(:,[ones(1,3) 2*ones(1,3) 3*ones(1,3)])';
  Jt                                     = Jt(:);
  Kt                                     = GG';
  Kt                                     = Kt(:);
  KtB                                    = GM';
  KtB                                    = KtB(:);
elseif(dim==3)
  It                                     = repmat(t(:,1:(dim+1)),1,dim+1)';
  It                                     = It(:);
  Jt                                     = t(:,[ones(1,4) 2*ones(1,4) 3*ones(1,4) 4*ones(1,4)])';
  Jt                                     = Jt(:);
  Kt                                     = GG';
  Kt                                     = Kt(:);
  KtB                                    = GM';
  KtB                                    = KtB(:);
end

K                                        = sparse(It,Jt,Kt,np,np);
K = K + sparse(Ibord,Ibord,1e15*ones(1,length(Ibord)),size(K,1),size(K,2));
if(nargout>1)
  M                                      = sparse(It,Jt,KtB,np,np);
end




return

max(abs(M(:)-MS(:)))
max(abs(K(:)-KS(:)))

return









%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Inversion manuelled un ensemble de matrices 4x4
function B = invmat4x4(A)


nm                 = size(A,1)/4;
B                  = zeros(size(A));
for k=0:(nm-1)
  At               = A((4*k+1):(4*k+4),:);
  a                = At(1);
  b                = At(2);
  c                = At(3);
  d                = At(4);
  e                = At(5);
  f                = At(6);
  g                = At(7);
  h                = At(8);
  i                = At(9);
  j                = At(10);
  k                = At(11);
  l                = At(12);
  m                = At(13);
  n                = At(14);
  o                = At(15);
  p                = At(16);
  Bt               = [ f*k*p - f*l*o - g*j*p + g*l*n + h*j*o - h*k*n, e*l*o - e*k*p + g*i*p - g*l*m - h*i*o + h*k*m, e*j*p - e*l*n - f*i*p + f*l*m + h*i*n - h*j*m, e*k*n - e*j*o + f*i*o - f*k*m - g*i*n + g*j*m;
 b*l*o - b*k*p + c*j*p - c*l*n - d*j*o + d*k*n, a*k*p - a*l*o - c*i*p + c*l*m + d*i*o - d*k*m, a*l*n - a*j*p + b*i*p - b*l*m - d*i*n + d*j*m, a*j*o - a*k*n - b*i*o + b*k*m + c*i*n - c*j*m;
 b*g*p - b*h*o - c*f*p + c*h*n + d*f*o - d*g*n, a*h*o - a*g*p + c*e*p - c*h*m - d*e*o + d*g*m, a*f*p - a*h*n - b*e*p + b*h*m + d*e*n - d*f*m, a*g*n - a*f*o + b*e*o - b*g*m - c*e*n + c*f*m;
 b*h*k - b*g*l + c*f*l - c*h*j - d*f*k + d*g*j, a*g*l - a*h*k - c*e*l + c*h*i + d*e*k - d*g*i, a*h*j - a*f*l + b*e*l - b*h*i - d*e*j + d*f*i, a*f*k - a*g*j - b*e*k + b*g*i + c*e*j - c*f*i];
  dt               = a*f*k*p - a*f*l*o - a*g*j*p + a*g*l*n + a*h*j*o - a*h*k*n - b*e*k*p + b*e*l*o + b*g*i*p - b*g*l*m - b*h*i*o + b*h*k*m + c*e*j*p - c*e*l*n - c*f*i*p + c*f*l*m + c*h*i*n - c*h*j*m - d*e*j*o + d*e*k*n + d*f*i*o - d*f*k*m - d*g*i*n + d*g*j*m;
  Bt                   = Bt/dt;
  B((4*k+1):(4*k+4),:) = Bt;
end



return




%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% version non vectorielle
It2                                      = zeros(nt*(dim+1)^2,1);
Jt2                                      = zeros(nt*(dim+1)^2,1);
Kt2                                      = zeros(nt*(dim+1)^2,1);
Kt2B                                     = zeros(nt*(dim+1)^2,1);




for j = 1:nt
 
  valloc                                 = (1+(j-1)*(dim+1)^2):((dim+1)^2+(j-1)*(dim+1)^2);
  if(dim==2)
    It2(valloc)                          = [t(j,1:3),t(j,1:3),t(j,1:3)];
    Jt2(valloc)                          = [t(j,1),t(j,1),t(j,1),t(j,2),t(j,2),t(j,2),t(j,3),t(j,3),t(j,3)];
  elseif(dim==3)
    It2(valloc)                          = [t(j,1:4),t(j,1:4),t(j,1:4),t(j,1:4)];
    Jt2(valloc)                          = [t(j,1),t(j,1),t(j,1),t(j,1),t(j,2),t(j,2),t(j,2),t(j,2),t(j,3),t(j,3),t(j,3),t(j,3),t(j,4),t(j,4),t(j,4),t(j,4)];
  end
  vt                                     = sb_stima3(p(t(j,1:(dim+1)),1:dim),ar,j);

  Kt2(valloc)                            = vt(:)'; 
  if(nargout>1)
    vt                                  = sb_mass3(p(t(j,1:(dim+1)),1:dim));
    Kt2B(valloc)                        = vt(:)'; 
  end
end

max(abs(It-It2)),max(abs(Jt-Jt2)),max(abs(Kt-Kt2)),max(abs(KtB-Kt2B))















%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function M = sb_stima3(vertices,ar,j)


d   = size(vertices,2);
G   = [ones(1,d+1);vertices'] \ [zeros(1,d);eye(d)];




if(d==2)
  M = 2*ar(j)* G * G' / prod(1:d);
elseif(d==3)
  M = 6*ar(j)* G * G' / prod(1:d);
end

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function M = sb_mass3(vertices)


%
%
% Illustration des formules : 2D + syms x y;int(int((-x-y+1)*y,y,0,1-x),x,0,1);int(int((-x-y+1)^2,y,0,1-x),x,0,1)
%                             3D + syms x y z;int(int(int((-x-y-z+1)*x,z,0,1-x-y),y,0,1-x),x,0,1);syms x y z;int(int(int(y*x,z,0,1-x-y),y,0,1-x),x,0,1)
%

d   = size(vertices,2);
if(d==2)
  M = 1/24 * det([ones(1,d+1);vertices']) * [2 1 1;1 2 1;1 1 2];
else
  M = 1/120 * det([ones(1,d+1);vertices']) * [2 1 1 1;1 2 1 1;1 1 2 1;1 1 1 2];
end
