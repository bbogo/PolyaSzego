function [K,M] = dir_assemKM_polya(p,t,opts)

%  dir_assemKM_polya(p,t,opts)
%
%  Assembly for P1 finite elements
%  the inputs p,t characterize a mesh for the regular 
%  polygon made of congruent triangles
%  p    : an array containing coordinate of points on rows
%  t    : an array containing indexes for triangle vertices on rows
%  opts : contains additional information
%
%  the output
%  K - rigidity matrix
%  M - mass matrix

nt                                       = size(t,1);
np                                       = size(p,1);


arM = opts.arM;   % area of a small triangle
n   = opts.n;     % number of vertices for the polygon

% computing constants: Proposition 4.1 of the article
% computations made using intervals in INTLAB
Pi = intval('pi');

val1 = -1/tan(2*Pi/n)/2;
val2 = -tan(2*pi/(2*n))/2;

% Explicit assembly for ridigity matrix

  it1                                    = t(:,1);
  it2                                    = t(:,2);
  it3                                    = t(:,3);

  c3                                     = val2;
  c1                                     = val1;
  c2                                     = val2;
  
  K                                      = sparse(it1,it2,c3,np,np);
  K                                      = K+sparse(it2,it3,c1,np,np);
  K                                      = K+sparse(it3,it1,c2,np,np);
  K                                      = K+K.';
  K                                      = K+sparse(it1,it1,-c2-c3,np,np);
  K                                      = K+sparse(it2,it2,-c3-c1,np,np);
  K                                      = K+sparse(it3,it3,-c1-c2,np,np);

% Explicit assembly for mass matrix

  aod                                    = arM/12; % Off diagonal element
  ad                                     = 2*aod; % Diagonal element
  M                                      = sparse(it1,it2,aod,np,np);
  M                                      = M+sparse(it2,it3,aod,np,np);
  M                                      = M+sparse(it3,it1,aod,np,np);
  M                                      = M+M.';
  M                                      = M+sparse(it1,it1,ad,np,np);
  M                                      = M+sparse(it2,it2,ad,np,np);
  M                                      = M+sparse(it3,it3,ad,np,np);

