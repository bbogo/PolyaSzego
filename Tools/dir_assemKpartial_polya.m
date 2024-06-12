function [Kxx,Kyy,Kxy,MU] = dir_assemKpartial_polya(p,tri,opts)
%
% dir_assemKpartial_polya(p,tri,opts)
%  
% assembles partial rigidity matrices 
% to compute integrals of the form: integral(slice)(d{x,y}(u) d{x,y}(U^i))
% 
% all computations are analytic: See section 4 in the paper
%
% Kxx is a structure containing the rigidity matrix on slice i in K{i}
% Kyy and Kxy are constructed in a similar way
%


np                                       = size(p,1);
dim                                      = size(p,2);

arM = opts.arM;
n   = opts.n;
reg = opts.reg;
Pi = intval('pi');
val = -sin(2*Pi/n)/2;
Theta = 2*Pi/n;

i_x = -intval(1);
i_y = -tan(Theta/2);

j_x = intval(1);
j_y = -1/tan(Theta);

k_x = intval(0);
k_y = 1/sin(Theta);

Ar  = sin(Theta)/2;

% loop on slices
for jj=0:n-1
  actr   = find(reg==jj);  % select triangles in current region
  t      = tri(actr,:);    % select triangles in current region
  nt     = size(t,1);

  ang    = jj*Theta;


  it1    = t(:,1);
  it2    = t(:,2);
  it3    = t(:,3);
  
  % matrix Kxx
  %c3 (i,j)
  c3  = cos(ang)^2*i_x*j_x-cos(ang)*sin(ang)*(i_x*j_y+i_y*j_x)+sin(ang)^2*i_y*j_y;
  
  %c1 (j,k)                                
  c1  = cos(ang)^2*k_x*j_x-cos(ang)*sin(ang)*(k_x*j_y+k_y*j_x)+sin(ang)^2*k_y*j_y;

  %c2  (i,k)                               
  c2  = cos(ang)^2*i_x*k_x-cos(ang)*sin(ang)*(i_x*k_y+i_y*k_x)+sin(ang)^2*i_y*k_y;
  
  
  K                                      = sparse(it1,it2,c3,np,np);
  K                                      = K+sparse(it2,it3,c1,np,np);
  K                                      = K+sparse(it3,it1,c2,np,np);
  K                                      = K+K.';
  K                                      = K+sparse(it1,it1,-c2-c3,np,np);
  K                                      = K+sparse(it2,it2,-c3-c1,np,np);
  K                                      = K+sparse(it3,it3,-c1-c2,np,np);

	Kxx{jj+1} = Ar*K; 
  
  % matrix Kyy
  %c3 (i,j)
  c3  = sin(ang)^2*i_x*j_x+cos(ang)*sin(ang)*(i_x*j_y+i_y*j_x)+cos(ang)^2*i_y*j_y;
  
  %c1 (j,k)
  c1  = sin(ang)^2*k_x*j_x+cos(ang)*sin(ang)*(k_x*j_y+k_y*j_x)+cos(ang)^2*k_y*j_y;

  %c2  (i,k)
  c2  = sin(ang)^2*i_x*k_x+cos(ang)*sin(ang)*(i_x*k_y+i_y*k_x)+cos(ang)^2*i_y*k_y;
  
  
  K                                      = sparse(it1,it2,c3,np,np);
  K                                      = K+sparse(it2,it3,c1,np,np);
  K                                      = K+sparse(it3,it1,c2,np,np);
  K                                      = K+K.';
  K                                      = K+sparse(it1,it1,-c2-c3,np,np);
  K                                      = K+sparse(it2,it2,-c3-c1,np,np);
  K                                      = K+sparse(it3,it3,-c1-c2,np,np);

  Kyy{jj+1} = Ar*K; 

  % matrix Kxy
  %c3 (i,j)
  c3  = 2*sin(ang)*cos(ang)*i_x*j_x+(cos(ang)^2-sin(ang)^2)*(i_x*j_y+i_y*j_x)-2*sin(ang)*cos(ang)*i_y*j_y;
  
  %c1 (j,k)
  c1  = 2*sin(ang)*cos(ang)*k_x*j_x+(cos(ang)^2-sin(ang)^2)*(k_x*j_y+k_y*j_x)-2*sin(ang)*cos(ang)*k_y*j_y;

  %c2  (i,k)
  c2  = 2*sin(ang)*cos(ang)*i_x*k_x+(cos(ang)^2-sin(ang)^2)*(i_x*k_y+i_y*k_x)-2*sin(ang)*cos(ang)*i_y*k_y;
  
  
  K                                      = sparse(it1,it2,c3,np,np);
  K                                      = K+sparse(it2,it3,c1,np,np);
  K                                      = K+sparse(it3,it1,c2,np,np);
  K                                      = K+K.';
  K                                      = K+sparse(it1,it1,-c2-c3,np,np);
  K                                      = K+sparse(it2,it2,-c3-c1,np,np);
  K                                      = K+sparse(it3,it3,-c1-c2,np,np);

  Kxy{jj+1} = Ar*K; 
  
  aod                                    = arM/12; % Off diagonal element
  ad                                     = 2*aod; % Diagonal element
  M                                      = sparse(it1,it2,aod,np,np);
  M                                      = M+sparse(it2,it3,aod,np,np);
  M                                      = M+sparse(it3,it1,aod,np,np);
  M                                      = M+M.';
  M                                      = M+sparse(it1,it1,ad,np,np);
  M                                      = M+sparse(it2,it2,ad,np,np);
  M                                      = M+sparse(it3,it3,ad,np,np);
  
  MU{jj+1}   = M;
end

