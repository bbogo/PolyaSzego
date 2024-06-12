function K = dir_assemKM(p,t,opts)

%  dir_assemKMxy_polya_minus(p,t,opts)
%
%  Assembly of rigidity matrix for integrals dx*dy+dy*dx on T_-
%  Analytical computations evaluated using intervals




nt                                       = size(t,1);
np                                       = size(p,1);


arM = opts.arM;  % Area of a triangle
n   = opts.n;
Pi = intval('pi');
Theta = 2*Pi/n;
val = (tan(Theta/2)*sin(Theta)-cos(Theta))/2;


  it1                                    = t(:,1);
  it2                                    = t(:,2);
  it3                                    = t(:,3);  
  
  c3                                     = intval(zeros(nt,1));
  c1                                     = intval(zeros(nt,1));
  c2                                     = intval(zeros(nt,1));
  
  c3                                     = 1/intval(2);
  c1                                     = -1/intval(2);
  c2                                     = val;
  
  
  
  K                                      = sparse(it1,it2,c3,np,np);
  K                                      = K+sparse(it2,it3,c1,np,np);
  K                                      = K+sparse(it3,it1,c2,np,np);
  K                                      = K+K.';
  K                                      = K+sparse(it1,it1,-c2-c3,np,np);
  K                                      = K+sparse(it2,it2,-c3-c1,np,np);
  K                                      = K+sparse(it3,it3,-c1-c2,np,np);



