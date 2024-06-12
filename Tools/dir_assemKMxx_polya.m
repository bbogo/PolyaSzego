function K = dir_assemKMxx_polya(p,t,opts)
%  dir_assemKMxx_polya(p,t,opts)
%
%  Assembly of rigidity matrix for integrals dx(phi_i)*dx(phi_j)
%  Analytical computations evaluated using intervals

nt                                       = size(t,1);
np                                       = size(p,1);


arM = opts.arM;
n   = opts.n;
Pi = intval('pi');
val = -sin(2*Pi/n)/2;

  it1                                    = t(:,1);
  it2                                    = t(:,2);
  it3                                    = t(:,3);
  
  nhalf = size(t,1)/2;
  c3                                     = intval(zeros(nt,1));
  c3(1:nhalf)                            = val;
  c1                                     = intval(zeros(nt,1));

  c2                                     = intval(zeros(nt,1));
  c2((nhalf+1):end)                      = val;
  
  K                                      = sparse(it1,it2,c3,np,np);
  K                                      = K+sparse(it2,it3,c1,np,np);
  K                                      = K+sparse(it3,it1,c2,np,np);
  K                                      = K+K.';
  K                                      = K+sparse(it1,it1,-c2-c3,np,np);
  K                                      = K+sparse(it2,it2,-c3-c1,np,np);
  K                                      = K+sparse(it3,it3,-c1-c2,np,np);

  



