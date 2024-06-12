function K = dir_assemKMyy_polya(p,t,opts)

%  dir_assemKMyy_polya(p,t,opts)
%
%  Assembly of rigidity matrix for integrals dy(phi_i)*dy(phi_j)
%  Analytical computations evaluated using intervals

nt                                       = size(t,1);
np                                       = size(p,1);


arM = opts.arM;
n   = opts.n;
Pi = intval('pi');
Theta = 2*Pi/n;
val1 = 0.5*cos(Theta)*tan(Theta/2);  % 0.1122...
val2 = 1/(2*tan(Theta));             % 0.16245...
val3 = 1/2*tan(Theta/2);             % 0.36327..

nhalf = size(t,1)/2;

  it1                                    = t(:,1);
  it2                                    = t(:,2);
  it3                                    = t(:,3);
  
  
  c3                                     = intval(zeros(nt,1));
  c3(1:nhalf)                            = val1;
  c3((nhalf+1):end)                      = val3;
  
  c1                                     = intval(zeros(nt,1));
  c1(1:nhalf)                            = -val2;
  c1((nhalf+1):end)                      = val2;
  
  c2                                     = intval(zeros(nt,1));
  c2(1:nhalf)                            = -val3;
  c2((nhalf+1):end)                      = -val1;
  
  K                                      = sparse(it1,it2,c3,np,np);
  K                                      = K+sparse(it2,it3,c1,np,np);
  K                                      = K+sparse(it3,it1,c2,np,np);
  K                                      = K+K.';
  K                                      = K+sparse(it1,it1,-c2-c3,np,np);
  K                                      = K+sparse(it2,it2,-c3-c1,np,np);
  K                                      = K+sparse(it3,it3,-c1-c2,np,np);

  

