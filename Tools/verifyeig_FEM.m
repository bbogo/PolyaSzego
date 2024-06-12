function [L,X] = verifyeig_FEM(A,lambda,xs,B,opts)

% Eigenvalue verification based on verifyeig
%
% opts should contain a lower bound for the smallest
% eigenvalue of B:   opts.minEig
% 
%
%
% calls a verified linear system sovler based on residual


rndold = getround;
  if rndold
    setround(0)
  end
  
  [n,k] = size(xs);

  [~,I] = sort(sum(abs(xs),2));
  u = I(1:n-k);
  v = I(n-k+1:n);

   midA = mid(A);
   midB = mid(B);

    % one floating point iteration
    R = midA - lambda*midB;
    R(:,v) = -midB*xs;
    spy(R)
    
    %
    disp("Linear System");
    
    %y = R\(midA*xs-lambda*midB*xs);
    setup = struct('type','ilutp','droptol',1e-6);
	[L,U] = ilu(R,setup);    disp("Done LU");
    
    y = tfqmr(R,(midA*xs-lambda*midB*xs),1e-12,2000,L,U);
    %y = tfqmr(R,(midA*xs-lambda*midB*xs),1e-12,2000);

    
    

    xs(u,:) = xs(u,:) - y(u,:);
    lambda = lambda - sum(diag(y(v,:)))/k;

    R = midA - lambda*midB;
    R(:,v) = -midB*xs;
    
    %R = inv( R );
    C = A - intval(lambda)*B;
    
    Z = - verifylss(R,( C * xs ));
    C(:,v) = -B*xs;
    %C = speye(n) - R * C;
    Y = Z;
    Eps = 0.1*mag(Y)*hull(-1,1) + midrad(0,realmin);
    m = 0;
    mmax = 15 * ( sum(sum(mag(Z(v,:))>.1)) + 1 );
    ready = 0;
    while ( ~ready ) && ( m<mmax ) && ( ~any(isnan(Y(:))) )
      m = m+1;
      X = Y + Eps;
      XX = X;
      XX(v,:) = 0;
      %Y = Z + C*X + R*((B*XX*X(v,:)));
      Y1 =  (R-C)*X+B*XX*X(v,:);
      Y = Z+verifylss(R,Y1);
      ready = all(all(in0(Y,X)));
    end



  if ready
    M = mag(Y(v,:));                         % eigenvalue correction
    if length(v)==1                          % one eigenvalue
      L = midrad(lambda,M);
    else                                     % eigenvalue cluster
      [Evec,Eval] = eig(M);
      [~,index] = max(abs(diag(Eval)));
      Perronx = abs(Evec(:,index));
      setround(1);
      rad = max( ( M*Perronx ) ./ Perronx );   % upper bound for Perron root
      setround(0)
      L = cintval(midrad(lambda,rad));
    end
    Y(v,:) = 0;
    X = xs + Y;
  else
    % disp('no inclusion achieved')
    X = intval(NaN*ones(size(xs)));
    L = intval(NaN);
  end

    
  if rndold
    setround(rndold)
  end
