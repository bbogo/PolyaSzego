function v=sb_simpvol(p,t)

%
%  v = sb_simpvol(p,t)
%
% SIMPVOL Simplex volume.
%   si t tetra    --> volume
%   si t triangle --> aire 2D/3D
%   
%

switch(size(t,2)-1)
 case 1
  d12   = p(t(:,2),:)-p(t(:,1),:);
  v     = d12;
 case 2
  d12   = p(t(:,2),:)-p(t(:,1),:);
  d13   = p(t(:,3),:)-p(t(:,1),:);
  if(size(p,2)==3)
    C   = cross(d12,d13,2);
    v   = 1/2*sqrt(sum(C.^2,2));
  else
    v   = (d12(:,1).*d13(:,2)-d12(:,2).*d13(:,1))/2;
  end
 case 3
  d12   = p(t(:,2),:)-p(t(:,1),:);
  d13   = p(t(:,3),:)-p(t(:,1),:);
  d14   = p(t(:,4),:)-p(t(:,1),:);
  v     = abs(dot(cross(d12,d13,2),d14,2))/6;
 otherwise 
end
