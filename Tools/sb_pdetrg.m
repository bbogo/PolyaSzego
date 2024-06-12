function [g1x,g1y,g1z,g2x,g2y,g2z,g3x,g3y,g3z]=sb_pdetrg(p,t,ar)

%
%
% PDETRG Triangle geometry data.
%
%       [g1x,g1y,g1z,g2x,g2y,g2z,g3x,g3y,g3z]=PDETRG(P,T,AR) returns the  the gradient
%       components of the triangle base functions eventually projected on the tangent vector space
%
% ATTENTION cas Beltrami bugge (DONC non utilise)  
%
dim            = size(p,2);

% Corner point indices
a1             = t(:,1);
a2             = t(:,2);
a3             = t(:,3);

% Triangle sides
r23x           = p(a3,1)-p(a2,1);
r23y           = p(a3,2)-p(a2,2);
r31x           = p(a1,1)-p(a3,1);
r31y           = p(a1,2)-p(a3,2);
r12x           = p(a2,1)-p(a1,1);
r12y           = p(a2,2)-p(a1,2);

if(dim>2) % pour Beltrami
  r12z         = p(a2,3)-p(a1,3);
  r23z         = p(a3,3)-p(a2,3);
  r31z         = p(a1,3)-p(a3,3);
end

if(dim==2)
  g1x          = -0.5*r23y./ar;
  g1y          =  0.5*r23x./ar;
  g2x          = -0.5*r31y./ar;
  g2y          =  0.5*r31x./ar;
  g3x          = -0.5*r12y./ar;
  g3y          =  0.5*r12x./ar;
  g1z          = zeros(size(g1x));
  g2z          = zeros(size(g2x));
  g3z          = zeros(size(g3x));
else         
  % calcul des normales
  normalsf     = sb_compute_normals(p,t(:,1:3));
  centresf     = (p(t(:,1),:) + p(t(:,2),:) + p(t(:,3),:))/3;
  cf           = centresf + 1e-1*normalsf;
 
  dets         = (p(a1,1)-cf(:,1)).*((p(a2,2)-cf(:,2)).*(p(a3,3)-cf(:,3)) - (p(a2,3)-cf(:,3)).*(p(a3,2)-cf(:,2))) - ...
                 (p(a1,2)-cf(:,2)).*((p(a2,1)-cf(:,1)).*(p(a3,3)-cf(:,3)) - (p(a2,3)-cf(:,3)).*(p(a3,1)-cf(:,1))) + ...
                 (p(a1,3)-cf(:,3)).*((p(a2,1)-cf(:,1)).*(p(a3,2)-cf(:,2)) - (p(a2,2)-cf(:,2)).*(p(a3,1)-cf(:,1)));


  % gradient associes a la fonction de base du premier sommet de chaque triangle
  grads1       = crossp(p(a2,:)-cf,p(a3,:)-cf);
  grads1       = bsxfun(@rdivide,grads1,dets);
  
  % gradient associes a la fonction de base du second sommet de chaque triangle
  grads2       = -crossp(p(a1,:)-cf,p(a3,:)-cf);
  grads2       = bsxfun(@rdivide,grads2,dets);
  
  % gradient associes a la fonction de base du troisieme sommet de chaque triangle
  grads3       = crossp(p(a1,:)-cf,p(a2,:)-cf);
  grads3       = bsxfun(@rdivide,grads3,dets);
  
  
  % projection des gradients sur les plans tangents
  scalt        = sum(grads1.*normalsf,2);
  grads1       = grads1 - bsxfun(@times,normalsf,scalt);
  scalt        = sum(grads2.*normalsf,2);
  grads2       = grads2 - bsxfun(@times,normalsf,scalt);
  scalt        = sum(grads3.*normalsf,2);
  grads3       = grads3 - bsxfun(@times,normalsf,scalt);
  
  g1x          = grads1(:,1);
  g1y          = grads1(:,2);
  g1z          = grads1(:,3);
  
  g2x          = grads2(:,1);
  g2y          = grads2(:,2);
  g2z          = grads2(:,3);
  
  g3x          = grads3(:,1);
  g3y          = grads3(:,2);
  g3z          = grads3(:,3);
  
 
  
  %return
  centresf     = (p(t(:,1),:) + p(t(:,2),:) + p(t(:,3),:))/3;
  plotsurf(p,t,'FaceColor',[0.1 0.7 0]);hold on
  %quiver3(centresf(:,1),centresf(:,2),centresf(:,3),normalsf(:,1), ...
  %        normalsf(:,2),normalsf(:,3));
  centresf = centresf + 1e-2*normalsf;
  h=quiver3(centresf(:,1),centresf(:,2),centresf(:,3),grads1(:,1), ...
            grads1(:,2),grads1(:,3));set(h,'Color','r');
  h=quiver3(centresf(:,1),centresf(:,2),centresf(:,3),grads2(:,1), ...
            grads2(:,2),grads2(:,3));set(h,'Color','y');
  h=quiver3(centresf(:,1),centresf(:,2),centresf(:,3),grads3(:,1), ...
            grads3(:,2),grads3(:,3));set(h,'Color','m');
  axis equal;camlight
  mkbug
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function z = crossp(x,y)
% x and y are (m,3) dimensional
z      = x;
z(:,1) = x(:,2).*y(:,3) - x(:,3).*y(:,2);
z(:,2) = x(:,3).*y(:,1) - x(:,1).*y(:,3);
z(:,3) = x(:,1).*y(:,2) - x(:,2).*y(:,1);
