function res = PolyaMesh(n,m)

% PolyaMesh(n,m)
%
% Construct symmetric mesh of the regular polygon from slices
% not efficient for large m
%
% the mesh is constructed using floating point arithmetic, but
% only the connectivity information is used in the computations
% assembly of finite elements is done explicitly
%
% n : number of vertices of the regular n-gon
% m : number of division points on a ray
%
% output : a structure containing mesh information

% plotting parameter: only use for small m
% code may take a lot of time if plotting is allowed
pl=0;

% n gives the regular polygon the slice comes from 
% m gives the number of points per ray for the division

% angle at center
theta = 2*pi/n; 



% construct the mesh for a slice
[pts,tri,Inside,ends] = MeshSlice(n,m);


M1.pts = pts;
M1.tri = tri;
M1.Inside = Inside;
M1.ends   = ends;
M1.Index0 = 1:((m+1)*(m+2)/2);

M0 = M1;

% construct new mesh incrementally
% by rotating M1 arount the origin
% and merging vertices that overlap

for j=1:(n-1)
	
	R = [cos(j*theta) -sin(j*theta);
    	 sin(j*theta)  cos(j*theta)];

	M2.pts = (R*M0.pts')';
	M2.tri = tri;
	M2.Inside = M0.Inside;
	M2.ends   = M0.ends;

	M1 = MeshMerge(M1,M2,m+1);
end

pts = M1.pts;
tri = M1.tri;

Inside = M1.Inside;
Index0 = M1.Index0;


% for vertices 1 to m+1 find duplicates

dupl = [];

for i=1:m+1
	pp = pts(i,:); 
	dd = sqrt((pts(:,1)-pp(1)).^2+(pts(:,2)-pp(2)).^2);
	ii = find(dd<1e-8);
	if(isempty(ii)==0)
		if length(ii)>1
			dupl = [dupl; ii'];
		end
	end
end


% remove duplicate vertices
% remove from pts

% find indices of points without duplicates
simInd = setdiff(1:size(pts,1),dupl(:,2));

% replace in tri
for i=1:size(dupl,1)
	tri(tri==dupl(i,2))= dupl(i,1);
end

% relabel tri: tri(old) = new
newInd = 1:length(simInd);
J = zeros(size(pts,1),1);
J(simInd) = 1:length(simInd);
J(dupl(:,2)) = dupl(:,1);

tri = J(tri);
pts = pts(simInd,:);
Index0 = Index0(simInd);

% Inside nodes
Inside = unique(J(Inside));

if pl==1
	clf
	patch('Faces',tri,'Vertices',pts,'FaceColor','none','Edgecolor','k');
	patch('Faces',tri,'Vertices',pts,'FaceVertexCData',Index0(:),'FaceColor','inter','Edgecolor','k');

	% uncomment if you want to see the indices for
	% the extension from the slice to the whole polygon
	% text(pts(:,1),pts(:,2),num2str(Index0(:)));

	axis equal
	colorbar
end 

res.pts = pts;
res.tri = tri;
% indices for nodes inside the mesh
res.Inside = Inside;
res.Index0 = Index0;

