function res =  MeshMerge(M1,M2,m);

% MeshMerge(M1,M2,m)
%
% function for mesh merging 
% (constructed to work for the particular meshes considered here)
%
% The two meshes overlap only on a line
% therefore, the new mesh will be a union of all triangles
% 
% one only needs to find duplicate points and merge them
%
% in addition: we know that the first m nodes in M2 are
% in the intersection

% plotting, only used for testing
pl = 0;

p1 = M1.pts;
% first m nodes in p2 are in p1
p2 = M2.pts;

Inside1 = M1.Inside;
Inside2 = M2.Inside;

ends1 = M1.ends;
ends2 = M2.ends;

Index0 = M1.Index0;

t1 = M1.tri;
t2 = M2.tri;

% now fix triangles
% for each one of the first m nodes in M2 search the 
% corresponding node in M1

for i=1:(m)
	pp = p2(i,:); 
	dd = sqrt((p1(:,1)-pp(1)).^2+(p1(:,2)-pp(2)).^2);
	ii = find(dd<1e-8);
	if(isempty(ii)==0)
		t2(t2==i) = m-size(p1,1)+ii;
	end
end


% shift nodes in t2 with size of p2-m
t2 = t2+size(p1,1)-m;

% build merged mesh
pts = [p1;
       p2(m+1:end,:)];


tri = [t1;t2];


res.pts = pts;
res.tri = tri;
ends = union(ends1,ends2+size(p1,1)-m);

res.ends = ends;
Inside = setdiff(1:size(pts,1),ends);

res.Inside = Inside;
res.Index0 = [Index0  (m+1):((m+1)*(m)/2)];
% =================

if pl==1
	clf


	patch('faces',tri,'vertices',pts,'facecolor','red');
	%patch('faces',tri2,'vertices',pts,'facecolor','green');

	hold on
	plot(pts(Inside,1),pts(Inside,2),'.g','MarkerSize',20);
	plot(pts(ends,1),pts(ends,2),'*m','MarkerSize',20);

	axis equal 
	axis off
	hold off
end

