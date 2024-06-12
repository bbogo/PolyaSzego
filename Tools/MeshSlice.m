function [pts,tri,Inside,ends] = MeshSlice(n,m)

% Construct the mesh of a triangular slice of the regular n-gon
% 
% n gives the regular polygon the slice comes from 
% m gives the number of points per ray for the division
%
% The mesh is constructed using floating point arithmetic for poitnts
% but the connectivity information is constructed manually
%
% Only the connectivity information is used in the validation process
% The assembly is done analytically. 

% plotting parameter: only use for testing, for small m
% leave it to 0 for efficiency
pl = 0; 


% angle at center
theta = 2*pi/n; 

% mesh parameters
h = 1/m;
hx = 1/m*cos(theta);
hy = 1/m*sin(theta);

% number of points
% (m+1), m, ..., 1
npts = (m+1)*(m+2)/2;
% number of triangles with first orientation
ntri = m*(m+1)/2;
% number of triangles with second orientation
ntri2 = m*(m-1)/2;

% arrays for points and triangle indices
pts = zeros(npts,2);
tri = zeros(ntri,3);
tri2 = zeros(ntri2,3);

% initial points on the horizontal ray
pts(1:(m+1),1) = (0:m)*h;

% start indices for horizontal lines
starts = 1+cumsum((m+1):-1:2);
starts = [1 starts];

% end indices for horizontal lines
ends = cumsum((m+1):-1:1);

% keep track of nodes inside the regular polygon
Inside = setdiff(1:npts,ends);

% iteratively construct the points
oldpos = 1;
pos    = m+2;
for i=1:m
	nx = m+1-i; % number of points per line
	pts(pos:pos+nx-1,1) = pts(oldpos:oldpos+nx-1,1)+hx;
	pts(pos:pos+nx-1,2) = pts(oldpos:oldpos+nx-1,2)+hy;
	
	pos = pos+nx;
	oldpos = oldpos+(nx+1);
end

% iteratively construct the triangles: first orientation
pos = 1;
for i=1:m
	nt = m+1-i; % number of tri per line
	tri(pos:pos+(nt-1),1) = starts(i):starts(i)+(nt-1);
	tri(pos:pos+(nt-1),2) = 1+(starts(i):starts(i)+(nt-1));
	tri(pos:pos+(nt-1),3) = nt+1+(starts(i):starts(i)+(nt-1));
	pos = pos+nt;
end

% iteratively construct the triangles: second orientation
pos = 1;
for i=2:m
	nt = m+1-i; % number of tri per line
	tri2(pos:pos+(nt-1),2) = starts(i):starts(i)+(nt-1);
	tri2(pos:pos+(nt-1),1) = 1+(starts(i):starts(i)+(nt-1));
	tri2(pos:pos+(nt-1),3) = -nt-1+(starts(i):starts(i)+(nt-1));
	pos = pos+nt;
end

% plotting, essentially only used for testing purposes
if pl==1
	clf
	hold on
	plot(pts(:,1),pts(:,2),'.','MarkerSize',10);
	plot(pts(Inside,1),pts(Inside,2),'.r','MarkerSize',25);
	plot(pts(ends,1),pts(ends,2),'.m','MarkerSize',20);


	%patch('faces',tri,'vertices',pts,'facecolor','red');
	%patch('faces',tri2,'vertices',pts,'facecolor','green');

	axis equal 
	axis off

	hold off

end

% regroup all triangles
tri = [tri; tri2];




