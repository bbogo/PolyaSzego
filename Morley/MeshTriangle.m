function [pts,tri,Inside,ends] = MeshSliceAngle(a,b,m)

% n gives the regular polygon the slice comes from 
% m gives the number of points per ray for the division

% angle at center

pl = 1; % plotting

h = 1/m;
hx = 1/m*a;
hy = 1/m*b;

% number of points
% (m+1), m, ..., 1
npts = (m+1)*(m+2)/2;
ntri = m*(m+1)/2;
ntri2 = m*(m-1)/2;

% points (X x 2)
pts = zeros(npts,2);
tri = zeros(ntri,3);
tri2 = zeros(ntri2,3);

pts(1:(m+1),1) = (0:m)*h;


starts = 1+cumsum((m+1):-1:2);
starts = [1 starts];

ends = cumsum((m+1):-1:1);

Inside = setdiff(1:npts,ends);

oldpos = 1;
pos    = m+2;
for i=1:m
	nx = m+1-i; % number of points per line
	pts(pos:pos+nx-1,1) = pts(oldpos:oldpos+nx-1,1)+hx;
	pts(pos:pos+nx-1,2) = pts(oldpos:oldpos+nx-1,2)+hy;
	
	pos = pos+nx;
	oldpos = oldpos+(nx+1);
end

pos = 1;
for i=1:m
	nt = m+1-i; % number of tri per line
	tri(pos:pos+(nt-1),1) = starts(i):starts(i)+(nt-1);
	tri(pos:pos+(nt-1),2) = 1+(starts(i):starts(i)+(nt-1));
	tri(pos:pos+(nt-1),3) = nt+1+(starts(i):starts(i)+(nt-1));
	pos = pos+nt;
end

pos = 1;
for i=2:m
	nt = m+1-i; % number of tri per line
	tri2(pos:pos+(nt-1),2) = starts(i):starts(i)+(nt-1);
	tri2(pos:pos+(nt-1),1) = 1+(starts(i):starts(i)+(nt-1));
	tri2(pos:pos+(nt-1),3) = -nt-1+(starts(i):starts(i)+(nt-1));
	pos = pos+nt;
end

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


tri = [tri; tri2];




