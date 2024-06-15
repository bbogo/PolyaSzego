function PolyaSzegoRunAll(n,ms)
% Runs all test cases
%
% n  - number of vertices
% ms - mesh parameters (tested for 100<=m<=600)
 

if nargin<2
	% the computations might take a while...
	% default ones for making the plots in the article
	% if you want the computation to finish quicker
	% set the upper bound to 400 in the line below
	ms = 100:50:600
end


for i=1:length(ms)
	m = ms(i);
	fprintf("Run computations for n=%d m=%d\n",n,m);
	if m<=400
		PolyaHessInterval(n,m);
	else
		% sparse version for m large
		PolyaHessInterval_Sparse(n,m);
	end
	PolyaHessIntervalU(n,m)
end
