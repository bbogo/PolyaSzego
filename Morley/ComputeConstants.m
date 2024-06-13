function ComputeConstants
%
% Computes upper bounds for P1 interpolation error
% on triangles with angles 2pi/n, pi/2-pi/n, pi/2-pi/n
% for n=3,4,5,6,7,8,9,10
%
% these triangles appear in the meshes of the regular
% n gon used in our paper
%
% B. Bogosel, D. Bucur
% Polygonal Faberk-Krahn inequality: local minimality via
% validated computing

Pi = intval('pi');
C  = intval(zeros(1,10));
m  = 10; % mesh parameter

for n=3:10
	a = cos(2*Pi/n);
	b = sin(2*Pi/n);
	
	C(n) = P1InterpolationConstant(a,b,m);
	fprintf("n = %d | Upper bound %.4f\n",n,C(n).sup);
end

OptimalC = C;

save('ConstantsP1.mat','OptimalC');
