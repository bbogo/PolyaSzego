function TestingFormulas(n)

theta = 2*pi/n;

for k=1:(n-1)
	for j=0:n-1
	
		MA = [-sin((2*j+1)*theta) cos((2*j+1)*theta); cos((2*j+1)*theta) sin((2*j+1)*theta)];
		
		MC = [-cos((2*j+1)*theta) -sin((2*j+1)*theta); -sin((2*j+1)*theta) cos((2*j+1)*theta)];
	
		M1 = (cos((k+1)*k*theta)-cos(j*k*theta) 
		                                                      
		[v,d] = eig(M1)
		
		pause
	end
end
