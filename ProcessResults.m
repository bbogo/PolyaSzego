function ProcessResults(n)
%
% Process results for a given n and a given sequence of mesh sizes
% PolyaHessIntervalU should be run before on every one of the cases

% Save images or not: uses 'export_fig'
saving = 1;

% list of mesh parameters m 
ms = 200:50:600;

% Various initialization
radFEM = zeros(size(ms)); % radii FEM eigenvalues
radFin = zeros(size(ms)); % radii final eigenvalues

Aest = zeros(size(ms));

LamDn = zeros(size(ms));
LamUp = zeros(size(ms));

Err   = zeros(size(ms));

AllEig = zeros(2*n-4,length(ms));
AllErr = zeros(2*n-4,length(ms));

% load all information from computed data
for ii = 1:length(ms)
	m = ms(ii)
	str = ['./Results/CompPolya_',num2str(n),'_',num2str(m),'.mat']
	load(str)
	radFEM(ii) = max(struc.AsFEM.rad)
	radFin(ii) = max(struc.AsFin.rad)
	Aest(ii)   = mean(struc.AsEstimate.sup)
	struc.l1FEM.mid
	struc.l2FEM.mid
	All =  [struc.l1Fin(3:end-1).mid;
	        struc.l2Fin(2:end).mid];
	AllE = [struc.l1Fin(3:end-1).rad;
	        struc.l2Fin(2:end).rad];            
	[All,I] = sort(All)
	AllE = AllE(I);            
	AllEig(:,ii) = All;
	%LamDn(ii) = 
	AllErr(:,ii)  = AllE;
	%pause
end

set(groot,'defaulttextinterpreter','latex');  
set(groot, 'defaultAxesTickLabelInterpreter','latex');  
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'defaultLegendFontSize',24);
set(groot, 'defaultLegendFontWeight','bold');




figure(1)
hh = 1./ms;

convrate(radFEM)
convrate(Aest)
%radFEM(1)/hh(1)^(-5)*hh.^(-5)
loglog(hh,[radFEM; Aest; 
           ],'LineWidth',2)
legend('FEM interval radius','A priori estimate','Location','best')
set(gca,'FontSize',16);
xlabel('$h=1/m$','interpreter','latex');
title('Estimates: interval size, \emph{a priori}','interpreter','latex');

axis tight
if(saving==1)
eval(['export_fig EstimatesFEM_Er_',num2str(n),' -r300']);
end

size(AllEig)
size(hh)



figure(3)
clf
hold on
for i=1:size(AllEig,1)
	errorbar(ms,AllEig(i,:),AllErr(i,:),'.','LineWidth',2,'MarkerSize',20);
end
plot(ms,[zeros(size(ms))],'b','LineWidth',2)
set(gca,'FontSize',16);
xlabel('$m$','interpreter','latex');
title('Certified intervals: non-zero eigenvalues')

if(saving==1)
eval(['export_fig Intervals_All_',num2str(n),' -r300'])
end


figure(4)
clf
hold on
for i=1:1
	errorbar(ms,AllEig(i,:),AllErr(i,:),'.','LineWidth',2,'MarkerSize',20);
end
plot(ms,[zeros(size(ms))],'b','LineWidth',2)
set(gca,'FontSize',16);
xlabel('$m$','interpreter','latex');
title('Certified intervals: smallest non-zero eigenvalue')
if(saving==1)
eval(['export_fig Intervals_Smallest_',num2str(n),' -r300'])
end

function q = convrate(v)

v = fliplr(v)

n = length(v);

dv = diff(v);
for i=1:n-3
	q(i) = log(abs((v(i+3)-v(i+2))/(v(i+2)-v(i+1))))...
	      /log(abs((v(i+2)-v(i+1))/(v(i+1)-v(i))));
end




