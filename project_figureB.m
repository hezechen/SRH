OVX = [8.26 11.21 9.08 10.03 9.845 15.49 19.23];
Bf_OVX = [0.4 0.41 0.448 0.455 0.45 0.633 0.732];

NFR = [9.4 17.08 20.16 14.8 21.03 24.21];
Bf_NFR = [0.54 0.62 0.65 0.665 0.7 0.73];

CON = [12.38 14.15 19.63 15.77 25.26 16.25];
Bf_CON = [0.658 0.66 0.74 0.75 0.77 0.775];

t = [OVX NFR CON];
bf = [Bf_OVX Bf_NFR Bf_CON];
trnd = fitlm(bf.^1.46,t,'Intercept',false);
cof=round(trnd.Coefficients.Estimate,2);
rsq = trnd.Rsquared.Ordinary;
Bfrac = 0.4:0.01:0.8;
tau = 30.44*(Bfrac.^1.46);

figure;
plot(Bf_OVX,OVX,'LineStyle','none','Marker','s','MarkerEdgeColor','[0.8500 0.3250 0.0980]','MarkerFaceColor','[0.8500 0.3250 0.0980]')
hold on
plot(Bf_NFR,NFR,'LineStyle','none','Marker','^','MarkerEdgeColor','blue','MarkerFaceColor','blue')
plot(Bf_CON,CON,'LineStyle','none','Marker','o','MarkerEdgeColor','[0.4660 0.6740 0.1880]','MarkerFaceColor','[0.4660 0.6740 0.1880]')
plot(Bfrac,tau,'k')

xlim([0.3 0.9])
ylim([0 30])

xlabel('Ct.BV/TV (mm^3/mm^3)')
ylabel('tau (MPa)')

legend('OVX', 'NFR', 'CON')

text(.325,25,['G=',num2str(cof),'*BVF^{1.46} ','R^2=',num2str(rsq)])