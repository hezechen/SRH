%% data points extrapolated from figure 
a=[...
1.05    0.4
1.4     0.41
1.2     0.448
1.25	0.455
1.3 	0.45
1.6     0.633
2.65    0.732
1.25	0.54
2.15	0.62
2.5     0.65
2.05	0.665
2.68	0.7
3.2     0.73
1.95    0.658
2.23	0.66
2.5 	0.74
2.4     0.75
3.25    0.77
2.48	0.775];
%% trend line deffinition
dlm = fitlm(a(:,2).^(1.49),a(:,1),'Intercept',false);
cof=round(dlm.Coefficients.Estimate,2);
rsq=round(dlm.Rsquared.Ordinary,2);
x=(.4:.01:.8);
calc=cof.*(x.^(1.49));
%% Remake figure a
figure
scatter(a(1:7,2),a(1:7,1),[],[1,.7,0],'filled','s');
axis([.3, .8, 0,4])
hold on
scatter(a(8:13,2),a(8:13,1),[],'b','filled','^');
scatter(a(14:end,2),a(14:end,1),[],'g','filled');
plot(x,calc,'k')
legend('OVX','NFR','CON',[num2str(cof),'*(BVF)^{1.49}'])
xlabel('Bone Volume Fraction (mm^3/mm^3)')
ylabel('Shear Modulus (GPa)')
title('Recreation of figure 3a')
text(.325,3.75,['G=',num2str(cof),'*BVF^{1.49} ','R^2=',num2str(rsq)])
hold off