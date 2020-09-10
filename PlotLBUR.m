
clc;
clear all;

tic();
    
[KAPPA,KAPPA_B,UR,URG,ITERELAPSED_BICG,TOL_1,TOL_2,TOL_3,TOL_4,UR_FOR_ERROR_BICG,UR_FOR_ERROR_GMRES,TOLG_1,TOLG_2,TOLG_3,TOLG_4] = UncerRatioPlotter(100);


toc();
%uncertainty plotter for bicg - uncertainty in rel_error and
%error_estimator in predicting relative error for BiCG

figure(1)
loglog(KAPPA,UR,'bo')
hold on
loglog(KAPPA,UR_FOR_ERROR_BICG,'ko');
xlabel('Forward condition number (K(A,x))')
ylabel('Uncertainty in estimation of relative error')
title('Uncertainty of relative residual and error estimator in estimation of relative error for BiCG')
legend('Uncertainty of relative residual','Uncertainty of error estimator')


%uncertainty plotter for gmres - uncertainty in rel_error and
%error_estimator in predicting relative error for BiCG

figure(2)
loglog(KAPPA_B,URG,'ro')
hold on
loglog(KAPPA_B,UR_FOR_ERROR_GMRES,'ko');
xlabel('Backward condition number (K(A,b))')
ylabel('Uncertainty in estimation of relative error')
title('Uncertainty of relative residual and error estimator in estimation of relative error for GMRES')
legend('Uncertainty of relative residual','Uncertainty of error estimator')

%sorting UR and iterations elapsed as keyed to each other

% [UR_sorted, UR_sorted_order] = sort(UR);
% [URG_sorted, URG_sorted_order] = sort(URG);
% ITERELAPSED_BICG_sorted = ITERELAPSED_BICG(UR_sorted_order);


[UR_sorted, UR_sorted_order] = sort(UR);
TOL_1_sorted = TOL_1(UR_sorted_order);
TOL_2_sorted = TOL_2(UR_sorted_order);
TOL_3_sorted = TOL_3(UR_sorted_order);
TOL_4_sorted = TOL_4(UR_sorted_order);
% TOL_5_sorted = TOL_5(UR_sorted_order);
%TOL_6_sorted = TOL_6(UR_sorted_order);

%removing elements in which algorithm didnt go desired tolerance 
keyOfNonZeroElementsBicg = find(TOL_1_sorted==0);
TOL_1_sorted(keyOfNonZeroElementsBicg)=[];TOL_2_sorted(keyOfNonZeroElementsBicg)=[];TOL_3_sorted(keyOfNonZeroElementsBicg)=[];TOL_4_sorted(keyOfNonZeroElementsBicg)=[];
UR_sorted(keyOfNonZeroElementsBicg)=[];

[URG_sorted, URG_sorted_order] = sort(URG);
TOLG_1_sorted = TOLG_1(URG_sorted_order);
TOLG_2_sorted = TOLG_2(URG_sorted_order);
TOLG_3_sorted = TOLG_3(URG_sorted_order);
TOLG_4_sorted = TOLG_4(URG_sorted_order);
% TOLG_5_sorted = TOLG_5(URG_sorted_order);
%TOLG_6_sorted = TOLG_6(URG_sorted_order);

keyOfNonZeroElementsGmres = find(TOLG_1_sorted==0);
TOLG_1_sorted(keyOfNonZeroElementsGmres)=[];TOLG_2_sorted(keyOfNonZeroElementsGmres)=[];TOLG_3_sorted(keyOfNonZeroElementsGmres)=[];TOLG_4_sorted(keyOfNonZeroElementsGmres)=[];
URG_sorted(keyOfNonZeroElementsGmres)=[];

%window = ceil(max(size(UR_sorted))/10);

points = 5;

UR_averaged = movMean(UR_sorted,points);
URG_averaged  = movMean(URG_sorted,points);

TOL_1_averaged = movMean(TOL_1_sorted,points);
TOL_2_averaged = movMean(TOL_2_sorted,points);
TOL_3_averaged = movMean(TOL_3_sorted,points);
TOL_4_averaged = movMean(TOL_4_sorted,points);
% TOL_5_averaged = movMean(TOL_5_sorted,window);
%TOL_6_averaged = movMean(TOL_6_sorted,window);

TOLG_1_averaged = movMean(TOLG_1_sorted,points);
TOLG_2_averaged = movMean(TOLG_2_sorted,points);
TOLG_3_averaged = movMean(TOLG_3_sorted,points);
TOLG_4_averaged = movMean(TOLG_4_sorted,points);
% TOLG_5_averaged = movMean(TOLG_5_sorted,window);
%TOLG_6_averaged = movMean(TOLG_6_sorted,window);

figure(3)
semilogx(UR_averaged,TOL_1_averaged);
hold on;
semilogx(UR_averaged,TOL_2_averaged,'g');
semilogx(UR_averaged,TOL_3_averaged,'r');
semilogx(UR_averaged,TOL_4_averaged,'black');
% semilogx(UR_averaged,TOL_5_averaged,'cyan');
%loglog(UR_averaged,TOL_6_averaged,'yellow');
xlabel('Uncertainty of relative residual')
ylabel('Iterations saved/lost')
title('Iterations saved/lost with uncertainty of relative residual in estimation of relative error (BiCG)')
legend('10^{-1} tolerance','10^{-2} tolerance','10^{-3} tolerance','10^{-4} tolerance')%'10^{-5} tolerance')

figure(4)
%semilogx(URG_averaged,TOLG_1_averaged);
%hold on
semilogx(URG_averaged,TOLG_2_averaged,'g');
hold on;
semilogx(URG_averaged,TOLG_3_averaged,'r');
semilogx(URG_averaged,TOLG_4_averaged,'black');
% semilogx(URG_averaged,TOLG_5_averaged,'cyan');
%loglog(URG_averaged,TOLG_6_averaged,'yellow');
xlabel('Uncertainty of relative residual')
ylabel('Iterations saved/lost')
title('Iterations saved/lost with uncertainty of relative residual in estimation of relative error (GMRES)')
%legend('10^{-1} tolerance','10^{-2} tolerance','10^{-3} tolerance','10^{-4} tolerance')%,'10^{-5} tolerance')
legend('10^{-2} tolerance','10^{-3} tolerance','10^{-4} tolerance')%,'10^{-5} tolerance')

% figure(5)
% plot(UR_averaged,TOL_1_averaged);
% hold on;
% plot(UR_averaged,TOL_2_averaged,'g');
% plot(UR_averaged,TOL_3_averaged,'r');
% plot(UR_averaged,TOL_4_averaged,'black');
% % semilogx(UR_averaged,TOL_5_averaged,'cyan');
% %loglog(UR_averaged,TOL_6_averaged,'yellow');
% xlabel('Uncertainty of relative residual')
% ylabel('Iterations saved/lost')
% title('Iterations saved/lost with uncertainty of relative residual in estimation of relative error (BiCG)')
% legend('10^{-1} tolerance','10^{-2} tolerance','10^{-3} tolerance','10^{-4} tolerance')%'10^{-5} tolerance')
















% figure(3)
% [urbins,tolbins] = findMeanForEachTolerance(UR_sorted,TOL_1_sorted);
% loglog(urbins,tolbins,'go');
% hold on;
% [urbins,tolbins] = findMeanForEachTolerance(UR_sorted,TOL_2_sorted);
% loglog(urbins,tolbins,'bo');
% hold on;
% [urbins,tolbins] = findMeanForEachTolerance(UR_sorted,TOL_3_sorted);
% loglog(urbins,tolbins,'yo');
% hold on;
% [urbins,tolbins] = findMeanForEachTolerance(UR_sorted,TOL_4_sorted);
% loglog(urbins,tolbins,'yo');
% hold on;

% loglog(UR_sorted,TOL_1_sorted,'go');
% hold on 
% loglog(UR_sorted,TOL_2_sorted,'bo');
% hold on 
% loglog(UR_sorted,TOL_3_sorted,'yo');
% hold on
% loglog(UR_sorted,TOL_4_sorted,'co');
% hold on
% loglog(UR_sorted,TOL_5_sorted,'ko');
% hold on 
% loglog(UR_sorted,TOL_6_sorted,'ro');
% 
% figure(4)
% [urbins,tolbins] = findMeanForEachTolerance(URG_sorted,TOLG_1_sorted);
% loglog(urbins,tolbins,'go');
% hold on
% [urbins,tolbins] = findMeanForEachTolerance(URG_sorted,TOLG_2_sorted);
% loglog(urbins,tolbins,'bo');
% hold on
% [urbins,tolbins] = findMeanForEachTolerance(URG_sorted,TOLG_3_sorted);
% loglog(urbins,tolbins,'yo');
% hold on
% [urbins,tolbins] = findMeanForEachTolerance(URG_sorted,TOLG_4_sorted);
% loglog(urbins,tolbins,'yo');
% hold on 


% loglog(URG_sorted,TOLG_1_sorted,'go');
% hold on 
% loglog(URG_sorted,TOLG_2_sorted,'bo');
% hold on
% loglog(URG_sorted,TOLG_3_sorted,'yo');
% hold on
% loglog(URG_sorted,TOLG_4_sorted,'co');
% hold on
% loglog(URG_sorted,TOLG_5_sorted,'ko');
% hold on 
% loglog(URG_sorted,TOLG_6_sorted,'ro');


% hold on 
% loglog(KAPPA,LUR,'go')