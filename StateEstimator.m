%State estimation
function [iter,V]=StateEstimator(Y,loadPower,VtruePrim,R,rdnPQ, rdnV, PMUnodes)

loadPowerNoisy(:,1)=loadPower(:,1);
loadPowerNoisy(:,[2,3])=loadPower(:,[2,3]).*(1+rdnPQ);
loadPowerNoisy(loadPowerNoisy<0)=0;
Vnoisy=VtruePrim.*(1+rdnV);
linePowerTransferNoisy=[];
%% Defining measurement data
% [z,zType]=getMeasurements(Vtrue, linePowerTransfer, loadPower);%for testing the accuracy of WLS method without any noise
[z,zType]=getMeasurements(Vnoisy, linePowerTransferNoisy, loadPowerNoisy, PMUnodes);%WLS with noisy data

%% performing SE
iter_max=10;
threshold=1e-7;%stopping criteria: WLS stopes when delta_x is smaller than threshold or iteration number>iter_max
[V,iter]=WLS_SE(Y,z,zType,R,iter_max,threshold,VtruePrim);