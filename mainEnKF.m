%running the 33 bus system in OpenDss
clear
clc

%% parameters for different runs
delay=1;%number of semples that PMU data is taking to be processed and communicated to us
currentPhasor = false;%if you want to consider the current phasor readings from PMU
resol=60;%data resolution in sec. we take a PMU reading every resol seconds
PMUnumber = 5; %number of installed PMUs
sigmaPMU=.01*.05/3;%micro-PMU accuracy

%% setting up OpenDSS
text = pwd;%current path
dss_path = [text, '\master33Full.dss'];
disp(dss_path);

[DSSObj, flag] = DSSStartup(dss_path);
% the OpenDSS Object alows interacting with the COM interface
% flag is set when the object is successfully instantiated
if ~flag
    error('Fatal error!\n Failed to create the COM object to interface with OpenDSS');
end

% DSSText=DSSObj.text;
% DSSCircuit=DSSObj.ActiveCircuit;
% DSSText.Command='Set mode=snapshot';
% DSSText.Command='solve';
% DSSText.command='show voltages';

%% Getting Ybus matrix
%manual Ybus
Y=manualYMatrix33ph;

%% Defining loads over simulation period
%base case loads
Pref=[100;90;120;60;60;200;200;60;60;45;60;60;120;60;60;60;90;90;90;90;90;90;420;420;60;60;60;120;200;150;210;60];
Qref=[60;40;80;30;20;100;100;20;20;30;35;35;80;10;20;20;40;40;40;40;40;50;200;200;25;25;20;70;600;70;100;40];

secNodes=55;%number of load nodes in the secondary power system
dur=3600;% the pseudo meas. will be available for every dur seconds
%Assigne load curves to nodes (using ADRES data set)
%getting the required load indices for connecting houses to secondary nodes
%and getting predicted loads for the next day
factor=10;
% [P,Q,delayedP,delayedQ,PpMean,PpStd,QpMean,QpStd,PdiffStd,QdiffStd,Ppsi,Qpsi]=loadProcessAssignHalf1ph1(Pref,secNodes,dur,resol,delay,factor);
% [P,Q]=loadProcessAssign1ph1Week(Pref,secNodes,dur,resol,delay,factor);
caseName = ['newLoadF10D1R' num2str(resol) 'Dur3600'];
load(caseName)
if delay == 0
    delayedP = P;
    delayedQ = Q;
end
% PpMean=PpMean/3;
% PpStd=PpStd/3;
% PdiffStd=PdiffStd/3;
% QpMean=QpMean/3;
% QpStd=QpStd/3;
% QdiffStd=QdiffStd/3;

%% perunit values
sBase=100e3/3; %single phase power base
vBasePri=12.66e3/sqrt(3); %single phase primary network voltage base
vBaseSec=.416/sqrt(3); %single phase secondary network voltage base
yBase=sBase/vBasePri^2;
priNodes=33*3; %all 3 phases
Y=Y/yBase;
Ibase=sBase/vBasePri;
% secNodes=32*110;

%% initialization
%assume the following nodes have PMUs
PMUmap=[1,33,32,31,18,17,30,16,29,15,14,13,28,12,11,10,9,8,27,26,7,6,25,24,5,4,23,3,22,21,20,19,2];
PMUnodes=sort(PMUmap(1:PMUnumber+1));
% PMUnodes=[1,33,18,22,25,12,6,3,29,9,15];
% PMUnodes=[1,33,18,22,25,6];
% PMUnodes=[1,33,18,22,25,6,3,31,29,27,20,26,14,12,10,8,4,24,23,4,2];
% PMUnodes=[1,33,18];
iter=zeros(1,length(P(1,:)));
V_WLS=zeros(length(P(1,:)),priNodes);
MSE_WLS=zeros(1,length(P(1,:)));
%ensemble Kalman filter
L=500;%ensemble size
X=zeros(2*length(PpMean(:,1)),L);
hX=zeros(2*3*length(PMUnodes),L);
Xave=zeros(2*length(PpMean(:,1)),length(P(1,:)));
V_EKF=zeros(length(P(1,:)),priNodes);
MSE_EKF=zeros(1,length(P(1,:)));
%PMUs at all nodes
MSE_PMU=zeros(1,length(P(1,:)));

%loading loadflow results
% load('LoadFlowResultsR60.mat');

%% main loop
for k=1:length(P(1,:))
    disp("hour: "+ k);
%     %% setting up OpenDSS
    text = pwd;%current path
    dss_path = [text, '\master33Full.dss'];
    [DSSObj, flag] = DSSStartup(dss_path);
    % the OpenDSS Object alows interacting with the COM interface
    % flag is set when the object is successfully instantiated
    if ~flag
        error('Fatal error!\n Failed to create the COM object to interface with OpenDSS');
    end
    %% Running the power flow
    [Vtrue, lineCurrent, ~, ~]=runPFfull(DSSObj,P(:,k),Q(:,k));
    [Vdelayed, lineCurrentDelayed, ~, ~]=runPFfull(DSSObj,delayedP(:,k),delayedQ(:,k));

    % calculating load power at each primary bus
    % method 1: adding the load powers together
    PseudoLoadPower=zeros(32*3,3);
    PseudoLoadPower(:,1)=4:33*3;
    PseudoLoadPower(:,2)=PpMean(:,ceil(k*resol/dur));%pseudo measurment Ppseudo
    PseudoLoadPower(:,3)=QpMean(:,ceil(k*resol/dur));%pseudo measurment Qpseudo
    
    %% converting to per-unit
    Vtrue(1:priNodes)=Vtrue(1:priNodes)/vBasePri;
    Vtrue(priNodes+1:end)=Vtrue(priNodes+1:end)/vBaseSec*sqrt(3);
%     linePowerTransfer(:,[3,4])=linePowerTransfer(:,[3,4])/sBase;
    PseudoLoadPower(:,[2,3])=PseudoLoadPower(:,[2,3])/sBase*1000;%PseudoLoadPower is in kW
    VtruePrim=Vtrue(1:priNodes);%getting the primary nodes voltages
    disp("initial error: "+ sqrt(1/99*sum(abs(VtruePrim-repmat([1,exp(-2*pi*1i/3),exp(-4*pi*1i/3)],1,33)).^2)))
    VdelayedPrim=Vdelayed(1:priNodes)/vBasePri;
    
    
    %% adding noise
    rdnV=sigmaPMU*randn(1,length(VtruePrim));%adding noise to magnitude only
    for rr=1:length(rdnV)
        alpha=rand*2*pi;
        rdnV(rr)=rdnV(rr)*(cos(alpha)+1i*sin(alpha));
    end
    % calculate the max angle difference
    Vnoisy=VdelayedPrim+rdnV;%adding total vector error
    VnoisyAngleDiff=max(angle(Vnoisy)-angle(VdelayedPrim));
    if currentPhasor == true
        %per unit
        lineCurrentDelayed(:,3:8)=lineCurrentDelayed(:,3:8)/Ibase;
        %adding noise to current readings
        rdnI=sigmaPMU/sqrt(2)*randn(2*32,6);
        lineCurrentNoisy=lineCurrentDelayed(:,1:2);
        lineCurrentNoisy(:,3:8)=0;
        lineCurrentNoisy(:,3:8)=lineCurrentDelayed(:,3:8).*(1+rdnI);
    else
        lineCurrentNoisy = [];
    end
%     %load flow results
%     Vnoisy=AllVnoisy(k,:);
%     lineCurrentNoisy=AllLineCurrentNoisy(:,:,k);
%     if delay==1
%         VtruePrim=AllVtruePrim(k,:);
%     else
%         VtruePrim=AllVdelayedPrim(k,:);
%     end

    %covariance matrix of the error
    if currentPhasor == true
        R=diag([PpStd(:,ceil(k*resol/dur)).'.^2,QpStd(:,ceil(k*resol/dur)).'.^2,...%Pseudo measurements
            sigmaPMU^2*ones(1,3*2*length(PMUnodes)-1),...%PMU voltage measurement
            sigmaPMU^2*ones(1,3*2*length(PMUnodes))]);%PMU current measurement
    else
        R=diag([PpStd(:,ceil(k*resol/dur)).'.^2,QpStd(:,ceil(k*resol/dur)).'.^2,...%Pseudo measurements
            sigmaPMU^2*ones(1,3*2*length(PMUnodes)-1)]);
    end
    
    %error if we had PMUs at all nodes
    MSE_PMU(k)=(1/99*sum(abs(VtruePrim-Vnoisy).^2));
    
    
    %% Defining measurement data for WLS SE
    % [z,zType]=getMeasurements(Vtrue, linePowerTransfer, PseudoLoadPower);%for testing the accuracy of WLS method without any noise
    [z,zType]=getMeasurements(Vnoisy, lineCurrentNoisy, PseudoLoadPower, PMUnodes, currentPhasor);%WLS with noisy data
    
    %% performing WLS SE
    iter_max=10;
    threshold=1e-7;%stopping criteria: WLS stopes when delta_x is smaller than threshold or iteration number>iter_max
    exit=0;
    while exit<3
        [V_WLS(k,:),iter(k)]=WLS_SE(Y,z,zType,R,iter_max,threshold,VtruePrim);
        if iter(k)>1
            exit=3;
        else
            exit=exit+1;
        end
    end
    disp("iteration number: "+ iter(k));
    %result error
    MSE_WLS(k)=(1/99*sum(abs(VtruePrim-V_WLS(k,:)).^2));
    disp("WLS: RMS error in voltage phasors: "+ sqrt(MSE_WLS(k)));
    
    %% setting up OpenDSS for EKF
    text = pwd;%current path
    dss_path = [text, '\master33.dss'];
    [DSSObj, flag] = DSSStartup(dss_path);
    % the OpenDSS Object alows interacting with the COM interface
    % flag is set when the object is successfully instantiated
    if ~flag
        error('Fatal error!\n Failed to create the COM object to interface with OpenDSS');
    end

    
    %% Ensemble Kalman Filter
    if k==1
        %initial ensemble
        X(:,:)=[PpMean(:,k);QpMean(:,k)]+[PpStd(:,k).*randn(length(PpMean(:,1)),L);QpStd(:,k).*randn(length(PpMean(:,1)),L)];
    else
        %state integration
        X(:,:)=X(:,:)+[PdiffStd(:,1).*randn(length(PdiffStd(:,1)),L);QdiffStd(:,1).*randn(length(QdiffStd(:,1)),L)];
    end
    %% Pseudo measurement assimilation
    pseudoStep=ceil(k*resol/dur);
    PSI=diag([Ppsi(:,pseudoStep);Qpsi(:,pseudoStep)]);
    QModel=diag([PdiffStd(:,1).^2;QdiffStd(:,1).^2]);
    R=diag([PpStd(:,pseudoStep).^2;QpStd(:,pseudoStep).^2]);
    H=eye(length(PSI));
    %ensemble
    d=[PpMean(:,pseudoStep);QpMean(:,pseudoStep)];%pseudo measurement
    D=d+[PpStd(:,pseudoStep).*randn(32*3,L);QpStd(:,pseudoStep).*randn(32*3,L)];
    %intermediary matrices
    Hdelta=H-PSI*H;
    C=QModel*H.'*PSI.';
    Ddelta=D-PSI*D;
    Rdelta=R-PSI*R*PSI.'+PSI*H*QModel*H.'*PSI.';
    %update equations
    Xp=X(:,:);
    E_X=mean(Xp,2);
    covXp=(Xp-E_X)*(Xp-E_X).';
    E=Hdelta*covXp*Hdelta.'+Rdelta+Hdelta*C+C.'*Hdelta.';
    K=(covXp*Hdelta.'+C)*E^-1;
    Xu=Xp+K*(Ddelta-Hdelta*Xp);
    
    %% PMU assimilation
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %PMU voltage phasor assimilation
    PMUnodes1ph=sort([3*PMUnodes-2,3*PMUnodes-1,3*PMUnodes]);
    z_pmu=[abs(Vnoisy(PMUnodes1ph))';angle(Vnoisy(PMUnodes1ph))';z];
    % PMU ensemble
%     Z=z_pmu+sigmaPMU*randn(length(z_pmu),L);
    rdnZV=sigmaPMU*randn(length(PMUnodes1ph),L);
    for ll=1:L
        for rr=1:length(PMUnodes1ph)
            alpha=rand*2*pi;
            rdnZV(rr,ll)=rdnZV(rr,ll)*(cos(alpha)+1i*sin(alpha));
        end
    end
    Zrand=Vnoisy(PMUnodes1ph).'+rdnZV;
    Z=[abs(Zrand);angle(Zrand)];
    % now we need to do a load flow for each ensemble of Xu
    for l=1:L
        %load flow
        [Vxl,Ixl]=runPF(DSSObj,Xu(1:length(PpMean(:,1)),l),Xu(length(PpMean(:,1))+1:end,l));
        Vxl=Vxl/vBasePri;
        hX(:,l)=[abs(Vxl(PMUnodes1ph))';angle(Vxl(PMUnodes1ph))'];
    end
    %Augmented ensemble X_hat
    X_hat=[Xu;hX];
    %defining H_hat
    H_hat=[zeros(2*length(PMUnodes1ph),2*length(PpMean(:,1))),diag(ones(1,2*length(PMUnodes1ph)))];
    %calculating Kalman Filter K
    E_X=mean(Xu,2);
    E_HX=mean(H_hat*X_hat,2);
    E_Z=mean(Z,2);
    K=(Xu-E_X)*(H_hat*X_hat-E_HX).'*...
        ((H_hat*X_hat-E_HX)*(H_hat*X_hat-E_HX).'+(Z-E_Z)*(Z-E_Z).')^-1;
    %Calculating Xa at iteration k
    Xa=Xu+K*(Z-H_hat*X_hat);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if currentPhasor == true
        %PMU Current phasor assimilation
        PMUCurrentMap=[1,2:2:64];
        PMU_Iidx=PMUCurrentMap(PMUnodes);
        ZI_PMU=[lineCurrentNoisy(PMU_Iidx,3:2:7);lineCurrentNoisy(PMU_Iidx,4:2:8)];
        ZI_PMU=reshape(ZI_PMU.',[],1);
        Z=ZI_PMU.*(1+sigmaPMU/sqrt(2)*randn(length(ZI_PMU),L));
        % now we need to do a load flow for each ensemble of Xa
        for l=1:L
            %load flow
            [Vxl,Ixl]=runPF(DSSObj,Xa(1:length(PpMean(:,1)),l),Xa(length(PpMean(:,1))+1:end,l));
            Ixl(:,3:8)=Ixl(:,3:8)/Ibase;
            Ixl_PMU=[Ixl(PMU_Iidx,3:2:7);Ixl(PMU_Iidx,4:2:8)];
            hX(:,l)=reshape(Ixl_PMU.',[],1);
        end
    end
    %Augmented ensemble X_hat
    X_hat=[Xa;hX];
    %defining H_hat
    H_hat=[zeros(2*length(PMUnodes1ph),2*length(PpMean(:,1))),diag(ones(1,2*length(PMUnodes1ph)))];
    %calculating Kalman Filter K
    E_X=mean(Xa,2);
    E_HX=mean(H_hat*X_hat,2);
    E_Z=mean(Z,2);
    K=(Xa-E_X)*(H_hat*X_hat-E_HX).'*...
        ((H_hat*X_hat-E_HX)*(H_hat*X_hat-E_HX).'+(Z-E_Z)*(Z-E_Z).')^-1;
    %Calculating Xa at iteration k
    Xa=Xa+K*(Z-H_hat*X_hat);
    
    %% replacing X_k with Xa_k
    X(:,:)=Xa;
    %averaging the ensemble
    Xave(:,k)=mean(X(:,:),2);
    %do a final load flow to get the node voltages
    V_EKF(k,:)=runPF(DSSObj,Xave(1:length(PpMean(:,1)),k),Xave(length(PpMean(:,1))+1:end,k));
    V_EKF(k,:)=V_EKF(k,:)/vBasePri;
    MSE_EKF(k)=(1/99*sum(abs(VtruePrim-V_EKF(k,:)).^2));
    disp("EKF: RMS error in voltage phasors: "+ sqrt(MSE_EKF(k)));
end

disp("Average RMS Error PMU: "+ sqrt(mean(MSE_PMU)));
disp("Average RMS Error WLS: "+ sqrt(mean(MSE_WLS)));
disp("Average RMS Error EKF: "+ sqrt(mean(MSE_EKF)));

runName = ['runF10VI' num2str(currentPhasor) 'D' num2str(delay) 'R' num2str(resol) 'P' num2str(PMUnumber) 'S' num2str(sigmaPMU) '.mat'];
save(runName, 'MSE_PMU', 'MSE_WLS', 'MSE_EKF')




