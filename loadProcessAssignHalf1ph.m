%redusing the time difference between samples to 0.5 seconds instead of 1s
%linear interapolation
function [P,Q,delayedP,delayedQ,PpMean,PpStd,QpMean,QpStd,PdiffStd,QdiffStd,Ppsi,Qpsi]=loadProcessAssignHalf1ph(Pref,secNodes,dur,resol,delay,factor)
%assigning households load data to each secondary node of the system
%the data will be extracted rom ADRES dataset
%we need to get instantanous PMU measurements and pseudo measurements
%which are forecasts for each period of time
dur=2*dur;
resol=2*resol;
% %loading ADRES data
% load('ADRES_Daten_120208.mat');
% %convert to total 3-phase PQ data
% %houseLoadData contains 420 load consumption data for 24 hours taken every
% %seconds
% for k=1:30%house index
%     for day=1:14
%         startIdx=(24*(day-1))*3600+1;
%         endIdx=(24*day)*3600;
%         data=[Data.PQ(startIdx:endIdx,6*k-1)+Data.PQ(startIdx:endIdx,6*k-3)+...
%               Data.PQ(startIdx:endIdx,6*k-5),Data.PQ(startIdx:endIdx,6*k)+...
%               Data.PQ(startIdx:endIdx,6*k-2)+Data.PQ(startIdx:endIdx,6*k-4)];
%         houseLoadDataOR.P(:,14*(k-1)+day)=data(:,1);
%         houseLoadDataOR.Q(:,14*(k-1)+day)=data(:,2);
%     end
% end
% 
% houseLoadDataOR.P=houseLoadDataOR.P.'/1000;%convert to kW- original time resolution
% houseLoadDataOR.Q=houseLoadDataOR.Q.'/1000;%convert to kVar- original time resolution
% NH=length(houseLoadDataOR.P(:,1));%number of original houses
% NS=length(houseLoadDataOR.P(1,:));%number of original samples
% for h=1:NH
%     houseLoadData.P(h,:)=reshape([houseLoadDataOR.P(h,:); zeros(1,NS)],[],1);
%     houseLoadData.Q(h,:)=reshape([houseLoadDataOR.Q(h,:); zeros(1,NS)],[],1);
%     for t=1:NS-1
%         houseLoadData.P(h,2*t)=(houseLoadData.P(h,2*t-1)+houseLoadData.P(h,2*t+1))/2;
%         houseLoadData.Q(h,2*t)=(houseLoadData.Q(h,2*t-1)+houseLoadData.Q(h,2*t+1))/2;
%     end
%     houseLoadData.P(h,2*NS)=(houseLoadData.P(h,2*NS-1)+houseLoadData.P(h,1))/2;
%     houseLoadData.Q(h,2*NS)=(houseLoadData.Q(h,2*NS-1)+houseLoadData.Q(h,1))/2;
% end
% 
% clear Data houseLoadDataOR
% 
% 
% 
% %creting new load traces
% % factor=2;%the new data set will be "factor" times biggr than the initial one
% 
% psize=size(houseLoadData.P);
% for k=1:factor
%     rnd=.1*randn(psize(1),24);
%     rnd=repmat(rnd,psize(2)/24,1);
%     rnd=reshape(rnd,psize(1),[]);
%     houseLoadData.P((k-1)*psize(1)+1:k*psize(1),:)=houseLoadData.P(1:psize(1),:).*(1+rnd);
%     houseLoadData.Q((k-1)*psize(1)+1:k*psize(1),:)=houseLoadData.Q(1:psize(1),:).*(1+rnd);
% end
%load all household data
load('HouseholdData.mat');

%now we need to extract the pseudo measurements for the given duration
allHouses=length(houseLoadData.P(:,1));% total number of household load profile we have
T=length(houseLoadData.P(1,:));% total time steps
housePmean=zeros(allHouses,T/dur);%mean
housePstd=zeros(allHouses,T/dur);%std
houseQmean=zeros(allHouses,T/dur);%mean
houseQstd=zeros(allHouses,T/dur);%std
for k=1:allHouses %house index
    for t=1:T/dur%chosen time period for doing the load assignment according to Pref
        housePmean(k,t)=mean(houseLoadData.P(k,(t-1)*dur+1:t*dur));
        housePstd(k,t)=std(houseLoadData.P(k,(t-1)*dur+1:t*dur));
        houseQmean(k,t)=mean(houseLoadData.Q(k,(t-1)*dur+1:t*dur));
        houseQstd(k,t)=std(houseLoadData.Q(k,(t-1)*dur+1:t*dur));
    end
end

houseIdx=zeros(32,1500);
for n=1:length(Pref)%node index
    Ptotal=0;
    count=1;
    while Ptotal<Pref(n)
        rnd=ceil(rand*allHouses);
        houseIdx(n,count)=rnd;
        count=count+1;
        Ptotal=Ptotal+housePmean(rnd,10);
    end
end

%variable aggregation level depending on the laod demand of the base case
loadIdx=cell(secNodes*32,1);%32 is the number of nodes with low voltage 
%network connected to them (aggregated primary loads)
%now we need to connect the loads to each low voltage secondary node
%we read the house number from houseIdx
PsecMean=zeros(length(Pref)*secNodes,length(housePmean(1,:)));%pseudo meas. mean
PsecStd=zeros(length(Pref)*secNodes,length(housePmean(1,:)));%pseudo meas. standard deviation
QsecMean=zeros(length(Pref)*secNodes,length(housePmean(1,:)));%pseudo meas. mean
QsecStd=zeros(length(Pref)*secNodes,length(housePmean(1,:)));%pseudo meas. standard deviation
for k=1:length(Pref)
    for h=1:length(housePmean(1,:))
        c=1;
        sn=1;%secondary node
        while houseIdx(k,c)>0
            if h==1
                temp=[loadIdx{(k-1)*secNodes+sn},houseIdx(k,c)];
                loadIdx{(k-1)*secNodes+sn}=temp;%storing the house index ..
                %connected to each secondary node
            end
            PsecMean((k-1)*secNodes+sn,h)=PsecMean((k-1)*secNodes+sn,h)+housePmean(houseIdx(k,c),h);
            PsecStd((k-1)*secNodes+sn,h)=PsecStd((k-1)*secNodes+sn,h)+housePstd(houseIdx(k,c),h)^2;
            QsecMean((k-1)*secNodes+sn,h)=QsecMean((k-1)*secNodes+sn,h)+houseQmean(houseIdx(k,c),h);
            QsecStd((k-1)*secNodes+sn,h)=QsecStd((k-1)*secNodes+sn,h)+houseQstd(houseIdx(k,c),h)^2;
            c=c+1;
            sn=sn+1;
            if sn>secNodes
                sn=1; 
            end
        end
        for sn=1:secNodes
            PsecStd((k-1)*secNodes+sn,h)=sqrt(PsecStd((k-1)*secNodes+sn,h));
            QsecStd((k-1)*secNodes+sn,h)=sqrt(QsecStd((k-1)*secNodes+sn,h));
        end
    end
end

%defining the phase of loads connected to each secondary node
loadPhase=[1,2,1,1,1,2,2,3,1,2,2,3,2,1,2,3,3,3,3,1,1,1,2,3,1,2,3,3,1,1,1,3,3,1,2,2,2,2,3,2,2,3,3,2,2,1,3,1,1,2,1,1,2,1,1];
% loadPhase=[repmat([1,2,3],1,18),1];


PpMean=zeros(32*3,length(PsecMean(1,:)));
PpStd=zeros(32*3,length(PsecMean(1,:)));
QpMean=zeros(32*3,length(PsecMean(1,:)));
QpStd=zeros(32*3,length(PsecMean(1,:)));



for n=1:32%primary nodes
    for t=1:length(PsecMean(1,:))
        for ph=1:3 %phase index
            phaseLoadIdx=(n-1)*secNodes+find(loadPhase==ph);
            PpMean(3*(n-1)+ph,t)=sum(PsecMean(phaseLoadIdx,t));
            PpStd(3*(n-1)+ph,t) =sqrt(sum(PsecMean(phaseLoadIdx,t).^2));
            QpMean(3*(n-1)+ph,t)=sum(QsecMean(phaseLoadIdx,t));
            QpStd(3*(n-1)+ph,t) =sqrt(sum(QsecMean(phaseLoadIdx,t).^2));
        end
    end
end

%extracting load data every resol seconds for load flow to get PMU reading
housesP=houseLoadData.P(:,1+delay:resol:end);
housesQ=houseLoadData.Q(:,1+delay:resol:end);

delayedHousesP=houseLoadData.P(:,1:resol:end);
delayedHousesQ=houseLoadData.Q(:,1:resol:end);

%now we need to get the aggregated loads at each secondary node to do load
%flow and get PMU readings
P=zeros(length(loadIdx),length(housesP(1,:)));
Q=zeros(length(loadIdx),length(housesP(1,:)));

for k=1:length(loadIdx)
    for t=1:length(housesP(1,:))
        P(k,t)=sum(housesP(loadIdx{k},t));
        Q(k,t)=sum(housesQ(loadIdx{k},t));
    end
end

delayedP=zeros(length(loadIdx),length(housesP(1,:)));
delayedQ=zeros(length(loadIdx),length(housesP(1,:)));

for k=1:length(loadIdx)
    for t=1:length(housesP(1,:))
        delayedP(k,t)=sum(delayedHousesP(loadIdx{k},t));
        delayedQ(k,t)=sum(delayedHousesQ(loadIdx{k},t));
    end
end

%now we need to find the load model which is an ARMA(1) model
%first we need to get the aggregate load at each primary node
PaggPri=zeros(length(Pref)*3,length(houseLoadData.P(1,:)));
QaggPri=zeros(length(Pref)*3,length(houseLoadData.P(1,:)));
for k=1:length(Pref)
    for ph=1:3 %phase index
        nodeIdx=loadIdx((k-1)*secNodes+1:k*secNodes);%all phases
        phaseLoadIdx=(nodeIdx(loadPhase==ph));
        phaseLoadIdx=[phaseLoadIdx{:}];
        PaggPri(3*(k-1)+ph,:)=sum(houseLoadData.P(phaseLoadIdx,:));
        QaggPri(3*(k-1)+ph,:)=sum(houseLoadData.Q(phaseLoadIdx,:));
    end
end

%mean and std deviation of all the samples
% PallMean=mean(PaggPri,2);
% PallStd=std(PaggPri,0,2);
% QallMean=mean(QaggPri,2);
% QallStd=std(QaggPri,0,2);

%now we calculate load variation for the given resol(time step)
Pdiff=PaggPri(:,resol+1:end)-PaggPri(:,1:end-resol);
Qdiff=QaggPri(:,resol+1:end)-QaggPri(:,1:end-resol);
%load variation standard deviation used for integrating the ensembels
PdiffStd=std(Pdiff,0,2);
QdiffStd=std(Qdiff,0,2);

%now we need to find the error time correlation coefficients
Ppsi=zeros(32*3,length(PaggPri(1,:))/dur);
Qpsi=zeros(32*3,length(PaggPri(1,:))/dur);
for pn=1:3*32%primary node index
    %finding the prediction error
    for k=1:length(PaggPri(1,:))/dur
        Perror=PaggPri(pn,(k-1)*dur+1:k*dur)-PpMean(pn,k);
        Qerror=QaggPri(pn,(k-1)*dur+1:k*dur)-QpMean(pn,k);
        %autocorrelation of the error
        Pacf=autocorr(Perror,'NumLags',resol);
        Ppsi(pn,k)=Pacf(end);
        Qacf=autocorr(Qerror,'NumLags',resol);
        Qpsi(pn,k)=Qacf(end);
    end
end






