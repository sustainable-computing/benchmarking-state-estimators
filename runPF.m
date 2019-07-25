function [V,lineCurrent] = runPF(DSSObj,P,Q)
DSSText=DSSObj.text;

%adjust loads
nNodes=33;
for ind=2:nNodes
    for ph=1:3 %phase index
        loadNumber=num2str(3*(ind-1)+ph);
        DSSText.Command = ['Edit Load.' loadNumber ' kW=' num2str(P(3*(ind-2)+ph)) ' kvar=' num2str(Q(3*(ind-2)+ph))];
    end
end

%solve power flow
DSSCircuit=DSSObj.ActiveCircuit;
DSSText.Command='Set mode=snapshot';
DSSText.Command='solve';

% DSSText.command='show voltages';
% DSSText.command='show Power';

% get bus voltages
vckt = DSSCircuit.YNodeVarray;
% vckt_len = length(vckt);
vpri_len=33*3*2;
V = vckt(1:2:vpri_len)+1i*vckt(2:2:vpri_len);

%get line currents
lineCurrent = zeros(32*2,8);
Lindex = DSSCircuit.Lines.First();
lineCounter = 0;
while Lindex ~=0
    lineCounter=lineCounter+1;
    currents = DSSCircuit.ActiveCktElement.Currents();
    lineCurrent(2*lineCounter-1,1:2) = [str2double(DSSCircuit.ActiveCktElement.BusNames{1}),str2double(DSSCircuit.ActiveCktElement.BusNames{2})];
    lineCurrent(2*lineCounter-1,3:8) = currents(1:6);
    lineCurrent(2*lineCounter,1:2) = lineCurrent(2*lineCounter-1,[2,1]);
    lineCurrent(2*lineCounter,3:8) = currents(7:12);
    Lindex = DSSCircuit.Lines.Next();
end

% % get power transfer
% linePowerTransfer = zeros(DSSCircuit.Lines.count*2,4);
% Lindex = DSSCircuit.Lines.First();
% lineCounter = 0;
% while Lindex ~= 0
%     lineCounter=lineCounter+1;
%     EPowers = DSSCircuit.ActiveCktElement.Powers();
%     linePowerTransfer(2*lineCounter-1,1:2) = [str2double(DSSCircuit.ActiveCktElement.BusNames{1}),str2double(DSSCircuit.ActiveCktElement.BusNames{2})];
%     linePowerTransfer(2*lineCounter-1,3:4) = 3e3*EPowers([1,2]);
%     linePowerTransfer(2*lineCounter,1:2) = linePowerTransfer(2*lineCounter-1,[2,1]);
%     linePowerTransfer(2*lineCounter,3:4) = 3e3*EPowers([7,8]);
%     Lindex = DSSCircuit.Lines.Next();
% end
% 
% loadPower = zeros(32,3);%bus number, Pload, Qload
% Lindex = DSSCircuit.Loads.First();
% loadCounter = 0;
% while Lindex ~= 0
%     loadCounter = loadCounter+1;
%     EPowers = DSSCircuit.ActiveCktElement.Powers();
%     loadPower(loadCounter,1) = str2double(DSSCircuit.ActiveCktElement.BusNames{1});
%     loadPower(loadCounter,2:3) = 3e3*EPowers(1:2);
%     Lindex = DSSCircuit.Loads.Next();
% end
