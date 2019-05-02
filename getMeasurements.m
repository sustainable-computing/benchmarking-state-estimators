function [z,zType]=getMeasurements(Vtrue, ~, loadPower, PMUnodes)
%measurements
%Load power
z=loadPower(:,2);%real powers of the nodes
zType(:,2)=loadPower(:,1);%defining the corresponding nodes 
%Types: 1:Pij 2:Pi 3:Qij 4:Qi 5:|Vi| 6:theta Vi 7:|Iij|
zType(:,1)=2;%defining measurement type
z=[z;loadPower(:,3)];%reactive powers of the nodes
zTypeAdd(:,2)=loadPower(:,1);%defining the corresponding nodes 
zTypeAdd(:,1)=4;%defining measurement type
zType=[zType;zTypeAdd];
zType(:,3)=0;
z=-z;%injected powers

%Adding 1st node with zero load
% z=[z;0;0];
% zTypeAdd=[2,1,0;4,1,0];
% zType=[zType;zTypeAdd];

%voltage magnitudes
PMUnodes1ph=sort([3*PMUnodes-2,3*PMUnodes-1,3*PMUnodes]);
z=[z;abs(Vtrue(PMUnodes1ph))'];
zTypeAdd=5*ones(length(PMUnodes1ph),1);%type 5 is voltage magnitude
zTypeAdd(:,2)=PMUnodes1ph';
zTypeAdd(:,3)=0;
zType=[zType;zTypeAdd];

%voltage angles
z=[z;(angle(Vtrue(PMUnodes1ph(2:end)))-angle(Vtrue(1)))'];
zTypeAdd=6*ones(length(PMUnodes1ph)-1,1);%type 6 is voltage angle
zTypeAdd(:,2)=(PMUnodes1ph(2:end))';
zTypeAdd(:,3)=0;
zType=[zType;zTypeAdd];
end