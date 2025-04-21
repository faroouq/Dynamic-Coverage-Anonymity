%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Main:  This simulation is for the proposed system. It keeps the nodes
%       (CHs) evidence in the file titled, "evidence.txt" which we put in
%       the belief c-program to get the beleif of each cluster.
%
%       The c-program is from Dr. Jon R. Ward and Dr. Rania El-Badry for 
%       her assistance in implementation for the evidence theory used in 
%       this paper.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function anonymity_beta()
global nodeGridx;      %BS is at the center of grid
global nodeGridy;
global clusterHeads;
global clusterState;
global droppedPackets;
global BSx;
global BSy;
global topoArea;
global coverThreshold;
global nodesPerCluster;


global energy;
global T_trans;
global P_tx;                                %transmission power in Watts
global P_rx;                                %reception power in watts

global Tx_matrix;                           %This matrix tells us who send to whom
%rng(2);
saveFile=true;

maxSensorPerNode=10;                        %No. of nodes per cluster
nodeGridx=100;                                %BS is at the center of grid (previous value is 9)

nodeGridy=100;                                %odd number (previous value is 9)
clusterHeads=zeros(nodeGridx, nodeGridy);   %this matrix stores traffic of the cluster head
clusterState=ones(nodeGridx, nodeGridy);    %This tells use whether cluster (including Head) are 
                                            %sleeping (0) or awake(1)
droppedPackets=zeros(nodeGridx, nodeGridy); %is node in path is sleepy drop packet
nodesPerCluster=round(rand(nodeGridx, nodeGridy)*maxSensorPerNode); %v How many clusters are there in each node

sleepWakeupThreshold=10;                 %what time to sleep or wakeup
highLevel=0;                                 %High traffic nodes reduce traffic according to 10 levels
lowLevel=10000;                                  %Low traffic nodes increase traffic according to 10 levels
%noSleepingNode=300;                         %How many regions can I make to sleep per Wakeup cycle 
simTimeMax=10000;                             %Unit in seconds

BSx=ceil(nodeGridx/2);                      %Location of sink node
BSy=ceil(nodeGridy/2);


nodesInRowGroup=2;                          %No. of nodes in a row
nodesInColGroup=2;                          %No. of nodes in a col for the block
coverThreshold=50;                          %percentage coverage beyound which you are not allowed to sleep, sorry :-)
topoArea=zeros(nodeGridx,nodeGridy);        %Topology group


energy=zeros(nodeGridx, nodeGridy);         %energy consumption of the CH

deploy_size=nodeGridx*nodeGridy;            %how many nodes are in the network?
Tx_matrix=zeros(deploy_size, deploy_size);	%who transmits to who :-)


dataRate=250e+03;                           %250Kbps is datarate of zigbee
pkt_size=128;                               %packet size: 128 bytes, including a variable payload of up to 104 bytes. 
T_trans=pkt_size/dataRate;

P_tx=90e-03;                                %transmission power in Watts
P_rx=72e-03;                                %reception power in watts



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Main()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
buildBlocks(nodesInRowGroup, nodesInColGroup);  %Group up the nodes

eCnt=1;
eTrack=[];                                      %Track the energy

bCnt=1;
bTrack=[];                                      %Track the beleif of path


%%%%%%%%%%%%%%%%%%%%%generating the path to monitor%%%%%%%%%%%%%%%%%%%%%%%%
 res=simTimeMax+1;                                         %get belief after every res
% 
% x=1;
% y=1;
% pathCnt=1;
% path1=[x y];
% path2=path1;
% while(or(x~=BSx,y~=BSy))
%     pathCnt=pathCnt+1;
%     if(x<BSx)
%         x=x+1;
%     elseif (y<BSy)
%         y=y+1;
%     end
%     path1(pathCnt,:)=[x,y];
% end
% 
% 
% BBsx=3;
% BBsy=2;
% pathCnt=1;
% x=1;
% y=1;
% while(or(x~=BBsx,y~=BBsy))
%     pathCnt=pathCnt+1;
%     if(x<BBsx)
%         x=x+1;
%     elseif (y<BBsy)
%         y=y+1;
%     end
%     path2(pathCnt,:)=[x,y];
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%tests=[1,9;9,1;9,9;1,1];       % Removed by Farouq on 31 December 2022 it
                                % was there to test the what happens if
                                % only on CH transmits. It was there for
                                % testing ony. So, I removed it afterwards
                                

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Simulation begins
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for t=1:1:simTimeMax                            %The start simulation
    if(random('exp',1)>0.5)                     %A cluster head wants to transmit
        %------------Removed by farouq on 31st December 2022---------------
        %clusterCoordx=tests(mod(t,4)+1,1);%round(random('uniform',1,nodeGridx));
        %clusterCoordy=tests(mod(t,4)+1,2);%round(random('uniform',1,nodeGridy));
        %--------------It is replaced wth the 2 lines below----------------
        
        clusterCoordx=round(random('uniform',1,nodeGridx));
        clusterCoordy=round(random('uniform',1,nodeGridy));
        
        
        while(or(clusterCoordx==BSx,clusterCoordy==BSy))
            clusterCoordx=round(random('uniform',1,nodeGridx));
            clusterCoordy=round(random('uniform',1,nodeGridy));
        end
        
        
        if not(and(and(clusterCoordx==BSx, clusterCoordy==BSy),clusterState(clusterCoordx, clusterCoordy)==0))      %BS does not sense and region not asleep
            nodesAwake = round((clusterState(clusterCoordx, clusterCoordy)) * nodesPerCluster(clusterCoordx, clusterCoordy));
            count = nodesAwake; %round(random('uniform',0,nodesAwake));                       %how many sensor are reporting data
            if(count>0)
                clusterHeads(clusterCoordx, clusterCoordy) = clusterHeads(clusterCoordx, clusterCoordy) + count;    %add to clusterHeads total traffic
                energy(clusterCoordx, clusterCoordy) = energy(clusterCoordx, clusterCoordy)+ P_rx*T_trans*count;    %energy due to traffic on CH
                route(clusterCoordx, clusterCoordy);      %aggregated data is then sent to the BS
            end
        end
    end
    
    if(mod(t, sleepWakeupThreshold)==0)
       clusterState=ones(nodeGridx, nodeGridy);     %all sleeping nodes should wakeup
       
       lowerTraffic(highLevel);                     %Lower traffic from areas with low traffic
       heightenTraffic(lowLevel);                   %Heighten traffic from areas with low traffic

    end
    
   	energy(BSx,BSy)=0;                                      %Lets remove the base station
	eTrack(1,eCnt)=sum(sum(energy))/(nodeGridx*nodeGridy-1);%average energy consumed
 	eTrack(2,eCnt)=t;
	eCnt=eCnt+1;  
    
    if(mod(t,res)==0)
        myStats=myEval;
        bCnt=bCnt+1;
    end
end

Traffic = sum(sum(clusterHeads));
PktDropped = sum(sum(droppedPackets));
BSTraffic = clusterHeads(BSx,BSy);
temporary=clusterHeads;
temporary(BSx,BSy)=0;
MaxTraffic = max(temporary(:));
[rr,cc]=find(temporary>=BSTraffic);
NoMaxTraffic=size(rr,1);

figure(simTimeMax);
subplot(2,1,1)
image(clusterHeads);
colorbar;
text=['Simlation time = ', num2str(simTimeMax)];
title(text);
xlabel('X Coordinates')
ylabel('Y Coordinates')


subplot(2,1,2)
surfc(clusterHeads);
title(text);
xlabel('X Coordinates')
ylabel('Y Coordinates')


figure(2);
hold on;
plot(eTrack(2,:), eTrack(1,:), '-r', 'linewidth',2);
title('Average Energy Consumption against Time');
xlabel('Time');
ylabel('Energy (J)');
grid on;
hold off;


% figure(3);
% hold on;
% plot(bTrack(3,:), bTrack(1,:), '-*r', 'linewidth',2);
% plot(bTrack(3,:), bTrack(2,:), '-db', 'linewidth',2);
% title('Evidence against Time');
% xlabel('Time');
% ylabel('Evidence (%)');
% str1=sprintf('Path to BS (%d,%d)',BSx,BSy);
% str2=sprintf('Path to a Node (%d,%d)', BBsx, BBsy);
% legend(str1,str2);
% grid on;
% hold off;

figure; imagesc(clusterHeads); colorbar;
text=['Simlation time = ', num2str(simTimeMax)];
title(text);
xlabel('X Coordinates')
ylabel('Y Coordinates')


figure; surfc(clusterHeads);
title(text);
xlabel('X Coordinates')
ylabel('Y Coordinates')
zlabel('Traffic')

    if(saveFile)
        fID=fopen('evidence10000.txt','w');
        for i=1:size(Tx_matrix,1)
            for j=1:size(Tx_matrix,1)
                fprintf(fID,"%.2f\t",double(Tx_matrix(i,j)));
            end

            fprintf(fID,"\n");
        end
        fclose(fID);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This function controll the level of areas with high traffic
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function lowerTraffic(level)
    global nodeGridx;       %BS is at the center of grid
    global nodeGridy;
    global clusterHeads;
    global clusterState;
    global droppedPackets;
    global BSx;
    global BSy;
    global topoArea;
    global coverThreshold;
    global nodesPerCluster;
    
    
    
    temp=clusterHeads;      %we cant afford to loose original data
    temp(BSx,BSy)=0;        %lets exclude the BS because without it there is no network
    
    
    while (level>0)                                 %(nds<=noSleepingNode)
       [coordX, coordY]=find(temp==max(temp(:)));   %step1 - Find the maximum

       %step2 - check percentage coverage
       if(temp(coordX(1), coordY(1))==-1)   
            break;  %if all nodes are assigned traffic level then leave
       end


        for ptr = 1:1:length(coordX)
            myBlock=topoArea(coordX(ptr), coordY(ptr));         %a)which blck r u 4rm
            [boundx,boundy]=find(topoArea==myBlock);            %b)get boundary of block
            myNodes=nodesPerCluster(coordX(ptr), coordY(ptr));  %c)get nodes managed by that clusterHead
                                                                %d)get total nodes covering the area:
                                                                    %i)get all nodes in block
                                                                    %ii)subtract those that are asleep (since they are nor covering at this time)

            allNodesInBlock=nodesPerCluster(boundx(1):boundx(size(boundx)),boundy(1):boundy(size(boundy)));
            nodesAsleep=clusterState(boundx(1):boundx(size(boundx)),boundy(1):boundy(size(boundy)));
            totalBlockNodes=sum(sum(allNodesInBlock.*nodesAsleep));
            percentageCoverage = (myNodes/totalBlockNodes)*100;

            %Lest's check if coverage and of the node has not exceed
            %threshold and that the node traffic is not zero (this happens
            %at the beggining of the simulation.
           if(not(and(and(and(coordX(ptr)<=0, coordY(ptr)<=0),max(temp(:)==0)), percentageCoverage>coverThreshold)))  
                clusterState(coordX(ptr), coordY(ptr))=1/level;   %make node with maximum traffic sleep
                temp(coordX(ptr), coordY(ptr))=-1;              %clear record so that we can get the next guy
           end
        end 

        level = level - 1;              %second wave
    end
end 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This function controll the level of areas with low traffic
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function heightenTraffic(level)
    global nodeGridx;       %BS is at the center of grid
    global nodeGridy;
    global clusterHeads;
    global clusterState;
    global droppedPackets;
    global BSx;
    global BSy;
    global topoArea;
    global coverThreshold;
    global nodesPerCluster;
    
    
    
    temp=clusterHeads;      %we cant afford to loose original data
    temp(BSx,BSy)=0;        %lets exclude the BS because without it there is no network
    
    
    while (level>0)                                 %(nds<=noSleepingNode)
       [coordX, coordY]=find(temp==min(temp(:)));   %step1 - Find the minimum
       
       %step2 - check percentage coverage
       if(temp(coordX(1), coordY(1))==-1)   
            break;  %if all nodes are assigned traffic level then leave
       end


        for ptr = 1:1:length(coordX)
            myBlock=topoArea(coordX(ptr), coordY(ptr));         %a)which blck r u 4rm
            [boundx,boundy]=find(topoArea==myBlock);            %b)get boundary of block
            myNodes=nodesPerCluster(coordX(ptr), coordY(ptr));  %c)get nodes managed by that clusterHead
                                                                %d)get total nodes covering the area:
                                                                    %i)get all nodes in block
                                                                    %ii)subtract those that are asleep (since they are nor covering at this time)

            allNodesInBlock=nodesPerCluster(boundx(1):boundx(size(boundx)),boundy(1):boundy(size(boundy)));
            nodesAsleep=clusterState(boundx(1):boundx(size(boundx)),boundy(1):boundy(size(boundy)));
            totalBlockNodes=sum(sum(allNodesInBlock.*nodesAsleep));
            percentageCoverage = (myNodes/totalBlockNodes)*100;

            %Lest's check if coverage and of the node has not exceed
            %threshold and that the node traffic is not zero (this happens
            %at the beggining of the simulation.
           if(not(and(and(and(coordX(ptr)<=0, coordY(ptr)<=0),max(temp(:)==0)), percentageCoverage>coverThreshold)))  
                clusterState(coordX(ptr), coordY(ptr))=level; %make node with minimum traffic sleep
                temp(coordX(ptr), coordY(ptr))=-1;              %clear record so that we can get the next guy
           end
        end 

        level = level - 1;              %second wave
    end
end 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This function route data to the BS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function route(nodex, nodey)
global nodeGridx;      %BS is at the center of grid
global nodeGridy;
global clusterHeads;
global clusterState;
global droppedPackets;
global BSx;
global BSy;

global energy;
global T_trans;
global P_tx;                                %transmission power in Watts
global P_rx;                                %reception power in watts

global Tx_matrix;                           %to record transmission


finalx=BSx;                                 %initialize destination
finaly=BSy;

    if not(and(nodex==BSx, nodey==BSy))     %BS does not sense
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       %First route
       if(nodex<BSx)
           for i=nodex:1:BSx-1
              if (clusterState(i, nodey)==0)   %if node is asleep then drop packet
                   droppedPackets(i, nodey)=droppedPackets(i, nodey)+1; %record dropped packet
                   return
              end
               
              clusterHeads(i, nodey) =  clusterHeads(i, nodey) + 1;
              energy(i+1,nodey)=energy(i+1,nodey)+P_rx*T_trans;     %record reception
              energy(i,nodey)=energy(i,nodey)+P_tx*T_trans;         %record transmission
              
              [txer]=getMyID(i,nodey);                             	%get the transmitters ID
              [rxer]=getMyID(i+1,nodey);                           	%get receivers ID
              Tx_matrix(txer,rxer)=Tx_matrix(txer,rxer)+1;         	%not the transmission
              finalx=i+1;
           end

       elseif(nodex>BSx)
           for i=nodex:-1:BSx+1
               if (clusterState(i, nodey)==0)   %if node is asleep then drop packet
                   droppedPackets(i, nodey)=droppedPackets(i, nodey)+1; %record dropped packet
                   return
               end
               
              clusterHeads(i, nodey) =  clusterHeads(i, nodey) + 1;
              energy(i-1,nodey)=energy(i-1,nodey)+P_rx*T_trans;     %record reception
              energy(i,nodey)=energy(i,nodey)+P_tx*T_trans;      	%record transmission
              
              [txer]=getMyID(i,nodey);                           	%get the transmitters ID
              [rxer]=getMyID(i-1,nodey);                           	%get receivers ID
              Tx_matrix(txer,rxer)=Tx_matrix(txer,rxer)+1;          %not the transmission
              finalx=i-1;
           end           
       end
       
       %Movement along x-axis
       if(nodey<BSy)
           for i=nodey:1:BSy-1                     %since we reached this point we dont have to trace it again
              if (clusterState(finalx, i)==0)   %if node is asleep then drop packet
                   droppedPackets(finalx, i)=droppedPackets(finalx, i)+1; %record dropped packet
                   return
              end
               
              clusterHeads(finalx, i) =  clusterHeads(finalx, i) + 1;
              energy(finalx,i)=energy(finalx,i)+P_tx*T_trans;           %record transmission
              energy(finalx,i+1)=energy(finalx,i+1)+P_rx*T_trans;       %record reception
              
              [txer]=getMyID(finalx,i);                                	%get the transmitters ID
              [rxer]=getMyID(finalx,i+1);                           	%get receivers ID
              Tx_matrix(txer,rxer)=Tx_matrix(txer,rxer)+1;             	%not the transmission
           end
           
           
       elseif(nodey>BSy)
           for i=nodey:-1:BSy+1
               if (clusterState(finalx, i)==0)   %if node is asleep then drop packet
                   droppedPackets(finalx, i)=droppedPackets(finalx, i)+1; %record dropped packet
                   return
               end
               
              clusterHeads(finalx, i) =  clusterHeads(finalx, i) + 1;
              energy(finalx,i)=energy(finalx,i)+P_tx*T_trans;           %record transmission
              energy(finalx,i-1)=energy(finalx,i-1)+P_rx*T_trans;       %record reception
              
              [txer]=getMyID(finalx,i);                               	%get the transmitters ID
              [rxer]=getMyID(finalx,i-1);                             	%get receivers ID
              Tx_matrix(txer,rxer)=Tx_matrix(txer,rxer)+1;              %not the transmission
           end           
       end
       
       clusterHeads(finalx, finaly) =  clusterHeads(finalx, finaly) + 1;
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       
       
       
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       %Second route       
       %Movement along y-axis
       if(nodey<BSy)
           for i=nodey:1:BSy-1
               if (clusterState(nodex, i)==0)   %if node is asleep then drop packet
                   droppedPackets(nodex, i)=droppedPackets(nodex, i)+1; %record dropped packet
                   return
               end
               
              clusterHeads(nodex, i) =  clusterHeads(nodex, i) + 1;
              energy(nodex,i)=energy(nodex,i)+P_rx*T_trans;         %record reception
              energy(nodex,i+1)=energy(nodex,i+1)+P_tx*T_trans;     %record transmission
              
              [rxer]=getMyID(nodex,i);                              %get the transmitters ID
              [txer]=getMyID(nodex,i+1);                            %get receivers ID
              Tx_matrix(txer,rxer)=Tx_matrix(txer,rxer)+1;          %not the transmission
              finaly=i+1;
           end
           
           
       elseif(nodey>BSy)
           for i=nodey:-1:BSy+1
               if (clusterState(nodex, i)==0)   %if node is asleep then drop packet
                   droppedPackets(nodex, i)=droppedPackets(nodex, i)+1; %record dropped packet
                   return
               end
               
              clusterHeads(nodex, i) =  clusterHeads(nodex, i) + 1;
              energy(nodex,i)=energy(nodex,i)+P_rx*T_trans;         %record reception
              energy(nodex,i-1)=energy(nodex,i-1)+P_tx*T_trans;     %record transmission
                                                                    %going reverse
              [rxer]=getMyID(nodex,i);                              %get the transmitters ID
              [txer]=getMyID(nodex,i-1);                            %get receivers ID
              Tx_matrix(txer,rxer)=Tx_matrix(txer,rxer)+1;          %not the transmission
              finaly=i-1;
           end           
       end
       
       %Movement along x-axis
       if(nodex<BSx)
           for i=nodex:1:BSx-1
               if (clusterState(i, finaly)==0)   %if node is asleep then drop packet
                   droppedPackets(i, finaly)=droppedPackets(i, finaly)+1; %record dropped packet
                   return
               end
              clusterHeads(i, finaly) =  clusterHeads(i, finaly) + 1;
              energy(i+1,finaly)=energy(i+1,finaly)+P_tx*T_trans;      	%record transmission
              energy(i,finaly)=energy(i,finaly)+P_rx*T_trans;           %record reception
              
              [rxer]=getMyID(i,finaly);                               	%get the transmitters ID
              [txer]=getMyID(i+1,finaly);                              	%get receivers ID
              Tx_matrix(txer,rxer)=Tx_matrix(txer,rxer)+1;             	%not the transmission
           end

       elseif(nodex>BSx)
           for i=nodex:-1:BSx+1
               if (clusterState(i, finaly)==0)   %if node is asleep then drop packet
                   droppedPackets(i, finaly)=droppedPackets(i, finaly)+1; %record dropped packet
                   return
               end
              clusterHeads(i, finaly) =  clusterHeads(i, finaly) + 1;
              energy(i-1,finaly)=energy(i-1,finaly)+P_tx*T_trans;       %record transmission
              energy(i,finaly)=energy(i,finaly)+P_rx*T_trans;           %record reception
              
              [rxer]=getMyID(i,finaly);                                 %get the transmitters ID
              [txer]=getMyID(i-1,finaly);                               %get receivers ID
              Tx_matrix(txer,rxer)=Tx_matrix(txer,rxer)+1;              %not the transmission
           end           
       end
       

       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This function takes the location of a node and returns its ID.
%Node ID is assigned incrementally from left-right and from top-down 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ur_id]=getMyID(nodeX, nodeY)
global nodeGridx;               %BS is at the center of grid


ur_id=nodeX+(nodeY-1)*nodeGridx;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This function takes the ID of a node and returns its location.
%Node ID is assigned incrementally from left-right and from top-down 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [nodeX, nodeY]=getMyCoords(ur_id)
    global nodeGridx;               %BS is at the center of grid

    nodeX=-1;
    nodeY=-1;
    if(ur_id>0)
        nodeX=ceil(ur_id/nodeGridx);
        nodeY=mod(ur_id,nodeGridx);
        if (nodeY==0)
            nodeY=nodeGridx;
        end
    end
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This function groups the sensor nodes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function buildBlocks(rowSize, colSize)
    global nodeGridx;      %BS is at the center of grid
    global nodeGridy;
    global topoArea;
    

    regionRowCount=0;
    regionColCount=1;
    for xCoords=1:1:nodeGridx
       for yCoords=1:1:nodeGridy
            topoArea(xCoords, yCoords) = regionColCount + regionRowCount*ceil(nodeGridy/colSize);
            if (mod(yCoords, colSize)==0)
                regionColCount=regionColCount+1;
            end
       end

       regionColCount=1;
       if(mod(xCoords, rowSize)==0)
           regionRowCount=regionRowCount+1;
       end
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This beleif calculator
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [bel]=belief(path)
    global Tx_matrix;                           %This matrix tells us who send to whom
    
    u=[];
    uCnt=1;
    
    v=[];
    vCnt=1;
    for i=1:size(path,1)
    	if(i<size(path,1))
            u(uCnt)=Tx_matrix(getMyID(path(i,1),path(i,2)), getMyID(path(i+1,1),path(i+1,2)));
            uCnt=uCnt+1;
        end
        
        for j=1:size(path,1)
           v(vCnt)=Tx_matrix(getMyID(path(i,1),path(i,2)), getMyID(path(j,1),path(j,2)));
           vCnt=vCnt+1;
        end
    end
    
    minU=min(u);
    v(vCnt)=minU;
    u(uCnt)=minU;
    
    denum=sum(v);
    bel=sum(u/denum);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Full facotrial search for the right path
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [bel, realPath]=myEval(pass)
    global nodeGridx;                       %BS is at the center of grid
    global nodeGridy;
    
    realPath=[];                            %initializing variables
    possiblePath=[];
    maxNodeID=nodeGridx*nodeGridy;          %ID of the last node
    myPath=zeros(1, maxNodeID);
    myPath(end)=0;                          %Initialize least significant digit
    
    cnt=1;                                  %keep track of the belief
    bel=zeros(101,1);
    
    [~,c]=find(myPath==maxNodeID);          %check whether we have reach the end
    while(size(c,2)~=maxNodeID)             %when we have not finish with all the nodes
        digit=maxNodeID;                    %current digit
        myPath(digit)=myPath(digit)+1;      %increment like a counter
        
        while(myPath(digit)>maxNodeID)      %rolling over, incrementing other digits
            myPath(digit)=1;                %Least Significant Digit (LSD) goes back to one
            digit=digit-1;                  %Next Significant Digit
            myPath(digit)=myPath(digit)+1;  %...increment
        end
        
        
        possiblePath=[];                            %initialize
        if(sum(myPath==0)<maxNodeID-1)          %if pathe consists of atleast two nodes
            for i=1:maxNodeID                   %lets get the coordinates of each node
                [x,y]=getMyCoords(myPath(i));   %Get nodes coordinates
                if(and(x>0, y>0))               %if node is valid
                    possiblePath=[possiblePath; [x y]];
                end
            end
        
        
            bel_per=floor(belief(possiblePath)*100);%sum of beliefs in percentile
            if(isnan(bel_per))
                bel_per=0;
                bel(bel_per+1)=bel(bel_per+1)+1;
            elseif(isinf(bel_per))
                bel_per=100;
                bel(bel_per+1)=bel(bel_per+1)+1;
            else
                bel(bel_per+1)=bel(bel_per+1)+1;
            end
            
            [~,c]=find(myPath==maxNodeID);      %check if we reach the end
            %fprintf('%d\t', myPath);
            %fprintf('\n');
        end
    end
    
    
end
