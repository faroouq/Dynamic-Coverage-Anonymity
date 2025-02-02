%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Topology discovery phase
%========================
%This code simulates how the nrnadmoly deployed nodes choose
%cluster heads (a.k.a fake base station). The code also shows
%how the network is setup.
%
%Assumption: Nodes have the ability to localize relaive to the BS
%
%Author:    Farouq Aliyu
%Date:      22th January, 2018
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [time, num_clusters]=topology_Discovery_final(n,z,r)
    %Global variables
    %========================
    global BSX                      %Base station is at the center
    global BSY
    global node_Xcoord              %Location of nodes
    global node_Ycoord 
    global radius                   %Radius of sensor nodes
    global cluster_ID             
    global cluster_HD
    global no_nodes
    global state_of_node
    global runs
    global zeta                     %number of 
    global time                     %simulation time
    
    %Delaraction of variables
    %========================
    zeta=z;                         %random number generator [1,z]
    display_sim=0;                  %a switch that shows the -- 
                                    %simulation state 1=On and 0=Off
    no_nodes=n;
    no_of_clusters=10;
    runs = 2;                       %Number of times election is re-run 
                                    %before a cluster head is actually
                                    %selected as a fake node in order to
                                    %avoid re-election in the respective
                                    %clusters.
    
    area_of_deployment=1600;        %unit (m^2)
    length=area_of_deployment^0.5;  %length=breadth of area
    BSX=round((length/2), 2);       %Base station is at the center
    BSY=round((length/2), 2);
    
    %radius = 5/no_of_clusters; %transmission range of nodes
    radius = r;

    
    %1. Deployment
    %====================================================================
    %   Nodes are uniformly deployed at random
    node_Xcoord = round(random('uniform', 0, length, 1, no_nodes), 2);
    node_Ycoord = round(random('uniform', 0, length, 1, no_nodes), 2);
    time=1;
    %display(display_sim);           %Display simultation
    
    
    %2. Selecting cluster heads
    %====================================================================
    vote_CH();
    
    
    
    %3. Fake node selection
    %====================================================================
    %vote_FBS();
    
    num_clusters=max(cluster_ID);
    display(display_sim);           %Display simultation
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This function is the selection process for cluster head
%Assumption: Nodes have the ability to localize relaive to the BS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function vote_FBS()
    %Global variables
    %========================
    global BSX                      %Base station is at the center
    global BSY
    global node_Xcoord              %Location of nodes
    global node_Ycoord 
    global radius                   %Radius of sensor nodes
    global cluster_ID             
    global cluster_HD
    global no_nodes
    global state_of_node
    global runs
    global time


    %====================================================================
    %backup_time:   each node select a random number n broadcast it to its
    %============   neighbours the one with the smallest number becomes the 
    %               the cluster head (or Fake Base Station)
    %
    %State_of_nodes:tells us whether node has join a cluster (1) or it is 
    %============== free (0) and FBS (3) and Cluster HD (2)
    %====================================================================
    node_in_dir=zeros(1,no_nodes); 	%worst case all nodes are to the left
    count_node_in_dir=0;            %count each side independently

    broadcast_Xnode=BSX;            %The node currently broadcasting
    broadcast_Ynode=BSY;
    dir=1;                          %Direction of propgation
    timeD=[0,0,0,0];                %Time for each direction (the max is 
                                    %the time since everything is happens
                                    %simultaneously)
                                    
    %First BS broadcast to its one hop neigbours
    while (dir <= 4)
        timeD(dir)=timeD(dir)+1;    %let's not time for each direction 
        for allnodes=1:no_nodes
            rSquareBS = (broadcast_Xnode-node_Xcoord(allnodes))^2 + (broadcast_Ynode-node_Ycoord(allnodes))^2;
            if not(rSquareBS==0)            %if is it not myself
                if (radius^2 >= rSquareBS)	%if node is one hop away
                    dx=broadcast_Xnode-node_Xcoord(allnodes);
                    dy=broadcast_Ynode-node_Ycoord(allnodes);
                    if(and(and((dx<0),(abs(dx)>abs(dy))), (dir==1)))%west
                        count_node_in_dir = count_node_in_dir + 1;	%ptr nodes
                        node_in_dir(count_node_in_dir) = allnodes;  %rec nodes
                    elseif(and(and((dx>0),(abs(dx)>abs(dy))), (dir==2))) %east
                        count_node_in_dir = count_node_in_dir + 1;	%ptr nodes
                        node_in_dir(count_node_in_dir) = allnodes;  %rec nodes
                    elseif(and(and((dy<0),(abs(dy)>abs(dx))), (dir==3))) %north
                        count_node_in_dir = count_node_in_dir + 1;	%ptr nodes
                        node_in_dir(count_node_in_dir) = allnodes;  %rec nodes
                    elseif(and(and((dy>0),(abs(dy)>abs(dx))), (dir==4))) %south
                        count_node_in_dir = count_node_in_dir + 1;	%ptr nodes
                        node_in_dir(count_node_in_dir) = allnodes;  %rec nodes
                    end
                end
            end
        end
        
        %if there are nodes one hop a way select FBS from the pool
        if (count_node_in_dir > 0) 
        %polling will be carried out by each side so that only one node is 
        %selected from nodes in the broadcast range of the BS at each side 
        backup_time=(round(random('uniform', 0, 100, 1, size(node_in_dir, 2)), 2)./node_in_dir).*node_in_dir;

            minTime=min(backup_time);
            [~,c]=find(backup_time==minTime);	%min time and who has it

            electionCount=0;    %if a cluster head is elected twice then their 
                                %a probability that it is alone in the
                                %broadcast range. This means their has to be
                                %re-election in that cluster

            %while there is a collision or it is a cluster HD and it is not
            %more than the number of re-runs then vote again.
            while (or((size(c,2)>1), and((state_of_node(node_in_dir(c)) == 2),(electionCount>runs))))	
                backup_time=(round(random('uniform', 0, 100, 1, size(node_in_dir, 2)), 2)./node_in_dir).*node_in_dir;
                minTime=min(backup_time);
                [~,c]=find(backup_time==minTime);	%min time and who has it

                if(cluster_HD(c) > 0)                   %if cluster HD wins
                    electionCount = electionCount + 1;  %count number of wins
                end
            end

            %After election, swearing-in ceremony follows
            %if winner is a Cluster HD
            if (state_of_node(node_in_dir(c)) == 2) 
                [~,cc] = find(cluster_HD==node_in_dir(c));
                cluster_HD(cc) = 0;              %remove it from cluster head
                state_of_node(node_in_dir(c))=3; %update status to FBS
                

                %get the cluster in ANARCHY :-p
                theClusterID = cluster_ID(node_in_dir(c));
                reElect(theClusterID);      %do a re-election for the cluster
                
                
                %reinitialize for the next iteration
                broadcast_Xnode=node_Xcoord(node_in_dir(c)); %next propagation
                broadcast_Ynode=node_Ycoord(node_in_dir(c));
                node_in_dir=zeros(1,no_nodes);
                count_node_in_dir=0;
            else
                state_of_node(node_in_dir(c))=3;  %update status to FBS
                
                %reinitialize for then next iteration
                broadcast_Xnode=node_Xcoord(node_in_dir(c)); %next propagation
                broadcast_Ynode=node_Ycoord(node_in_dir(c));
                node_in_dir=zeros(1,no_nodes);
                count_node_in_dir=0;
            end
        else
            dir = dir + 1;          %let's move to the next direction
            
            %reinitialize for new direction
            broadcast_Xnode=BSX;    %start all over
            broadcast_Ynode=BSY;
            node_in_dir=zeros(1,no_nodes);
            count_node_in_dir=0;    %just to be safe
        end  
    end
    
    time = time + max(timeD);       %The maximum is convergence time
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This function is the selection process for cluster head
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function reElect(myclusterID)
    %Global variables
    %========================
    %global BSX                      %Base station is at the center
    %global BSY
    global node_Xcoord              %Location of nodes
    global node_Ycoord 
    global radius                   %Radius of sensor nodes
    global cluster_ID             
    global cluster_HD
    global no_nodes
    global state_of_node
    global zeta                     %number of
    global time                     %Simulation time
    
    count = 1;
    clusterCount = max(cluster_ID);             %the last cluster_ID so 
                                                %that we can continue from
                                                %where we stop. 
                                                
    [~, c] = find(cluster_ID == myclusterID);   %get all the nodes in the 
                                                %cluster
                                                
    while(size(c,2)>0)    %foreach node in the given cluster
        time = time + 1;    
        
        %to the poll!!!!
        reBackupTime = round(random('uniform', 1, zeta, 1, size(c,2)));
        winner = min(reBackupTime);
        [~, cc] = find(reBackupTime==winner);
        
        while(size(cc,2)>1)     %redo while there is a tie
            time = time + 1;
            reBackupTime = round(random('uniform', 1, zeta, 1, size(c,2)));
            winner = min(reBackupTime);             %who's the winner
            [~, cc] = find(reBackupTime==winner);   %it better not be a tie
        end
        
        time = time + 1;
        for i=c                         %foreach node that's not the winner
            winningNode = c(cc);        %note c and cc form a link-list 
                                        %cc as the index.
                                        %registring new Cluster Head
            cluster_ID(winningNode) = count + clusterCount;
            state_of_node(winningNode) = 2;
            cluster_HD(count + clusterCount) = winningNode;
            
            if not(i==winningNode)      %if I'm not talking to myself
                range = (node_Xcoord(winningNode) - node_Xcoord(i))^2 - (node_Ycoord(winningNode) - node_Ycoord(i))^2;
                
                %check if node is with range
                if and((radius^2>=range), (state_of_node(i)<3))    
                    cluster_ID(i) = count + clusterCount;   %update cluster
                end
            end
        end
        
        [~, c] = find(cluster_ID == myclusterID);   %get all the nodes in the 
    end                                             %cluster
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This function is the selection process for cluster head
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function vote_CH()
    %Global variables
    %========================
    %global BSX                      %Base station is at the center
    %global BSY
    global node_Xcoord              %Location of nodes
    global node_Ycoord 
    global radius                   %Radius of sensor nodes
    global cluster_ID             
    global cluster_HD
    global no_nodes
    global state_of_node
    global zeta                     %number of
    global time
    

    %====================================================================
    %backup_time:   each node select a random number n broadcast it to its
    %============   neighbours the one with the smallest number becomes the 
    %               the cluster head (or Fake Base Station)
    %
    %State_of_nodes:tells us whether node has join a cluster (1) or it is 
    %============== free (0) and cluster head (2) and FBS(3)
    %====================================================================
    backup_time = round(random('uniform', 1, zeta, 1, no_nodes), 2);
    state_of_node = zeros(1, no_nodes);
    cluster_ID = zeros(1, no_nodes);
    cluster_HD = zeros(1, no_nodes);
    
    count = 1;      %Count is the number of times polling has been done
    cluster_Count = 0;
    temp_State_of_node = state_of_node;
    
    time = time + 1;
    while (count<no_nodes)
        min_num=min(backup_time);           %what is the minimum backup  
        [~,c]=find(backup_time==min_num);   %time and who has it
        for i = c
            cluster_Count = cluster_Count + 1;  %lets give each cluster ID
            temp_cluster_ID = cluster_ID;       
            temp_cluster_ID(i) = cluster_Count;
            temp_cluster_HD = cluster_HD;
            temp_State_of_node(i) = 2;
            
            %foreach node
            for j=1:no_nodes
                %if the node has not join a cluster and it is not myself
                if(and(not(i==j), (temp_cluster_ID(j)==0)))
                    range=(node_Xcoord(i) - node_Xcoord(j))^2 + (node_Ycoord(i) - node_Ycoord(j))^2;
                    %if node is within the range and backup time is not the
                    %same as that of the min, then there is no collission
                    %and node joins
                    if (and((radius^2 >= range), (backup_time(j))>min_num))
                        temp_State_of_node(j) = 1;
                        temp_cluster_ID(j) = cluster_Count;
                    elseif(and((radius^2 >= range), (backup_time(j))==min_num))
                    %else if they have the same backup time, then there is
                    %collission and only both nodes should get new backup
                    %time and the polling will be redone (re-Election ;-))
                        backup_time(i)=round(random('uniform', 0, 100));
                        backup_time(j)=round(random('uniform', 0, 100));
                        
                        %abort Election
                        temp_State_of_node = state_of_node; 
                        temp_cluster_ID = cluster_ID;                        
                        cluster_Count = cluster_Count - 1;
                        break;
                    end
                end
                
                temp_cluster_HD(cluster_Count) = i;     %note cluster head
            end
            %update state of node
            state_of_node = temp_State_of_node;
            cluster_ID = temp_cluster_ID;
            cluster_HD = temp_cluster_HD;
            
            %lets make all those nodes that have a cluster have no
            %backuptime this can be done by making the backup time bigger
            %then the max. expected value (100).
            backup_time = 200*state_of_node + backup_time;  
        end
        
        %if all nodes have joined a cluster
        [~,not_Join] = find(state_of_node==0);
        total = size(not_Join);
        if(total(2)==0)
            break;
        end
        
        %increase iteration count
        count = count + 1;
        time = time + 1;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This function dsiplays the state of the simulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function display(display_state)
    %Global variables
    %========================
    global BSX          %Base station is at the center
    global BSY
    global node_Xcoord	%Location of nodes
    global node_Ycoord
    global radius
    global cluster_ID
    global cluster_HD
    global no_nodes
    global state_of_node
    
    
    radius_Scale=130;    %number fits my screen calibration
    
    for i = 1:no_nodes
       if (mod(cluster_ID(i),3) == 0)
           node_Xcoord_0(i)=node_Xcoord(i);
           node_Ycoord_0(i)=node_Ycoord(i);
       elseif(mod(cluster_ID(i),3) == 1)
           node_Xcoord_1(i)=node_Xcoord(i);
           node_Ycoord_1(i)=node_Ycoord(i);
       else
           node_Xcoord_2(i)=node_Xcoord(i);
           node_Ycoord_2(i)=node_Ycoord(i);
       end
    end
    
    
    %display deployed nodes as blue fat dots
    if (display_state==1)
        max(cluster_ID)
        figure;
        hold on;
        
        plot(node_Xcoord_0, node_Ycoord_0, 'bo', 'MarkerFaceColor', 'b');
        %plot(node_Xcoord_0, node_Ycoord_0, 'bo', 'MarkerSize', radius*radius_Scale);
        
        plot(node_Xcoord_1, node_Ycoord_1, 'go', 'MarkerFaceColor', 'g');
        %plot(node_Xcoord_1, node_Ycoord_1, 'go', 'MarkerSize', radius*radius_Scale);

        plot(node_Xcoord_2, node_Ycoord_2, 'ko', 'MarkerFaceColor', 'k');
        %plot(node_Xcoord_2, node_Ycoord_2, 'ko', 'MarkerSize', radius*radius_Scale);
        

        label=zeros([1,no_nodes]);
        for i=1:no_nodes
            if not(any(cluster_ID(i)==label))
                text(node_Xcoord(i), node_Ycoord(i), num2str(cluster_ID(i)), 'FontSize', 14);
                label(i)=cluster_ID(i);
            end
        end
        
        
        for i=cluster_HD
            if not(i==0)
                text(node_Xcoord(i), node_Ycoord(i), '*', 'FontSize', 14, 'color', 'r');
            end
        end
        
        [~,c]=find(state_of_node==3);
        for i=c
            text(node_Xcoord(i), node_Ycoord(i), '_{FBS}', 'FontSize', 14, 'color', 'r');
        end
        
        plot(BSX, BSY, 'ro', 'MarkerFaceColor', 'r');
        plot(BSX, BSY, 'ro', 'MarkerSize', 2*radius*radius_Scale);
        grid on;
        
        mystr = '(* are cluster heads and numbers are cluster ID: No. Nodes =';
        mystr = strcat(mystr, num2str(no_nodes));
        mystr = strcat(mystr, ')');
        
        title({'Node Deployment', mystr});
    end
end