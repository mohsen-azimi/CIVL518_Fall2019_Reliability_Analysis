%% Reliability analysys
%% 
clear;clc; close all; 
%%
global paramMRdamper  paramStructure paramSSI ij myStructure Controlled Nv MRdamper mytry earthquake

tic
for Nv   = [5000] % This is number of Random Structures

g    = 9.806;   
dt   = 0.01;
Tend = 40;
eartquakes = {'ElCentro';'Northridge';'Kobe';'ChiChi';'SanFernando';'LomaPrieta';'Erzincan'};

ctrlMethod = {'UnControlled'; 'Active';'PassiveMax'};%;'FLC';'ClippedOptLQR'};


% Random Parameters 
myStructure.ns = 15;   ns = myStructure.ns; % 
myStructure.r  = [1,1,1,1,1,1,1,1,1,1,1,1,1,1,1]'; % which stories has MR-damper? 1=true; 0=false


paramStructure.m1 = f_Rnd(0.10,37500  * 2);            
paramStructure.m2 = f_Rnd(0.10,37500  * 1); 
paramStructure.k1 = f_Rnd(0.10,600000 * 2);            
paramStructure.k2 = f_Rnd(0.10,600000 * 1);            

paramMRdamper.c0A    = f_Rnd(0.05,50.30*1000);     % N.s/m
paramMRdamper.c0B    = f_Rnd(0.05,48.70*1000);     % N.s/m
paramMRdamper.c1A    = f_Rnd(0.05,8106.20*1000);   % N.s/m
paramMRdamper.c1B    = f_Rnd(0.05,7807.90*1000);   % N.s/m V
paramMRdamper.k0     = f_Rnd(0.05,0.0054*1000);    % N/m
paramMRdamper.k1     = f_Rnd(0.05,0.0087*1000);    % N/m
paramMRdamper.alphaA = f_Rnd(0.05,8.70*1000);      % N
paramMRdamper.alphaB = f_Rnd(0.05,6.40*1000);      % N
paramMRdamper.gamma  = f_Rnd(0.05,496);            % m^-2
paramMRdamper.beta   = f_Rnd(0.05,496);            % m^-2 
paramMRdamper.A      = f_Rnd(0.05,810.50);         
paramMRdamper.n      = f_Rnd(0.05,2);              
paramMRdamper.eta    = f_Rnd(0.05,190);            % s^-1
paramMRdamper.x0     = f_Rnd(0.05,0.0);            % m 
paramMRdamper.Vmax   = f_Rnd(0.05,5);                %vmax ~ 9     
paramMRdamper.Fmax   = f_Rnd(0.05,500*10^3);         %Fmax ~500



% %   SSI Rnd Prameters
paramSSI.ksx = f_Rnd(0.08, 9.32e7); 
paramSSI.ksy = f_Rnd(0.08, 9.32e7);

paramSSI.kt  = f_Rnd(0.08, 1.08e5); 

paramSSI.krx = f_Rnd(0.08, 1.06e5); 
paramSSI.kry = f_Rnd(0.08, 1.06e5); 

paramSSI.csx = f_Rnd(0.08, 5.24e4); 
paramSSI.csy = f_Rnd(0.08, 5.24e4); 

paramSSI.ct  = f_Rnd(0.08, 1.65e4); 

paramSSI.crx = f_Rnd(0.08, 1.68e4);  
paramSSI.cry = f_Rnd(0.08, 1.68e4); 
...___________________________________________




for myctrl=[1:3]
    for SSI = [true]
        MR  = true;  %dont change this (I guess!)
        for myEQ = 4

            earthquake = f_loadEarthquake(eartquakes{myEQ},g,dt,Tend);
            
            myDir=['D:\Results\',ctrlMethod{myctrl},'_',earthquake.name,'\'];
            mkdir(myDir);
            
            for ij=1:Nv
                placeMRdampers
                f_myStructure(SSI, MR)
                f_controller_Linear(ctrlMethod{myctrl})
                save([myDir,myStructure.SoilType,'_',ctrlMethod{myctrl},'_',earthquake.name,'_Nv_',num2str(ij),'_',num2str(Nv),'.mat'],'Controlled')
            end
            
        end
    end
end

%%
end
% Reliability

elpst=toc; 
    
    
    
    
