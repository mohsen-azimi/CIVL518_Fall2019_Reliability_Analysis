%% Reliability analisys
% clear; clc; 
close all; 
%%
Nv   = 5000; % This is number of Random Structures


ns = 15;
eqs = {'ElCentro'; 'ChiChi'; 'Kobe';'LomaPrieta';'Northridge';'SanFernando'; 'Erzincan'};
ctrls = {'UnControlled';'PassiveMax';'Active'};

for ctrl = [1 2 3]
for eq = 5
for nv=1:Nv
    nv
    load(['D:\Results\',ctrls{ctrl},'_',eqs{eq},'\SoftSoil_',ctrls{ctrl},'_',eqs{eq},'_Nv_',num2str(nv),'_3000.mat'],'Controlled')
%     Response{ij}=Controlled;
    
    
    for story=1:ns
        Max.displX(nv, story) = max(Controlled.displ(3*story-2,:));
        Max.displY(nv, story) = max(Controlled.displ(3*story-1,:));
        
        Max.intrstDrftX(nv, story) = max(Controlled.interstoryDrift(3*story-2,:))/3.2;  % ratio
        Max.intrstDrftY(nv, story) = max(Controlled.interstoryDrift(3*story-1,:))/3.2;  % ratio
        
    end

end



Nv = nv; 

%%
% Reliability

LSF.intrstDrft     = 0.02;
margSF = 1.2; %just for ploting


for story=1:ns
    ....______________________________________________________
    Dir = 0;
    for lsf = 0 :(LSF.intrstDrft)/1000: margSF*LSF.intrstDrft
        Dir = Dir+1;
        Rel.intrstDrftX(Dir,story)=sum(Max.intrstDrftX(:,story)<lsf)/Nv;
        Rel.intrstDrftY(Dir,story)=sum(Max.intrstDrftY(:,story)<lsf)/Nv;

        LSF.intrstDrft_lsf(Dir,1) = lsf;
    end
    ....______________________________________________________
    ....______________________________________________________
    
end

%% Plots
colors = [[0.00,0.45,0.74];  
          [0.85,0.33,0.10];
          [0.49,0.18,0.56];
          [0.47,0.67,0.19];
          [0.00,0.00,1.00]];

styles  = {'-';'-';'-';'-';'-';
           '-.';'-.';'-.';'-.';'-.';
           '--';'--';'--';'--';'--';};

       figure(ctrl);
       for Dir=1:2
           subplot(1,2,Dir)
           
           for story = 1:ns
               
               if Dir ==1
                   plot(100*LSF.intrstDrft_lsf,Rel.intrstDrftX(:,story),...
                       'linestyle',styles{story},'color',colors(5-rem(story,5),:),...
                       'linewidth',1.0,...
                       'DisplayName',['Story ', num2str(story)]);
               else
                   plot(100*LSF.intrstDrft_lsf,Rel.intrstDrftY(:,story),...
                       'linestyle',styles{story},'color',colors(5-rem(story,5),:),...
                       'linewidth',1.0,...
                       'DisplayName',['Story ', num2str(story)]);
                   
               end
               hold on; grid off; box on;
               
           end
           
           set(gca, 'LineWidth',1, 'FontWeight','normal', 'FontName','Times New Roman', 'FontSize',10)
           ylabel('Reliability (1-P_f)', 'fontsize',12,'fontname','Times New Roman','FontWeight','Bold')
           xlabel('LSF, Inter-story drift ratio (%)', 'fontsize',12,'fontname','Times New Roman','FontWeight','Bold')
           
           if Dir==1
               title({[ctrls{ctrl},', ',eqs{eq}];' (Inter-story Drift in {\itx}-dir)'}, 'fontsize',12,'fontname','Times New Roman','FontWeight','Bold');
               legend('show','location','southeast','numcolumns',3);

           else
               title({[ctrls{ctrl},', ',eqs{eq}];' (Inter-story Drift in {\ity}-dir)'}, 'fontsize',12,'fontname','Times New Roman','FontWeight','Bold');
               
           end
%            xlim([0 4])
       end


end

end
