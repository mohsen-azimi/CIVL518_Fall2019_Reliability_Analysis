%% Reliability analisys
% clear; clc; 
close all; 
%%
Nv   = 5000; % This is number of Random Structures
mytry = 1;


ns = 15;



Nv = nv; 

%%
% Reliability

LSF.intrstDrft     = 0.02;
margSF = 1.1; %just for ploting


for story=1:ns
    ....______________________________________________________
    i = 0;
    for lsf = 0 :(LSF.intrstDrft)/1000: margSF*LSF.intrstDrft
        i = i+1;
        Rel.intrstDrftX(i,story)=sum(Max.intrstDrftX(:,story)<lsf)/Nv;
        Rel.intrstDrftY(i,story)=sum(Max.intrstDrftY(:,story)<lsf)/Nv;

        LSF.intrstDrft_lsf(i,1) = lsf;
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

       close all; figure(1);
       for i=1:2
           subplot(1,2,i)
           
           for story = 1:ns
               
               if i ==1
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
           
           if i==1
               title({'Reliabiltiy Analysis';' (Inter-story Drift in {\itx}-dir)'}, 'fontsize',12,'fontname','Times New Roman','FontWeight','Bold');
               legend('show','location','southeast','numcolumns',3);

           else
               title({'Reliabiltiy Analysis';' (Inter-story Drift in {\ity}-dir)'}, 'fontsize',12,'fontname','Times New Roman','FontWeight','Bold');
               
           end
           
       end

