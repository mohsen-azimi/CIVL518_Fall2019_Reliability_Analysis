% Based on Chen's book, page ~193; The procedure is repeated for both
% elastic and plastic parts of the base isolation (bilinear)
function f_controller_Linear(myControlCase)

global MRdamper ij Controlled earthquake myStructure Nv mytry Z
%% Load UnControlled results for FLC Controller   
% % % if strcmp(myControlCase,'FLC')== 1
% % %     
% % %     FIS = f_FISCreate2(); % Create and Save FIS for FLC  control
% % %     % load uncontrolled data
% % %     unControlled = load([myDirectory,'\Results & Figures\',...
% % %         myStructure.SoilType,'_','UnControlled','_',earthquake.name,'_111111111111111.mat']);
% % %     
% % %     MaxU.displX  = max(  abs(  unControlled.Controlled.Z (28,:) )  );
% % %     MaxU.displY  = max(  abs(  unControlled.Controlled.Z (29,:) )  );
% % %     if strcmp(myStructure.SoilType,'SoftSoil')== 1
% % %         MaxU.velX    = max(  abs(  unControlled.Controlled.Z (63,:) )  );
% % %         MaxU.velY    = max(  abs(  unControlled.Controlled.Z (64,:) )  );
% % %     else
% % %         MaxU.velX    = max(  abs(  unControlled.Controlled.Z (58,:) )  );
% % %         MaxU.velY    = max(  abs(  unControlled.Controlled.Z (59,:) )  );
% % %         
% % %     end
% % % end

%% Recall some data

xddot_gX = earthquake.xddot_gX;
xddot_gY = earthquake.xddot_gY;
t        = earthquake.t;
dt       = earthquake.dt;

n        = myStructure.n;    % number of stories

......................................................................
% state-space parameters (for Elastic part of bilinear-base-isolation)

A   = myStructure.A;
BrX = myStructure.BrX;
BrY = myStructure.BrY;
Bu  = myStructure.Bu;    
......................................................................

%%

R  = 10^-6*eye(size(Bu,2));              % LQR, R
Q  = 10^6*eye(2*n);                    % LQR, Q

......................................................................    
%      P  = care(A,Bu,Q,R);           % P in Ricatti equation, (Eqn. 4-99)
%      G0 = R^-1*Bu'*P;               % Control Gain Matrix from (Eqn. 4-101)

G_lqr = lqr(A, Bu, Q, R);

%         G = 0*G_lqr;
.......................................................................

%%  Calc T which is the eig vectror of A (See page 504, 512 for examples)
      
[v,~] = eig(A);

for ii = 1:n
    vv(:,ii) = v(:,2*ii-1);
end
.......................................................................  
    
for ii = 1:n
    T_help(:,ii)      =  vv(:,n+1-ii);
    T(:,(2*(ii-1)+1)) = real(T_help(:,ii));
    T(:,(2*(ii-1)+2)) = imag(T_help(:,ii));
end

Phi = T^-1*A*T;

for ii = 1:n
    Phi(ii,n+1:end) = 0;
    Phi(n+1:end,ii) = 0;
end
 .......................................................................           
 
  %% i = 1, t = 0: Initial Conditions  ['(1)' is at t = 0]
  .........................................................
  % 1- Prealloactions: 

i = 1;
u =zeros(size(Bu,2),max(size(xddot_gX)));
Z =zeros(2*n,max(size(xddot_gX)));
    
Psi=zeros(2*n,max(size(xddot_gX)));
Gamma = zeros(2*n,max(size(xddot_gX)));
Lambda =zeros(2*n,max(size(xddot_gX)));

% 2- Assignment (not needed in this case!):
xddot_gX(i) =0;
xddot_gY(i) =0;
    
% From Equations 4-28, 4-33, and 4-39, we have:
Gamma(:,i) = T^-1*Bu*u(:,i) + T^-1*(BrX*xddot_gX(i)+BrY*xddot_gY(i));
Lambda(:,i)= expm(zeros(2*n))*Gamma(:,i);
Psi(:,i)   = T^-1*Z(:,i);

%% i = 2, t = 0+dt = 0.01 sec.             (u = 0 between t = 0  t = dt; u = u(:,1) still!)
......................................................... 
    
i = 2;
Gamma(:,i) = ( T^-1*Bu*u(:,i))+( T^-1*((BrX*xddot_gX(i)+BrY*xddot_gY(i))));   % Eqn 4-33 (Chen's book, pge193)
Psi(:,i) =   Lambda(:,i-1)+Gamma(:,i)*dt/2;                   % Eqn 4-40

Z(:,i) = T*Psi(:,i);                                        % Eqn 4-42
u(:,i+1) = -G_lqr*Z(:,i)*0;     % Eqn 4-23 For LQR METHOD
    
%% i = 3+, t = 0+2dt = 0.02+ sec.
 
for i = 3:max(size(t))
    
    clc;
    disp([myControlCase,'  ',earthquake.name,' Nv=',num2str(ij),'/',num2str(Nv),'    try=',num2str(mytry),'    Time=',num2str(i*earthquake.dt)])
    
    Gamma(:,i)   =( T^-1*Bu*u(:,i) )+( T^-1*(BrX*xddot_gX(i)+BrY*xddot_gY(i)) ); % Eqn 4-33
    Lambda(:,i-1)=expm(Phi*dt)*(Lambda(:,i-2)+Gamma(:,i-1)*dt);                  % Eqn 4-41a??
    
    %Lambda(:,i-1) = expm(Phi*dt)*Gamma(:,i-1)*dt;                               % Example 4.4.2
    
    Psi(:,i) = Lambda(:,i-1)+Gamma(:,i)*dt/2;                                    % Eqn 4-40
    
    % Structural Response
    Z(:,i) = T*Psi(:,i);                                                         % Eqn 4-42
    %% Add Controlling Forces
    %% %%%%%%%%%%% Update MR-damper displ & velocity %%%%%%%%%%%%%%%%%%%%%%%%
    
    if strcmp(myControlCase,'UnControlled')== 1 
        u(:,i+1) = -G_lqr*Z(:,i)*0;                 % Eqn 4-23 For LQR METHOD
    elseif strcmp(myControlCase,'Active')== 1       
        u(:,i+1) = -G_lqr*Z(:,i);                   % Eqn 4-23 For LQR METHOD
    else
        f_updateMRdamperDispl()                     % fisrt update MRdamper displacement and velocity
        switch myControlCase
% % % %             case 'ClippedOptLQR'
% % % %                 f_lqr = -G_lqr*Z(:,i);              % desired control force vector (=ActiveLQR)
% % % %                 
% % % %                 % first, calc the current force for each MR damper. based on the previous voltage                    
% % % %                 for story = 1:myStructure.ns        %  MRdampers on all stories
% % % %                     for MR = 1:4                    %  4 Dampers: in [X X  Y Y]-dir
% % % %                         MRdamper{story,MR}.V(i) = MRdamper{story,MR}.V(i-1);
% % % %                     end
% % % %                 end
% % % %                 
% % % %                 f_callMRdamper_mBW();               % Update all MR-damper forces
% % % %                 % second, use the clippedOpt algorithm Ref: "An Experimental Study of MR Dampers for Seismic Protection" by S. Dyke
% % % %                 
% % % %                 for story = 1:myStructure.ns       %  MRdampers on all stories
% % % %                     for MR = 1:4                    %  4 Dampers: in [X X  Y Y]-dir
% % % %                         Fmr=MRdamper{story,MR}.F(i);
% % % %                         Flqr=f_lqr(4*story-4+MR);
% % % %                         MRdamper{story,MR}.V(i) = heaviside((Flqr-Fmr)*Fmr)* MRdamper{story,MR}.Vmax;
% % % %                     end
% % % %                 end
                
                %%%%%%%%%%%%%%%%%%
% % % %             case 'PassiveMin'
% % % %                 for story = 1:myStructure.ns;  %  MRdampers on all stories
% % % %                     for MR = 1:4               %  4 Dampers: in [X X  Y Y]-dir
% % % %                         MRdamper{story,MR}.V(i) = 0;
% % % %                     end
% % % %                 end
                %%%%%%%%%%%%%%%%%%
            case 'PassiveMax'
                for story = 1:myStructure.ns;  %  MRdampers on all stories
                    for MR = 1:4               %  4 Dampers: in [X X  Y Y]-dir
                        MRdamper{story,MR}.V(i) = MRdamper{story,MR}.Vmax;
                    end
                end
                
                %%%%%%%%%%%%%%%%%%
% % % %             case 'FLC'
% % % %                 %   Calc MRdamper Force using Fuzzy Controller based on x, xdot 
% % % %                 x    = 1*Z(28,i)/MaxU.displX;
% % % %                 y    = 1*Z(29,i)/MaxU.displY;
% % % %                 
% % % %                 if strcmp(myStructure.SoilType,'SoftSoil')== 1         
% % % %                     xdot = 1*Z(63,i)/MaxU.velX;
% % % %                     ydot = 1*Z(64,i)/MaxU.velY; 
% % % %                 else xdot = 1*Z(58,i)/MaxU.velX; 
% % % %                     ydot = 1*Z(59,i)/MaxU.velY;
% % % %                 end
% % % %                 FISinputsX = [x, xdot];         % input range:[-1,1]; [-1,1]
% % % %                 V_x = max((evalfis(FISinputsX,FIS)),0); %                  
% % % %                 FISinputsY = [y, ydot];         % input range:[-1,1]; [-1,1]
% % % %                 V_y = max((evalfis(FISinputsY,FIS)),0); % 
% % % %                         
% % % %             ..................................................
% % % %             for story = 1:myStructure.ns;       %  MRdampers on all stories  
% % % %                 for MR = 1:2                    %  we have 4 Dampers: in [X X  Y Y]-dir
% % % %                     MRdamper{story,MR}.V(i) = 3*V_x*MRdamper{story,MR}.Vmax;% MRdamper{story,MR}.Vmax; % 
% % % %                     
% % % %                 end
% % % %             .................................................
% % % %             for MR = 3:4                        %  we have 4 Dampers: in [X X  Y Y]-dir
% % % %                 MRdamper{story,MR}.V(i) = 3*V_y*MRdamper{story,MR}.Vmax;% MRdamper{story,MR}.Vmax; % 
% % % %             end
% % % %             end
            ......................
            otherwise
            pause
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        f_callMRdamper_mBW();                   % Update all MR-damper forces
        f_MR = zeros(myStructure.ns,1); jj=0;
        for story = 1:myStructure.ns            % place MRdamper on all stories
            for MR = 1:4                        %  4 Dampers: in [X X  Y Y]-dir
                jj=jj+1;
                f_MR(jj) = MRdamper{story,MR}.F(i)*myStructure.r(story);
            end
        end                                     % Apply Control force from all  sources
        u(:,i+1) = -G_lqr*0*Z(:,i);             % Eqn 4-23 For LQR METHOD
        u(:,i+1)   = u(:,i+1)- f_MR(:);         % add MRdamperForce; SSI=True/False
    end
end
 
%%

u(:,i+1) = [];  % for dimension match for  plot

for i = 1:n
    
    % D_max(i) =  max (abs(Z(i,:)))*100;
    % end
    % Export Data corrected
    
    Controlled.Case=myControlCase;
    Controlled.Z=Z;
    
    % Inter-Story Drift; Verlosity     
    Z1 = Z;
    Z2 = Z;
    Z1 = [zeros(3,size(Z2,2));Z1 ];
    Z1(n+1:n+3,:) = 0;
    Z1(end-2:end,:)=[];
    dZ = Z2-Z1;         
    Controlled.displ           = Z(1:n,:);      % absolute    displ
    Controlled.interstoryDrift = dZ(1:n,:);     % inter-story displ
    
    Controlled.vel             = Z(n+1:end,:);  % absolute    vel
    Controlled.interstoryVel   = dZ(n+1:end,:); % inter-story vel
    
% % % Accelarion 
% % clear Z1 Z2 dZ
% %          Z1 = Z;
% %          Z2 = Z;
% %          Z2(:,1) = [];
% %          Z2 = [Z2 zeros(2*n,1)];
% %          Z1(:,end) = zeros(2*n,1);
% %          dZ = Z2-Z1;
% %          dZ(:,end-5:end) = zeros(2*n,6); %  remove the last results
% %     !THIS IN NOT ABSOLUTE ACCEL Controlled.Accel           = dZ(n+1:end,:)/earthquake.dt;         

    Controlled.u = u;
end





end

