% ... based on Ebrahim's Paper (SSI)
function f_myStructure(SSI,MR)
% %   % defaults are from Ebrahim :2016
% % %   http://www.sciencedirect.com/science/article/pii/S0267726116300835
% %      % main Equation  ==>  [M]xddot(t)+[C]xdot(t)+[K]x(t) = [gamma]u(t)+[delta]xddot_g(t)
% % disp('Unit:           N, m, kg')
% % disp('Model source:   Ebrahim Nazarimofrad 2016 (15-Story/15*3+5-DOF)') 
global  MRdamper  paramStructure paramSSI ij myStructure 

%% Mass Matrix Assembly ||||||||| Unit ??????

% Mass of each Stories : [ Mass, Xdim, Ydim, Xcoord, Ycoord]

m1=paramStructure.m1(ij);   % large mass for light slabs
m2=paramStructure.m2(ij);   % light mass for heavy slabs

... here, we tried to put more weight on right side of the plan
... and more rigidity on the left-down corner...
ns = myStructure.ns;
for story=1:ns
    mm(:,:,story) = [ m1/9.806     5  	5	   2.5     2.5
                      m1/9.806     5  	5	   4.5     2.5
                      m2/9.806     5  	5	   12.5	   2.5
                      m2/9.806     5  	5	   17.5    2.5
                  
                      m1/9.806     5	5	   2.5     4.5
                      m1/9.806     5	5	   4.5     4.5
                      m2/9.806     5  	5	   12.5	   4.5
                      m2/9.806     5  	5	   17.5	   4.5
                     
                      m2/9.806     5  	5	   2.5     12.5
                      m2/9.806     5	5	   4.5     12.5
                      0            5   	5	   12.5    12.5
                      0            5	5      17.5    12.5
                  
                      m2/9.806     5	5	   2.5     17.5
                      m2/9.806     5	5	   4.5     17.5
                      0            5  	5      12.5	   17.5
                      0            5  	5	   17.5	   17.5];         
end 

m=zeros(ns,1);     
for story=1:ns
    m(story)  = sum(mm(:,1,story));   % mass of each floor
end

for i = 1:ns
    COM(i,1) = sum(mm(:,1,i).*mm(:,4,i))/sum(mm(:,1,i)); % Center of Mass (X-dir) (eqn. 21)
    COM(i,2) = sum(mm(:,1,i).*mm(:,5,i))/sum(mm(:,1,i)); % Center of Mass (Y-dir)
end

%--------------- Moment of Inertia bout x,y,z  -------------------------------------------------------
% Iz == I0 in the paper Eqn. 20

for i=1:ns  % (Correct on the other file)
    Ix (i) = sum(  mm(:,1,i)/4.*(mm(:,3,i).^2)...
            + mm(:,1,i).*(( mm(:,4,i)-COM(i,1) ).^2 + ( mm(:,5,i)-COM(i,2) ).^2)); % ??? Why divided by 4?

    Iy (i) = sum(  mm(:,1,i)/4.*(mm(:,2,i).^2)...
            + mm(:,1,i).*(( mm(:,4,i)-COM(i,1) ).^2 + ( mm(:,5,i)-COM(i,2) ).^2)); % ??? Why divided by 4?

    Iz (i) = sum(  mm(:,1,i)/4.*(mm(:,2,i).^2+mm(:,3,i).^2)...
            + mm(:,1,1).*(( mm(:,4,i)-COM(i,1) ).^2 + ( mm(:,5,i)-COM(i,2) ).^2)); % ??? Why divided by 4?
end

Ix = Ix(:);     % Mass Moment of Inertia vector about-x 
Iy = Iy(:);     % Mass Moment of Inertia vector about-y 
Iz = Iz(:);     % Mass Moment of Inertia vector about-Z 

Ixb = Ix(1);    % Foundation Mass Moment of Inertia  about-x 
Iyb = Iy(1);    % Foundation Mass Moment of Inertia  about-y 
Izb = Iz(1);    % Foundation Mass Moment of Inertia  about-z
    
%-------------total mass matrix-----------------------------------
% = M*   Eqn. 19

Ms=[];
for i=1:ns
    Mm = blkdiag(m(i), m(i), Iz(i));  % {Correct on the other file}
    Ms = blkdiag(Ms,Mm);
end

%% Stiffness Matrix Assembly ||||||||| Unit ??????

% Stiffness of each Stories : [ Mass, Xdim, Ydim, Xcoord, Ycoord]

%--------------stiffness of x direction ----------------------
%         [ ID,   Kxj,  Xcoord, Ycoord] (EQn 14) for each column 

Kx = repmat(zeros(21,4),1,1,ns);    % 21 columns; populize the cells (dimention match)
Ky = repmat(zeros(21,4),1,1,ns);    % 21 columns; populize the cells (dimention match)

.........................................
% Story: 1-15

k1 = paramStructure.k1(ij);   
k2 = paramStructure.k2(ij);  
for story=1:ns
    Kx(:,:,story)=[ 1	k1     0     	0
                    2	k1     5	0
                    3	k1     10	0
                    4	k2     15	0
                    5	k2     20	0
                    
                    6	k1     0	    5
                    7	k1     5	5
                    8	k1     10	5
                    9	k2     15	5
                    10	k2     20	5
                
                    11	k1     0	    10
                    12	k1     5	10
                    13	k1     10	10
                    14	k2     15	10
                    15	k2     20	10
                
                    16	k2     0	    15
                    17	k2     5	15
                    18	k2     10	15
                    
                    19	k2     0	    20
                    20	k2     5	20
                    21	k2     10	20];    % (EQn 14)  [ ID, Kxj, Xcoord, Ycoord]

    Ky(:,:,story)=Kx(:,:,story);
end
    
%--------------stiffness of y direction ----------------------
%         [ ID, Kxj,        Xcoord, Ycoord] (EQn 14) 

for i = 1:ns
    Kxx(i) = sum(Kx(:,2,i));
    Kyy(i) = sum(Ky(:,2,i));
end

Kxx = Kxx(:); % Vertical vector
Kyy = Kyy(:); % Vertical vector

%---------Center Of Rigidity of story j ----------------

COR = zeros(ns,2);
for i = 1:ns
    COR(i,1) = sum(Ky(:,2,i).*Ky(:,3,i))/sum(Ky(:,2,i));
    COR(i,2) = sum(Kx(:,2,i).*Kx(:,4,i))/sum(Kx(:,2,i));
end

%-------------------Eccentric------------------------

Ecc=COM-COR;
Ki = repmat(zeros(3),1,1,ns);

for i =1:ns
    Kyt(i)=sum((COM(i,1)-Ky(:,3,i)).*Ky(:,2,i));
    Kxt(i)=sum((COM(i,2)-Kx(:,4,i)).*Kx(:,2,i));
    Ktt(i)=sum((Ky(:,2,i).*(Ky(:,4,i)-COM(i,1)).^2)...
               +(Kx(:,2,i).*(Kx(:,4,i)-COM(i,2)).^2));
    
    Ki(:,:,i)=[Kxx(i)    0     Kxt(i)
                 0     Kyy(i)  Kyt(i)
               Kxt(i)  Kyt(i)  Ktt(i)];
end

% K => K* Eqn 16

Ks  = zeros(3*ns);
Ks(1:3,1:3) =  Ki(:,:,1) + Ki(:,:,2); % corrected!

Ks(end-2:end,end-2:end)  =  Ki(:,:,ns);
Ks(end-2:end,end-5:end-3)= -Ki(:,:,ns);
Ks(end-5:end-3,end-2:end)= -Ki(:,:,ns); 

for i = 2:ns-1
    LL = i*3-2 : i*3;
    Ks(LL,LL-3) = -Ki(:,:,i);
    Ks(LL,LL)   =  Ki(:,:,i) + Ki(:,:,(i+1)) ;
    Ks(LL-3,LL) = -Ki(:,:,i);
end

%%  Damping [C]

Meig = Ms(1:3*(ns),1:3*(ns));
Keig = Ks(1:3*(ns),1:3*(ns));
   
% ...........

[Phi,val] = eig(Keig,Meig);
wn = sqrt(diag(val)); 
Tn = 2*pi()./wn;   

% %%%%%%%%%%%%  Calculation of effective modal mass ratios (m_eff)
% % %     Mhat = Phi'*Meig*Phi; % http://www.vibrationdata.com/tutorials2/ModalMass.pdf
% % %     rhat = ones(3*(ns),1);
% % %     Lbar = Phi'*Meig*rhat;
% % %     Gammai = Lbar./diag(Mhat);
% % %     m_eff  = Lbar.^2./diag(Mhat)/sum(diag(Meig));


xi1 = 0.05; xi2 = 0.05;             % daming ratios of first and second natural frequencies
n1 =1;      n2  = 2;

Ralpha = 2*wn(n1)*wn(n2)*((xi1*wn(n2)-xi2*wn(n1))/(wn(n2)^2-wn(n1)^2));
Rbeta  = 2*(xi2*wn(n2)-xi1*wn(n1))               /(wn(n2)^2-wn(n1)^2);
Cs = Ralpha*Meig+Rbeta*Keig;   % Cs =  Reighly Damping Matrix
            
%%    

if SSI == false
    myStructure.SoilType='FixedBase';
    n = 3*ns;
    %  gamma and delta
    % [gamma]*u(t) (example: [1 -1 0;0 1 -1;0 0 1]) 
    %%%%%%%% gamma %%%%%%%%%%%:
    % gamma=zeros(3*ns,2*ns); % two possible inputs in each floor;
    
    if MR ==true
        % nInput = 4*ns;
        gamma=zeros(n,4*ns); % four possible inputs in each floor;
        
        for story=1:ns  
            
            distuy = [MRdamper{story,1}.position(1,2),MRdamper{story,2}.position(1,2)];
            distux = [MRdamper{story,3}.position(1,1),MRdamper{story,4}.position(1,1)];
            
            % ---------------------------
            By = abs([distux(2)-COM(story,2)    COM(story,2)-distux(1)]);  % Correct
            Bx = abs([distuy(2)-COM(story,1)    COM(story,1)-distuy(1)]);
            
            Di =       [1       1       0           0;
                        0       0       1           1;
                        By(1)  -By(2)   -Bx(1)  Bx(2)]...
                       *myStructure.r(story); % Eqn 11 (*r==> does a story has MR dampers?)
            % ---------------------------  
            
            if story==1
                [rD, cD] = size(Di);
                gamma(1:rD,1:cD) = Di;
            else
                gamma(story*rD-2*rD+1:story*rD, (story-1)*cD+1:story*cD)=[-Di;Di];
            end
        end
    end
    
    %% delta:
    
    for j=1:ns
        mX(3*j-2:3*j,1) = [sum(mm(:,1,j));0;0]; % effective mass in X-dir  Corrected
        mY(3*j-2:3*j,1) = [0;sum(mm(:,1,j));0]; % effective mass in X-dir
    end
    
    deltaX = -mX;
    deltaY = -mY;
    
    %% %%%%%%%%%%%%%%%%%%%%% State Space Parameters  %%%%%%%%%%%%%%%%%%%%%%% 
    % ******************** Input A,B,C,D Matrixes *********************
    % State Vector: {Z} = {x;xdot}
    %               {Zdot} = [A]Z+Bu{u}+{Br}xddot_g
    
    O = zeros(n);
    I = eye(n);
    A = [O            I
        (-Ms^-1)*Ks   (-Ms^-1)*Cs];
    Bu = [ zeros(size(gamma))
            (Ms^-1)*gamma  ];    % shows size(mm,1)mber of inputs
    BrX = [zeros(n,1)
          (Ms^-1)*deltaX];
    BrY = [zeros(n,1)
          (Ms^-1)*deltaY];
    % Br_x, Br_y, Br_theta (Be careful!; be recalled from Controlled_LQR; Br*Xddot (summ up these))
    B = [   O
         (Ms^-1) ];          % Refer to Nikoo's files
    C = [I  O
         O  I];   % 2nx2n==> 2 outputs (=x, x_dot)
    D = zeros(2*n,n) ;  % out put again
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
if SSI==true
    n = 3*ns+5;
    %% Soil soft Properties
    hi = [1:ns]'*3.2; % Story heights
    
    % mb=25000;
    %             mb  = mm(1,1,1);  % same as the 1st floor
    myStructure.SoilType = 'SoftSoil';
    
    %           Swaying     Twisting	Rocking
    ...---------------------------------------------
        % r         9.77        8.64        10.39
    % k         9.32E+08	6.18E+10	1.06E+11
    % Damping	5.24E+07	8.01E+08	1.65E+09
    
    %          ksx = 6.15e7;  kt = 1.68e9;  krx = 2.19e9;
    %          ksy = 6.15e7;                kry = 2.19e9;
    % ....................................................
    %          csx = 1.15e7;  ct  = 5.31e7; crx = 7.99e7;       
    %          csy = 1.15e7;                cry = 7.99e7;       
    
    mb  = mm(1,1,1);
    
%     ksx = 9.32e7;
    ksx = paramSSI.ksx(ij);

%     ksy = 9.32e7;
    ksy = paramSSI.ksy(ij);
%     
%     kt  = 1.08e5;
    kt  = paramSSI.kt(ij);
%     
%     krx = 1.06e5;
    krx = paramSSI.krx(ij);
%     kry = 1.06e5;
    kry = paramSSI.kry(ij);
%     
%     csx = 5.24e4;
    csx = paramSSI.csx(ij);
%     csy = 5.24e4;
    csy = paramSSI.csy(ij);
%     
%     ct  = 1.65e4;
    ct  = paramSSI.ct(ij);
%     
%     crx = 1.68e4;
    crx = paramSSI.crx(ij);
%     cry = 1.68e4;
    cry = paramSSI.cry(ij);
    
    % mb=25000;ksx=6.15e7;ksy=6.15e7;kt=1.68e9;krx=2.19e9;kry=2.19e9;csx=1.15e7;csy=1.15e7;ct=5.31e7;crx=7.99e7;cry=7.99e7;
    % % myStructure.SoilType = 'DenseSoil';
    % % ksx = @@;
    % % ksy = @@;
    % % kt  = @@;
    % % krx = @@;
    % % kry = @@;
    % % csx = @@;
    % % csy = @@;
    % % ct  = @@;
    % % crx = @@;
    % % cry = @@;
    % % mb=25000;
    
    Ksoil = diag([ksx;ksy;kt;krx;kry]);
    Csoil = diag([csx;csy;ct;crx;cry]);
    
    %% M ==> M+Soil ; K ==> K+Soil  (with SSI)
    
    for i = 1:ns
        Ms(3*ns+1 , i*3-2  ) = m(i)*0.1;
        Ms(3*ns+2 , i*3-1  ) = m(i)*0.1;
        Ms(3*ns+3 , i*3    ) = Iz(i)*0.1;
        
        Ms(i*3-2  , 3*ns+1 ) = m(i)*0.1;
        Ms(i*3-1  , 3*ns+2 ) = m(i)*0.1;
        Ms(i*3    , 3*ns+3 ) = Iz(i)*0.1;
        
        Ms(i*3-2  , 3*ns+4 ) = m(i)*hi(i)*0.1;
        Ms(i*3-1  , 3*ns+5 ) = m(i)*hi(i)*0.1;
        Ms(3*ns+4 , i*3-2  ) = m(i)*hi(i)*0.1;
        Ms(3*ns+5 , i*3-1  ) = m(i)*hi(i)*0.1;
    end
    
    Ks = blkdiag(Ks,Ksoil);
    Cs = blkdiag(Cs,Csoil);
    
    Ms(3*ns + 1, 3*ns + 1) = sum(m)*0.1 + mb; 
    Ms(3*ns + 2, 3*ns + 2) = sum(m)*0.1 + mb;
    
    Ms(3*ns + 4, 3*ns + 1) = sum(m.*hi)*0.1;
    Ms(3*ns + 5, 3*ns + 2) = sum(m.*hi)*0.1;
    
    Ms(3*ns + 1, 3*ns + 4) = sum(m.*hi)*0.1;
    Ms(3*ns + 2, 3*ns + 5) = sum(m.*hi)*0.1;
    
    Ms(3*ns + 3, 3*ns + 3) = sum(Iz) + Izb;
    Ms(3*ns + 4, 3*ns + 4) = sum(Iy + m.*hi.^2*0.1) + Iyb;
    Ms(3*ns + 5, 3*ns + 5) = sum(Ix + m.*hi.^2*0.1) + Ixb;
    
    
    %%  gamma and delta (with SSI)
    % [gamma]*u(t) (example: [1 -1 0;0 1 -1;0 0 1]) 
    ...............%%%%%%%%%%% gamma:
    gamma=zeros(3*ns+5,4*ns); % Eqn 9; 10  (5 is for 5-DOFs of the soil)
    
    for story=1:ns   
        distuy = [MRdamper{story,1}.position(1,2),MRdamper{story,2}.position(1,2)];
        distux = [MRdamper{story,3}.position(1,1),MRdamper{story,4}.position(1,1)];
        
        % ---------------------------
        By = abs([distux(2)-COM(story,2)    COM(story,2)-distux(1)]);  % Correct
        Bx = abs([distuy(2)-COM(story,1)    COM(story,1)-distuy(1)]);
        
        Di =       [1       1       0           0;
                    0       0       1           1;
                    By(1)  -By(2)   -Bx(1)  Bx(2)]...
                   *myStructure.r(story); % Eqn 11 (*r==> does a story has MR dampers?)
        % ---------------------------    
        
        if story==1
            [rD, cD] = size(Di);
            gamma(1:rD,1:cD) = Di;
            gamma(3*ns+1:3*ns+3,1:cD) = -Di;  % effect of the controllers on the soil DOFs, only sway and twisting is considered
        else
            gamma(story*rD-2*rD+1:story*rD, (story-1)*cD+1:story*cD)=[-Di;Di];
        end
    end
    
    % gamma = blkdiag(gamma, zeros(5,5));  % Soil DOFs effect
    ................%%%%%%%%%%% delta:
    for j=1:ns
        mX(3*j-2:3*j,1) = [sum(mm(:,1,j));0;0]; % effective mass in X-dir
        mY(3*j-2:3*j,1) = [0;sum(mm(:,1,j));0]; % effective mass in X-dir
    end
    
    deltaX = [-mX;-[Ms(3*ns+1,3*ns+1);0;0;sum(m.*hi);0]];  % for Soild DOFs==>zeros(5,1) (Correct this on...)
    deltaY = [-mY;-[0;Ms(3*ns+2,3*ns+2);0;0;sum(m.*hi)]];  % corrected
    
    %% %%%%%%%%%%%%%%%%%%% State Space Parameters (with SSI)%%%%%%%%%%%%%%%%%% 
    % ******************** Input A,B,C,D Matrixes ****************************
    %                      State Vector: {Z} = {x;xdot}
    %                     {Zdot} = [A]Z+Bu{u}+{Br}xddot_g
    
    O = zeros(n);
    I = eye(n);
    
    A = [O            I
        (-Ms^-1)*Ks   (-Ms^-1)*Cs];
    Bu = [zeros(size(gamma))
         (Ms^-1)*gamma];    % shows size(mm,1)mber of inputs
    BrX = [zeros(n,1)
          (Ms^-1)*deltaX];  %maybe wrone! need correction Br must =-1
    %                 BrX = [zeros(n,1)
    %                       (diag(diag(Ms)))^-1*deltaX];  % Wrong! need correction Br must =-1
    BrY = [zeros(n,1)
          (Ms^-1)*deltaY];
    B = [   O
         (Ms^-1) ];          % Refer to Nikoo's files
    C = [I  O
         O  I];   % 2nx2n==> 2 outputs (=x, x_dot)
    D = zeros(2*n,n+ns) ;  % out put again
    
end    

%%

myStructure.ns = ns;
myStructure.n  = n;
myStructure.hasSSI  = SSI;
myStructure.hasMR   = MR;

% myStructure.nInput  = nInput;

myStructure.m  = m;   

% myStructure.m_eff  = m_eff;    

myStructure.Ms = Ms;
myStructure.Ks = Ks;
myStructure.Cs = Cs;
myStructure.COM = COM;
myStructure.COR = COR;
myStructure.Ecc = Ecc;

myStructure.Tn = Tn;
myStructure.wn = wn;

myStructure.gamma = gamma;
myStructure.deltaX = deltaX;
myStructure.deltaY = deltaY;

myStructure.A   = A;
myStructure.Bu  = Bu;
myStructure.BrX = BrX;
myStructure.BrY = BrY;    
myStructure.B   = B;
myStructure.C   = C;
myStructure.D   = D;

end


