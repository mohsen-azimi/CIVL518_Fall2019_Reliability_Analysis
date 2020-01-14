function f_updateMRdamperDispl()

global MRdamper  Z i myStructure

%       -------X========X----------
%       |          MR2            |
%       |                         |
%       Y                         Y
%       \            COM          \
%   MR3 \           o             \ MR4
%       \              O          \
%       Y              COR        Y
%       |                         |
%       |          MR1            |
%       -------X========X----------

%   (y)
%    |
%    |
%    |
%    -------> (x)
%% Update new Locations MR-dampers
for story = 1:myStructure.ns;  % place MRdamper on all stories
    for MR = 1:4               %  4 Dampers: in [X X  Y Y]-dir
% MRdamper{story,MR}.position =  [12,0;18,0];% Place MR-damper in X-dir & lower one [xi,yi; xj,yj]

X     = Z(3*story-2,i);  % displacement in X-dir for Story(i)
Y     = Z(3*story-1,i);  % displacement in X-dir for Story(i)
Theta = Z(3*story  ,i);  % rotation  for Story(i)

vel_X     = Z(myStructure.n + 3*story-2,  i);  % velocity in X-dir for Story(i)
vel_Y     = Z(myStructure.n + 3*story-1,  i);  % velocity in X-dir for Story(i)
vel_Theta = Z(myStructure.n + 3*story  ,  i);  % rotation-velocity  for Story(i)



dx_myStructure.CORi = abs(MRdamper{story,MR}.position(1)-myStructure.COR(story,1));  % delta-X from Center of stiffness (i)
dy_myStructure.CORi = abs( MRdamper{story,MR}.position(2)-myStructure.COR(story,1));  % delta-Y from Center of stiffness (i)

dx_myStructure.CORj = abs( MRdamper{story,MR}.position(3)-myStructure.COR(story,1));  % delta-X from Center of stiffness (j)
dy_myStructure.CORj = abs( MRdamper{story,MR}.position(4)-myStructure.COR(story,1));  % delta-X from Center of stiffness (j)


    end
    
 
    MRdamper{story,1}.displT(i) = (X+Theta*(dy_myStructure.CORi)+X+Theta*(dy_myStructure.CORi))/2;     % total displ; in X-dir
    MRdamper{story,2}.displT(i) = (X+Theta*(dy_myStructure.CORj)+X+Theta*(dy_myStructure.CORj))/2;     % total displ; in X-dir   
    MRdamper{story,3}.displT(i) = (Y+Theta*(dx_myStructure.CORi)+Y+Theta*(dx_myStructure.CORi))/2;     % total displ; in Y-dir   
    MRdamper{story,4}.displT(i) = (Y+Theta*(dx_myStructure.CORj)+Y+Theta*(dx_myStructure.CORj))/2;     % total displ; in Y-dir  

    MRdamper{story,1}.velT(i) = (vel_X+vel_Theta*(dy_myStructure.CORi)+vel_X+vel_Theta*(dy_myStructure.CORi))/2;     % total vel; in X-dir
    MRdamper{story,2}.velT(i) = (vel_X+vel_Theta*(dy_myStructure.CORj)+vel_X+vel_Theta*(dy_myStructure.CORj))/2;     % total vel; in X-dir   
    MRdamper{story,3}.velT(i) = (vel_Y+vel_Theta*(dx_myStructure.CORi)+vel_Y+vel_Theta*(dx_myStructure.CORi))/2;     % total vel; in Y-dir   
    MRdamper{story,4}.velT(i) = (vel_Y+vel_Theta*(dx_myStructure.CORj)+vel_Y+vel_Theta*(dx_myStructure.CORj))/2;     % total vel; in Y-dir  

    
   if story==1 
    MRdamper{story,1}.displ(i) = MRdamper{story,1}.displT(i);     % inter-story displ;in X-dir
    MRdamper{story,2}.displ(i) = MRdamper{story,2}.displT(i);     % inter-story displ;in X-dir   
    MRdamper{story,3}.displ(i) = MRdamper{story,3}.displT(i);     % inter-story displ;in Y-dir   
    MRdamper{story,4}.displ(i) = MRdamper{story,4}.displT(i);     % inter-story displ;in Y-dir  
    
    MRdamper{story,1}.vel(i) = MRdamper{story,1}.velT(i);     % inter-vel displ;in X-dir
    MRdamper{story,2}.vel(i) = MRdamper{story,2}.velT(i);     % inter-vel displ;in X-dir   
    MRdamper{story,3}.vel(i) = MRdamper{story,3}.velT(i);     % inter-vel displ;in Y-dir   
    MRdamper{story,4}.vel(i) = MRdamper{story,4}.velT(i);     % inter-vel displ;in Y-dir  

   else
    
    
    MRdamper{story,1}.displ(i) = MRdamper{story,1}.displT(i)-MRdamper{story-1,1}.displT(i);     % in X-dir
    MRdamper{story,2}.displ(i) = MRdamper{story,2}.displT(i)-MRdamper{story-1,2}.displT(i);     % in X-dir   
    MRdamper{story,3}.displ(i) = MRdamper{story,3}.displT(i)-MRdamper{story-1,3}.displT(i);     % in Y-dir   
    MRdamper{story,4}.displ(i) = MRdamper{story,4}.displT(i)-MRdamper{story-1,4}.displT(i);     % in Y-dir  

    MRdamper{story,1}.vel(i) = MRdamper{story,1}.velT(i) - MRdamper{story-1,1}.velT(i);     % in X-dir
    MRdamper{story,2}.vel(i) = MRdamper{story,2}.velT(i) - MRdamper{story-1,2}.velT(i);     % in X-dir   
    MRdamper{story,3}.vel(i) = MRdamper{story,3}.velT(i) - MRdamper{story-1,3}.velT(i);     % in Y-dir   
    MRdamper{story,4}.vel(i) = MRdamper{story,4}.velT(i) - MRdamper{story-1,4}.velT(i);     % in Y-dir  

   end
   
end
    











end