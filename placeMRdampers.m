% global    myStructure MRdamper ...
global MRdamper myStructure ij 
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
............................................................................................................

MRdamper=[];
for story = 1:myStructure.ns  % place MRdamper on all stories
    for Damper = 1:4               %  4 Dampers: in [X X  Y Y]-dir
        MRdamper{story,Damper} = f_loadMRdamper_mBW();
    end
    
    MRdamper{story,1}.position = [5,  0;   10 ,0 ]; MRdamper{story,1}.position = [5,  0;   10 ,0 ];   % Place MR-damper in X-dir & lower one [xi,yi; xj,yj]
    MRdamper{story,2}.position = [5, 15;   10 ,15]; MRdamper{story,2}.position = [5, 15;   10 ,15];   % Place MR-damper in X-dir & upper one [xi,yi; xj,yj]
    MRdamper{story,3}.position = [0,  5;   0  ,10]; MRdamper{story,3}.position = [0,  5;   0  ,10];   % Place MR-damper in Y-dir & left  one [xi,yi; xj,yj]
    MRdamper{story,4}.position = [15, 5;   15 ,10]; MRdamper{story,4}.position = [15, 5;   15 ,10];   % Place MR-damper in Y-dir & right one [xi,yi; xj,yj]
    
end


