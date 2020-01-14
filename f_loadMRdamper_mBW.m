%% Parameters are borrowed from "Modeling and control of magnetorheological dampers for seismic response reduction"
%  disp('Device:   Mehdi"s paper      " (N m kg)')
function MR = f_loadMRdamper_mBW()
%% MRdamper Constant Parameters

global paramMRdamper ij earthquake
% Random Parameters put in its place
MR.c0A    = paramMRdamper.c0A(ij);     % N.s/m 
MR.c0B    = paramMRdamper.c0B(ij);     % N.s/m
MR.c1A    = paramMRdamper.c1A(ij);     % N.s/m
MR.c1B    = paramMRdamper.c1B(ij);     % N.s/m V
MR.k0     = paramMRdamper.k0(ij);      % N/m
MR.k1     = paramMRdamper.k1(ij);      % N/m
MR.alphaA = paramMRdamper.alphaA(ij);  % N
MR.alphaB = paramMRdamper.alphaB(ij);  % N
MR.gamma  = paramMRdamper.gamma(ij);   % m^-2
MR.beta   = paramMRdamper.beta(ij);    % m^-2 
MR.A      = paramMRdamper.A(ij);
MR.n      = paramMRdamper.n(ij);
MR.eta    = paramMRdamper.eta(ij);     % s^-1
MR.x0     = paramMRdamper.x0(ij);      % m 

%% capacity

MR.Vmax = paramMRdamper.Vmax(ij);      % volt 
MR.Fmax = paramMRdamper.Fmax(ij);      % N 

%% other Pre-allocations

MR.c0     = [MR.c0A,    zeros(1,max(size(earthquake.t))-1)];               % Preallocations
MR.c1     = [MR.c1A,    zeros(1,max(size(earthquake.t))-1)];               % Preallocations
MR.alpha  = [MR.alphaA, zeros(1,max(size(earthquake.t))-1)];               % Preallocations

MR.z      = zeros(1,max(size(earthquake.t)));                           % Preallocations
MR.y      = zeros(1,max(size(earthquake.t)));                           % Preallocations
MR.yd     = zeros(1,max(size(earthquake.t)));                           % Preallocations

MR.F      = zeros(1,max(size(earthquake.t)));                          % Preallocations
MR.V      = zeros(1,max(size(earthquake.t)));                          % Preallocations
MR.u      = zeros(1,max(size(earthquake.t)));                          % Preallocations



end
