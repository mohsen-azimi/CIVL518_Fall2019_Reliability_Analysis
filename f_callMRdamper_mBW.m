
 % 1=previous step; 2=current step; 3=next step
 % Simple Bouc Wen Model of MR damper using Runge Kutta 4th order method
function f_callMRdamper_mBW()
%% [f,y3,yd3,z3] = ModBouc_F_fn2(x1,y2,  v1,v2,  z2,y2, yd2,U, dt,kd,b)
global  i MRdamper myStructure earthquake 

%% Update from the previous state (call from the saved global variables)

for story = 1:myStructure.ns;  % place MRdamper on all stories
    for MR = 1:4               %  4 Dampers: in [X X  Y Y]-dir
       x1 = MRdamper{story,MR}.displ(i-1); % velocity from previous time-step
       x2 = MRdamper{story,MR}.displ(i);   % velocity from current  time-step

        
       v1 = MRdamper{story,MR}.vel(i-1); % velocity from previous time-step
       v2 = MRdamper{story,MR}.vel(i);   % velocity from current  time-step

    


................................................
dt = earthquake.dt;

c0A    = MRdamper{story,MR}.c0A;     % N.s/m
c0B    = MRdamper{story,MR}.c0B;     % N.s/m . V

c1A    = MRdamper{story,MR}.c1A;     % N.s/m
c1B    = MRdamper{story,MR}.c1B;     % N.s/m .V

 k0     = MRdamper{story,MR}.k0;      % N.s/m
 k1     = MRdamper{story,MR}.k1;      % N/m

alphaA = MRdamper{story,MR}.alphaA;  % N/m
alphaB = MRdamper{story,MR}.alphaB;  % N/m . V

gamma  = MRdamper{story,MR}.gamma;   % m^-2
beta   = MRdamper{story,MR}.beta;    % m^-2 
A      = MRdamper{story,MR}.A; 
n      = MRdamper{story,MR}.n;
eta    = MRdamper{story,MR}.eta;     % s^-1
x0     = MRdamper{story,MR}.x0;      % m 
     
................................................
z2  = MRdamper{story,MR}.z(i);                   
y2  = MRdamper{story,MR}.y(i);                   
yd2  = MRdamper{story,MR}.yd(i);                   
V   = MRdamper{story,MR}.V(i);
u      = MRdamper{story,MR}.u(i-1); 
uD     = -eta*( u  -  V );
u      =  uD*dt  +  u;



%%
      
% Damper constants computation                       
   c0    = c0A    + c0B.*u;
   c1    = c1A    + c1B.*u;
   alpha = alphaA + alphaB.*u;
   

   %***********Determination of slopes using RK4***************
 
% y(i+1) = y(i)+dy
% dy     = 1/6 (k1+2K2+2k3+k4)h

...   k1=f(x(i),y(i))
...   k2=f(x(i)+h/2, y(i)+k1*h/2)
...   k3=f(x(i)+h/2, y(i)+k2*h/2)
...   k4=f(x(i)=h, y(i)+k3*h)
%     h=dt
%***********Determination of slopes using RK4 alt***************
 
kz1=dt*(  -gamma.*abs(v1-yd2).*z2.*abs(z2).^(n-1)-beta.*(v1-yd2).*abs(z2).^n + A.*(v1-yd2) );
ky1=dt*(  (1./(c0+c1)).*(alpha.*z2 + c0.*v1 +  k0.*(x1-y2)) );

kz2=dt*(  -gamma.*abs((v1+v2)/2-(yd2+ky1/2)).*(z2+kz1/2).*abs(z2+kz1/2).^(n-1)-beta.*( (v1+v2)/2-(yd2+ky1/2) ).*abs(z2+kz1/2).^n + A.*( (v1+v2)/2- (yd2+ky1/2) ) );
ky2=dt*(  (1./(c0+c1)).*(alpha.*(z2+kz1/2) + c0.*(v1+v2)/2 +  k0.*( (x1+x2)/2 - (y2+ky1/2) )));
 

kz3=dt*(  -gamma.*abs((v1+v2)/2-(yd2+ky2/2)).*(z2+kz2/2).*abs(z2+kz2/2).^(n-1)-beta.*( (v1+v2)/2-(yd2+ky2/2) ).*abs(z2+kz2/2).^n + A.*( (v1+v2)/2- (yd2+ky2/2) ) );
ky3=dt*(  (1./(c0+c1)).*(alpha.*(z2+kz2/2) + c0.*(v1+v2)/2 +  k0.*( (x1+x2)/2 - (y2+ky2/2))) );
 


kz4=dt*(  -gamma.*abs(v2)-(yd2+ky3).*(z2+kz3).*abs(z2+kz3/2).^(n-1)-beta.* v2-(yd2+ky3).*abs(z2+kz3).^n + A.* v2 - (yd2+ky3)  );
ky4=dt*(  (1./(c0+c1)).*(alpha.*(z2+kz3) + c0.*v2 +  k0.* x2 - (y2+ky3))  );

                                                  
%Next value of  variable (z)and variable (y)

        
        dz=(kz1+2*kz2+2*kz3+kz4)/6;
        dy=(ky1+2*ky2+2*ky3+ky4)/6;
        z3=z2 + dz;  
        y3=y2 + dy; 
      

%************Damper force using Modified Bouc-Wen model****************
   
  yd3=(1./(c0+c1)).*(alpha.*z3 + c0.*v2 +  k0.*(x2-y3));
  F = c1.*yd3 +  k1.*(x2-x0);

  %******************end********************************************
  F = min(max(-MRdamper{story,MR}.Fmax,F),MRdamper{story,MR}.Fmax);
  
  
  %% Save the current state for the next step
MRdamper{story,MR}.c0(i)     = c0; 
MRdamper{story,MR}.c1(i)     = c1; 
MRdamper{story,MR}.alpha(i)  = alpha; 

MRdamper{story,MR}.z(i+1)    = z3;                  

MRdamper{story,MR}.y(i+1)    = y3;                  
MRdamper{story,MR}.yd(i+1)    = yd3;                  


MRdamper{story,MR}.F(i)      = F; 

MRdamper{story,MR}.u(i)      = u;                  


    end 
    
end
     

end