%% Introduce normaly Distributed random numbers

function R_z=f_Rnd(fcov,fIn)

global Nv

[fm,fn]=size(fIn);
R_z1=zeros(fm*fn,Nv);

for fi=1:fm
    for fj=1:fn
        fk=(fi-1)*fn+fj;
        fpar=fIn(fi,fj);
        stdpar=fcov*fpar;
        if fpar==0
            R_z1(fk,:)=0;
        else 
            R_z1(fk,:)=fpar+stdpar*randn(1,Nv);
        end
    end
end

R_z=zeros(Nv*fm,fn);
for fi=1:Nv
    fb=fi*fm;
    fa=fb-fm+1;
    R_z(fa:fb,:)=reshape(R_z1(:,fi),fn,fm)';
end


end







