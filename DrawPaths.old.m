function DrawPaths(cptin)
global FID
global lenchrom
global Localminima
global CTN

Patients=[2,4,2,4,0,5,6,4,1;7,6,4,3,5,5,3,0,1];
CTPatient=struct('NVP',[],'ListVP',[]);
X=[];
Y=[];
hold off
fprintf(FID,'Cptin=%d',cptin);
figure('Name','Team Care Traveled paths')

for k=1:cptin
    for l=1:CTN
    CTPatient.NVP(l)=0;
    end
    disp(Localminima.BDecChrom(k));
    
    for l=1:CTN
        indexPatient=0;
        for m=1:lenchrom
              if (Localminima.BDecChrom(k,m)==l)
                indexPatient=indexPatient+1;
                CTPatient.ListVP(l,indexPatient)=m;
              end
        end
        CTPatient.NVP(l)=indexPatient;
    end
    for l=1:CTN
       if (CTPatient.NVP(l)~=0)
            X(1)=0;
            Y(1)=0;
            for k=2:CTPatient.NVP(l)
                X(k)=Patients(1,CTPatient.ListVP(l,k));
                Y(k)=Patients(2,CTPatient.ListVP(l,k));
            end
          %  X(CTPatient.NVP(l)+1)=0;
          %  Y(CTPatient.NVP(l)+1)=0;     
       end
        hold on
        LWidths=k;
        plot(X,Y,'-o','LineWidth',LWidths);
    end 
end  

