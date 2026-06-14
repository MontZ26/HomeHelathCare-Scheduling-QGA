function DrawPaths()
global lenchrom
global Globalminima
global QGlobalminima
global CTN

Patients=[2,4,2,4,0,5,6,4,1;7,6,4,3,5,5,3,0,1];
CTPatient=struct('NVP',[],'ListVP',[]);
QCTPatient=struct('NVP',[],'ListVP',[]);

X=[];
Y=[];
X1=[];
Y1=[];
hold off
figure('Name','Team Care Traveled paths')

 for l=1:CTN
    CTPatient.NVP(l)=0;
    QCTPatient.NVP(l)=0;
 end
 fprintf('\n------------------------------------------\n');
 fprintf('\n---------Global Minima------\n');
 fprintf('\n------------------------------------------\n');

 fprintf('fitness=%1.3f',Globalminima.fitness);
 for m=1:lenchrom
     fprintf('%d  ',Globalminima.DecChrom(m));
 end
 fprintf('\n------------------------------------------\n');
  fprintf('\n---------QGlobal Minima------\n');
 fprintf('\n------------------------------------------\n');
 fprintf('fitness=%1.3f',QGlobalminima.fitness);
 for m=1:lenchrom
     fprintf('%d  ',QGlobalminima.DecChrom(m));
 end

 for l=1:CTN
        indexPatient=0;
        QindexPatient=0;
        for m=1:lenchrom
              if (Globalminima.DecChrom(m)==l)
                indexPatient=indexPatient+1;
                CTPatient.ListVP(l,indexPatient)=m;
              end
              if (QGlobalminima.DecChrom(m)==l)
                QindexPatient=QindexPatient+1;
                QCTPatient.ListVP(l,QindexPatient)=m;
              end
        end
        CTPatient.NVP(l)=indexPatient;
        QCTPatient.NVP(l)=QindexPatient;
  end
  for l=1:CTN
     if (CTPatient.NVP(l)~=0)
            X(1)=0;
            Y(1)=0;
            for k=2:CTPatient.NVP(l)
                X(k)=Patients(1,CTPatient.ListVP(l,k));
                Y(k)=Patients(2,CTPatient.ListVP(l,k));
            end
            X(CTPatient.NVP(l)+1)=0;
            Y(CTPatient.NVP(l)+1)=0;
     end
     if (QCTPatient.NVP(l)~=0)
            X1(1)=0;
            Y1(1)=0;
            for k=2:QCTPatient.NVP(l)
                X1(k)=Patients(1,QCTPatient.ListVP(l,k));
                Y1(k)=Patients(2,QCTPatient.ListVP(l,k));
            end
            X1(QCTPatient.NVP(l)+1)=0;
            Y1(QCTPatient.NVP(l)+1)=0;  
     end
      fprintf('\n==================================');
      for k=1:CTPatient.NVP(l)
                fprintf('X(%d)=%d',k,X(k));
      end
      fprintf('\n==================================');
      for k=1:CTPatient.NVP(l)
                fprintf('Y(%d)=%d',k,Y(k));
      end
      fprintf('\n==================================');
      for k=1:QCTPatient.NVP(l)
                fprintf('X1(%d)=%d',k,X1(k));
      end
      fprintf('\n==================================');
      for k=1:QCTPatient.NVP(l)
                fprintf('Y1(%d)=%d',k,Y1(k));
      end
      %plot(X,Y,'b',X1,Y1,'g--o');
      %plot(X,Y,'*');
      plot(X,Y,X1,Y1);
      grid on

  end 
end  

