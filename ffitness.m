function ret=ffitness(index)
global FID
global PDist
global sizepop
global lenchrom
global individuals
global CTN
global PReq
global CTask
global CTaskN
global Tgene

%DecChrom=zeros(1,lenchrom);
% The array of binary into decimal group

fprintf(FID,'I am in fitness lenchrom=%d \n',lenchrom);

%for i=1:10
%    for j=1:10

%    fprintf(FID,'I am in fitness PDist(%d,%d)=%2.2f \n',i,j,PDist(i,j));

%    end
%end
fprintf('===========In Fitness==============')
 
  fprintf('\n individuals(%d)   ',index);   
  fprintf('%1.3f',individuals.fitness(index));
  for j=1:lenchrom
   fprintf(' %d',individuals.DecChrom(index,j));
  end
  fprintf('###');
   for j=1:lenchrom*Tgene
    fprintf(' %d',individuals.BinChrom(index,j));
   end
   fprintf('###');
  for j=1:lenchrom*Tgene
    fprintf(' %1.2f',individuals.Alphachrom(index,j));
  end
     fprintf('###');
  for j=1:lenchrom*Tgene
    fprintf(' %1.2f',individuals.Betachrom(index,j));
  end
  fprintf('\n');

 fprintf('================================================================\n');
fprintf(FID,'List of patients and care team to vsisit them \n');

for k=1:lenchrom
    
fprintf(FID,'Care team %d should visit patient %d \n',individuals.DecChrom(index,k),k);
    
%drop bad solutions ---------check on PReq and CTask tables
    
Found=0;
    T=1;
    
CT=individuals.DecChrom(index,k);
    
while (T<=CTaskN && Found==0)
        
if(PReq(T,k)==1 && CTask(CT,T)==1)
            
Found=1;
        
end
        
T=T+1;
    
end
      
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (Found==1)

DistCGT=zeros(1,CTN);

for k=1:CTN
Prec=0;
    
for lc=1:lenchrom
    
     
if (individuals.DecChrom(index,lc)==k)
        DistCGT(k)=DistCGT(k)+PDist(Prec+1,lc+1);
        Prec=lc;
end
end
lc=1;
DistCGT(k)=DistCGT(k)+PDist(lc,Prec+1);
fprintf(FID,'I am in fitness lc=%d Prec=%d\n',lc,Prec);
%fprintf(FID,'I am in fitness DistCGT(%d)=%2.2f distance from patient %d to HHC is PDist(lc,Prec+1)=%2.2f\n',k,DistCGT(k),Prec,PDist(lc,Prec+1));
end
% The summation of all care teams distances
TotDist=0.0;
for k=1:CTN 
    TotDist=TotDist+DistCGT(k);
    fprintf(FID,'Total distance traveled by team %d is %2.2f',k,DistCGT(k));
end

for k=1:sizepop
    for l=1:lenchrom
         fprintf(FID,'%d ',individuals.DecChrom(k,l)); 
    end
    fprintf(FID,'\n');
end
 ret=1/(1+TotDist);
else
 ret=0;
end