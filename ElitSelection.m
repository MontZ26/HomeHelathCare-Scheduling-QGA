function ElitSelection()
global sizepop
global lenchrom
global Nindividuals
global Oindividuals
global newfitnessT
global oldfitnessT
global Tgene
NT=struct('NIF',zeros(1,sizepop),'NEWF',zeros(1,sizepop));
OT=struct('OIF',zeros(1,sizepop),'OLDF',zeros(1,sizepop));
%INDEXT=zeros(1,sizepop);
fprintf('==========================================\n');
fprintf('======SNFT Old Fitness Table before sorting==============\n');
fprintf('==========================================\n');
for k=1:sizepop
    fprintf('k=%d oldfitnessT(k)=%1.3f\n',k,oldfitnessT(k));
end
fprintf('==========================================\n');
fprintf('======SOFT NEW Fitness Table before sorting==============\n');
fprintf('==========================================\n');

for k=1:sizepop
    fprintf('k=%d  newfitnessT(k)=%1.3f \n',k,newfitnessT(k));
end
% suppose 'oldfitnessT' is the struct array. 'oldfitnessT' is the field that contains date and time.
[SofT,INDEXT] = sort(oldfitnessT,'descend'); % sort the table by 'oldfitness'in descending order
fprintf('==================================================\n');
    %fprintf('OI=%d ',OI);
        fprintf('OI=%d ',INDEXT);

 fprintf('==================================================\n');
 %INDEXT=OI;
for k=1:sizepop
    fprintf('k=%d SofT(k)=%1.3f\n',k,SofT(k));
    OT.OIF(k)=INDEXT(k);
    OT.OLDF(k)=SofT(k); 
    fprintf('k=%d  OT.OIF(k)=%d OT.OLDF(k)=%1.3f \n',k,OT.OIF(k),OT.OLDF(k));
end
% suppose 'newfitnessT' is the struct array. 'newfitnessT' is the field that contains date and time.
[SnfT,NI] = sort(newfitnessT); % sort the table by 'newfitness'
fprintf('==========================================\n');
fprintf('==================================================\n');
    fprintf('NI=%d ',NI);
 fprintf('==================================================\n');
  %INDEXT=OI;
for k=1:sizepop
    fprintf('k=%d  SnfT(k)=%1.3f \n',k,SnfT(k));
    NT.NIF(k)=INDEXT(k);
    NT.NEWF(k)=SnfT(k); 
    fprintf('k=%d  NT.NIF(k)=%d NT.NEWF(k)=%1.3f \n',k,NT.NIF(k),NT.NEWF(k));
end
    fprintf('=================END OF RWS==========\n');
fprintf('==========================================');
 %copy of 20% of old population in new population

    fprintf('==========================================')
    s=1;
    for k=1:ceil(0.2*sizepop):ceil(0.8*sizepop)
    SofTindex=OT.OIF(k);
    SnfTindex=NT.NIF(s);
   % fprintf('k=%d Softindex=%d Snftindex=%d\n',k,SofTindex,SnfTindex);
    Nindividuals.fitness(SnfTindex)=Oindividuals.fitness(SofTindex);
    for j=1:lenchrom
            Nindividuals.DecChrom(SnfTindex,j)=Oindividuals.DecChrom(SofTindex,j);
     end
     for j=1:lenchrom*Tgene
            Nindividuals.BinChrom(SnfTindex,j)=Oindividuals.BinChrom(SofTindex,j);
            Nindividuals.Alphachrom(SnfTindex,j)=Oindividuals.Alphachrom(SofTindex,j);
            Nindividuals.Betachrom(SnfTindex,j)=Oindividuals.Betachrom(SofTindex,j);
 
     end
     s=s+1;
    end
    %======================================================================
%Stochastic Universal Selection algorithm is used to copy the remaining 60%
Sumfitness=0;
for k=1:sizepop
     Sumfitness=Sumfitness+oldfitnessT(k);
end
AverageP=Sumfitness/sizepop;
Alpha=rand(0,1);
Sumfitness=oldfitnessT(1);
Delta=Alpha*AverageP;
k=1;
% SU selection of 80% of individuls to join Nindividuals
while 1 
        if (Delta<Sumfitness)
            Delta=Delta+Sumfitness;
            SofTindex=OT.OIF(k);
            Nindividuals.fitness(ceil(0.2*sizepop)+k)=Oindividuals.fitness(SofTindex);
            for j=1:lenchrom 
                Nindividuals.DecChrom(ceil(0.1*sizepop)+k,j)=Oindividuals.DecChrom(SofTindex,j);
            end
            for j=1:lenchrom*Tgene
          
            Nindividuals.BinChrom(ceil(0.2*sizepop)+k,j)=Oindividuals.BinChrom(SofTindex,j);  
            Nindividuals.Alphachrom(ceil(0.2*sizepop)+k,j)=Oindividuals.Alphachrom(SofTindex,j);
            Nindividuals.Betachrom(ceil(0.2*sizepop)+k,j)=Oindividuals.Betachrom(SofTindex,j);

            end
        else
            k=k+1;
            Sumfitness=Sumfitness+oldfitnessT(k);
        end
if (k<ceil(0.8*sizepop))
    break;
end
end %end function ElitSelection