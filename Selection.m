function Selection()
global sizepop
global lenchrom
global OGAindividuals
global NGAindividuals
global newfitnessT
global oldfitnessT
global Tgene
global individuals
global GAindividuals
global Oindividuals
global Nindividuals

fprintf('\n======================================================');
fprintf('\n============Starting Selection: I am in Selection====================');
fprintf('\n============Old individuals...Table OGAindividuals====================');
fprintf('\n======================================================');
 for i=1:sizepop
     NGAindividuals.fitness(i)=GAindividuals.fitness(i);
     fprintf('\n OGAindividuals(%d)   ',i);   
     fprintf('%1.3f',OGAindividuals.fitness(i));
  for j=1:lenchrom
   fprintf(' %d',OGAindividuals.DecChrom(i,j));
  end
  fprintf('###');
  for j=1:lenchrom*Tgene
    fprintf(' %d',OGAindividuals.BinChrom(i,j));
  end
   fprintf('###');
   fprintf('\n');
 end

fprintf('\n======================================================');
fprintf('\n============I am in Roulette Wheel====================');
fprintf('\n============New individuals Table GAindividuals====================');
fprintf('\n======================================================');
 for i=1:sizepop
     fprintf('\n GAindividuals(%d)   ',i);   
     fprintf('%1.3f',GAindividuals.fitness(i));
  for j=1:lenchrom
   fprintf(' %d',GAindividuals.DecChrom(i,j));
  end
  fprintf('###');
  for j=1:lenchrom*Tgene
    fprintf(' %d',GAindividuals.BinChrom(i,j));
  end
   fprintf('###');
   fprintf('\n');
 end

NT=struct('NIF',zeros(1,sizepop),'NEWF',zeros(1,sizepop));
OT=struct('OIF',zeros(1,sizepop),'OLDF',zeros(1,sizepop));

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
% suppose 'oldfitnessT' is the struct array. 
[SofT,OI] = sort(oldfitnessT,'descend'); % sort the table by 'oldfitness'in descending order
fprintf('==================================================\n');
for k=1:sizepop
    OT.OIF(k)=OI(k);
    OT.OLDF(k)=SofT(k); 
end
% suppose 'newfitnessT' is the struct array. 
[SnfT,NI] = sort(newfitnessT); % sort the table by 'newfitness'
fprintf('==================================================\n');
for k=1:sizepop
    NT.NIF(k)=NI(k);
    NT.NEWF(k)=SnfT(k); 
end
    fprintf('=================END OF RWS==========\n');
   for k=1:sizepop
      Nindividuals.fitness(k)=GAindividuals.fitness(k);
    for j=1:lenchrom
            Nindividuals.DecChrom(k,j)=GAindividuals.DecChrom(k,j);
     end
     for j=1:lenchrom*Tgene
            Nindividuals.BinChrom(k,j)=GAindividuals.BinChrom(k,j);
     end
    end
fprintf('==========================================');
   %copy of 20% of old population innew population
for k=1:ceil(0.2*sizepop)
    fprintf('k=%d OT.OLDF(%d)=%1.3f',k,OT.OIF(k),OT.OLDF(k));
    fprintf('NT.NEWF(%d)=%1.3f\n',NT.NIF(k),NT.NEWF(k));
    Nindividuals.fitness(NT.NIF(k))=Oindividuals.fitness(OT.OIF(k));
    for j=1:lenchrom
            Nindividuals.DecChrom(NT.NIF(k),j)=Oindividuals.DecChrom(OT.OIF(k),j);
     end
     for j=1:lenchrom*Tgene
            Nindividuals.BinChrom(NT.NIF(k),j)=Oindividuals.BinChrom(OT.OIF(k),j);
            Nindividuals.Alphachrom(NT.NIF(k),j)=Oindividuals.Alphachrom(OT.OIF(k),j);
            Nindividuals.Betachrom(NT.NIF(k),j)=Oindividuals.Betachrom(OT.OIF(k),j);
     end
end
fprintf('\n======================================================');
fprintf('\n============Situation After filling 30 percent in Roulette Wheel===');
fprintf('\n======================================================');
 for i=1:sizepop
     fprintf('\n Nindividuals(%d)   ',i);   
     fprintf('%1.3f',Nindividuals.fitness(i));
  for j=1:lenchrom
   fprintf(' %d',Nindividuals.DecChrom(i,j));
  end
  fprintf('###');
   for j=1:lenchrom*Tgene
    fprintf(' %d',Nindividuals.BinChrom(i,j));
   end
   fprintf('###');
   fprintf('\n');
 end
for i=1:sizepop
       GAindividuals.fitness(i)=Nindividuals.fitness(i);
       individuals.fitness(i)=Nindividuals.fitness(i);
    for j=1:lenchrom
            GAindividuals.DecChrom(i,j)=Nindividuals.DecChrom(i,j);
            individuals.DecChrom(i,j)=Nindividuals.DecChrom(i,j);
     end
     for j=1:lenchrom*Tgene
            GAindividuals.BinChrom(i,j)=Nindividuals.BinChrom(i,j);
            individuals.BinChrom(i,j)=Nindividuals.BinChrom(i,j);
            individuals.Alphachrom(i,j)=Nindividuals.Alphachrom(i,j);
            individuals.Betachrom(i,j)=Nindividuals.Betachrom(i,j);
     end
end

%======================================================================
%Stochastic Universal Selection algorithm is used to copy the remaining 60%
Sumfitness=0;
for k=1:sizepop
     Sumfitness=Sumfitness+newfitnessT(k);
     NGAindividuals.PS(k)=0;
end
AverageP=Sumfitness/sizepop;
Alpha=rand;
Sumfitness=newfitnessT(1);
Delta=Alpha*AverageP;
k=1;
% SU selection of individuls to join Nindividuals
Select_counter=0;
    while 1 
        %fprintf('\n k=%d Select_counter=%d  Delta=%1.3f  Sumfitness=%1.3f',k,Select_counter,Delta,Sumfitness);
        if (Delta<Sumfitness) && (NGAindividuals.PS(k)==0) && (Select_counter <=sizepop)
            Delta=Delta+Sumfitness;
            NGAindividuals.PS(k)=1;
            NGAindividuals.fitness(k)=GAindividuals.fitness(k);
            for j=1:lenchrom 
            NGAindividuals.DecChrom(k,j)=GAindividuals.DecChrom(k,j);
            end
            for j=1:lenchrom*2
            NGAindividuals.BinChrom(k,j)=GAindividuals.BinChrom(k,j);
            end
            Select_counter=Select_counter+1;
        else
            k=k+1;
        end
        if (k>sizepop)&& (Select_counter==sizepop)
            break
        elseif(k>sizepop)&& (Select_counter<sizepop)
            k=1;
            Alpha=rand;
            Sumfitness=newfitnessT(1);
            Delta=Alpha*AverageP;
        end
    end
end %end function Selection