function Roulette_Wheel_Selection()
global sizepop
global lenchrom
global GAindividuals
global OGAindividuals
global NGAindividuals
global newfitnessT
global oldfitnessT
global Tgene
fprintf('\n============I just started Roulette Wheel=================');
fprintf('\n============Table GAindividuals====================');
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
for i=1:sizepop
    fprintf('\n oldfitnessT(%d)=%1.3f',i,oldfitnessT(i));
end
fprintf('==========================================\n');
fprintf('======SOFT NEW Fitness Table before sorting==============\n');
fprintf('==========================================\n');

for i=1:sizepop
    fprintf('\n newfitnessT(%d)=%1.3f ',i,newfitnessT(i));
end
fprintf('####################################################################################\n');
% suppose 'oldfitnessT' is the struct array. 
[SofT,OI] = sort(oldfitnessT,'descend'); % sort the table by 'oldfitness'in descending order
fprintf('==================================================\n');
for i=1:sizepop
    OT.OIF(i)=OI(i);
    OT.OLDF(i)=SofT(i); 
end
% suppose 'newfitnessT' is the struct array. 
[SnfT,NI] = sort(newfitnessT,'descend'); % sort the table by 'newfitness'in descending order
fprintf('==================================================\n');
for i=1:sizepop
    NT.NIF(i)=NI(i);
    NT.NEWF(i)=SnfT(i); 
end
fprintf('\n==========================================\n');
fprintf('\n======SNFT Old Fitness Table after sorting==============\n');
fprintf('\n==========================================\n');
for i=1:sizepop
  
    fprintf('\n OldSortedFitnessT(%d)=%1.3f',i,SofT(i));
end
for i=1:sizepop
       fprintf('\n NewSortedFitnessT(%d)=%1.3f ',i,SnfT(i));
end
fprintf('==========================================');
   %copy of old sorted ascendent population in new population
   
%for i=1:ceil(0.2*sizepop):ceil(0.8*sizepop)
for i=1:ceil(0.2*sizepop)
   fprintf('k=%d OT.OLDF(%d)=%1.3f',i,OT.OIF(i),OT.OLDF(i));
   fprintf('NT.NEWF(%d)=%1.3f\n',NT.NIF(i),NT.NEWF(i));
   NGAindividuals.fitness(NT.NIF(sizepop-(i-1)))= OGAindividuals.fitness(OT.OIF(i));
   for j=1:lenchrom
       NGAindividuals.DecChrom(NT.NIF(sizepop-(i-1)),j)=OGAindividuals.DecChrom(OT.OIF(i),j);
   end
   for j=1:lenchrom*Tgene
       NGAindividuals.BinChrom(NT.NIF(sizepop-(i-1)),j)=OGAindividuals.BinChrom(OT.OIF(i),j);
   end
end
fprintf('\n======================================================');
fprintf('\n====== Situation After filling 20 percent in Roulette Wheel ===');
fprintf('\n======================================================');
 for i=1:ceil(0.8*sizepop)
    NGAindividuals.fitness(i)= GAindividuals.fitness(NT.NIF(i));
   for j=1:lenchrom
    NGAindividuals.DecChrom(i,j)=GAindividuals.DecChrom(NT.NIF(i),j);
    end
    for j=1:lenchrom*Tgene
    NGAindividuals.BinChrom(i,j)=GAindividuals.BinChrom(NT.NIF(i),j);
    end
 end
fprintf('\n============================================================');
fprintf('\n==== Transfer from NGAindividuals table to GAindividuals ===');
fprintf('\n============================================================');
 for i=1:sizepop
    GAindividuals.fitness(i)= NGAindividuals.fitness(i);
   for j=1:lenchrom
    GAindividuals.DecChrom(i,j)=NGAindividuals.DecChrom(i,j);
    end
    for j=1:lenchrom*Tgene
    GAindividuals.BinChrom(i,j)=NGAindividuals.BinChrom(i,j);
    end
 end
fprintf('\n===================================================================================');
fprintf('\n====Table GAindividuals After transfer of the best 20 per cent of old population===');
fprintf('\n===================================================================================');
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
 fprintf('\n======================================================================================');

%===============================================================
%Stochastic Universal Selection algorithm is used to copy the remaining 70%
Sumfitness=0;
for i=1:sizepop
     Sumfitness=Sumfitness+newfitnessT(i);
end
AverageP=Sumfitness/sizepop;
Alpha=rand;
Sumfitness=newfitnessT(1);
Delta=Alpha*AverageP;
fprintf('Average=%1.3f Delta=%1.3f  Sumfitness=%1.3f',AverageP, Delta,Sumfitness);

% SU selection of 80% of individuls to join Nindividuals
fprintf('\n ============================================================');
fprintf('\n ceil(0.2*sizepop)=(%d) ceil(0.8*sizepop)=(%d)',ceil(0.2*sizepop),ceil(0.8*sizepop));
fprintf('\n ============================================================');
    cptSelect=0;
    k=1;
    while 1 
            if(Delta<Sumfitness)
                Delta=Delta+Sumfitness;
                fprintf('\n ======================================================================================');
                fprintf('\n Parent ceil(0.2*sizepop)=%d k=%d \n',ceil(0.2*sizepop),k);
                fprintf('\n of fitness=%1.3f is selected\n',GAindividuals.fitness(k));
                fprintf('\n ======================================================================================');
                NGAindividuals.fitness(k)=GAindividuals.fitness(k);
                for j=1:lenchrom 
                NGAindividuals.DecChrom(k,j)=GAindividuals.DecChrom(k,j);
                end
                for j=1:lenchrom*Tgene
                NGAindividuals.BinChrom(k,j)=GAindividuals.BinChrom(k,j);
                end
                cptSelect=cptSelect+1;
            else
                if k~=0
                Sumfitness=Sumfitness+newfitnessT(k);
                end
                k=k+1;
                fprintf('\nk=%d Sumfitness=%1.3f\n',k,Sumfitness);
            end
            if ((k>(ceil(0.8*sizepop))) && (cptSelect<ceil(0.8*sizepop))) 
                Sumfitness=newfitnessT(1);
                Delta=Alpha*AverageP;
                fprintf('\n ============================================================');
                fprintf('CptSelect=%d',cptSelect);
                fprintf('\n ============================================================\n');
                k=1;
            end
        if (k>(ceil(0.8*sizepop)))
            break;
        end
        fprintf('\n ============================================================\n');
        fprintf('\n ========End of iteration of probabilistic selection  ========\n');
        fprintf('\n ============================================================\n');
    end
fprintf('\n============================================================');
fprintf('\n==== Transfer from NGAindividuals table to GAindividuals ===');
fprintf('\n============================================================');
 for i=1:sizepop
    GAindividuals.fitness(i)= NGAindividuals.fitness(i);
   for j=1:lenchrom
    GAindividuals.DecChrom(i,j)=NGAindividuals.DecChrom(i,j);
    end
    for j=1:lenchrom*Tgene
    GAindividuals.BinChrom(i,j)=NGAindividuals.BinChrom(i,j);
    end
 end
fprintf('\n=================================================================');
fprintf('\n====Table GAindividuals at the end of Roulette Wheel Selection===');
fprintf('\n=================================================================');
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
 fprintf('\n======================================================================================');
end %end function Roulette Wheel Selection