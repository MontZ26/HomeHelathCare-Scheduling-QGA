function Roulette_Wheel_Selection()
global sizepop
global lenchrom
global GAindividuals
global OGAindividuals
global NGAindividuals
global Nindividuals
global newfitnessT
global oldfitnessT
global Tgene


fprintf('\n======================================================');
fprintf('\n============I am in Roulette Wheel====================');
fprintf('\n============Table OGAindividuals====================');
fprintf('\n======================================================');
 for i=1:sizepop
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
    fprintf('i=%d oldfitnessT(i)=%1.3f\n',i,oldfitnessT(i));
end
fprintf('==========================================\n');
fprintf('======SOFT NEW Fitness Table before sorting==============\n');
fprintf('==========================================\n');

for i=1:sizepop
    fprintf('k=%d  newfitnessT(k)=%1.3f \n',i,newfitnessT(i));
end
% suppose 'oldfitnessT' is the struct array. 
[SofT,OI] = sort(oldfitnessT,'descend'); % sort the table by 'oldfitness'in descending order
fprintf('==================================================\n');
for i=1:sizepop
    OT.OIF(i)=OI(i);
    OT.OLDF(i)=SofT(i); 
end
% suppose 'newfitnessT' is the struct array. 
[SnfT,NI] = sort(newfitnessT); % sort the table by 'newfitness'
fprintf('==================================================\n');
for i=1:sizepop
    NT.NIF(i)=NI(i);
    NT.NEWF(i)=SnfT(i); 
end
fprintf('=================END OF RWS==========\n');
for i=1:sizepop
     NGAindividuals.fitness(i)=GAindividuals.fitness(i);
     for j=1:lenchrom
            NGAindividuals.DecChrom(i,j)=GAindividuals.DecChrom(i,j);
     end
     for j=1:lenchrom*Tgene
            NGAindividuals.BinChrom(i,j)=GAindividuals.BinChrom(i,j);
     end
end
fprintf('==========================================');
   %copy of 20% of old population innew population
for i=1:ceil(0.2*sizepop)
   fprintf('k=%d OT.OLDF(%d)=%1.3f',i,OT.OIF(i),OT.OLDF(i));
   fprintf('NT.NEWF(%d)=%1.3f\n',NT.NIF(i),NT.NEWF(i));

   NGAindividuals.fitness(NT.NIF(i))=OGAindividuals.fitness(OT.OIF(i));
    for j=1:lenchrom
            NGAindividuals.DecChrom(NT.NIF(i),j)=OGAindividuals.DecChrom(OT.OIF(i),j);
    end
    for j=1:lenchrom*Tgene
            NGAindividuals.BinChrom(NT.NIF(i),j)=OGAindividuals.BinChrom(OT.OIF(i),j);
    end
end
fprintf('\n======================================================');
fprintf('\n====== Situation After filling 20 percent in Roulette Wheel ===');
fprintf('\n======================================================');
for i=1:sizepop
     fprintf('\n Nindividuals(%d)=',i);   
     fprintf('%1.3f',NGAindividuals.fitness(i));
     for j=1:lenchrom
        fprintf(' %d',NGAindividuals.DecChrom(i,j));
     end
     fprintf('###');
     for j=1:lenchrom*Tgene
        fprintf(' %d',NGAindividuals.BinChrom(i,j));
     end
     fprintf('###');
     fprintf('\n');
 end

%======================================================================
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

% SU selection of 780% of individuls to join Nindividuals
fprintf('\n ============================================================');
fprintf('\n ceil(0.2*sizepop)=(%d) ceil(0.8*sizepop)=(%d)',ceil(0.2*sizepop),ceil(0.8*sizepop));
fprintf('\n ============================================================');
cptSelect=0;
k=1;
while 1 
        if ceil(0.2*sizepop)+k >sizepop
                    break;
        end
        if(Delta<Sumfitness)
            Delta=Delta+Sumfitness;
            fprintf('\n ======================================================================================');
            fprintf('\n Parent ceil(0.2*sizepop)=%d k=%d \n',ceil(0.2*sizepop),k);
            fprintf('\n of fitness=%1.3f is selected\n',GAindividuals.fitness(k));
            fprintf('\n ======================================================================================');

            for j=1:lenchrom 
                NGAindividuals.DecChrom(ceil(0.2*sizepop)+k,j)=GAindividuals.DecChrom(k,j);
            end
            for j=1:lenchrom*Tgene
                NGAindividuals.BinChrom(ceil(0.2*sizepop)+k,j)=GAindividuals.BinChrom(k,j);
            end
            cptSelect=cptSelect+1;
       else
            if k~=0
                Sumfitness=Sumfitness+newfitnessT(k);
            end
            k=k+1;
            fprintf('\nk=%d Sumfitness=%1.3f\n',k,Sumfitness);
       end
       if ((k>(ceil(0.9*sizepop)-2)) && (cptSelect<ceil(0.8*sizepop)-1)) 
            Sumfitness=newfitnessT(1);
            Delta=Alpha*AverageP;
            fprintf('\n ============================================================');
            fprintf('CptSelect=%d',cptSelect);
            fprintf('\n ============================================================\n');
            k=1;
       end
       if (k>(ceil(0.9*sizepop)))
            break;
       end
end
for i=1:sizepop
    Nindividuals.fitness(i)=NGAindividuals.fitness(i);
    for j=1:lenchrom
        Nindividuals.DecChrom(i,j)=NGAindividuals.DecChrom(i,j);
    end
    for j=1:lenchrom*Tgene
        Nindividuals.BinChrom(i,j)=NGAindividuals.BinChrom(i,j);
    end
end
fprintf('==========end function Roulette Wheel Selection================================');
end