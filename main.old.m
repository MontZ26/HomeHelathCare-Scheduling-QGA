% This is quantum genetic algorithm(QGA)
% In this function ,it fulfils quantum genetic algorithm
clc;
clear;
global FID
global sizepop 
global lenchrom
global judge
global maxgen
global individuals
global Nindividuals
global GAindividuals
global NGAindividuals
global OGAindividuals
global Oindividuals
global PDist
global CTN
global pmutation
global pcross
global PReq
global CTask
global CTaskN
global gen
global best
global oldfitnessT
global newfitnessT


CTask=[1,0,1,0,1;0,1,1,1,0;1,1,0,1,1;1,1,1,0,0];
PReq=[1,0,1,0,1,1,0,1,0;1,1,0,1,0,1,1,0,1;0,1,1,0,0,1,1,0,1;0,1,1,0,0,1,0,1,1;1,0,1,0,1,0,1,1,0;0,1,1,0,1,1,0,1,1];
CTaskN=5;
%--------10-------These can be modified as you like-----------------------
maxgen=200;                         % maximum generation
sizepop=27;                        % size of population 
CTN=4; % Care Team Number

lenchrom=9; % length of bit of every chromosom varible
%BinChrom=[1,lenchrom*2]; % Bianry chromosom of patients and care givers team visiting them
%DecChrom=[1 lenchrom]; % matrix of patients and care givers team visiting them
%PDist matrix of distances between patients DistCGT sum distance of team
%GlobalDist global distance sum= DistCGT[1]+DistCGT[2]+DistCGT[3] 
%ChromDist=0;
PDist=[0 7.28 7.81 4.47 5 5 7.07 6.71 4 1.41;7.28 0 2.2 3 4 3 3.6 5.7 7.1 6.1;7.81 2.24 0 4 3 4 1.4 3.6 6 5.8;4.47 3 4 0 2 3 4.1 6 4.5 2.2;5 4.47 3 2.2 0 4 2.2 2 4 3.6;5 2.83 4.5 3.2 4 0 5 6.3 6.4 4.1;7.07 3.61 1.4 4.1 2 5 0 2.2 5.1 5.4;6.71 5.66 3.6 6 2 6 2.2 0 3.6 5.4;4 7.14 6 4.5 4 6 5.1 3.6 0 3.2;1.41 6.08 5.8 2.2 4 4 5.4 5.4 3.2 0];
pcross=0.8;
pmutation=0.008;

QGAbestindex=0;
QGAbestfitness=0;
QGAbestfitnessT=zeros(sizepop,1);
QGAGlobalbestfitness=0;

judge=0;                 % disaster counter

individuals=struct('fitness',zeros(1,sizepop*2),'DecChrom',zeros(sizepop*2,lenchrom),'BinChrom',zeros(sizepop*2,lenchrom*2),'chrom',zeros(sizepop*2,lenchrom*2));   % structure of population
Oindividuals=struct('fitness',zeros(1,sizepop*2),'DecChrom',zeros(sizepop*2,lenchrom),'BinChrom',zeros(sizepop*2,lenchrom*2),'chrom',zeros(sizepop*2,lenchrom*2));   % structure of population
Nindividuals=struct('fitness',zeros(1,sizepop*2),'DecChrom',zeros(sizepop*2,lenchrom),'BinChrom',zeros(sizepop*2,lenchrom*2),'chrom',zeros(sizepop*2,lenchrom*2));   % New structure of population
GAindividuals=struct('fitness',zeros(1,sizepop),'DecChrom',zeros(sizepop,lenchrom),'BinChrom',zeros(sizepop,lenchrom*2));   % structure of GA population
NGAindividuals=struct('fitness',zeros(1,sizepop),'DecChrom',zeros(sizepop,lenchrom),'BinChrom',zeros(sizepop,lenchrom*2));   % structure of GA population
OGAindividuals=struct('fitness',zeros(1,sizepop),'DecChrom',zeros(sizepop,lenchrom),'BinChrom',zeros(sizepop,lenchrom*2));   % structure of GA population

best=struct('Bfitness',0,'BDecChrom',zeros(1,lenchrom),'Bbinary',zeros(1,lenchrom*2),'Bchrom',zeros(2,lenchrom*2));   % best individual 
oldfitnessT=zeros(1,sizepop);
newfitnessT=zeros(1,sizepop);
FileName='results';
FID = fopen(FileName, 'w');
if FID < 0, error('Cannot open file'); end
%++++++++++++++++++++++++++++++++Quantum Genetic Algoritm++++++++++++++++++++++++++++++++
%Initialization
fprintf('**************************************************************\n');
fprintf('****************Quantum Genetic Algorithm*********************\n');
fprintf('**************************************************************\n');

for i=1:sizepop*2
    for j=1:lenchrom
     individuals.DecChrom(i,j)=0;
    end
    for j=1:lenchrom*2
     individuals.chrom(i,j)= 1/sqrt(2);
     individuals.BinChrom(i,j)=0;
   end
    individuals.fitness(i)=0;
end
QGAGlobalbestfitness=0;
 Qcollapse();
%Transfer from new population to old population

  for k=1:sizepop*2
        Oindividuals.fitness(k)=individuals.fitness(k);
        for j=1:lenchrom
        Oindividuals.DecChrom(k,j)=individuals.DecChrom(k,j);
        end
        for j=1:lenchrom*2
        Oindividuals.chrom(k,j)=individuals.chrom(k,j);
        Oindividuals.BinChrom(k,j)=individuals.BinChrom(k,j);
        end
   end
%================================================================
% evolution starts of Quantum Genetic Algorithm
%================================================================    
for gen=1:maxgen
fprintf('****Starting Quantum Collapse at the Beginning of iteration');

 %{
  % avoid disaster in quantum GA
  % disaster
    if judge>0.2*maxgen
       Qdisaster();
   end
 %}
    fprintf('===========NEW GENERATION==============')
 for i=1:sizepop*2
     fprintf('\n individuals(%d)   ',i);   
     fprintf('%1.3f',individuals.fitness(i));
  for j=1:lenchrom
   fprintf(' %d',individuals.DecChrom(i,j));
  end
  fprintf('###');
   for j=1:lenchrom*2
    fprintf(' %d',individuals.BinChrom(i,j));
   end
   fprintf('###');
  for j=1:lenchrom*2
    fprintf(' %1.2f',individuals.chrom(i,j));
  end
  fprintf('\n');
 end
fprintf('===========NEWFITNESS==============')
% calculate fitness 
% ffitness(x) is the function for calculating fitnes
    for k=1:2:2*sizepop-1
        fprintf('limitsize=%d\n',2*sizepop-1);
        % calculate fitness 
        % ffitness(x) is the function for calculating fitness 
        individuals.fitness(k)=ffitness(k);        % ffitness(x) is the function for calculating fitness
        fprintf('===========NEWFITNESS-individuals==============');
        fprintf('\n individuals(%d)   ',k);   
        fprintf('%1.3f',individuals.fitness(k));
        if gen==1
          Oindividuals.fitness(k)= individuals.fitness(k);
        end
        %newfitnessT.ni(ceil(k/2))=k;
        newfitnessT(ceil(k/2))=individuals.fitness(k);
        individuals.fitness(k+1)=0;
        %oldfitnessT.oi(ceil(k/2))=k;
        oldfitnessT(ceil(k/2))=Oindividuals.fitness(k);
        fprintf('k=%d  newfitness(k)=%1.3f\n',k,newfitnessT(ceil(k/2)));
        fprintf('k=%d oldfitness=%1.3f\n',k,oldfitnessT(ceil(k/2)));

    end
    %======================= record the best individual to=================
    [QGAbestfitness,QGAbestindex]=max(individuals.fitness);     % find maximum value which is the best
    %======================= record the best individual to=================
    best.Bfitness=QGAbestfitness;
    for j=1:lenchrom
    best.BDecChrom(j)=individuals.DecChrom(QGAbestindex,j);
    end
    for j=1:lenchrom*2
    best.Bchrom(1,j)=individuals.chrom(QGAbestindex,j);
    best.Bchrom(2,j)=individuals.chrom(QGAbestindex+1,j);
    best.Bbinary(1,j)=individuals.BinChrom(QGAbestindex,j);
    best.Bbinary(2,j)=individuals.BinChrom(QGAbestindex+1,j);
    end
    fprintf('===========================================')
    %Transfer from new population to old population
    for k=1:sizepop*2
        Oindividuals.fitness(k)=individuals.fitness(k);
        for j=1:lenchrom
            Oindividuals.DecChrom(k,j)=individuals.DecChrom(k,j);
        end
        for j=1:lenchrom*2
            Oindividuals.chrom(k,j)=individuals.chrom(k,j);
            Oindividuals.BinChrom(k,j)=individuals.BinChrom(k,j);
        end
    end

    %fprintf('****Starting Quantum Crossover******\n');
    % quantum cross
    Qcross();
    %fprintf('****Starting Quantum mutation******\n');
    % quantum mutation 
    Qmutation();
    fprintf('****Starting Quantum Gate******\n');
    % quantum gate
    Hgate();

    % record the QGA best individual to "best"
    fprintf(FID,'**********************************************************\n');
    QGAbestfitnessT(gen)=QGAbestfitness;
    if QGAbestfitness>QGAGlobalbestfitness
          QGAGlobalbestfitness=QGAbestfitness;
    end
fprintf(FID,'********************************************************************************************************************\n');
fprintf('********************************************************************************************************************\n');
fprintf('Generation %d  QGAbestIndex=%d QGABestfitnessT(%d)=%2.3f QGAbestfitness=%2.3f    QGAGlobalbestfitness=%2.3f\n',gen,gen, QGAbestindex,QGAbestfitnessT(gen),QGAbestfitness,QGAGlobalbestfitness);
fprintf(FID,'Generation %d  QGAbestIndex=%d QGABestfitnessT(%d)=%2.3f QGAbestfitness=%2.23f    QGAGlobalbestfitness=%2.3f\n',gen,gen, QGAbestindex,QGAbestfitnessT(gen),QGAbestfitness,QGAGlobalbestfitness);
fprintf('********************************************************************************************************************\n');
fprintf(FID,'*******************************************************************************************************************\n');
 fprintf('****Starting Quantum Collapse at the Beginning of iteration');
 % The quantum collapse.
 Qcollapse();
 end
fprintf('**************************************************************\n');
fprintf('*********************Genetic Algorithm************************\n');
fprintf('**************************************************************\n');
fprintf('**************************************************************\n');
fprintf('*********************Genetic Algorithm************************\n');
fprintf('**************************************************************\n');
fprintf('**************************************************************\n');
fprintf('*********************Genetic Algorithm************************\n');
fprintf('**************************************************************\n');
fprintf('**************************************************************\n');
fprintf('*********************Genetic Algorithm************************\n');
fprintf('**************************************************************\n');
fprintf('======================================================================');
fprintf('===========evolution starts of Simple Elit Genetic Algorithm==========');
fprintf('======================================================================');

%Initialization
pcross=0.8;
pmutation=0.2;

bestindex=0;
bestfitness=0;
bestfitnessT=zeros(sizepop,1);
Globalbestfitness=0;
for i=1:sizepop*2
    for j=1:lenchrom
     individuals.DecChrom(i,j)=0;
    end
    for j=1:lenchrom*2
        individuals.chrom(i,j)= rand;
     %individuals.chrom(i,j)= 1/sqrt(2);
     individuals.BinChrom(i,j)=0;
   end
    individuals.fitness(i)=0;
end
%{

  for k=1:sizepop*2
        individuals.fitness(k)=individualsT.fitness(k);
        Oindividuals.fitness(k)=individualsT.fitness(k);
        for j=1:lenchrom
        individuals.DecChrom(k,j)=individualsT.DecChrom(k,j);
        Oindividuals.DecChrom(k,j)=individualsT.DecChrom(k,j);
        end
        for j=1:lenchrom*2
        individuals.chrom(k,j)=individualsT.chrom(k,j);
        Oindividuals.chrom(k,j)=individualsT.chrom(k,j);
        Oindividuals.BinChrom(k,j)=individualsT.BinChrom(k,j);
        individuals.BinChrom(k,j)=individualsT.BinChrom(k,j);
        end
   end
%}
 for i=1:sizepop*2
     fprintf('\n individuals(%d)   ',i);   
     fprintf('%1.3f',individuals.fitness(i));
  for j=1:lenchrom
   fprintf(' %d',individuals.DecChrom(i,j));
  end
  fprintf('###');
   for j=1:lenchrom*2
    fprintf(' %d',individuals.BinChrom(i,j));
   end
   fprintf('###');
   for j=1:lenchrom*2
    fprintf(' %1.3f',individuals.chrom(i,j));
   end
   fprintf('\n');
 end
 fprintf('\n======================================================');
 fprintf('\n===========COLLAPSE2==================================');
 fprintf('\n======================================================');
 Qcollapse2();
%Transfer from new population to old population

  for k=1:sizepop
        for j=1:lenchrom
        OGAindividuals.DecChrom(k,j)=GAindividuals.DecChrom(k,j);
        GAindividuals.DecChrom(k,j)=GAindividuals.DecChrom(k,j);
        end
        for j=1:lenchrom*2
        OGAindividuals.BinChrom(k,j)=GAindividuals.BinChrom(k,j);
        GAindividuals.BinChrom(k,j)=GAindividuals.BinChrom(k,j);
        end
  end
for gen=1:maxgen
     fprintf('===========NEW GENERATION==============')
 for i=1:sizepop
     fprintf('\n GAindividuals(%d)   ',i);   
     fprintf('%1.3f',GAindividuals.fitness(i));
  for j=1:lenchrom
   fprintf(' %d',GAindividuals.DecChrom(i,j));
  end
  fprintf('###');
   for j=1:lenchrom*2
    fprintf(' %d',GAindividuals.BinChrom(i,j));
   end
   fprintf('###');
   fprintf('\n');
 end
fprintf('===========NEWFITNESS==============');
% calculate GA fitness 
% GAfitness(x) is the function for calculating fitnes
    for k=1:sizepop
        % calculate fitness 
        % GAfitness(x) is the function for calculating fitness 
        OGAindividuals.fitness(k)= GAindividuals.fitness(k);
        GAindividuals.fitness(k)=GAfitness(k);        % ffitness(x) is the function for calculating fitness
        fprintf('===========NEWFITNESS-individuals==============');
        fprintf('\n individuals(%d)   ',k);   
        fprintf('%1.3f',GAindividuals.fitness(k));
        newfitnessT(k)=GAindividuals.fitness(k);
        oldfitnessT(k)=OGAindividuals.fitness(k);
        fprintf('k=%d  newfitness(k)=%1.3f\n',k,newfitnessT(k));
        fprintf('k=%d oldfitness=%1.3f\n',k,oldfitnessT(k));

    end
    %======================= record the best individual to=================
    [bestfitness,bestindex]=max(GAindividuals.fitness);     % find maximum value which is the best
    %======================= record the best individual to=================
    best.Bfitness=bestfitness;
    for j=1:lenchrom
    best.BDecChrom(j)=GAindividuals.DecChrom(bestindex,j);
    end
    for j=1:lenchrom*2
    best.Bbinary(j)=GAindividuals.BinChrom(bestindex,j);
    end
    fprintf('===========================================')
    %Transfer from new population to old population
   if gen==1 
   
    for k=1:sizepop
        OGAindividuals.fitness(k)=GAindividuals.fitness(k);
        for j=1:lenchrom
            OGAindividuals.DecChrom(k,j)=GAindividuals.DecChrom(k,j);
        end
        for j=1:lenchrom*2
            OGAindividuals.BinChrom(k,j)=GAindividuals.BinChrom(k,j);
        end
    end
   end

     fprintf('\n======================================================');
     fprintf('\n======================================================');
     % Selection using Stochastic Universal Selection
    fprintf('****Starting  Roulette Wheel Selection******\n');
   %ElitSelection();
   %{-----------------------%}
   Roulette_Wheel_Selection();
   %Selection();
   
    %Transfer from new population to old population

  for k=1:sizepop
        GAindividuals.fitness(k)=NGAindividuals.fitness(k);
        for j=1:lenchrom
            GAindividuals.DecChrom(k,j)=NGAindividuals.DecChrom(k,j);
        end
        for j=1:lenchrom*2
            GAindividuals.BinChrom(k,j)=NGAindividuals.BinChrom(k,j);
        end
   end
   
   
    fprintf('\n====================================================');

    fprintf('****Starting GA Crossover******\n');
    % GA crossover
  cross();
     fprintf('===========NEW GENERATION==============')
 for i=1:sizepop
      fprintf('===========GAindividuals==============')
     fprintf('\n individuals(%d)   ',i);   
     fprintf('%1.3f',GAindividuals.fitness(i));
  for j=1:lenchrom
   fprintf(' %d',GAindividuals.DecChrom(i,j));
  end
  fprintf('###');
   for j=1:lenchrom*2
    fprintf(' %d',GAindividuals.BinChrom(i,j));
   end
   fprintf('###');
   fprintf('===========OGAindividuals==============')
   fprintf('\n OGAindividuals(%d)   ',i);   
   fprintf('%1.3f',OGAindividuals.fitness(i));
   for j=1:lenchrom
   fprintf(' %d',OGAindividuals.DecChrom(i,j));
  end
     fprintf('###');
   for j=1:lenchrom*2
    fprintf(' %d',OGAindividuals.BinChrom(i,j));
   end
   fprintf('###');
   fprintf('\n');
 end
    fprintf('****Starting GA mutation******\n');
    % GA mutation 
   mutation();
    % record the QGA best individual to "best"
    fprintf(FID,'**********************************************************\n');
    bestfitnessT(gen)=bestfitness;
    if bestfitness>Globalbestfitness
          Globalbestfitness=bestfitness;
    end
fprintf(FID,'********************************************************************************************************************\n');
fprintf('********************************************************************************************************************\n');
fprintf('Generation %d  bestIndex=%d BestfitnessT(%d)=%2.3f bestfitness=%2.3f    Globalbestfitness=%2.3f\n',gen,gen, bestindex,bestfitnessT(gen),bestfitness,Globalbestfitness);
fprintf(FID,'Generation %d  bestIndex=%d BestfitnessT(%d)=%2.3f bestfitness=%2.23f    Globalbestfitness=%2.3f\n',gen,gen, bestindex,bestfitnessT(gen),bestfitness,Globalbestfitness);
fprintf('********************************************************************************************************************\n');
fprintf(FID,'*******************************************************************************************************************\n');
end

fprintf(FID,'****************************************************************************************\n');
    fprintf('-------------------------------------------------------------\n');
    fprintf('----gen--bestfitnessT(gen))----------------------QGAbestfitnessT(gen)\n');
    fprintf('-------------------------------------------------------------\n');
for gen=1:maxgen
     fprintf('----|%d-| |%1.5f|                         |%1.5f|\n',gen,bestfitnessT(gen),QGAbestfitnessT(gen));
end

%DrawPaths(cpt);

hold off

X=1:maxgen;
%figure
figure('Name','Comparison of Quantum Elit Genetic Algorithm & Elit Genetic Algorithm')
if (QGAbestfitnessT ~= 0)
    Y1=QGAbestfitnessT;
end
if (bestfitnessT ~= 0)
    Y2=bestfitnessT;
end
plot(X,Y1,X,Y2,'--');
fclose(FID);

