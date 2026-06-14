function Roulette_Wheel_Selection()
global sizepop
global lenchrom
global individuals
global Nindividuals
global oldfitnessT
global newfitnessT
global Tgene
SofT=struct('oi',zeros(1,sizepop),'oldfitness',zeros(1,sizepop));
SnfT=struct('ni',zeros(1,sizepop),'newfitness',zeros(1,sizepop));
Exchange=struct('oi',0,'oldfitness',0);

 fprintf('===========NEWFITNESS==============')

for k=1:2:2*sizepop-1
    fprintf('limitsize=%d\n',2*sizepop-1);
    % calculate fitness 
    % ffitness(x) is the function for calculating fitness 
    Nindividuals.fitness(k)=ffitness(k);        % ffitness(x) is the function for calculating fitness
    newfitnessT.ni(ceil(k/2))=k;
    newfitnessT.newfitness(ceil(k/2))=Nindividuals.fitness(k);
    Nindividuals.fitness(k+1)=0;
    oldfitnessT.oi(ceil(k/2))=k;
    oldfitnessT.oldfitness(ceil(k/2))=individuals.fitness(k);
   
    fprintf('k=%d ffitness(%d)=%1.3f newfitnessT.ni(ceil(k/2)=%d)=%d',k,k,ffitness(k),ceil(k/2),newfitnessT.ni(ceil(k/2)));
        fprintf('k=%d newfitness.ni(k)=%d newfitness.newfitness(k)=%1.3f\n',k,newfitnessT.ni(ceil(k/2)),newfitnessT.newfitness(ceil(k/2)));
        fprintf('k=%d oldfitness.oi(k)=%d oldfitness.oldfitness(k)=%1.3f\n',k,oldfitnessT.oi(ceil(k/2)),oldfitnessT.oldfitness(ceil(k/2)));

end
 fprintf('===========================================')

% suppose 'oldfitnessT' is the struct array. 'oldfitnessT' is the field that contains date and time.
T = struct2table(oldfitnessT); % convert the struct array to a table
sortedT = sortrows(T, 'oldfitness'); % sort the table by 'oldfitness'
SofT = table2struct(sortedT); % change it back to struct array if necessary

% suppose 'newfitnessT' is the struct array. 'newfitnessT' is the field that contains date and time.
T = struct2table(newfitnessT); % convert the struct array to a table
sortedT = sortrows(T, 'newfitness'); % sort the table by 'newfitness'
SnfT = table2struct(sortedT); % change it back to struct array if necessary
fprintf('==========================================')
for k=1:ceil(sizepop/2)-1
    Exchange.oi=SofT.oi(k);
    Exchange.oldfitness=SofT.oldfitness(k);
    SofT.oi(k)=SofT.oi(sizepop-k+1);
    SofT.oldfitness(k)=SofT.oldfitness(sizepop-k+1);
    SofT.oi(sizepop-k+1)=Exchange.oi;
    SofT.oldfitness(sizepop-k+1)=Exchange.oldfitness;
    fprintf('k=%d SofT.oi(k)=%d SofT.oldfitness(k)=%1.3\n',k,SofT.oi(k),SofT.oldfitness(k));
end
    fprintf('==========================================')
    %{
for k=1:ceil(0.1*sizepop)
    SofTindex=SofT.oi(k);
    SnfTindex=SnfT.ni(k);
    fprintf('k=%d Softindex=%d Snftindex=%d\n',k,SofTindex,SnfTindex);
    Nindividuals.fitness(SnfTindex)=individuals.fitness(SofTindex);
    for j=1:lenchrom
            SofTindex=SofT.oi(k);
            SnfTindex=SnfT.ni(k);
            Nindividuals.DecChrom(SnfTindex,j)=individuals.DecChrom(SofTindex,j);
     end
     for j=1:lenchrom*2
            SofTindex=SofT.oi(k);
            SnfTindex=SnfT.ni(k);
            Nindividuals.chrom(SnfTindex,j)=individuals.chrom(SofTindex,j);
            Nindividuals.BinChrom(SnfTindex,j)=individuals.BinChrom(SofTindex,j);
     end
end
   %} 
    s=1;
    for k=1:ceil(0.1*sizepop):ceil(sizepop/2)
    SofTindex=SofT.oi(k);
    SnfTindex=SnfT.ni(s);
   % fprintf('k=%d Softindex=%d Snftindex=%d\n',k,SofTindex,SnfTindex);
    Nindividuals.fitness(SnfTindex)=individuals.fitness(SofTindex);
    for j=1:lenchrom
            SofTindex=SofT.oi(k);
            SnfTindex=SnfT.ni(s);
            Nindividuals.DecChrom(SnfTindex,j)=individuals.DecChrom(SofTindex,j);
     end
     for j=1:lenchrom*Tgene
            SofTindex=SofT.oi(k);
            SnfTindex=SnfT.ni(s);
            Nindividuals.chrom(SnfTindex,j)=individuals.chrom(SofTindex,j);
            Nindividuals.BinChrom(SnfTindex,j)=individuals.BinChrom(SofTindex,j);
     end
     s=s+1;
    end
end %end function ElitSelection