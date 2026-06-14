function mutation()
global sizepop
global pmutation
global lenchrom
global GAindividuals
global Tgene
fprintf('\n I AM IN MUTATION   I AM IN MUTATION   I AM   IN MUTATION');
% quantum mutation  , Quantum Not Gate
for i=1:sizepop
    
    % select child at random
    pick=rand;
    while pick==0
        pick=rand;
    end
    % probability of mutation
    pick=rand;
    if pick>pmutation
      continue; % go to next "for" loop 
    end
    %===================================================================
    %============== Mutation possible ==================================
    %===================================================================
    % mutation  index (position of chromosome in population) 
    pick=rand;
    while pick==0
         pick=rand;
    end
    index=ceil(pick*(sizepop));
    % exchange o in 1 and 1 in 0 ? and ?
    pick=rand;
    while pick==0
         pick=rand;
    end
    pos=ceil(pick*lenchrom*Tgene);
    fprintf('%d %1.3f ',index,GAindividuals.BinChrom(index,pos));
    GAindividuals.BinChrom(index,pos)=1-GAindividuals.BinChrom(index,pos);
    fprintf('\n=================End of Mutation===================');
end
end