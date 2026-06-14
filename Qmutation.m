function Qmutation()

global sizepop
global pmutation
global lenchrom
global individuals
global Tgene

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
        fprintf('=====================================================\n');

    fprintf('I AM IN Quantum MUTATION\n');
       fprintf('index=%d Alpha=%1.2f Beta=%1.2f Alpha^2+Beta^2=%1.2f\n',index,individuals.Alphachrom(index,pos),individuals.Betachrom(index,pos), individuals.Alphachrom(index,pos).^2+individuals.Betachrom(index,pos).^2);

    % v=individuals.chrom(2*index-1,pos);
    X=individuals.Alphachrom(index,pos);
    individuals.Alphachrom(index,pos)= individuals.Betachrom(index,pos);
    individuals.Betachrom(index,pos)=X;
    fprintf('=====================================================\n');
         fprintf('I AM IN Quantum MUTATION\n');
       fprintf('index=%d Alpha=%1.2f Beta=%1.2f Alpha^2+Beta^2=%1.2f\n',index,individuals.Alphachrom(index,pos),individuals.Betachrom(index,pos), individuals.Alphachrom(index,pos).^2+individuals.Betachrom(index,pos).^2);

    % individuals.chrom(2*index,pos)=v;
end