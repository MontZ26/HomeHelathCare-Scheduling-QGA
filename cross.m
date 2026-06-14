function cross()
global lenchrom
global GAindividuals
global sizepop
global pcross
global Tgene
% The crossover
%pcross the probability of crossover is the probability that crossover will occur at a particular mating; 
%that is, not all matings must reproduce by crossover, but one could choose Pc=1.0.
% Single Point Crossover
  fprintf('\n======================================================');
fprintf('\n============Starting  CROSSOVER====================');
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
%{
 crosscount=0;
  for i=1:sizepop
      GAindividuals.PS(i)=0;
  end
%}
for i=1:sizepop
     % select child at random
    pick=rand;
    while pick==0
        pick=rand;
    end
    % probability of crossover
    pick=rand;
    if pick>pcross
          continue;   % go to next "for" loop 
    end
  % select two individuals at random
       pick=rand(1,2);
       index=ceil(pick.*(sizepop));
       fprintf('I AM IN CROSSOVER\n');
       fprintf('index(1)=%d index(2)=%d \n',index(1),index(2));
%{
pass=1;
        while (pass==1) || (prod(pick)==0) || (index(1)==index(2)) 
           pick=rand(1,2);
          index=ceil(pick.*(sizepop));
          if (GAindividuals.PS(index(1))==0 && GAindividuals.PS(index(2))==0)
              pass=0;
          end
        end
%}
% random position of crossover
%if (pass==0)&& (prod(pick)~=0)&&(index(1)~=index(2))
        pick=rand;
        while pick==0
            pick=rand;
        end
         pos=ceil(pick*(lenchrom*Tgene));
       fprintf('=====================================================\n');
       fprintf('index(1)=%d index(2)=%d POS_CROSSOVER=%d\n',index(1),index(2), pos);

% exchang tail
%        GAindividuals.PS(index(1))=1;
%        GAindividuals.PS(index(2))=1;
%        crosscount=crosscount+1;
        for j=pos:(lenchrom*Tgene)
             tail1=GAindividuals.BinChrom(index(1),j);
             GAindividuals.BinChrom(index(1),j)= GAindividuals.BinChrom(index(2),j);
             GAindividuals.BinChrom(index(2),j)=tail1;
        end
%        if (crosscount==round(sizepop/2))
%            break
%        end
end
%end
fprintf('\n============End of CROSSOVER===================='); 
end


