function Qcross()
global lenchrom
global individuals
global sizepop
global pcross
global Tgene

% The crossover
%pcross the probability of crossover is the probability that crossover will occur at a particular mating; 
%that is, not all matings must reproduce by crossover, but one could choose Pc=1.0.
% Single Point Crossover
%crosscount=0;
for i=1:sizepop
      individuals.PS(i)=0;
end
for i=1:sizepop
       % probability of crossover
       pick=rand;
       if pick>pcross
          continue;   % go to next "for" loop 
       end
  % select two individuals at random
       pick=rand(1,2);
       index=ceil(pick.*sizepop);
       fprintf('=====================================================\n');
         fprintf('I AM IN Quantum CROSSOVER\n');
       fprintf('index(1)=%d index(2)=%d \n',index(1),index(2));
    
        pick=rand;
        while pick==0
            pick=rand;
        end
        pos=ceil(pick.*(lenchrom*Tgene));
        fprintf('=====================================================\n');
        fprintf('I AM IN Quantum CROSSOVER\n');
        fprintf('index(1)=%d index(2)=%d  POS_CROSSOVER=%d\n',index(1),index(2),pos);
%{
       pass=1;
        while (pass==1) || (prod(pick)==0) || (index(1)==index(2))
           pick=rand(1,2);
          index=ceil(pick.*(sizepop));
          if (individuals.PS(index(1))==0 && individuals.PS(index(2))==0)
              pass=0;
          end
        end
if (pass==0)&& (prod(pick)~=0)&&(index(1)~=index(2))
% random position of crossover
        pick=rand;
        while pick==0
            pick=rand;
        end
         pos=ceil(pick.*(lenchrom*Tgene));
   fprintf('=====================================================\n');
         fprintf('I AM IN Quantum CROSSOVER\n');
       fprintf('index(1)=%d index(2)=%d  POS_CROSSOVER=%d\n',index(1),index(2),pos);
%}
% exchang tail 
        individuals.PS(index(1))=1;
        individuals.PS(index(2))=1;
%        crosscount=crosscount+1;
         for j=pos: (lenchrom*Tgene)
         %Crossover of Alpha
             tail1=individuals.Alphachrom(index(1),j);
             individuals.Alphachrom(index(1),j)= individuals.Alphachrom(index(2),j);
             individuals.Alphachrom(index(2),j)=tail1;
         %Crossover of Beta
             tail2=individuals.Betachrom(index(1),j);
             individuals.Betachrom(index(1),j)= individuals.Betachrom(index(2),j);
             individuals.Betachrom(index(2),j)=tail2;
         end
%        if (crosscount==round(sizepop/2))
%            break
%        end
%end
end

end