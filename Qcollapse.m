function Qcollapse()
global individuals
global sizepop
global lenchrom
global Tgene
% The quantum collapse.

  for i=1:sizepop
     for j=1:lenchrom*Tgene
        pick=rand;
        if (pick<(individuals.Alphachrom(i,j))^2)
            individuals.BinChrom(i,j)=1;
        else
            individuals.BinChrom(i,j)=0;
        end
      end
  end
end
