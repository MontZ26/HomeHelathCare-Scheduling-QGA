function Qdisaster()

%  quantum disaster
global best
global individuals
global sizepop
global judge 
global Tgene

judge=0;        % Catastrophes Counter to zero

% Keep only the best individuals, the rest of the individual to generate

individuals.chrom(1,:)=best.Bchrom(1,:);
individuals.chrom(2,:)=best.Bchrom(2,:);
%{
for i=3:sizepop*2
    for j=1:2*lenchrom
        individuals.chrom(i,j)=1/sqrt(2);
    end 
end
%}
for i=3:2:sizepop*2-1
    for j=1:lenchrom*Tgene
     Alpha=rand;
     individuals.chrom(i,j)=sqrt(Alpha);
     individuals.chrom(i+1,j)=sqrt(1-Alpha);
    end
end

