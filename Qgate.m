function Qgate()
global FID
% On the basis of preparation
%  quantum gate ( another )
global lenchrom
global individuals
global best
global sizepop
global Tgene

%global maxgen gen
fprintf(FID,'best binary   best binary  best binary best binary ');
for i=1:lenchrom*Tgene
    fprintf(FID,'%d ',best.Bbinary(i));
end
for i=1:sizepop
    for j=1:Tgene*lenchrom
        Alpha=individuals.Alphachrom(i,j);    % ?
        Beta=individuals.Betachrom(i,j);      % ?
        x=individuals.BinChrom(i,j);
        b=best.Bbinary(j);
        deltae=0;                % The size of the rotation angle
        s=0;                     % s The symbol for the rotation angle, that is, the direction of rotation
        if (x==0 && b==1 &&(individuals.fitness(i)< best.Bfitness))
            deltae=0.5*pi; %0.025
            if Alpha*Beta>0   
                s=-1;
            elseif Alpha*Beta<0
                s=+1;
            elseif Alpha==0
                s=sign(randn);
            elseif Beta==0
                s=0;
            end
        elseif (x==1 && b==0 && (individuals.fitness(i)< best.Bfitness))
            deltae=0.5*pi; %0.025
            if Alpha*Beta>0   
                s=1;
            elseif Alpha*Beta<0
                s=-1;
            elseif Alpha==0
                s=0;
            elseif Beta==0
                s=sign(randn);
            end
        elseif (x==1 && b==0 &&(individuals.fitness(i)>= best.Bfitness))
            deltae=0.01*pi; %0.005 
            if Alpha*Beta>0   
                s=-1;
            elseif Alpha*Beta<0
                s=+1;
            elseif Alpha==0
                s=sign(randn);
            elseif Beta==0
                s=0;
            end
        elseif (x==1 && b==1 && (individuals.fitness(i)< best.Bfitness))
            deltae=0.5*pi; %0.025
            if Alpha*Beta>0   
                s=1;
            elseif Alpha*Beta<0
                s=-1;
            elseif Alpha==0
                s=0;
            elseif Beta==0
                s=sign(randn);
            end
        elseif (x==1 && b==1 && (individuals.fitness(i)>= best.Bfitness))
            deltae=0.01*pi; %0.005
            if Alpha*Beta>0   
                s=1;
            elseif Alpha*Beta<0
                s=-1;
            elseif Alpha==0
                s=0;
            elseif Beta==0
                s=sign(randn);
            end
            
        end
        
        e=s*deltae;
        u=[cos(e) -sin(e);sin(e) cos(e)];
        y=u*[Alpha Beta]';
        individuals.Alphachrom(i,j)=y(1);
        individuals.Betachrom(i,j)=y(2);
    end
end
