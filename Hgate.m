function Hgate()
% On the basis of preparation
%  H gate 
global individuals
global Qbest
global sizepop
global lenchrom
global Tgene

for i=1:sizepop
    for j=1:Tgene*lenchrom
        A=individuals.chrom(2*i-1,j);    % ?
        B=individuals.chrom(2*i,j);      % ?
        x=individuals.BinChrom(2*i-1,j);
        b=Qbest.Bbinary(j);
        deltae=0;                % The size of the rotation angle
        s=0;                     % s is the rotation angle, that is, the direction of rotation
        
        if (x==0 && b==1 && (individuals.fitness(2*i-1)>=Qbest.Bfitness))
            deltae=0.005*pi;
            if A*B>0   
                s=1;
            elseif A*B<0
                s=-1;
            elseif A==0
                s=0;
            elseif B==0
                s=sign(randn);
            end
        elseif (x==0 && b==1 && (individuals.fitness(2*i-1)< Qbest.Bfitness))
             deltae=0.1 *pi;
            if A*B>0   
                s=-1;
            elseif A*B<0
                s=1;
            elseif A==0
                 s=sign(randn);
            elseif B==0
               s=0;
            end
        elseif (x==1 && b==0 && (individuals.fitness(2*i-1)>= Qbest.Bfitness))
            deltae=0.005*pi;
            if A*B>0   
                s=-1;
            elseif A*B<0
                s=1;
            elseif A==0
                s=sign(randn);
            elseif B==0
                s=0;
            end
        elseif (x==1 && b==0 &&(individuals.fitness(2*i-1)< Qbest.Bfitness))
            deltae=0.1*pi;
            if A*B>0   
                s=1;
            elseif A*B<0
                s=-1;
            elseif A==0
                 s=0;
            elseif B==0
               s=sign(randn);
            end
        elseif (x==1 && b==1 && (individuals.fitness(2*i-1)>= Qbest.Bfitness))
            deltae=0.005*pi;
            if A*B>0   
                s=1;
            elseif A*B<0
                s=-1;
            elseif A==0
                s=0;
            elseif B==0
                s=sign(randn);
            end
        elseif (x==1 && b==1 && (individuals.fitness(2*i-1)< Qbest.Bfitness))
            deltae=0.1*pi;
            if A*B>0   
                s=1;
            elseif A*B<0
                s=-1;
            elseif A==0
                s=0;
            elseif B==0
                s=sign(randn);
            end
            elseif (x==0 && b==0 && (individuals.fitness(2*i-1)>= Qbest.Bfitness))
            deltae=0.005*pi;
            if A*B>0   
                s=1;
            elseif A*B<0
                s=-1;
            elseif A==0
                s=0;
            elseif B==0
                s=sign(randn);
            end
        elseif (x==0 && b==0 && (individuals.fitness(2*i-1)< Qbest.Bfitness))
            deltae=0.1*pi;
            if A*B>0   
                s=1;
            elseif A*B<0
                s=-1;
            elseif A==0
                s=0;
            elseif B==0
                s=sign(randn);
            end
            
        end
        e=s*deltae;             % e as the rotation angle
        u=[cos(e) -sin(e);sin(e) cos(e)];         % Quantum revolving door
        y=u*[A B]';             % y an updated quantum bit

    E=0.005;
        
        if y(1)^2<=E
            y=[sqrt(E) sqrt(1-E)]';
        elseif y(1)^2>=1-E
                y=[sqrt(1-E) sqrt(E)]';
        end
        individuals.chrom(2*i-1,j)=y(1);
        individuals.chrom(2*i,j)=y(2);
    end
end
