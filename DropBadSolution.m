function DropBadSolution()
global sizepop
global GAindividuals
global lenchrom
for i=1:sizepop
   Double=0;
   l=1;
   while(l<=lenchrom && Double==0)
       li=l;
       prec=GAindividuals.DecChrom(i,l);
       while(li<=lenchrom && Double==0)
           if GAindividuals.DecChrom(i,l)==prec
               Double=1;
           else
               li=li+1;
           end
           if (Double==0)
               l=l+1;
           end
       end
       if (Double==1)
          GAindividuals.PS(i)
       end
   end
end