function kCAFigPlotter(cells,mutant1,mutant2,mutant3,genNum)


grids = 0;

[M,N,z] = size(cells);

mutantType = zeros(M,N);

for a1=1:M
    for a2=1:N

        if (cells(a1,a2,genNum) == 1)
            mutantType(a1,a2) = 1 + 4*mutant1(a1,a2,genNum) + ...
                2*mutant2(a1,a2,genNum) + mutant3(a1,a2,genNum);
        end

    end
end




mutantTypeTemp = zeros(M,N);

for a1=1:M
    mutantTypeTemp(a1,:) = mutantType(M+1-a1,:);
end

mutantType = mutantTypeTemp;

for k1 = 1:M
    for k2 = 1:N

        m = mutantType(k1,k2);

        if (m == 1)

            mutantType(k1,k2) = 1/5;

        elseif ((m == 2)|(m == 3)|(m == 4)|(m == 6))

            mutantType(k1,k2) = 2/5;

        elseif(m == 5)

            mutantType(k1,k2) = 3/5;

        elseif(m == 7)

            mutantType(k1,k2) = 4/5;

        elseif(m == 8)

            mutantType(k1,k2) = 1;

        end

    end
end

% 0 vacant - white
% 0.2 normal - grey
% 0.4 other - black
% 0.6 hyper - pink
% 0.8 glyc - green
% 1 resist - yellow

kmap = [1,1,1; .5,.5,.5; 0,0,0; 1,.5,.5; 0,1,0; 1,1,0];

close all
colormap(kmap)

imagesc(mutantType,[0,1])
if (grids == 1)
    set(gca,'XTick',.5:1:(N+.5),'XTickLabel',[],'YTick',0.5:1:(M+.5),'YTickLabel',[],...
        'Xgrid','on','Ygrid','on','TickLength',[0,0],'GridLineStyle',':',...
        'PlotBoxAspectRatio',[N,M,1])
else
    set(gca,'XTick',[],'YTick',[],'PlotBoxAspectRatio',[N,M,1])
end