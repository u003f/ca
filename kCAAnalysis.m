function [m1,m2,m3] = kCAAnalysis(cells,mutant1,mutant2,mutant3,P)

% EDIT, 19 MAY 2014: ignore cell=-1 (dead, not yet removed)
cells = max(cells,0);

if (nargin < 5)
    P = 95;
end

[M,N,z] = size(cells);

c = reshape(sum(sum(cells),2),z,1);
m1 = reshape(sum(sum(mutant1)),z,1)./c;
m2 = reshape(sum(sum(mutant2)),z,1)./c;
m3 = reshape(sum(sum(mutant3)),z,1)./c;

close all
figure
hold on

% plot(1:z,m1,'b')
% plot(1:z,m2,'g')
% plot(1:z,m3,'r')
% legend('Hyperplastic','Glycolytic','Acid-resistant')
% xlabel('Generation number','FontSize',12)
% ylabel('Proportion of cells','FontSize',12)

plot(1:z,m1,'k--')
plot(1:z,m2,'k:')
plot(1:z,m3,'k-')
xlabel('Generation number','FontSize',12)
ylabel('Proportion of cells','FontSize',12)


temp = find(m1 >P/100);
if (isempty(temp))
    m1_P = 0;
else
    m1_P = temp(1);
end

temp = find(m2 >P/100);
if (isempty(temp))
    m2_P = 0;
else
    m2_P = temp(1);
end

temp = find(m3 >P/100);
if (isempty(temp))
    m3_P = 0;
else
    m3_P = temp(1);
end

axis([0,z,0,1])

disp(['Maximum system size: M=',int2str(M),' cells deep.']);
disp([int2str(P),'% cells with mutation 1 after ',int2str(m1_P),' generations.']);
disp([int2str(P),'% cells with mutation 2 after ',int2str(m2_P),' generations.']);
disp([int2str(P),'% cells with mutation 3 after ',int2str(m3_P),' generations.']);
box