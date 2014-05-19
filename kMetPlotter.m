function kMetPlotter(C,c0,h0,g0)

c = c0(:,:,end);
h = h0(:,:,end);
g = g0(:,:,end);

[a,t1,t2] = size(C);
X = 0:a;

C = mean(c,2); C = [1,C'];
H = mean(h,2); H = [0,H'];
G = mean(g,2); G = [1,G'];


close all
plot(X,G,'b-');
hold on; 
plot(X,C,'g-');
plot(X,H/max(H),'r-');
hold off
xlabel('Distance from basement membrane','FontSize',12)
ylabel('Metabolite level','FontSize',12)
axis([0,a,0,1.05])
legend('g','c','h')

% close all
% plot(X,G,'k--');
% hold on; 
% plot(X,C,'k-');
% plot(X,H/max(H),'k:');
% hold off
% xlabel('Distance from basement membrane','FontSize',12)
% ylabel('Metabolite level','FontSize',12)
% axis([0,a,0,1.05])
box