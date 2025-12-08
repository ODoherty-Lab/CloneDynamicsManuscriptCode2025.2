dataposstat5b = [1.42 2.38 8.19 10 10.4 10.99 11.33 12 1.71 10.09 12.61 13.2 13.74 14.24; 19 23 25 18 99 51 268 120 30 18 16 70 23 63; 255 556 516 222 1139 791 3601 1439 510 414 317 600 256 620];  
datanegstat5b = [1.42 2.38 8.19 10 10.4 10.99 11.33 12 1.71 10.09 12.61 13.2 13.74 14.24; 9 23 15 7 41 26 148 60 19 13 7 18 15 23; 255 556 516 222 1139 791 3601 1439 510 414 317 600 256 620]; 

data = dataposstat5b; % making a copy to edit later.
data1 = datanegstat5b;

paraset1 = [0.0002 0.0363];
paraset = [0.0033 0.0398];

paraset2 = [-0.0001 0.0359];
paraset3 = [0.0038 0.04];

m = paraset(1,1);b=paraset(1,2);m1 = paraset1(1,1);b1=paraset1(1,2);
m2 = paraset2(1,1);b2=paraset2(1,2);m3 = paraset3(1,1);b3=paraset3(1,2);

cost1 = 0;
cost = 0;
cost2 = 0;
cost3 = 0;
N=1e6;

for i = 1:8
    y1 = m1*data1(1,i)+b1;
    y = m*data(1,i)+b;
    cost1 =cost1+ log(pdf('Binomial',data1(2,i),data1(3,i),y1));
    cost =cost+ log(pdf('Binomial',data(2,i),data(3,i),y));
    end
likelihood1 = cost1;
likelihood = cost;

for i = 9:14
    y1 = m3*data1(1,i)+b3;
    y = m2*data(1,i)+b2;
    cost2 =cost2+ log(pdf('Binomial',data1(2,i),data1(3,i),y));
    cost3 =cost3+ log(pdf('Binomial',data(2,i),data(3,i),y1));
end
likelihood2 = cost2;
likelihood3 = cost3;


for k = 2:N
    m = random('Uniform',-0.002,0.008);
    b = random('Uniform',-0.01,0.1);
    m1 = random('Uniform',-0.005,0.005);
    b1 = random('Uniform',-0.01,0.1);
    m2 = random('Uniform',-0.005,0.005);
    b2 = random('Uniform',-0.01,0.1);
    m3 = random('Uniform',-0.005,0.01);
    b3 = random('Uniform',-0.01,0.1);
cost1 = 0;
cost = 0;
cost2 = 0;
cost3 = 0;
for i = 1:8
    y1 = m1*data1(1,i)+b1;
    y = m*data(1,i)+b;
    cost1 =cost1+ log(pdf('Binomial',data1(2,i),data1(3,i),y1));
    cost =cost+ log(pdf('Binomial',data(2,i),data(3,i),y));
end

for i = 9:14
    y1 = m3*data1(1,i)+b3;
    y = m2*data(1,i)+b2;
    cost2 =cost2+ max(log(pdf('Binomial',data1(2,i),data1(3,i),y)),-1e6);
    cost3 =cost3+ max(log(pdf('Binomial',data(2,i),data(3,i),y1)),-1e6);
end

if cost1>likelihood1(k-1)+log(rand(1))
    likelihood1 = [likelihood1;cost1];
    paraset1 = [paraset1;m1 b1];
else
    likelihood1 = [likelihood1;likelihood1(k-1)];
    paraset1 = [paraset1;paraset1(k-1,:)];
end

if cost>likelihood(k-1)+log(rand(1))
    likelihood = [likelihood;cost];
    paraset = [paraset;m b];
else
    likelihood = [likelihood;likelihood(k-1)];
    paraset = [paraset;paraset(k-1,:)];
end

if cost2>likelihood2(k-1)+log(rand(1))
    likelihood2 = [likelihood2;cost2];
    paraset2 = [paraset2;m2 b2];
else
    likelihood2 = [likelihood2;likelihood2(k-1)];
    paraset2 = [paraset2;paraset2(k-1,:)];
end

if cost3>likelihood3(k-1)+log(rand(1))
    likelihood3 = [likelihood3;cost3];
    paraset3 = [paraset3;m3 b3];
else
    likelihood3 = [likelihood3;likelihood3(k-1)];
    paraset3 = [paraset3;paraset3(k-1,:)];
end
end

    % m = random('Uniform',-0.002,0.008);
    % b = random('Uniform',-0.01,0.1);
    % m1 = random('Uniform',-0.005,0.005);
    % b1 = random('Uniform',-0.01,0.1);
    % m2 = random('Uniform',-0.005,0.005);
    % b2 = random('Uniform',-0.01,0.1);
    % m3 = random('Uniform',-0.005,0.01);
    % b3 = random('Uniform',-0.01,0.1);
figure;%sgtitle('Growth Related Genes','Fontsize',12,'FontWeight','bold');
subplot(4,2,1);h11=histogram(([paraset(:,1)]),'Normalization','probability','facecolor','0.64,0.84,0.61','edgecolor','none');hold on;line([-0.002;0.008],[0;0],'Color','r','LineWidth',2,'Marker','+','MarkerSize',10);title('Slope CT1 Sense');xlabel('Rate of Change (fraction of total genes per year)');legend('Posterior','Uniform Prior');
subplot(4,2,2);h1=histogram(([paraset(:,2)]),'Normalization','probability','facecolor','0.64,0.84,0.61','edgecolor','none');hold on;line(([-0.01;0.1]),[0;0],'Color','r','LineWidth',2,'Marker','+','MarkerSize',10);title('Intercept CT1 Sense');xlabel('Initial fraction of total genes');legend('Posterior','Uniform Prior');
subplot(4,2,3);h21=histogram(([paraset1(:,1)]),'Normalization','probability','facecolor','0.28,0.57,0.77','edgecolor','none');hold on;line([-0.005;0.005],[0;0],'Color','r','LineWidth',2,'Marker','+','MarkerSize',10);title('Slope CT1 Antisense');xlabel('Rate of Change (fraction of total genes per year)');legend('Posterior','Uniform Prior');
subplot(4,2,4);h2=histogram(([paraset1(:,2)]),'Normalization','probability','facecolor','0.28,0.57,0.77','edgecolor','none');hold on;line(([-0.01;0.1]),[0;0],'Color','r','LineWidth',2,'Marker','+','MarkerSize',10);title('Intercept CT1 Antisense');xlabel('Initial fraction of total genes');legend('Posterior','Uniform Prior');
subplot(4,2,5);h31=histogram(([paraset3(:,1)]),'Normalization','probability','facecolor','1,1,0','edgecolor','none');hold on;line([-0.005;0.005],[0;0],'Color','r','LineWidth',2,'Marker','+','MarkerSize',10);title('Slope CT2 Sense');xlabel('Rate of Change (fraction of total genes per year)');legend('Posterior','Uniform Prior');
subplot(4,2,6);h3=histogram(([paraset3(:,2)]),'Normalization','probability','facecolor','1,1,0','edgecolor','none');hold on;line(([-0.01;0.1]),[0;0],'Color','r','LineWidth',2,'Marker','+','MarkerSize',10);title('Intercept CT2 Sense');xlabel('Initial fraction of total genes');legend('Posterior','Uniform Prior');
subplot(4,2,7);h41=histogram(([paraset2(:,1)]),'Normalization','probability','facecolor','1,0,1','edgecolor','none');hold on;line([-0.005;0.01],[0;0],'Color','r','LineWidth',2,'Marker','+','MarkerSize',10);title('Slope CT2 Antisense');xlabel('Rate of Change (fraction of total genes per year)');legend('Posterior','Uniform Prior');
subplot(4,2,8);h4=histogram(([paraset2(:,2)]),'Normalization','probability','facecolor','1,0,1','edgecolor','none');hold on;line(([-0.01;0.1]),[0;0],'Color','r','LineWidth',2,'Marker','+','MarkerSize',10);title('Intercept CT2 Antisense');xlabel('Initial fraction of total genes');legend('Posterior','Uniform Prior');


figure;h11=histogram(([paraset(:,1)]),'Normalization','probability','facecolor','0.64,0.84,0.61','edgecolor','none');hold on;line([-0.008;0.008],[0;0],'Color','r','LineWidth',2,'Marker','+','MarkerSize',10);title('Growth Related Genes P1 Sense');xlabel('Rate of Change (fraction of total genes per year)');legend('Posterior','Uniform Prior');
figure;h1=histogram(([paraset(:,2)]),'Normalization','probability','facecolor','0.64,0.84,0.61','edgecolor','none');hold on;line(([-0.01;0.1]),[0;0],'Color','r','LineWidth',2,'Marker','+','MarkerSize',10);title('Growth Related Genes P1 Sense');xlabel('Initial fraction of total genes');legend('Posterior','Uniform Prior');
figure;h21=histogram(([paraset1(:,1)]),'Normalization','probability','facecolor','0.28,0.57,0.77','edgecolor','none');hold on;line([-0.005;0.005],[0;0],'Color','r','LineWidth',2,'Marker','+','MarkerSize',10);title('Growth Related Genes P1 Antisense');xlabel('Rate of Change (fraction of total genes per year)');legend('Posterior','Uniform Prior');
figure;h2=histogram(([paraset1(:,2)]),'Normalization','probability','facecolor','0.28,0.57,0.77','edgecolor','none');hold on;line(([-0.01;0.15]),[0;0],'Color','r','LineWidth',2,'Marker','+','MarkerSize',10);title('Growth Related Genes P1 Antisense');xlabel('Initial fraction of total genes');legend('Posterior','Uniform Prior');
figure;h31=histogram(([paraset2(:,1)]),'Normalization','probability','facecolor','1,1,0','edgecolor','none');hold on;line([-0.008;0.005],[0;0],'Color','r','LineWidth',2,'Marker','+','MarkerSize',10);title('Growth Related Genes P2 Antisense');xlabel('Rate of Change (fraction of total genes per year)');legend('Posterior','Uniform Prior');
figure;h3=histogram(([paraset2(:,2)]),'Normalization','probability','facecolor','1,1,0','edgecolor','none');hold on;line(([-0.1;0.2]),[0;0],'Color','r','LineWidth',2,'Marker','+','MarkerSize',10);title('Growth Related Genes P2 Antisense');xlabel('Initial fraction of total genes');legend('Posterior','Uniform Prior');
figure;h41=histogram(([paraset3(:,1)]),'Normalization','probability','facecolor','1,0,1','edgecolor','none');hold on;line([-0.005;0.01],[0;0],'Color','r','LineWidth',2,'Marker','+','MarkerSize',10);title('Growth Related Genes P2 Sense');xlabel('Rate of Change (fraction of total genes per year)');legend('Posterior','Uniform Prior');
figure;h4=histogram(([paraset3(:,2)]),'Normalization','probability','facecolor','1,0,1','edgecolor','none');hold on;line(([-0.01;0.1]),[0;0],'Color','r','LineWidth',2,'Marker','+','MarkerSize',10);title('Growth Related Genes P2 Sense');xlabel('Initial fraction of total genes');legend('Posterior','Uniform Prior');
figure;histogram(([paraset1(:,1)]),'Normalization','probability','facecolor','0.28,0.57,0.77','edgecolor','none');hold on;histogram(([paraset(:,1)]),'Normalization','probability','facecolor','0.64,0.84,0.61','edgecolor','none');histogram(([paraset2(:,1)]),'Normalization','probability','facecolor','1,0,1','edgecolor','none');histogram(([paraset3(:,1)]),'Normalization','probability','facecolor','1,1,0','edgecolor','none');
legend('P1 Antisense','P1 Sense','P2 Antisense','P2 Sense');title('Growth-Related Genes');xlabel('Rate of Change (fraction of total genes per year)');

[like1sort,like1ind]=sort(likelihood1);
[likesort,likeind]=sort(likelihood);
[like2sort,like2ind]=sort(likelihood2);
[like3sort,like3ind]=sort(likelihood3);
time = linspace(-1,15);
m = paraset(likeind(end),1);b = paraset(likeind(end),2);
y = m*time+b;ymax = y; ymin = y;
m1 = paraset1(like1ind(end),1);b1 = paraset1(like1ind(end),2);
y1 = m1*time+b1;y1max = y1; y1min = y1;
m2 = paraset2(like2ind(end),1);b2 = paraset2(like2ind(end),2);
y2 = m2*time+b2;y2max = y2; y2min = y2;
m3 = paraset3(like3ind(end),1);b3 = paraset3(like3ind(end),2);
y3 = m3*time+b3;y3max = y3; y3min = y3;

for i = round(N/20)+1:1:N
    m = paraset(likeind(i),1);b = paraset(likeind(i),2);
    ytest = m*time+b;ymax = max(ymax,ytest); ymin = min(ymin,ytest);
    m1 = paraset1(like1ind(i),1);b1 = paraset1(like1ind(i),2);
    y1test = m1*time+b1;y1max = max(y1max,y1test); y1min = min(y1min,y1test);
    m2 = paraset2(like2ind(i),1);b2 = paraset2(like2ind(i),2);
    y2test = m2*time+b2;y2max = max(y2max,y2test); y2min = min(y2min,y2test);
    m3 = paraset3(like3ind(i),1);b3 = paraset3(like3ind(i),2);
    y3test = m3*time+b3;y3max = max(y3max,y3test); y3min = min(y3min,y3test);
end
figure;plot(time,y1,'Color','0.28,0.57,0.77','LineWidth',2,'LineStyle','-');
hold on;
plot(time,y,'Color','0.64,0.84,0.61','LineWidth',2,'LineStyle','-');
plot(time,y2,'Color','1,0,1','LineWidth',2,'LineStyle','-');
plot(time,y3,'Color','1,1,0','LineWidth',2,'LineStyle','-');
fill([time, fliplr(time)],[y2min fliplr(y2max)],1,'facecolor','1,0,1','edgecolor','none','facealpha',0.2)
fill([time, fliplr(time)],[y3min fliplr(y3max)],1,'facecolor','1,1,0','edgecolor','none','facealpha',0.2)
fill([time, fliplr(time)],[y1min fliplr(y1max)],1,'facecolor','0.28,0.57,0.77','edgecolor','none','facealpha',0.2)
fill([time, fliplr(time)],[ymin fliplr(ymax)],1,'facecolor','0.64,0.84,0.61','edgecolor','none','facealpha',0.2)

plot(data(1,1:8),data(2,1:8)./data(3,1:8),'LineStyle','none','Marker','o','MarkerFaceColor','0.64,0.84,0.61','MarkerEdgeColor','0,0,0','MarkerSize',10);
plot(data(1,9:14),data(2,9:14)./data(3,9:14),'LineStyle','none','Marker','^','MarkerFaceColor','1,1,0','MarkerEdgeColor','0,0,0','MarkerSize',10);
plot(data1(1,1:8),data1(2,1:8)./data1(3,1:8),'LineStyle','none','Marker','o','MarkerFaceColor','0.28,0.57,0.77','MarkerEdgeColor','0,0,0','MarkerSize',10);
plot(data1(1,9:14),data1(2,9:14)./data1(3,9:14),'LineStyle','none','Marker','^','MarkerFaceColor','1,0,1','MarkerEdgeColor','0,0,0','MarkerSize',10);
hold off;
axis([0 14.5 -0.01 0.1])
legend('P1 Antisense','P1 Sense','P2 Antisense','P2 Sense');xlabel('Time since treatment initiation');ylabel('Fraction of integrated HIV genomes');title('Growth-Related Genes');

size(unique(paraset),1)/N*100
size(unique(paraset1),1)/N*100
size(unique(paraset2),1)/N*100
size(unique(paraset3),1)/N*100


prctile(paraset(:,1),[0 2.5 50 97.5 100])
prctile((paraset(:,2)),[0 2.5 50 97.5 100])
prctile(paraset1(:,1),[0 2.5 50 97.5 100])
prctile((paraset1(:,2)),[0 2.5 50 97.5 100])

prctile(paraset2(:,1),[0 2.5 50 97.5 100])
prctile((paraset2(:,2)),[0 2.5 50 97.5 100])
prctile(paraset3(:,1),[0 2.5 50 97.5 100])
prctile((paraset3(:,2)),[0 2.5 50 97.5 100])

%std(log(paraset(:,2)))
% To stop a line from showing in the legend, select in plot browser and
% type
% set(get(get(gco,'Annotation'),'LegendInformation'),'IconDisplayStyle','off')


time = [2.38 1.42 8.19 10.99 10 11.33 12 10.4];
data = [19 23 25 18 99 51 268 120]./[556 255 516 791 222 3601 1439 1139];
m = paraset(likeind(end),1);b = paraset(likeind(end),2);
y = m*time+b;ymax = y; ymin = y;
res = y-data;
data1 = [9 23 15 7 41 26 148 60]./[556 255 516 791 222 3601 1439 1139];
m1 = paraset1(like1ind(end),1);b1 = paraset1(like1ind(end),2);
y1 = m1*time+b1;y1max = y1; y1min = y1;
res1 = y1-data1;

time = [1.71 10.09 12.61 13.74 13.2 14.24];
data2 = [19 13 7 18 15 23]./[510 414 317 256 600 620];

m2 = paraset2(like2ind(end),1);b2 = paraset2(like2ind(end),2);
y2 = m2*time+b2;y2max = y2; y2min = y2;
res2 = y2-data2;
data3 = [30 18 16 70 23 63]./[510 414 317 256 600 620];
m3 = paraset3(like3ind(end),1);b3 = paraset3(like3ind(end),2);
y3 = m3*time+b3;y3max = y3; y3min = y3;
res3 = y3-data3;

figure;qqplot(res);title('CT1 Sense');figure;qqplot(res1);title('CT1 AntiSense');figure;qqplot(res2);title('CT2 AntiSense');figure;qqplot(res3);title('CT2 Sense')
