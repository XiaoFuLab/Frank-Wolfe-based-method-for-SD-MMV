clear

load results_Jan21.mat

p = [6,7,1:5];
q = [4,1:3];

x = categorical({'MAG1','MAG2','DBLP_{-1}','DBLP_{-2}','DBLP_{-3}','DBLP_{-4}','DBLP_{-5}'});
x = reordercats(x,{'MAG1','MAG2','DBLP_{-1}','DBLP_{-2}','DBLP_{-3}','DBLP_{-4}','DBLP_{-5}'});

src = diag(1./[3,3,13,16,16,16,15])*src(q,p)';
ttt = ttt(q,p)';

subplot(2,1,1)
bar(x,src)
ylabel('SRC_{avg}')

legend('CD-MVSI','GeoNMF','SPOC','tensor CPD','Orientation','horizontal')

subplot(2,1,2)
bar(x,ttt)
set(gca,'yscale','log')
ylabel('run time')
