%% Controls
dataForMatlabAltEstTick = importfile('dataForMatlabAltEstTick.csv')
years=dataForMatlabAltEstTick.year;
index=(dataForMatlabAltEstTick.FOMCused==1);
years(index==1)=[];
for i=1:9
yearcountC(i,:)=years==2003+i;
howmuchC(i)=round(mean(yearcountC(i,:)),7);
end
datesC{1}=linspace(2004,2005,howmuchC(1)*length(years));
datesC{2}=linspace(2005,2006,howmuchC(2)*length(years));
datesC{3}=linspace(2006,2007,howmuchC(3)*length(years));
datesC{4}=linspace(2007,2008,howmuchC(4)*length(years));
datesC{5}=linspace(2008,2009,howmuchC(5)*length(years));
datesC{6}=linspace(2009,2010,howmuchC(6)*length(years));
datesC{7}=linspace(2010,2011,howmuchC(7)*length(years));
datesC{8}=linspace(2011,2012,howmuchC(8)*length(years));
datesC{9}=linspace(2012,2013,howmuchC(8)*length(years));
timeC=[datesC{1} datesC{2} datesC{3} datesC{4} datesC{5} datesC{6} datesC{7} datesC{8} datesC{9}];
%% Treatments
years=dataForMatlabAltEstTick.year;
index=(dataForMatlabAltEstTick.FOMCused==0);
years(index==1)=[];
for i=1:11
yearcountP(i,:)=years==2003+i;
howmuchP(i)=round(mean(yearcountP(i,:)),50);
end
datesP{1}=linspace(2004,2005,howmuchP(1)*length(years));
datesP{2}=linspace(2005,2006,howmuchP(2)*length(years));
datesP{3}=linspace(2006,2007,howmuchP(3)*length(years));
datesP{4}=linspace(2007,2008,howmuchP(4)*length(years));
datesP{5}=linspace(2008,2009,howmuchP(5)*length(years));
datesP{6}=linspace(2009,2010,howmuchP(6)*length(years));
datesP{7}=linspace(2010,2011,howmuchP(7)*length(years));
datesP{8}=linspace(2011,2012,howmuchP(8)*length(years));
datesP{9}=linspace(2012,2013,howmuchP(9)*length(years));
datesP{10}=linspace(2013,2014,howmuchP(10)*length(years));
datesP{11}=linspace(2014,2015,howmuchP(11)*length(years));
timeP=[datesP{1} datesP{2} datesP{3} datesP{4} datesP{5} datesP{6} datesP{7} datesP{8} datesP{9} datesP{10} datesP{11}];

save('time.mat', 'timeC', 'timeP')