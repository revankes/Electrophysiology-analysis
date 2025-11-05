%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%START UP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Upload your file, this is the nwb file. The NWB folder and the abfanalysis
%need to be in your path, which are added using addpath(genpath(folder)).
%If place is changed, this needs to be edited. 

clear
close all   

addpath(genpath('')) ; %add location for file path
addpath(genpath('')) %add location for data to analyse

%here you chose the file, which gets assigned to fn by concatinate strings (the file and path)
[file,path] = uigetfile({'*.nwb'},'Select a file');
%[file,path] = uigetfile('*.nwb','Select a file');

fn = strcat(path,file) ;

%if you want to analyse gap free,       chose 1
%if you want to analyse ramps,          chose 2
%if you want to analyse input output,   chose 3

IWantToAnalyse=3;

if IWantToAnalyse==1
   stimsetfilter={'Gap'}; %gapfree recording
end

if IWantToAnalyse==2
   stimsetfilter={'Ramp'}; %ramp recording
end

if IWantToAnalyse==3
   stimsetfilter={'Input'}; %input output recording 
end

%Chose the sweeps you want to analyse. If you want to analyse all sweeps,
%enter 0:100. 
sweepfilter=(1:100); %normal 1:35

%here you use the nwb file to get the sweeps for the filter you chose
%before. It also starts the analysis, for every sweep and makes a variable
%for the nrofaps and the epochs
nwb = NWBfile(fn, stimsetfilter, sweepfilter);

swps = nwb.getstimset().getnwbchannel(1).getsweep;
analysis = swps().analysesweep;

% Plots all figures
for n=1:numel(swps)
    figure(n)
    swps(n).plot
end 

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% RHEOBASE/ %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% the rheobase determines the minimum input threshold at which the cell fires.
% input is in the form of a ramp, so this means that it increases in a
% linear fashion, according to y = ax+b, where y (input) is the input in pA
% and a (t) being the constant and b the intersction. This differs per
% ramp, make sure to fill in the correct numbers fitting to the ramp you
% used!!!!


%Fill in the sweep of interest (=swpofinterest) from which you see the
%first amp. The function will then return the type of ramp (=ramptype) and
%calculate the input depending on the type of ramp.

% Threshold time is the time until the first fire. This will be calculated
% for the swp of interest you want. This you need to fill in.
% will calculate for the swp of interest you want, this is something you
% need to fill in. This will be prompted. 

prompt = "What is your sweep of interest?" ;
z = inputdlg (prompt) ;
swpofinterest = str2double(z) ;

thresholdtime = analysis(1,swpofinterest).aps(1,1).thresh_time ;
mvthreshold = analysis(1,swpofinterest).aps(1,1).thresh ;
checktheramptype = analysis(1,swpofinterest).stimwavename ; 

if analysis(1,swpofinterest).stimwavename == 'Ramp_100_DA_0'
    input = (0.1 * thresholdtime)-25 ;
    Ramptype = '{100 dA}' ;  
end

if analysis(1,swpofinterest).stimwavename == 'Ramp_250_DA_0'
    input = (0.25 * thresholdtime)-62.5 ;
    Ramptype = '{250 dA}' ;
end

if analysis(1,swpofinterest).stimwavename == 'Ramp_450_DA_0'
    input = (0.45 * thresholdtime)-112.5 ;
    Ramptype = '{450 dA}' ;
end

if analysis(1,swpofinterest).stimwavename == 'Ramp_400_DA_0'
    input = (0.40 * thresholdtime)-100 ;
    Ramptype = '{400 dA}' ;
end

if analysis(1,swpofinterest).stimwavename == 'Ramp_600_DA_0'
    input = (0.6 * thresholdtime)-150;
    Ramptype = '{600 dA}' ;
end

if analysis(1,swpofinterest).stimwavename == 'Ramp_900_DA_0' 
    input = (1.5 * thresholdtime)-375 ;
    Ramptype = '{900 dA}' ;
end

%table with the results
rheobaseresults = table(convertCharsToStrings(file), input, thresholdtime, mvthreshold, cellstr(Ramptype)) ;

%%
% INPUT OUTPUT 1/4
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%% INPUT OUTPUT CURVE%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%In the input output curve, you want to plot the amount of aps against the
%input (pA). This you can then do by going through the sweeps with 'Input'.
%The function below gives the ap frequency (frequency = 1/nr of aps). per
%input. 

%your first ap has a freq of zero, because you cannot determine an ap with one
%point. you do not take this into account with your average. So: 
% column 1: the sweep number
% column 2: gives the input, stepwise increase by 25 pA
% column 3: the amount of ap per sweep
% column 4: the average of freq per sweep without the first freq (this is always zero).
% column 5: the average of freq per sweep     

%Gives the ap frequency per sweep, for ALL sweeps
aps_frequency=NaN(numel(analysis),5) ;
for i=1:numel(analysis)
    numaps=analysis(i).nrofaps;
     if numaps>0 
        freq_Swps=NaN(analysis(i).nrofaps,1);
        for j=1:analysis(i).nrofaps
            freq_Swps(j)=analysis(i).aps(j).freq ;
        end
        aps_frequency(i,4)= mean(freq_Swps([2:end],1)) ; % gives av freq per sweep, without first freq
        aps_frequency(i,5)= mean(freq_Swps);  %gives av freq per sweep
     end
     aps_frequency(i,1)= i ;%sweep nr
     aps_frequency(i,2)= -100 + 25*i ;%input in pA  depends on which protocol you use!!!!!!!!!!!!!!!!!!!!!!!!!!!!! now it starts at -75 and increases with 25 pA (09/24)
     aps_frequency(i,3)= numaps ; %gives nr of ap
end

%%
% INPUT OUTPUT 2/4
%1/3%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%  time constant %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The membrane time capacitance says something about the rise and fall of the action potential
% It is calculated from V(t)=Vmax=(1-e^-t/timeconstant), while the fall is
% described by V(t) = Vmax*e^t/timeconstant. Fill in the smallest negative
% stimulation amplitude. This is at the input of -25 pA. 
%if the swps start at a different number than 1, change this in here.
swpstart = 1;

% here the datapoints and the timepoints are put in a table and assigned to
% x and y. 
% Check where the slope is and adjust accordingly
t = analysis(1,(swpstart+2)).epochs(1,2).Time;
d = analysis(1,(swpstart+2)).epochs(1,2).Data;
A=[];
A(:,1)= t;
A(:,2)= abs(d);
x = A(:,1)-A(1,1); %so it starts at zero 
y = A(:,2);
xcut=x(x>149)-149; % make sure your data starts at 0, otherwise your fit goes haywire (-149)
ycut=y(x>149);

%
close all
% There is crosstalk between d and tau, so better to decrease range.
% a-b*exp(-(x)/c)-d*(x)

% here the data is fitted to an exponential function
 myFitType1 = fittype('a-b*exp(-(x)/c)','indep','x'); %

[myModel1,gof1] = fit(xcut,ycut,myFitType1);
%[myModel1,gof1] = fit(x,y,myFitType1);
myTimeConst1 = myModel1.c;
rSquare1 = gof1.rsquare;

hold on
plot(xcut,ycut,'x','color','cyan')
plot(myModel1, 'r')
hold off

% If the tail messes up the gof, we can remove it to improve 
% not necessary for new data file, check how your curve fits
% however, we do not want the 'tail' of the data, as this gives a distorted
% image. here we only take the data 10 times from the time constant.
timeMultiplier = 4;  %marginalizing the tau for later. This is either 10 or 4 (more stringent)
cutofftail=myTimeConst1*timeMultiplier;
% 
newA(:,1)= xcut;
newA(:,2)= ycut;
newx = xcut(xcut<cutofftail);
newy = ycut(xcut<cutofftail);

hold on
plot(newx,newy,'x','color','black')
hold off

%Now fit again and look at the differences. nu halen we ook de d weg 
myFitType2 = fittype('a-b*exp(-(x)/c)','indep','x');

[myModel2,gof2] = fit(newx,newy,myFitType2) ;
myTimeConst2 = myModel2.c;
rSquare2 = gof2.rsquare;

hold on
plot(newx,newy,'x','color','cyan')
plot(myModel2, 'b')
hold off

resultstimeconstant = table(myTimeConst1,rSquare1, myTimeConst2,rSquare2);

%2/3%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%% Membrane capacitance %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%The membrane capacitance is the amount of energy stored/seperated for a
%potential. It is a function of the properties of the lipid bilayer and the
%resistance (i.e. number of open ion channels). Capacitance = time constant
%/ membrane resistance. Check which sweep you want to use!

deltaV3= (analysis(1,(swpstart+2)).epochs(1,2).meansignal) - (analysis(1,(swpstart+2)).epochs(1,3).meansignal); %of the 3rd sweep

deltaI3 = (abs(analysis(1,(swpstart+2)).epochs(1,3).amplitude)); %of the 3rd sweep, -25

membraneresistance = deltaV3/deltaI3;

%Resistance/Capacitance
membranecapacitance = myTimeConst1 / membraneresistance; % check which timeconstant you add in 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% INPUT RESISTANCE %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%The input resistance reflects to which extent the membrane channels are open, 
%and therefore the size of the cell also influences the input resistance.
%You calculate the input resistance using the change in voltage (Baseline - Vss) in response to the change
%in I. Check for how many sweeps you want to do this.


V = NaN(5,1);
I = NaN(5,1);

for i=(swpstart:(swpstart+4))
    Vss = swps(i).getsampleusingtime(900,990).median;
    Baseline = swps(i).getsampleusingtime(160,230).median;
    V(i) = abs(Vss - Baseline);
    I(i) = abs(swps(i).epochs(3).amplitude);
end

% The data (v en I) are put in a matrix
[xData, yData] = prepareCurveData( I, V );

% Fit model to data.
[fitresult, gof] = fit( xData, yData, 'poly1' ); %poly1 is a linear polynomial curve, now you will get the y=p1x+p2, where p1 is your input resistance

InputResistance = fitresult.p1*1000;
GoodnessOfFit = gof.adjrsquare;
resultsinputresistance = table(InputResistance, GoodnessOfFit);

%3/3%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Sag ratio %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The sag ratio is the ratio between steady state decrease in voltage and
% the largest decrease in voltage following a hyperpolarizing current step.
% This is of approximately -7.5 mV. Check at which sweep and epoch
% 
% sag ratio
% This means Vss (steady state voltage) - Vmin (min value after current
% injection) / Vmin - Vrmp (resting membrane potential).

sagratio_persweep=NaN(numel(analysis),2);
for i=1:3
    Vmin = min(analysis(1,(i)).epochs(1,2).Data) ; 
    Vss= analysis (1,(i)).epochs(1,2).steadystate ; 
    Vrmp = analysis(1,(i)).epochs(1,1).meansignal ; 
    sagratio_persweep(i,1)= abs((Vss-Vmin)/(Vmin - Vrmp)) ;
    sagratio_persweep(i,2)= Vss-Vrmp;
end

%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%fill in here!%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%linear interpolatie if your change in voltage is too far from -7.5. 

sv = [] ; %values that needs to be interprolated or extrapolated

%sv1 and sv3 =change in sagratio
%sv2 and sv4 =change in voltage

s1 = sv(1) ;
s2 = sv(3) ;
v1 = sv(2);
v2 = sv(4) ;
sagq = -7.5; %voltage change wanted

sagratio = [s1, s2] ;
voltagechange = [v1, v2] ;
points = [1, 2] ;

sagratio = [s1, s2] ;
voltagechange = [v1, v2] ;
points = [1, 2] ;

vquerysag = interp1(voltagechange, sagratio, sagq, 'linear'); % calculates the sagratio between two points
vquerysagex = interp1(voltagechange, sagratio, sagq, 'linear', 'extrap'); %when extrapolation is needed

restresults = table(myTimeConst1,rSquare1,myTimeConst2, rSquare2, membraneresistance,membranecapacitance,InputResistance,GoodnessOfFit,vquerysag,vquerysagex);


%% INPUT OUTPUT 3/4
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% AP characteristics %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%in this part the amplitude, the afterhyperpolarization voltage (ahp), the
%ahp time to peak as well as the speed of the depolariazation (dvdtmax) and
%repolarization (dvdtmin) are given. 

%The data is in the analysis file, in the particular sweep and
%action potential. The particular sweep you need to fill in, this is where
%the input generates about 3-5 action potentials. Fill this in at
%swpofinterest.

prompt = "What is your sweep of interest?" ;
z = inputdlg (prompt) ;
swpofinterest = str2double(z) ;

numberofaps = numel(analysis(1,swpofinterest).aps);
amplitudes = NaN(numberofaps,1);         

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%AP amplitude of the action potentials in the sweep of interest 
%all amps AMPLITUDE

for x=1:numberofaps
    amplitudes (x) = analysis(1,swpofinterest).aps(1,x).amp;
end
ampaverage = mean(amplitudes);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%here the afterhyperpolarization is calculated, which is the
%hyperpolarizing fase of the ap. This is where it falls below the normal
%resting potential. this is the undershoot phase. 
%the choice is either a ahp_slow or ahp
%%%%% this you need to choose depending on your cell type

%all ap AFTER HYPERPOLARIZATION relative
for x=1:numberofaps
    afterhyperp (x) = analysis(1,swpofinterest).aps(1,x).relahp;
end
afterhyperpaverage = mean(afterhyperp, 'omitnan');

%all ap AFTER HYPERPOLARIZATION SLOW relative
for x=1:numberofaps
    afterhyperpslow (x) = analysis(1,swpofinterest).aps(1,x).relahp_slow;
end
afterhyperpaverage_slow = mean(afterhyperpslow, 'omitnan');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%here the ahp time is given, this is the time from the peak to the ahp peak
%the choice is either a ahp_slow_time or ahp_time. This you need to choose
%depending on your cell type

%all amps AFTER HYPERPOLARIZATION TIME
for x=1:numberofaps
    afterhyperptime (x) = (analysis(1,swpofinterest).aps(1,x).ahp_time)-(analysis(1,swpofinterest).aps(1,x).peak_time);
end
afterhyperptimeaverage = mean(afterhyperptime, 'omitnan');

%all amps AFTER HYPERPOLARIZATION TIME SLOW

for x=1:numberofaps
    afterhyperptimeslow (x) = (analysis(1,swpofinterest).aps(1,x).ahp_slow_time)-(analysis(1,swpofinterest).aps(1,x).peak_time);
end
afterhyperptimeaverage_slow = mean(afterhyperptimeslow, 'omitnan');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% maxdvdt is the max speed of the depolarization
for x=1:numberofaps
    dvdtmax (x) = analysis(1,swpofinterest).aps(1,x).maxdvdt;
end
dvdtmaxaverage = mean(dvdtmax);

% mindvdt is the speed of the repolarization
for x=1:numberofaps
    dvdtmin (x) = analysis(1,swpofinterest).aps(1,x).mindvdt;
end
dvdtminaverage = mean(dvdtmin);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%AP halfwidth
for x=1:numberofaps
    aphalfwidth (x) = analysis(1,swpofinterest).aps(1,x).halfwidth;
end
aphalfwidthaverage = mean(aphalfwidth);

%AP characteristics table
APcharacteristics = table (swpofinterest, ampaverage, afterhyperpaverage, afterhyperpaverage_slow, afterhyperptimeaverage, afterhyperptimeaverage_slow, dvdtmaxaverage, dvdtminaverage, aphalfwidthaverage);

%% INPUT OUTPUT 4/4
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% mean ISI interval - adaptation ratio %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Mean ISI interval and adaptation ratio (1st ISI/9th ISI) will tell the
% time between the first and second AP. Therefore, the first ISI is
% skipped, which is why it is now x+1. 

% Type in the sweep of interest.
% Type in the sweep that has at least 10 ap.

prompt = "What is your sweep of interest?" ;
z = inputdlg (prompt) ;
swpofinterest2 = str2double(z) ;

numberofaps = numel(analysis(1,swpofinterest2).aps);
amplitudes = NaN(numberofaps,1);

%Mean ISI. 
for x=1:(numberofaps-1)
    interspikein (x) = analysis(1,swpofinterest2).aps(1,(x+1)).isi;
end
interspikeinaverage = mean(interspikein);
interspikeinaverage2 = mean(interspikein(1:2));

%adaptation ratio.
%Type in the sweep that has at least 10 ap.

adaptationratio = (analysis(1,swpofinterest2).aps(1,(2)).isi)/(analysis(1,swpofinterest2).aps(1,(9)).isi);

%AP characteristics table2
APcharacteristics2 = table (swpofinterest2, interspikeinaverage,interspikeinaverage2,adaptationratio);

%AP characteristics table total
APchar = table (swpofinterest, ampaverage, afterhyperpaverage, afterhyperpaverage_slow, afterhyperptimeaverage, afterhyperptimeaverage_slow, dvdtmaxaverage, dvdtminaverage, aphalfwidthaverage,swpofinterest2, interspikeinaverage,interspikeinaverage2,adaptationratio) ;

%% RMP 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% RESTING MEMBRANE POTENTIAL %%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%determines your average and median resting membrane potential ,
%per determined in a certain sweep (time_meansignal) Using the ginput
%function, you are able to get the begin and endtime by clicking on the
%graph. These are returned in X (time) and put in new variables (begintime
%and endtime). Then these are used to calculate the rmp. 

%Swp number you still have to fill in. 
prompt = "What is your sweep of interest?" ;
z = inputdlg (prompt) ;
swpnumber = str2double(z) ; % str to double, otherwise can't use it in the function

[x,y] = ginput (2) ;
begintime = x(1) ;
endtime = x(2) ;

swps(swpnumber).getsampleusingtime(begintime,endtime).plot
meansignal_time=mean(swps(swpnumber).getsampleusingtime(begintime,endtime).Data);
mediansignal_time=median(swps(swpnumber).getsampleusingtime(begintime,endtime).Data);

%now putting everything in a table
restingmembraneresults=table(convertCharsToStrings(file), swpnumber, begintime, endtime, meansignal_time, mediansignal_time);


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% STOP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% PROTOCOL PER SWEEPS %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% gives the protocol per swp and the original swp number
stimprotocols= cell(6,3)

for n=1:numel(swps)
    stimprotocols(n,1)= {n};
    stimprotocols(n,2)= {swps(n).number};
    stimprotocols(n,3)= {swps(n).stimwavename};
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% PLOTS ALL SWEEPS %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for n=1:numel(swps)
    figure(n)
    swps(n).plot
end 

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% PLOTS part of SWEEPS %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%plost part of sweep, enter begin x value and end x value. To see the x
%value, use the plot of the whole sweep.

swps(x).getsampleusingtime(begin,end).plot

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% PLOTS ANY SWEEP %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plots a graph of sweep of choice, the number you fill in is the stimset nr, not total sweeps

swps(3).plot

%!!!important!!!!
%if it is the old version -> up to 08/21, use the filespath old. here in
%line 68 you will need two groups. in line 86      swps={info.Groups(1).Groups(2).Groups.Name};
%if it is in the new version -> from 09/21, use the filespath. Here you
%only need one group.  in line 86      swps={info.Groups(1).Groups().Groups.Name};

