%Variables = filename('MaestroTest137Cs.Spe'), NumPeaks(1) , xmin (408),xmax (649)
% centerx (540), centery(900), FWHM(5)
function prints = rip_me(filename, NumPeak, xmin, xmax, centerx, centery, FWHMguess)
% rip_me Returns the parameters after the fitting function runs
%   filename = String, location of the file
%   NumPeaks = int, number of peaks
%   xmin, xmax = Floats, the minimum and maximum possible x coordinates of
%   where the peaks are
%   centerx, centery = Floats, a guess of the x and y coordinates of the peak 
%   FWHMguess = Float, guess of how large the FWHM is 
%
%returns prints 
%A one dimensional array with 6 elements
% of the format [Centroid, dCentroid, FWHM, dFWHM, Area, dArea]







%filename='MaestroTest137Cs.Spe';
%NumPeak = 1;
%xmin = 408;
%xmax = 649;
%centerx = 540;
%centery = 900;
%FWHMguess = 5;

% CUT DATA %
Dat=fileread(filename);
FileSize=size(Dat);
for i=1:FileSize(2)
    if((Dat(i)=='D') && (Dat(i+1)=='A') && (Dat(i+2)=='T') && (Dat(i+3)=='A'))
        StartChar=i+13;
    end
    if((Dat(i)=='R') && (Dat(i+1)=='O') && (Dat(i+2)=='I'))
        EndChar=i-2;
    end
end

y=Dat(StartChar:EndChar);
y=str2num(y);
DataPoints=size(y);
x=zeros(DataPoints(1),1);

for i=1:DataPoints(1)
    x(i)=i;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% FITTING %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% USER INPUT%

% Number of peaks
NumPeaks=NumPeak;
Peaks=int64(NumPeaks);

% Initialise arrays to store the initial peak guesses
amps=zeros(Peaks);
cents=zeros(Peaks);
FWHMs=zeros(Peaks);

% Set the fit range




if(xmax < xmin)
    temp=xmin;
    xmin=xmax;
    xmax=temp;
end

ytrunc=y(xmin:xmax);
LineHeight=max(ytrunc);
LineHeight;

ylow=ytrunc(1);
yhigh=ytrunc(xmax-xmin);

if(ylow == 0)
    ylow=1;
end

if (yhigh == 0)
    yhigh=1;
end

xlinemin=[xmin xmin];
xlinemax=[xmax xmax];
yline=[0 LineHeight];

% Click on the peak centres to get their approximate centres and amplitudes

for i=1:Peaks
    amps(i)=centery;
    cents(i)=centerx;
    xlinepeak=[cents(i) cents(i)];
end

% Enter the FWHM of the peaks

for i=1:Peaks
    FWHMs(i)=FWHMguess;
end

% Calculate the initial guess for the linear background
grad=double((yhigh-ylow)/(xmax-xmin));
offset=double(yhigh-grad*xmax);

% Set peak amplitudes starting guesses and correct for the background
for i=1:Peaks
    centsapprox(i)=int64(cents(i));
    amps(i)=y(centsapprox(i));
    amps(i)=amps(i) - (grad*cents(i) + offset);
end
% Fill Peak starting guesses
% [Amp1 Centroid1 FWHM1 Amp2 Centroid2 FWHM 2 .... Bkd grad Bkd Offset]
params=zeros(1,Peaks*3 + 2);
for i=1:Peaks
    params((i-1)*3 + 1)=amps(i);
    params((i-1)*3 + 2)=cents(i);
    params((i-1)*3 + 3)=FWHMs(i);
end
params(3*Peaks + 1)=grad;
params(3*Peaks + 2)=offset;

%% FITTING PROCESS %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
NewX = double(x(xmin:xmax));
NewY = double(y(xmin:xmax));
NumParams=size(params);
NumParams=NumParams(2);
NumPeaks=(NumParams-2)/3;

init=params;

options=optimset('lsqcurvefit');
options.MaxFunEvals=50000;
options.MaxIter=5000;
options.TolFun=1e-9;
options.Display = 'none';

[NewParameters,error,residual,exitflag,output,lambda,Jac]=lsqcurvefit(@MultiGaussEqnLinearBkd,init,NewX,NewY);

new=MultiGaussEqnLinearBkd(NewParameters,NewX);

% Calculate the standard errors on the fit parameters
resid=NewY-new;
[Q,R] = qr(Jac,0);
mse = sum(abs(resid).^2)/(size(Jac,1)-size(Jac,2));
Rinv = inv(R);
Sigma = Rinv*Rinv'*mse;
se = sqrt(diag(Sigma));

% Plot the resulting fit over the data
prints = zeros(1,6);

for i=1:NumPeaks
    % Output fit parameters and their errors to the screen
    disp(' ');
    PeakLabel=['Peak ' num2str(i)];
    disp(PeakLabel);
    PeakCentroid=['Centroid ' num2str(abs(NewParameters(3*(i-1)+2))) ' +/- ' num2str(abs(se(3*(i-1)+2)))];
    disp(PeakCentroid);
    PeakWidth=['FWHM ' num2str(abs(NewParameters(3*(i-1)+3))) ' +/- ' num2str(abs(se(3*(i-1)+3)))];
    disp(PeakWidth);
    
    prints(1) = abs(NewParameters(3*(i-1)+2));
    prints(2) = abs(se(3*(i-1)+2));
    prints(3) = abs(NewParameters(3*(i-1)+3));
    prints(4) = abs(se(3*(i-1)+3));
    
    Cent=NewParameters(3*(i-1)+2);
    CentMin=NewParameters(3*(i-1)+2)-se(3*(i-1)+2);
    CentMax=NewParameters(3*(i-1)+2)+se(3*(i-1)+2);
    FWHM=NewParameters(3*(i-1)+3);
    FWHMMin=NewParameters(3*(i-1)+3)-se(3*(i-1)+3);
    FWHMMax=NewParameters(3*(i-1)+3)+se(3*(i-1)+3);
    sigma=FWHM/2.35;
    sigmaMin=FWHMMin/2.35;
    sigmaMax=FWHMMax/2.35;
    Amp=NewParameters(3*(i-1)+1);
    AmpMin=NewParameters(3*(i-1)+1)-se(3*(i-1)+1);
    AmpMax=NewParameters(3*(i-1)+1)+se(3*(i-1)+1);
    Area=Amp*sigma*sqrt(2*3.1415);
    AreaMin=AmpMin*sigmaMin*sqrt(2*3.1415);
    AreaMax=AmpMax*sigmaMax*sqrt(2*3.1415);
    AreaError=abs(AreaMax-AreaMin);
    

    
    PeakArea=['Area ' num2str(abs(Area)) ' +/- ' num2str(abs(AreaError))];
    disp(PeakArea);
    disp(' ');
    % Plot the individual peaks
    SinglePeakParams=[abs(Amp) abs(Cent) abs(FWHM) 0 0];
    yGauss=MultiGaussEqnLinearBkd(SinglePeakParams,NewX);
    
   
    prints(5) = abs(Area); 
    prints(6) = abs(AreaError); 
end

% Plot the background
% grad=NewParameters(NumParams-1);
% grad=NewParameters(NumParams);
% SinglePeakParams=[0 1 1 grad offset];
% yGauss=MultiGaussEqnLinearBkd(SinglePeakParams,NewX);
% plot(NewX,yGauss,'black','LineWidth',2);

%lsqcurvefit(fun,x0,x,y)
