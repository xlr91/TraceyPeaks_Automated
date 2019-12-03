clc; clear; clear all; close all;

%rips the filenames of each of the thingys
folderstring = '~/Desktop/Nuclear_Lab_files/Tracey_Peaker_V0_01_17/source/';
folderpath = fullfile(folderstring, '**');
filelist   = dir(folderpath);
name       = {filelist.name};
name       = name(~strncmp(name, '.', 1));   % No files starting with '.'
sizee = size(name);
parameters = zeros(sizee(2), 6);


%setting up variables for the peak
%prints = rip_me(filename, NumPeak, xmin, xmax, centerx, centery, FWHMguess)

x_min = 1330;
x_max = 1340;
center_x = 1335;
center_y = 15000;
FWHM_guess = 10;

%runs the traceypeakerV0 program for all files
for k=1:sizee(2)
    %disp(name(k))
    stringg = string(append(folderstring, name(k)));
    disp(stringg);
    print = rip_me(stringg, 1, x_min, x_max, center_x, center_y, FWHM_guess);
    parameters(k,:) = print;
end

fprintf('finito')
