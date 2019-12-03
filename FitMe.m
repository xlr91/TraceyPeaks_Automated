clc; clear; clear all; close all;


%Lines to edit depending on the peaks/folder location
%All files need to be put in a single folder
folderstring = '~/Desktop/Nuclear_Lab_files/Tracey_Peaker_V0_01_17/source/';
x_min = 1330;
x_max = 1340;
center_x = 1335;
center_y = 15000;
FWHM_guess = 10;


%rips the filenames of each of the thingys
folderpath = fullfile(folderstring, '**');
filelist   = dir(folderpath);
name       = {filelist.name};
name       = name(~strncmp(name, '.', 1));   % No files starting with '.'
sizee = size(name);
parameters = zeros(sizee(2), 6);




%runs the traceypeakerV0 program for all files
for k=1:sizee(2)
    %disp(name(k))
    stringg = string(append(folderstring, name(k)));
    disp(stringg);
    print = rip_me(stringg, 1, x_min, x_max, center_x, center_y, FWHM_guess);
    parameters(k,:) = print;
end

fprintf('finito')
