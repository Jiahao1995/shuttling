% % Constants
sparams.ee = 1.602E-19; % Charge of electron
sparams.me = 9.11E-31*0.191; % Effective mass of electron
sparams.hbar = 6.626E-34/(2*pi); % Reduced Planck's constant
sparams.c = 2.998E8; % Speed of light

% % Miscellaneous simulation parameters
sparams.updateWaitbar = 5000; % How many frames to wait before updating waitbar
sparams.updateFigure = 2500; % How many frames to wait before updating figure
sparams.updateFidelity = 4000; % How many frames to wait before calculating fidelity
sparams.updateInterpPot = 40000; % Interpolate the potential at 40000 points at a time
sparams.nFigureFrames = 75; % How many frames of the simulation to save for outputted gif

% % Main simulation parameters
% sparams.potFile = 'Shuttling3gates_121steps.xlsx'; % Files to load potentials from for simulation
sparams.potDir = 'simulatedPotentials/simulation_0328/'; % Directory where the pot files are
sparams.interpPotentials = 1; % Whether or not to interpolate the potentials in time domain for simulation
sparams.interpType = 'linear'; % What type of interpolation of the potentials to do in time domain

sparams.dt = 5E-17; % Time between each simulation frame [sec]
% sparams.totalTime = [5E-11,logspace(-10,-9,10),2E-9,5E-9,10E-9]; % Total time of simulted pulse [sec]
% sparams.totalTime = [1E-11,3E-11];
sparams.totalTime = 5E-10;
