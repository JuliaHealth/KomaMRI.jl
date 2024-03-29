seq=mr.Sequence();              % Create a new sequence object
fov=230e-3; Nx=100;             % Define FOV and resolution
alpha=30;                       % flip angle
sliceThickness=3e-3;            % slice
TE=10e-3;                        % TE; give a vector here to have multiple TEs (e.g. for field mapping)
TR=20e-3;                      % only a single value for now
Nr=fix(pi/2 * Nx);              % number of radial spokes
Ndummy=20;                      % number of dummy scans
delta=pi/ Nr;                  % angular increment; try golden angle pi*(3-5^0.5) or 0.5 of it

% more in-depth parameters
rfSpoilingInc=117;              % RF spoiling increment

% set system limits
sys = mr.opts('MaxGrad', 28, 'GradUnit', 'mT/m', ...
    'MaxSlew', 80, 'SlewUnit', 'T/m/s', 'rfRingdownTime', 20e-6, ...
    'rfDeadTime', 100e-6, 'adcDeadTime', 10e-6);

% Create alpha-degree slice selection pulse and gradient
[rf, gz] = mr.makeSincPulse(alpha*pi/180,'Duration',4e-3,...
    'SliceThickness',sliceThickness,'apodization',0.5,'timeBwProduct',4,'system',sys);

% Define other gradients and ADC events
deltak=1/fov;
gx = mr.makeTrapezoid('x','FlatArea',Nx*deltak,'FlatTime',6.4e-3,'system',sys);
adc = mr.makeAdc(Nx,'Duration',gx.flatTime,'Delay',gx.riseTime,'system',sys);
adc.delay = adc.delay;% - 0.5*adc.dwell; % compensate for the 0.5 samples shift
gxPre = mr.makeTrapezoid('x','Area',-gx.area/2,'Duration',2e-3,'system',sys);
gzReph = mr.makeTrapezoid('z','Area',-gz.area/2,'Duration',2e-3,'system',sys);

% gradient spoiling
gxSpoil=mr.makeTrapezoid('x','Area',0.5*Nx*deltak,'system',sys);
gzSpoil=mr.makeTrapezoid('z','Area',4/sliceThickness,'system',sys);

% Calculate timing
delayTE=ceil((TE - mr.calcDuration(gxPre) - gz.fallTime - gz.flatTime/2 ...
    - mr.calcDuration(gx)/2)/seq.gradRasterTime)*seq.gradRasterTime;
delayTR=ceil((TR - mr.calcDuration(gxPre) - mr.calcDuration(gz) ...
    - mr.calcDuration(gx) - delayTE)/seq.gradRasterTime)*seq.gradRasterTime;
assert(all(delayTR>=mr.calcDuration(gxSpoil,gzSpoil)));

rf_phase=0;
rf_inc=0;

for i=(-Ndummy):Nr
    for c=1:length(TE)
        rf.phaseOffset=rf_phase/180*pi;
        adc.phaseOffset=rf_phase/180*pi;
        rf_inc=mod(rf_inc+rfSpoilingInc, 360.0);
        rf_phase=mod(rf_phase+rf_inc, 360.0);
        seq.addBlock(rf,gz);
        phi=delta*(i-1);
        seq.addBlock(mr.rotate('z',phi,gxPre,gzReph));
        seq.addBlock(mr.makeDelay(delayTE(c)));
        if (i>0)
            seq.addBlock(mr.rotate('z',phi,gx,adc));
        else
            seq.addBlock(mr.rotate('z',phi,gx));
        end
        seq.addBlock(mr.rotate('z',phi,gzSpoil,mr.makeDelay(delayTR))); %Removed the gradient spoiler in the x direction
    end
end

seq.plot();

%
seq.setDefinition('FOV', [fov fov sliceThickness]);
seq.setDefinition('Name', 'gre_rad');

seq.write('gre_rad.seq')       % Write to pulseq file
%%
[ktraj_adc, t_adc, ktraj, t_ktraj, t_excitation, t_refocusing, ~, ~] = seq.calculateKspacePP();

time_axis=(1:(size(ktraj,2)))*sys.gradRasterTime;

figure; plot(t_ktraj*1e3, ktraj'); % plot the entire k-space trajectory
hold; plot(t_adc*1e3,ktraj_adc(1,:),'.'); % and sampling points on the kx-axis

figure; plot(ktraj(1,:),ktraj(2,:),'b',...
             ktraj_adc(1,:),ktraj_adc(2,:),'r.'); % a 2D plot
axis('equal'); % enforce aspect ratio for the correct trajectory display