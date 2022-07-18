cd("/home/ccp/Documents/pulseq-1.4.0/matlab/")
%%
% this is a demo low-performance EPI sequence
% which doesn"t use ramp-samping. It is only good for educational purposes

seq=mr.Sequence();              % Create a new sequence object
fov=230e-3; Nx=101; Ny=100;     % Define FOV and resolution
thickness=  2e-2; %1/3*1e-2;           % slice thinckness
slice_gap = 2/3*1e-2;           % gap between slices
Nslices=1;

% Set system limits
lims = mr.opts('MaxGrad',32,'GradUnit','mT/m',...
               'MaxSlew',130,'SlewUnit','T/m/s', ...
               'rfRingdownTime', 30e-6, 'rfDeadTime', 100e-6);


% Create 90 degree slice selection pulse and gradient
[rf, gz] = mr.makeSincPulse(pi/2,'system',lims,'Duration',3e-3,...
    'SliceThickness',thickness,'apodization',0.5,'timeBwProduct',4);

% Define other gradients and ADC events
deltak=1/fov;
kWidth = (Nx-1)*deltak;
dwellTime = 4e-6; % I want it to be divisible by 2
readoutTime = (Nx-1)*dwellTime;
flatTime=ceil(readoutTime*1e5)*1e-5; % round-up to the gradient raster
gx = mr.makeTrapezoid('x',lims,'Amplitude',kWidth/readoutTime,'FlatTime',flatTime);
adc = mr.makeAdc(Nx,'Duration',readoutTime,'Delay',gx.riseTime+flatTime/2-(readoutTime)/2);

% Pre-phasing gradients
preTime=8e-4;
gxPre = mr.makeTrapezoid('x',lims,'Area',-gx.area/2,'Duration',preTime); % removed -deltak/2 to aligh the echo between the samples
gzReph = mr.makeTrapezoid('z',lims,'Area',-gz.area/2,'Duration',preTime);
gyPre = mr.makeTrapezoid('y',lims,'Area',-Ny/2*deltak,'Duration',preTime);
gxReph = gxPre;
gxReph.amplitude = -gxReph.amplitude;
gyReph = gyPre;
gyReph.amplitude = gyReph.amplitude * 49/50;
% Phase blip in shortest possible time
dur = ceil(2*sqrt(deltak/lims.maxSlew)/10e-6)*10e-6;
gy = mr.makeTrapezoid('y',lims,'Area',deltak,'Duration',dur);

% Define sequence blocks
for s=1:Nslices
    rf.freqOffset=gz.amplitude*slice_gap*(s-1-(Nslices-1)/2);
    seq.addBlock(rf,gz);
    seq.addBlock(gxPre,gyPre,gzReph);
    seq.addBlock(mr.makeDelay(2e-3));
    for i=1:Ny
        seq.addBlock(gx,adc);           % Read one line of k-space
        if (i~=Ny) 
            seq.addBlock(gy);           % Phase blip
        end
        gx.amplitude = -gx.amplitude;   % Reverse polarity of read gradient
    end
    seq.addBlock(gxReph,gyReph);
    seq.addBlock(mr.makeDelay(20e-3));
end

seq.definitions("Name") = 'epi'; % FOR RECON

close all
seq.plot();             % Plot sequence waveforms
seq.write('epi.seq');   % Output sequence for scanner
%%
close all
[ktraj_adc, t_adc, ktraj, t_ktraj, t_excitation, t_refocusing, ~, ~] = seq.calculateKspacePP();

% time_axis=(1:(size(ktraj,2)))*seq.sys.gradRasterTime;
% 
% figure; plot(t_ktraj*1e3, ktraj'); % plot the entire k-space trajectory
% hold; plot(t_adc*1e3,ktraj_adc(1,:),'.'); % and sampling points on the kx-axis

figure; plot(ktraj(1,:),ktraj(2,:),'b',...
             ktraj_adc(1,:),ktraj_adc(2,:),'r.'); % a 2D plot
axis('equal'); % enforce aspect ratio for the correct trajectory display