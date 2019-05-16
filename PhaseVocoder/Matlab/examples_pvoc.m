%%Load the .NET library into Matlab
folder = pwd;
path=[folder '\..\PhaseVocoder\bin\debug\PhaseVocoder.dll'];
asm = NET.addAssembly(path);

% Create a test signal: sinusoidal, 8192 samples long, frequency 440 Hz, sampling rate 8192, = duration 1 s
t = 0:1/8192:1-1/8192;
x = sin(2*pi*440*t);

% TEST1. Test 2x tempo change on a 1s 440 Hz sine wave = 0.5s 440 Hz sine wave 
r = PhaseVocoder.PhaseVocoder(2, 1, 1, 1, 128);
y = r.apply(x);
subplot(2,2,1), plot(double(y));

% TEST2. Test 0.5x tempo change on a 1s 440 Hz sine wave = 2s 440 Hz sine wave
r = PhaseVocoder.PhaseVocoder(1, 2, 1, 1, 128);
y = r.apply(x);
subplot(2,2,2), plot(double(y));

% TEST3. Test 0.5x pitch shift on a 1s 440 Hz sine wave = 1s 220 Hz sine wave
r = PhaseVocoder.PhaseVocoder(1, 1, 1, 2, 128);
y = r.apply(x);
subplot(2,2,3), plot(double(y));

% TEST4. Test 2x pitch shift on a 1s 440 Hz sine wave = 1s 880 Hz sine wave
r = PhaseVocoder.PhaseVocoder(1, 1, 2, 1, 128);
y = r.apply(x);
subplot(2,2,4), plot(double(y));

% Continuous tests
% TEST5. Test 2x tempo change on a 10s 440Hz sine wave (block by block) = 5s 440 Hz sine wave in continuous blocks */
% Working here, changing block size in continuous mode does not work ...
r = PhaseVocoder.PhaseVocoder(2, 1, 1, 1, 128);
xb = zeros(1024,1);
y = [];
for b = 1:8;
    for n = 1:1024;
        xb(n) = x((b-1) * 1024 + n);
    end;
    
    y = [y double(r.applyContinuous(xb))];

end;

y = [y double(r.applyContinuous([])];


% TEST6. Test 0.5x tempo change on a 1s 440 Hz sine wave = 2s 440 Hz sine wave in continuous blocks */
r = PhaseVocoder.PhaseVocoder(1, 2, 1, 1, 128);
xb = zeros(1024,1);
y = [];
for b = 1:8;
    for n = 1:1024;
        xb(n) = x((b-1) * 1024 + n);
    end;
    
    y = [y double(r.applyContinuous(xb))];   
end;

y = [y double(r.applyContinuous([]))];

% TEST7. Test 2x pitch shift on a 1s 440 Hz sine wave = 1s 880 Hz sine wave in continuous blocks */
r = PhaseVocoder.PhaseVocoder(1, 1, 2, 1, 128);
xb = zeros(1024,1);
y = [];
for b = 1:8;
    for n = 1:1024;
        xb(n) = x((b-1) * 1024 + n);
    end;
    
    y = [y double(r.applyContinuous(xb))];   
end;

y = [y double(r.applyContinuous([]))];

% TEST8. Test 0.5x pitch shift on a 1s 440 Hz sine wave = 1s 220 Hz sine wave in continuous blocks */
r = PhaseVocoder.PhaseVocoder(1, 1, 1, 2, 128);
xb = zeros(1024,1);
y = [];
for b = 1:8;
    for n = 1:1024;
        xb(n) = x((b-1) * 1024 + n);
    end;
    
    y = [y double(r.applyContinuous(xb))];   
end;

y = [y double(r.applyContinuous([]))];

