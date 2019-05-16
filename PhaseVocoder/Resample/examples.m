%%Load the .NET library into Matlab
folder = pwd;
path=[folder '\Resample.dll'];
asm = NET.addAssembly(path);
%%
p = 11;
q = 4;
%Create a resampler object with default settings
r = Resample.Resampler(p, q);
%Create an example signal: Sine at 1/10 of Nyquist
x=sin((0:99)*2*pi*0.05);
%Resample by p/q
y= r.Resample(x);
%The result is wrapped in some .NET stuff, convert to double
y = double(y);
%Plot original vs resampled
hold off
stem(0:(length(x) - 1), x);
hold on
plot((0:(length(y) - 1))*q/p, y);
%%
p = 22; %P and Q can have a common denominator
q = 4;
delay = 20; % Filter length used internally in the resampler will be delay*2*max(p, q)+1
%Create a resampler object with  n = delay
r = Resample.Resampler(p, q, uint32(delay));%need to convert to uint32 for matlab to understand that its not a double
blocksize = 41;
nblocks = 5;
%Create an example signal: Cosine at 0.22 of Nyquist
x=cos((0:(nblocks*blocksize))*2*pi*0.11);
y = [];
r.ResetFilter(); %Not necessary the first time
%Resample by p/q blockwise
for i = 1:nblocks
    t = r.ResampleContinuous(x(((i-1)*blocksize +1):(i*blocksize)));
    y= [y double(t)];
end
%Plot original vs resampled
hold off
stem(0:(length(x) - 1), x);
hold on
%%Compensate for the delay in the plot
plot(((0:(length(y) - 1))*q/p)-delay, y);
%%
%Downsample - should work even if q and p have a GCD different from 1
q = 10;
p = 4;
r = Resample.Resampler(p, q);%need to convert to uint32 for matlab to understand that its not a double
blocksize = 77;
nblocks = 8;
%Create an example signal: Cosine at 0.11 of Nyquist
x=cos((0:(nblocks*blocksize))*2*pi*0.035);
y = [];
r.ResetFilter(); %Not necessary the first time
%Resample by p/q blockwise
for i = 1:nblocks
    t = r.ResampleContinuous(x(((i-1)*blocksize +1):(i*blocksize)));
    y= [y double(t)];
end
%Plot original vs resampled
hold off
plot(0:(length(x) - 1), x);
hold on
%%Compensate for the default delay - 10 - in the plot
stem((((0:(length(y) - 1))-10)*q/p), y);