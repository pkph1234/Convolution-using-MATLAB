%% Project1.2 Hemal Dave
clc
close all
clear all
disp 'selection for convolution'
disp 'Press 1:F(t)=Sinosudial Function & h(t)= Step Function '
disp 'Press 2:F(t)=Step Function & h(t)= Step Function'
disp 'Press 3:F(t)=Exponential Function & H(t)= Step Function'
disp 'Press 4:F(t)=Ramp Function & H(t)=Step Function'
disp 'Press 5:F(t)=Step Function & H(t)=Sinosudial Function'
disp 'Press 6:F(t)=Sinosudial Function & H(t)=Sinosudial Function'
disp 'Press 7:F(t)=Exponential Function & H(t)=Sinosudial Function'
disp 'Press 8:F(t)=Ramp Function & H(t)=Sinosudial Function'
disp 'Press 9:F(t)=Step Function & H(t)=Exponential Function'
disp 'Press 10:F(t)=Sinosudial Function & H(t)=Exponential Function'
disp 'Press 11:F(t)=Exponential Function & H(t)=Exponential Function'
disp 'Press 12:F(t)=Ramp Function & H(t)=Exponential Function'
disp 'Press 13:F(t)=Step Function & H(t)=Ramp Function'
disp 'Press 14:F(t)=Sinosudial Function & H(t)=Ramp Function'
disp 'Press 15:F(t)=Ramp Function & H(t)=Ramp Function '
disp 'Press 16:F(t)=Exponential Function & H(t)=Ramp Function'

n=input('Enter the number');
switch(n) 
    case 1
        clc;
close all;
clear all;
%% Step Function h(t):
disp('This is Step function:h(t)')
min='what is your min range?';
% Min=input(min);
Min=-10;
max='What is your max range?';
Max=input(max);
% Max=10;
sample='how many sample you want?';
Sample=input(sample);
% Sample=500;
a=Sample/(Max-Min);
delay='How much time shift do you need?';
Delay=input(delay);
% Delay=-5;
def=Delay-Min;
Stepwidth=input('How much is the step width?');
% Stepwidth=3;
n=Min:(Max-Min)/Sample:Max;% Scalling of the x axis according to sample
Amp=input('What is your Amplitude?');
% Amp=5;
w=(0-Min)*Sample;
%% Sin function f(t)
disp('This is sin function f(t)');
Minrange=input('What is your minimum range of the plot?');
% Minrange=-10;
Maxrange=input('What is your maximum range?');
% Maxrange=10;
sample=('How much samples do you need? ');
Sample=input(sample);
% Sample=500;
s=Sample/(Maxrange-Minrange);
t=Minrange:1/s:Maxrange;
frequency=('What is your frequency?');
Frequency=input(frequency);
% Frequency=1;
% frequency=2;
Amp=('what is your Amplitude?');
amp=input(Amp);
% amp=5;
% amp=4;
phase=('What is your phase shift? (in rad)');
Phase=input(phase);
% Phase=0;
delay=input('What is your delay? ');
% delay=0;
def=Delay-Minrange;
Def=delay-Min;
% close
Delay=Delay-(def-Def);
array_part=0;
 %% Intilization of variables for convolution
L1 = length(0:(1/a):Max);
L2 = length(1:1:length(n));
matrix = zeros(L1,L2);
yfinal=zeros(L1,L2);
zfinal=zeros(L1,L2);
l=zeros(L1,L2);
it=0;i=1;
yy = 0;
qw=100;
% asd=1:1:100;
as=zeros(L1,L2);
input=input('Choose from cos and sin function ? 1 for sin & 0 for cos.');
er=input;
%% Calculation of F(t) & h(T-t)
for v=0:(1/a):(Max+500)
   def=Delay-(v)-Minrange;
   Def=delay-Min;
   
   
   for x=1:1:length(n)

    if x>=def*a && x<def*a+(Stepwidth)*a% for selecting the Rabge where Step is heppening according to the step width
        
        y(x)=Amp;% take the value pf the amplitude till the step width
        
    else 
        
        y(x)=0;
        
    end
    if er==1
    
         if x>=Def*a && x<Def*a+(Stepwidth)*a
              z(x)=Amp*sin(2*pi*Frequency*x*(1/s)+Phase);
         else
              z(x)=0;
         end
         if er==0
             if x>=Def*a && x<Def*a+(Stepwidth)*a
              z(x)=Amp*cos(2*pi*Frequency*x*(1/s)+Phase);
         else
              z(x)=0;
             end
        end
    end
   end
L1 = length(0:(1/a):Max);
L2 = length(1:1:length(n));
% close
M=Min:(Max-Min)/Sample:Stepwidth;

% title(['Sin Function: ','Minimum range = ', num2str(Minrange),', Maximum range = ', num2str(Maxrange),', Sample = ',num2str(Sample),', Frequecy =', num2str(Frequency),', Phase Shift = ',num2str(Phase),', Amplitude = ',num2str(Amp),' .']);


% subplot(2,1,1)
%% Calculation for convolution:
% figure
% eval(sprintf('%v = y(x)', v))
% save('file.mat', sprintf('%v',v)) 
it = it +1;
cv=fliplr(y);
Yfinal(it,:)=cv; 
Zfinal(it,:)=z;
% g=Yfinal(it,it).*Zfinal(it,it)
% l(1:it,:)=Yfinal.*Zfinal; 
% v=Yfinal(i,round(def*a)+1).*Zfinal(i,round(Def*a)+1);
% for i=1:1:qw
% v=Yfinal(i,round(def*a)+1).*Zfinal(i,round(Def*a)+1)
v=Yfinal(it,:).*Zfinal(it,:);
as(it,:)=sum(v)+as(it,:);
% end
%% Plotting the convolution:
subplot(2,2,1)
  plot(as,'linewidth',3)
  xlabel('time');
  ylabel('Convolution function');
  title('Convolution');
  subplot(2,2,2)
  plot(n,cv,'linewidth',3)% plot the graph
  xlabel('time');
  ylabel('function h(T-t)');
  title('h(T-t)');
  pause(0.1)
% hold on
% plot(-n,y,'linewidth',3)
  subplot(2,2,3)
  plot(n,z,'r')
  xlabel('time');
  ylabel('function f(t)');
  title('f(t)');
  pause(0.1);
  subplot(2,2,4)
  plot(n,cv,'linewidth',3)% plot the graph
  hold on
  plot(n,z,'r')
  hold off
  xlabel('time');
  ylabel('Convolution function');
  title('f(t)*h(T-t)');
  pause(0.1)
% asdf=1:1:100;
% subplot(3,2,4)
% hold on
% plot(as)
% pause(0.1);
% hold off
% xlabel('Time (in sec)');
% ylabel('Amplitude')
% hold off
% pause(0.1)
% o=conv(y,z)
% close
% figure()
% plot(t,l)
% pause(1);
% clear(sprintf('%v',v))
% grid on;
% xlabel('time');
% ylabel('Amplitude');
  p=repmat(y,L1,L2);
end
% load('file.mat');
% chu=y();
% eval(sprintf('array = %v', v));
% v=1:1:10;
% bn=sum(Yfinal(v,:).*Zfinal(1,:))
    
    case 2
        clc;
close all;
clear all;
%% Step function H(t)
disp('Convolution of Step function with Step Function')
disp('');
disp('This is Step function:h(t)')
min='what is your min range?';
Min=input(min);
Min=-10;
max='What is your max range?';
Max=input(max);
% Max=10;
sample='how many sample you want?';
Sample=input(sample);
% Sample=1000;
a1=Sample/(Max-Min);
Delay='How much time shift do you need?';
delay=input(Delay);
% delay=0;
def=delay-Min;
Stepwidth1=input('How much is the step width?');
% Stepwidth1=5;
n=Min:(Max-Min)/Sample:Max;% Scalling of the x axis according to sample
% Amp1=10;
Amp1=input('What is your Amplitude?');
%% Step function F(t):
disp('This is Step function:f(t)')
min='what is your min range?';
Min=input(min);
% Min=-10;
max='What is your max range?';
Max=input(max);
% Max=10;
sample='how many sample you want?';
Sample=input(sample);
% Sample=1000;
a=Sample/(Max-Min);
Delay='How much time shift do you need?';
delay2=input(Delay);
% delay2=0;
Def=delay2-Min;
def=def-(def-Def);
Stepwidth2=input('How much is the step width?');
% Stepwidth2=5;
n=Min:(Max-Min)/Sample:Max;% Scalling of the x axis according to sample
% Amp2=10;
Amp2=input('What is your amplitue?');
%% Intilization of variables for convolution
array_part=0;
L1 = length(0:(1/a):Max);
L2 = length(1:1:length(n));
matrix = zeros(L1,L2);
yfinal=zeros(L1,L2);
zfinal=zeros(L1,L2);
l=zeros(L1,L2);
it=0;i=1;
yy = 0;
qw=100;
% asd=1:1:100;
as=zeros(L1,L2);
%% Calculation of F(t) & h(T-t)
for v=0:(1/a):(Max+500)
   def=delay-(v)-Min;
   Def=delay2-Min;
   for x=1:1:length(n)

    if x>=def*a && x<def*a+(Stepwidth1)*a% for selecting the Rabge where Step is heppening according to the step width
        
        y(x)=Amp1;% take the value pf the amplitude till the step width
        
    else 
        
        y(x)=0;
        
    end
   end
    for x=1:1:length(n)

    if x>=Def*a && x<Def*a+(Stepwidth2)*a% for selecting the Rabge where Step is heppening according to the step width
        
        z(x)=Amp2;% take the value pf the amplitude till the step width
        
    else 
        
        z(x)=0;
        
    end
    end
    it = it +1;
    cv=fliplr(y);
    Yfinal(it,:)=cv; 
    Zfinal(it,:)=z;
    v=Yfinal(it,:).*Zfinal(it,:);
   as(it,:)=sum(v)+as(it,:);
  % end
  ed=Min:Sample:Max;
   subplot(2,2,1)
  plot(as,'linewidth',3)
  xlabel('time');
  ylabel('Convolution function');
  title('Convolution');
  subplot(2,2,2)
  plot(n,cv,'linewidth',3)% plot the graph
  xlabel('time');
  ylabel('function h(T-t)');
  title('h(T-t)');
  pause(0.1)
% hold on
% plot(-n,y,'linewidth',3)
  subplot(2,2,3)
  plot(n,z,'r')
  xlabel('time');
  ylabel('function f(t)');
  title('f(t)');
  pause(0.1);
  subplot(2,2,4)
  plot(n,cv,'linewidth',3)% plot the graph
  hold on
  plot(n,z,'r')
  hold off
  xlabel('time');
  ylabel('Convolution function');
  title('f(t)*h(T-t)');
  pause(0.1)
% asdf=1:1:100;
% subplot(3,2,4)
% hold on
% plot(as)
% pause(0.1);
% hold off
% xlabel('Time (in sec)');
% ylabel('Amplitude')
% hold off
% pause(0.1)
% o=conv(y,z)
% close
% figure()
% plot(t,l)
% pause(1);
% clear(sprintf('%v',v))
% grid on;
% xlabel('time');
% ylabel('Amplitude');
  p=repmat(y,L1,L2);
     end
    case 3 
        clc;
close all;
clear all;
%% Step Function H(t)
disp('This is Step function:h(t)')
min='what is your min range?';
Min=input(min);
% Min=-5;
max='What is your max range?';
Max=input(max);
% Max=5;
sample='how many sample you want?';
Sample=input(sample);
% Sample=1000;
a=Sample/(Max-Min);
delay='How much time shift do you need?';
Delay=input(delay);
% Delay=0;
def=Delay-Min;
Stepwidth=input('How much is the step width?');
% Stepwidth=5;
n=Min:(Max-Min)/Sample:Max;% Scalling of the x axis according to sample
% Amp1=5;
Amp1=input('What is your amplitude?');
%% Exponential Function f(t)
disp('This is exponential function f(t)')
f=0;
% tstart=1;
% tstop=10;
min='What is the start point of the graph? Enter your Start point?';
Min=input(minrange);
% Min=-5;
max='what is the end point of the graph?';
Max=input(maxrange);
% Max=5;
Sample=('what is your sample size?');
sample=input(Sample);
% sample=1000;
amp='What is your amplitude?';
Amp2=input(amp);
% Amp2=1;
% amp=1;
decay='What is your decay rate?';
Decay=input(decay);
% Decay=1;
delay=Input('What is your delay?');
% delay=0;
stop=input('From where your graph will stop?');
% stop=1;
Def=delay-Min;
%% initialization of variables for convolution:
Delayx=Delay-(def-Def);
% delay=input('What is your delay?');
% decay=1;
% include time shiftsamp
% sample=input('Enter your sample size?');
a=sample/(Max-Min);
% t=Min:a:Max;
array_part=0;
L1 = length(0:(1/a):Max);
L2 = length(1:1:length(n));
matrix = zeros(L1,L2);
yfinal=zeros(L1,L2);
zfinal=zeros(L1,L2);
l=zeros(L1,L2);
it=0;i=1;
yy = 0;
qw=100;
% asd=1:1:100;
as=zeros(L1,L2);
%% plotting of 2 signals:
for v=0:(1/a):(Max+500)
   def=Delay-(v)-Min;
   Def=delay-Min;
   for x=1:1:length(n)

    if x>=def*a && x<def*a+(Stepwidth)*a% for selecting the Rabge where Step is heppening according to the step width
        
        y(x)=Amp1;% take the value pf the amplitude till the step width
        
    else 
        
        y(x)=0;
    end
   end
   for x=1:1:length(n)
   if x>Def*a && x<=Def*a+(stop)*a
       z(x)=Amp2*exp(Decay*x/a);
   else
    z(x)=0;
   end
   end
   %% Calculation of Convolution:
   it = it +1;
cv=fliplr(y);
Yfinal(it,:)=cv; 
Zfinal(it,:)=z;
v=Yfinal(it,:).*Zfinal(it,:);
as(it,:)=sum(v)+as(it,:);
% end
%% plotting of final Convolution:
subplot(2,2,1)
  plot(as,'linewidth',3)
  xlabel('time');
  ylabel('Convolution function');
  title('Convolution');
  subplot(2,2,2)
  plot(n,cv,'linewidth',3)% plot the graph
  xlabel('time');
  ylabel('function h(T-t)');
  title('h(T-t)');
  pause(0.1)
% hold on
% plot(-n,y,'linewidth',3)
  subplot(2,2,3)
  plot(n,z,'r','linewidth',3)
  xlabel('time');
  ylabel('function f(t)');
  title('f(t)');
  pause(0.1);
  subplot(2,2,4)
  plot(n,cv,'linewidth',3)% plot the graph
  hold on
  plot(n,z,'r')
  hold off
  xlabel('time');
  ylabel('Convolution function');
  title('f(t)*h(T-t)');
  pause(0.1)
% asdf=1:1:100;
% subplot(3,2,4)
% hold on
% plot(as)
% pause(0.1);
% hold off
% xlabel('Time (in sec)');
% ylabel('Amplitude')
% hold off
% pause(0.1)
% o=conv(y,z)
% close
% figure()
% plot(t,l)
% pause(1);
% clear(sprintf('%v',v))
% grid on;
% xlabel('time');
% ylabel('Amplitude');
  p=repmat(y,L1,L2);
end
    case 4
        clc;
close all;
clear all;
%% Step function H(t)
disp('Convolution of Step function with Step Function')
disp('');
disp('This is Step function:h(t)')
min='what is your min range?';
Min=input(min);
% Min=-10;
max='What is your max range?';
Max=input(max);
% Max=10;
sample='how many sample you want?';
Sample=input(sample);
% Sample=100;
a1=Sample/(Max-Min);
Delay='How much time shift do you need?';
delay=input(Delay);
% delay=0;
% def=delay-Min;
Stepwidth1=input('How much is the step width?');
% Stepwidth1=5;
n=Min:(Max-Min)/Sample:Max;% Scalling of the x axis according to sample
% Amp1=10;
Amp1=input('What is your Amplitude?');
%% Ramp Function:
%% Ramp function:
Min=input('What is your Minimum value of your plot?')
% Min=-10;
Max=input('What is your Maximum value of your plot?')
% Max=10;
% def=Max-Min;
sample=input('How mnay samples do you need?');
% sample=100;
a=sample/(Max-Min);
n=Min:1/a:Max;
start=0;
End=input('From where your ramp should end?')
% End=5;
Delay=input('What is your Dealy?');
% Delay=0;
decay=input('What is ypur decay rate? Should be 1 or -1');
% decay=1;
Def=0-Min;
% Delay=-1;
%% Intilization of variables for convolution
array_part=0;
L1 = length(0:(1/a):Max);
L2 = length(1:1:length(n));
matrix = zeros(L1,L2);
yfinal=zeros(L1,L2);
zfinal=zeros(L1,L2);
l=zeros(L1,L2);
it=0;i=1;
yy = 0;
qw=100;
% asd=1:1:100;
as=zeros(L1,L2);
%% Calculation of F(t) & h(T-t)
for v=0:(1/a):(Max+500)
   def=Delay-(v)-Min;
   Def=0-Min;
   for x=1:1:length(n)

    if x>=def*a && x<def*a+(Stepwidth1)*a% for selecting the Rabge where Step is heppening according to the step width
        
        y(x)=Amp1;% take the value pf the amplitude till the step width
        
    else 
        
        y(x)=0;
        
    end
    if x>(Def+delay)*a && x<=(Def+End+delay)*a
        z(x)=decay*n(x-(delay*a));
    else
       z(x)=0;
    end
   end
       it = it +1;
 
   
    cv=fliplr(y);
    Yfinal(it,:)=cv; 
    Zfinal(it,:)=z;
    v=Yfinal(it,:).*Zfinal(it,:);
     as(it,:)=sum(v)+as(it,:);
  % end
%   ed=Min:Sample:Max;
   subplot(2,2,1)
  plot(as,'linewidth',3)
  xlabel('time');
  ylabel('Convolution function');
  title('Convolution');
  subplot(2,2,2)
  plot(n,cv,'linewidth',3)% plot the graph
  xlabel('time');
  ylabel('function h(T-t)');
  title('h(T-t)');
  pause(0.1)
% hold on
% plot(-n,y,'linewidth',3)
  subplot(2,2,3)
  plot(n,z,'r')
  xlabel('time');
  ylabel('function f(t)');
  title('f(t)');
  pause(0.1);
  subplot(2,2,4)
  plot(n,cv,'linewidth',3)% plot the graph
  hold on
  plot(n,z,'r')
  hold off
  xlabel('time');
  ylabel('Convolution function');
  title('f(t)*h(T-t)');
  pause(0.1)
% asdf=1:1:100;
% subplot(3,2,4)
% hold on
% plot(as)
% pause(0.1);
% hold off
% xlabel('Time (in sec)');
% ylabel('Amplitude')
% hold off
% pause(0.1)
% o=conv(y,z)
% close
% figure()
% plot(t,l)
% pause(1);
% clear(sprintf('%v',v))
% grid on;
% xlabel('time');
% ylabel('Amplitude');
%   p=repmat(y,L1,L2);
    end
    case 5
        clc
close all
clear all
%% sin function: h(t)
disp('This is Sionosudial function h(t)')
% Min=-5;
Min=input('What is your minimum range?');
Max=input(maxrange);
% Max=5;
Sample=('How much samples do you need? ');
sample=input(Sample);
% sample=1000;
s=sample/(Max-Min);
t=Min:1/s:Max;
frequency=('What is your frequency?');
Frequency=input(frequency);
% Frequency=1;
% frequency=2;
amp=('what is your Amplitude?');
Amp=input(amp);
% Amp=5;
% amp=4;
phase=('What is your phase shift? (in rad)');
Phase=input(phase);
% Phase=0;
Delay=('What is your delay?');
% Delay=0;
Stepwidth=('From where your sinosudial function will stop plotting sinosuidial function?');
% Stepwidth=3;
def=Delay-Min;

% Def=delay-Min;
%% Step function f(t)
disp('This is Step function: f(t)')
min='what is your min range?';
Min=input(min);
% Min=-5;
max='What is your max range?';
Max=input(max);
% Max=5;
Sample='how many sample you want?';
sample=input(Sample);
% sample=1000;
a=sample/(Max-Min);
Delay='How much time shift do you need?';
delay=input(Delay);
% delay=0;
% Def=delay-Min;
Stepwidth1=input('How much is the step width?');
% Stepwidth1=1;
n=Min:(Max-Min)/sample:Max;% Scalling of the x axis according to sample
Amp1=input('What is your amplitude?');
% Amp1=5;
Def=delay-Min;
%% variables 
Delay=Delay-(def-Def);
n=Min:(Max-Min)/sample:Max;
% Delay=Delay-(def-Def);
% delay=input('What is your delay?');
% decay=1;
% include time shiftsamp
% sample=input('Enter your sample size?');
a=sample/(Max-Min);
% t=Min:a:Max;
% array_part=0;
L1 = length(0:(1/a):Max);
L2 = length(1:1:length(n));
matrix = zeros(L1,L2);
yfinal=zeros(L1,L2);
zfinal=zeros(L1,L2);
l=zeros(L1,L2);
it=0;i=1; 
yy = 0;
qw=100;
as=zeros(L1,L2);
input=input('For Sin function press 1 for cos function press 0');
hem=input;
%% plot the both function:
for v=0:(1/a):(Max+500)
   def=(Delay - (v)-Min);
   Def=delay-Min;
   input=0;
   for x=1:1:length(n)
       if hem==1
         if x>=(def*a) && x<(def*a+(Stepwidth)*a)
              y(x)=Amp1*sin(2*pi*Frequency*(x+v*a)*(1/s)+Phase);
         else
              y(x)=0;
         end
       end
    if hem==0
        if x>=(def*a) && x<(def*a+(Stepwidth)*a)
              y(x)=Amp1*cos(2*pi*Frequency*(x+v*a)*(1/s)+Phase);
         else
              y(x)=0;
         end
       end
   
    
    if x>=Def*a && x<=Def*a+(Stepwidth1)*a% for selecting the Rabge where Step is heppening according to the step width
        
        z(x)=Amp;% take the value pf the amplitude till the step width
 
    else 
        
        z(x)=0;
    end  
   end
        %% convolution
        it = it +1;
     cv=fliplr(y);
     Yfinal(it,:)=cv; 
    Zfinal(it,:)=z;
     v=Yfinal(it,:).*Zfinal(it,:);
     as(it,:)=sum(v)+as(it,:);
    %% plotting ther graph
    subplot(2,2,1)
  plot(as,'linewidth',3)
  xlabel('time');
  ylabel('Convolution function');
  title('Convolution');
  subplot(2,2,2)
  plot(n,cv,'linewidth',3)% plot the graph
  xlabel('time');
  ylabel('function h(T-t)');
  title('h(T-t)');
  pause(0.1)
% hold on
% plot(-n,y,'linewidth',3)
  subplot(2,2,3)
  plot(n,z,'r','linewidth',3)
  xlabel('time');
  ylabel('function f(t)');
  title('f(t)');
  pause(0.1);
  subplot(2,2,4)
  plot(n,cv,'linewidth',3)% plot the graph
  hold on
  plot(n,z,'r')
  hold off
  xlabel('time');
  ylabel('Convolution function');
  title('f(t)*h(T-t)');
  pause(0.1)
% asdf=1:1:100;
% subplot(3,2,4)
% hold on
% plot(as)
% pause(0.1);
% hold off
% xlabel('Time (in sec)');
% ylabel('Amplitude')
% hold off
% pause(0.1)
% o=conv(y,z)
% close
% figure()
% plot(t,l)
% pause(1);
% clear(sprintf('%v',v))
% grid on;
% xlabel('time');
% ylabel('Amplitude');
  p=repmat(y,L1,L2);
 end
    case 6
        clc
clear all
close all
% clear all
%% sin function:f(t)
% Min=-5;
disp ('This is sionosudial function h(t)')
Min=input('What is upur Minimum Range?');
max=('What is your Max range?');
Max=input(max);
% Max=5;
Sample=('How much samples do you need? ');
sample=input(Sample);
% sample=1000;
s=sample/(Max-Min);
t=Min:1/s:Max;
frequency=('What is your frequency?');
Frequency1=input(frequency);
% Frequency1=1;
% frequency=2;
amp=('what is your Amplitude?');
Amp1=input(amp);
% Amp1=5;
% amp=4;
phase=('What is your phase shift? (in rad)');
Phase1=input(phase);
% Phase1=0;
% Delay=0;
Delay=input('What is your delay?');
Stepwidth1=input('From where your graph should stop?');
% Stepwidth1=3;
def=Delay-Min;
%% Sin Function f(t):
% Min=-5;
disp ('This is sinosudial function f(t)')
Min=input('What is upur Minimum Range?');
max=('What is your Max range?');
Max=input(max);

% Max=5;
Sample=('How much samples do you need? ');
sample=input(Sample);
% sample=1000;
s=sample/(Max-Min);
t=Min:1/s:Max;
frequency=('What is your frequency?');
Frequency2=input(frequency);
% Frequency2=1;
% frequency=2;
amp=('what is your Amplitude?');
Amp2=input(amp);
% Amp2=5;
% amp=4;
phase=('What is your phase shift? (in rad)');
Phase2=input(phase);
% Phase2=0;
% delay=0;
delay=input('What is your delay?');
% Stepwidth2=3;
Stepwidth2=input('From where your graph should stop?');
Def=delay-Min;
%% initialization of variables:
Delay=Delay-(def-Def);
n=Min:(Max-Min)/sample:Max;
% Delay=Delay-(def-Def);
% delay=input('What is your delay?');
% decay=1;
% include time shiftsamp
% sample=input('Enter your sample size?');
a=sample/(Max-Min);
% t=Min:a:Max;
% array_part=0;
L1 = length(0:(1/a):Max);
L2 = length(1:1:length(n));
matrix = zeros(L1,L2);
yfinal=zeros(L1,L2);
zfinal=zeros(L1,L2);
l=zeros(L1,L2);
it=0;i=1; 
yy = 0;
qw=100;
input=input('Press 1 for sin function & press 0 for cosine function:');
input1=input('Press 1 fr sin function & press 0 for cosine function');

% asd=1:1:100;
as=zeros(L1,L2);
%% plot both functions:
for v=0:(1/a):(Max+500)
   def=(Delay - (v)-Min);
   Def=delay-Min;
%    input=0;
   for x=1:1:length(n)
       if input==1
         if x>=(def*a) && x<(def*a+(Stepwidth1)*a)
              y(x)=Amp1*sin(2*pi*Frequency1*(x+(v*a))*(1/s)+Phase1);
         else
              y(x)=0;
         end
       end
    if input==0
        if x>=(def*a) && x<(def*a+(Stepwidth1)*a)
              y(x)=Amp1*cos(2*pi*Frequency1*(x+(v*a))*(1/s)+Phase1);
         else
              y(x)=0;
         end
    end
       
%    input1=1;
 
       if input1==1
         if x>=(Def*a) && x<(Def*a+(Stepwidth2)*a)
              z(x)=Amp2*sin(2*pi*Frequency2*(x)*(1/s)+Phase2);
         else
              z(x)=0;
         end
       end
    if input1==0
        if x>=(Def*a) && x<(Def*a+(Stepwidth2)*a)
              z(x)=Amp2*cos(2*pi*Frequency2*(x)*(1/s)+Phase2);
         else
              z(x)=0;
         end
    end
   end
    %% calculating convolution:
    it = it +1;
     cv=fliplr(y);
     Yfinal(it,:)=cv; 
    Zfinal(it,:)=z;
     v=Yfinal(it,:).*Zfinal(it,:);
     as(it,:)=sum(v)+as(it,:);
     %% plotting convolution
     subplot(2,2,1)
  plot(as,'linewidth',3)
  xlabel('time');
  ylabel('Convolution function');
  title('Convolution');
  subplot(2,2,2)
  plot(n,cv,'linewidth',3)% plot the graph
  xlabel('time');
  ylabel('function h(T-t)');
  title('h(T-t)');
  pause(0.1)
% hold on
% plot(-n,y,'linewidth',3)
  subplot(2,2,3)
  plot(n,z,'r','linewidth',3)
  xlabel('time');
  ylabel('function f(t)');
  title('f(t)');
  pause(0.1);
  subplot(2,2,4)
  plot(n,cv,'linewidth',3)% plot the graph
  hold on
  plot(n,z,'r')
  hold off
  xlabel('time');
  ylabel('Convolution function');
  title('f(t)*h(T-t)');
  pause(0.1)
% asdf=1:1:100;
% subplot(3,2,4)
% hold on
% plot(as)
% pause(0.1);
% hold off
% xlabel('Time (in sec)');
% ylabel('Amplitude')
% hold off
% pause(0.1)
% o=conv(y,z)
% close
% figure()
% plot(t,l)
% pause(1);
% clear(sprintf('%v',v))
% grid on;
% xlabel('time');
% ylabel('Amplitude');
  p=repmat(y,L1,L2);
end
    case 7
        %% sino function:
clc
clear all
close all
disp('This is sinosudial function h(t)')
Min=input('What is yoyr minimum range?');
% Min=-5;
Max=input('What is your max range?');
% Max=5;
sample=('How much samples do you need? ');
Sample=input(sample);
% Sample=1000;
s=Sample/(Max-Min);
t=Min:1/s:Max;
frequency=('What is your frequency?');
Frequency=input(frequency);
% Frequency=1;
% frequency=2;
amp=('what is your Amplitude?');
Amp1=input(amp);
% Amp1=5;
% amp=4;
phase=('What is your phase shift? (in rad)');
Phase=input(phase);
% Phase=0;
% Delay=0;
Delay=input('Enter your delay');
Stepwidth=input('Freom where your sin fucntion stop?');
% Stepwidth=3;
def=Delay-Min;

% Def=delay-Min;
%% exponential function:
f=0;
disp 'This is exponential Function f(t)'
% tstart=1;
% tstop=10;
minrange='What is the start point of the graph? Enter your Start point?';
Min=input(minrange);
% Min=-5;
maxrange='what is the end point of the graph?';
Max=input(maxrange);
% Max=5;
Sample=('what is your sample size?');
sample=input(Sample);
sampTle=1000;
amp='What is your amplitude?';
Amp2=input(amp);
% Amp2=1;
% amp=1;
decay='What is your decay rate?';
Decay=input(decay);
% Decay=1;
delay=input('Enter your delay');
stop=input('Freom where your fucntion stop?');
% delay=0;
% stop=1;
Def=delay-Min;
%% initilization of variables
Delay=Delay-(def-Def);
n=Min:(Max-Min)/Sample:Max;
% Delay=Delay-(def-Def);
% delay=input('What is your delay?');
% decay=1;
% include time shiftsamp
% sample=input('Enter your sample size?');
a=sample/(Max-Min);
% t=Min:a:Max;
% array_part=0;
L1 = length(0:(1/a):Max);
L2 = length(1:1:length(n));
matrix = zeros(L1,L2);
yfinal=zeros(L1,L2);
zfinal=zeros(L1,L2);
l=zeros(L1,L2);
it=0;i=1; 
yy = 0;
qw=100;
input=('Press 1 for sin function & press 0 for cosine function');
% asd=1:1:100;
as=zeros(L1,L2);
%% plotting graphs:
for v=0:(1/a):(Max+500)
   def=(Delay - (v)-Min);
   Def=delay-Min;
%    input=0;
   for x=1:1:length(n)
       if input==1
         if x>=(def*a) && x<(def*a+(Stepwidth)*a)
              y(x)=Amp1*sin(2*pi*Frequency*(x+v*a)*(1/s)+Phase);
         else
              y(x)=0;
         end
       end
    if input==0
        if x>=(def*a) && x<(def*a+(Stepwidth)*a)
              y(x)=Amp1*cos(2*pi*Frequency*(x+v*a)*(1/s)+Phase);
         else
              y(x)=0;
         end
       end
        
 end
    for x=1:1:length(n)
   if x>Def*a && x<=Def*a+(stop)*a
       z(x)=Amp2*exp(Decay*x/a);
   else
    z(x)=0;
   end
   end
    
    %% Convolution:
       it = it +1;
     cv=fliplr(y);
     Yfinal(it,:)=cv; 
    Zfinal(it,:)=z;
     v=Yfinal(it,:).*Zfinal(it,:);
     as(it,:)=sum(v)+as(it,:);
    %% plotting ther graph
    subplot(2,2,1)
  plot(as,'linewidth',3)
  xlabel('time');
  ylabel('Convolution function');
  title('Convolution');
  subplot(2,2,2)
  plot(n,cv,'linewidth',3)% plot the graph
  xlabel('time');
  ylabel('function h(T-t)');
  title('h(T-t)');
  pause(0.1)
% hold on
% plot(-n,y,'linewidth',3)
  subplot(2,2,3)
  plot(n,z,'r','linewidth',3)
  xlabel('time');
  ylabel('function f(t)');
  title('f(t)');
  pause(0.1);
  subplot(2,2,4)
  plot(n,cv,'linewidth',3)% plot the graph
  hold on
  plot(n,z,'r')
  hold off
  xlabel('time');
  ylabel('Convolution function');
  title('f(t)*h(T-t)');
  pause(0.1)
% asdf=1:1:100;
% subplot(3,2,4)
% hold on
% plot(as)
% pause(0.1);
% hold off
% xlabel('Time (in sec)');
% ylabel('Amplitude')
% hold off
% pause(0.1)
% o=conv(y,z)
% close
% figure()
% plot(t,l)
% pause(1);
% clear(sprintf('%v',v))
% grid on;
% xlabel('time');
% ylabel('Amplitude');
  p=repmat(y,L1,L2);
end
    case 8
        clc
clear all;
close all;
disp('This is Sionosudial function h(t)')
% Min=-5;
Min=input('What is your minimum range?');
Max=input(maxrange);
% Max=5;
Sample=('How much samples do you need? ');
sample=input(Sample);
% sample=1000;
s=sample/(Max-Min);
t=Min:1/s:Max;
frequency=('What is your frequency?');
Frequency=input(frequency);
% Frequency=1;
% frequency=2;
amp=('what is your Amplitude?');
Amp1=input(amp);
% Amp1=5;
% amp=4;
phase=('What is your phase shift? (in rad)');
Phase=input(phase);
% Phase=0;
Delay=('What is your delay?');
% Delay=0;
Stepwidth=('From where your sinosudial function will stop plotting sinosuidial function?');
% Stepwidth=2;
def=Delay-Min;

% Def=delay-Min;
%% Ramp Function:
%% Ramp function:
disp 'This is Ramp Function f(t)'
Min=input('What is your Minimum value of your plot?')
% Min=-5;
Max=input('What is your Maximum value of your plot?')
% Max=5;
% def=Max-Min;
sample=input('How mnay samples do you need?');
% sample=1000;
a=sample/(Max-Min);
n=Min:1/a:Max;
length(n)
start=0;
End=input('From where your ramp should end?');
% End=1;
delay=input('What is your Dealy?');
% delay=0;
decay=input('What is ypur decay rate? Should be 1 or -1');
% decay=1;
% Def=0-Min;
% Delay=-1;
%% Intilization of variables for convolution
array_part=0;
L1 = length(0:(1/a):Max);
L2 = length(1:1:length(n));
matrix = zeros(L1,L2);
yfinal=zeros(L1,L2);
zfinal=zeros(L1,L2);
l=zeros(L1,L2);
it=0;i=1;
yy = 0;
qw=100;
% asd=1:1:100;
as=zeros(L1,L2);
% input=input('For Sin function press 1 for cos function press 0');
input=input('Press 1 for sinfunction & press 0 for cosine function');
hem=input;
%% plot the both function:
for v=0:(1/a):(Max+500)
   def=(Delay-(v)-Min);
   Def=0-Min;
%    input=0;
   for x=1:1:length(n)
       if hem==1
         if x>=(def*a) && x<(def*a+(Stepwidth)*a)
              y(x)=Amp1*sin(2*pi*Frequency*(x+v*a)*(1/s)+Phase);
         else
              y(x)=0;
         end
       end
    if hem==0
        if x>=(def*a) && x<(def*a+(Stepwidth)*a)
              y(x)=Amp1*cos(2*pi*Frequency*(x+v*a)*(1/s)+Phase);
         else
              y(x)=0;
        end
    end
       if x>(Def+delay)*a && x<=(Def+End+delay)*a
           z(x)=decay*n(x-(delay*a));
       else
          z(x)=0;
    end
   end
    %% convolution
        it = it +1;
        cv=fliplr(y);
        yfinal(it,:)=cv; 
       zfinal(it,:)=z;
     v=yfinal(it,:).*zfinal(it,:);
     as(it,:)=sum(v)+as(it,:);
    %% plotting ther graph
    subplot(2,2,1)
  plot(as,'linewidth',3)
  xlabel('time');
  ylabel('Convolution function');
  title('Convolution');
  subplot(2,2,2)
  plot(n,cv,'linewidth',3)% plot the graph
  xlabel('time');
  ylabel('function h(T-t)');
  title('h(T-t)');
  pause(0.1)
% hold on
% plot(-n,y,'linewidth',3)
  subplot(2,2,3)
  plot(n,z,'r','linewidth',3)
  xlabel('time');
  ylabel('function f(t)');
  title('f(t)');
  pause(0.1);
  subplot(2,2,4)
  plot(n,cv,'linewidth',3)% plot the graph
  hold on
  plot(n,z,'r')
  hold off
  xlabel('time');
  ylabel('Convolution function');
  title('f(t)*h(T-t)');
  pause(0.1)
% asdf=1:1:100;
% subplot(3,2,4)
% hold on
% plot(as)
% pause(0.1);
% hold off
% xlabel('Time (in sec)');
% ylabel('Amplitude')
% hold off
% pause(0.1)
% o=conv(y,z)
% close
% figure()
% plot(t,l)
% pause(1);
% clear(sprintf('%v',v))
% grid on;
% xlabel('time');
% ylabel('Amplitude');
%   p=repmat(y,L1,L2);
 end
   
    case 9
clc
clear all
close all
%% exponential function:
f=0;
disp 'This is exponential function h(t)'
% tstart=1;
% tstop=10;
minrange='What is the start point of the graph? Enter your Start point?';
Min=input(minrange);
% Min=-5;
maxrange='what is the end point of the graph?';
Max=input(maxrange);
% Max=5;
Sample=('what is your sample size?');
sample=input(Sample);
% sample=1000;
amp='What is your amplitude?';
Amp=input(amp);
% Amp=1;
% amp=1;
decay='What is your decay rate?';
Decay=input(decay);
% Decay=1;
Delay=input('What is your delay?');
% Delay=0;
stop1=input('From where your function will stop?');
% stop1=1;
def=Delay-Min;
%% Step Function
disp('This is Step function:f(t)')
min='what is your min range?';
Min=input(min);
% Min=-5;
max='What is your max range?';
Max=input(max);
% Max=5;
Sample='how many sample you want?';
sample=input(Sample);
% sample=1000;
a=sample/(Max-Min);
Delay1='How much time shift do you need?';
delay=input(Delay1);
% delay=0;
% Def=delay-Min;
Stepwidth1=input('How much is the step width?');
% Stepwidth1=5;
n=Min:(Max-Min)/sample:Max;% Scalling of the x axis according to sample
Amp1=5;
Amp1=input('What is your Amplitude?');
Def=delay-Min;
%% variables 
Delay=Delay-(def-Def);
n=Min:(Max-Min)/sample:Max;
% Delay=Delay-(def-Def);
% delay=input('What is your delay?');
% decay=1;
% include time shiftsamp
% sample=input('Enter your sample size?');
a=sample/(Max-Min);
% t=Min:a:Max;
% array_part=0;
L1 = length(0:(1/a):Max);
L2 = length(1:1:length(n));
matrix = zeros(L1,L2);
yfinal=zeros(L1,L2);
zfinal=zeros(L1,L2);
l=zeros(L1,L2);
it=0;i=1; 
yy = 0;
qw=100;
as=zeros(L1,L2);
%% plotting the graph
for v=0:(1/a):(Max+500)
   def=(Delay - (v)-Min);
   Def=delay-Min;
   input=1;
   for x=1:1:length(n)
   if x>def*a && x<=def*a+(stop1)*a
       y(x)=Amp*exp(Decay*(x+v*a)*a/sample);
   else
    y(x)=0;
   end
    if x>=Def*a && x<=Def*a+(Stepwidth1)*a% for selecting the Rabge where Step is heppening according to the step width
        
          z(x)=Amp;% take the value pf the amplitude till the step width
 
           else 
        
        z(x)=0;
    end  
      end
        %% convolution
        it = it +1;
     cv=fliplr(y);
     Yfinal(it,:)=cv; 
    Zfinal(it,:)=z;
     v=Yfinal(it,:).*Zfinal(it,:);
     as(it,:)=sum(v)+as(it,:);
    %% plotting ther graph
    subplot(2,2,1)
  plot(as,'linewidth',3)
  xlabel('time');
  ylabel('Convolution function');
  title('Convolution');
  subplot(2,2,2)
  plot(n,cv,'linewidth',3)% plot the graph
  xlabel('time');
  ylabel('function h(T-t)');
  title('h(T-t)');
  pause(0.1)
% hold on
% plot(-n,y,'linewidth',3)
  subplot(2,2,3)
  plot(n,z,'r','linewidth',3)
  xlabel('time');
  ylabel('function f(t)');
  title('f(t)');
  pause(0.1);
  subplot(2,2,4)
  plot(n,cv,'linewidth',3)% plot the graph
  hold on
  plot(n,z,'r')
  hold off
  xlabel('time');
  ylabel('Convolution function');
  title('f(t)*h(T-t)');
  pause(0.1)
% asdf=1:1:100;
% subplot(3,2,4)
% hold on
% plot(as)
% pause(0.1);
% hold off
% xlabel('Time (in sec)');
% ylabel('Amplitude')
% hold off
% pause(0.1)
% o=conv(y,z)
% close
% figure()
% plot(t,l)
% pause(1);
% clear(sprintf('%v',v))
% grid on;
% xlabel('time');
% ylabel('Amplitude');
  p=repmat(y,L1,L2);
end
    case 10
        clc
clear all
close all
%% exponential function:
disp('Exponential function h(t)')
f=0;
% tstart=1;
% tstop=10;
minrange='What is the start point of the graph? Enter your Start point?';
Min=input(minrange);
% Min=-5;
maxrange='what is the end point of the graph?';
Max=input(maxrange);
% Max=5;
Sample=('what is your sample size?');
sample=input(Sample);
% sample=1000;
amp='What is your amplitude?';
Amp=input(amp);
% Amp=1;
% amp=1;
decay='What is your decay rate?';
Decay=input(decay);
% Decay=1;
Delay=input('What is your delay?');
% Delay=0;
stop1('From where ypur function will stop?');
% stop1=1;
def=Delay-Min;
%% sin function
disp('Sinosudial function f(t)')
% Min=-10;
Min=input('what is the minimum range?');
Max=input(maxrange);
% Max=10;
Sample=('How much samples do you need? ');
sample=input(Sample);
% sample=500;
s=sample/(Max-Min);
t=Min:1/s:Max;
frequency=('What is your frequency?');
Frequency=input(frequency);
% Frequency=1;
% frequency=2;
Amp=('what is your Amplitude?');
amp=input(Amp);
% amp=5;
% amp=4;
phase=('What is your phase shift? (in rad)');
Phase=input(phase);
% Phase=0;
delay=input('What is your delay?');
% delay=0;
Stepwidth=input('From where your graph will stop?');
% Stepwidth=2;
def=Delay-Min;
Def=delay-Min;
% close
Delay=Delay-(def-Def);
%% varibles
Delay=Delay-(def-Def);
n=Min:(Max-Min)/sample:Max;
% Delay=Delay-(def-Def);
% delay=input('What is your delay?');
% decay=1;
% include time shiftsamp
% sample=input('Enter your sample size?');
a=sample/(Max-Min);
% t=Min:a:Max;
% array_part=0;
L1 = length(0:(1/a):Max);
L2 = length(1:1:length(n));
matrix = zeros(L1,L2);
yfinal=zeros(L1,L2);
zfinal=zeros(L1,L2);
l=zeros(L1,L2);
it=0;i=1; 
yy = 0;
qw=100;
input=input('Press1 for sin function & press 0 for cos function');
% asd=1:1:100;
as=zeros(L1,L2);
%% plotting of 2 graphs
for v=0:(1/a):(Max+500)
   def=(Delay - (v)-Min);
   Def=delay-Min;
%    input=1;
   for x=1:1:length(n)
   if x>def*a && x<=def*a+(stop1)*a
       y(x)=Amp*exp(Decay*(x+v*a)*a/sample);
   else
    y(x)=0;
   end
   if input==1
         if x>=(Def*a) && x<(Def*a+(Stepwidth)*a)
              z(x)=Amp*sin(2*pi*Frequency*(x)*(1/s)+Phase);
         else
              z(x)=0;
         end
       end
    if input==0
        if x>=(Def*a) && x<(Def*a+(Stepwidth)*a)
              z(x)=Amp1cos(2*pi*Frequency*(x)*(1/s)+Phase);
         else
              z(x)=0;
         end
       end
        
   end
 % Convolution:
       it = it +1;
     cv=fliplr(y);
     Yfinal(it,:)=cv; 
    Zfinal(it,:)=z;
     v=Yfinal(it,:).*Zfinal(it,:);
     as(it,:)=sum(v)+as(it,:);
    %% plotting ther graph
    subplot(2,2,1)
  plot(as,'linewidth',3)
  xlabel('time');
  ylabel('Convolution function');
  title('Convolution');
  subplot(2,2,2)
  plot(n,cv,'linewidth',3)% plot the graph
  xlabel('time');
  ylabel('function h(T-t)');
  title('h(T-t)');
  pause(0.1)
% hold on
% plot(-n,y,'linewidth',3)
  subplot(2,2,3)
  plot(n,z,'r','linewidth',3)
  xlabel('time');
  ylabel('function f(t)');
  title('f(t)');
  pause(0.1);
  subplot(2,2,4)
  plot(n,cv,'linewidth',3)% plot the graph
  hold on
  plot(n,z,'r')
  hold off
  xlabel('time');
  ylabel('Convolution function');
  title('f(t)*h(T-t)');
  pause(0.1)
% asdf=1:1:100;
% subplot(3,2,4)
% hold on
% plot(as)
% pause(0.1);
% hold off
% xlabel('Time (in sec)');
% ylabel('Amplitude')
% hold off
% pause(0.1)
% o=conv(y,z)
% close
% figure()
% plot(t,l)
% pause(1);
% clear(sprintf('%v',v))
% grid on;
% xlabel('time');
% ylabel('Amplitude');
  p=repmat(y,L1,L2);
end
case 11
%% exponential function h(t)
disp 'This is exponential function h(t)'
clc
clear all
close all
disp 'This is exponential function h(t)'
f=0;
% tstart=1;
% tstop=10;
minrange='What is the start point of the graph? Enter your Start point?';
Min=input(minrange);
% Min=-5;
maxrange='what is the end point of the graph?';
Max=input(maxrange);
% Max=5;
Sample=('what is your sample size?');
sample=input(Sample);
% sample=1000;
amp='What is your amplitude?';
Amp1=input(amp);
% Amp1=1;
% amp=1;
decay='What is your decay rate?';
Decay1=input(decay);
% Decay1=1;
Delay=input('What is your delay?');
% Delay=0;
stop1=input('From where your function should stop?');
% stop1=1;
def=Delay-Min;
%% exponential function f(t)
disp 'This is exponential function f(t)'
f=0;
% tstart=1;
% tstop=10;
minrange='What is the start point of the graph? Enter your Start point?';
Min=input(minrange);
% Min=-5;
maxrange='what is the end point of the graph?';
Max=input(maxrange);
% Max=5;
Sample=('what is your sample size?');
sample=input(Sample);
% sample=1000;
amp='What is your amplitude?';
Amp2=input(amp);
% Amp2=2;
% amp=1;
decay='What is your decay rate?';
Decay2=input(decay);
% Decay2=1;
% delay=0;
delay=input('What is your delay?');
% stop2=1;
stop2=input('From where your function should stop?');
Def=delay-Min;
%% initilization of variables
Delay=Delay-(def-Def);
n=Min:(Max-Min)/sample:Max;
% Delay=Delay-(def-Def);
% delay=input('What is your delay?');
% decay=1;
% include time shiftsamp
% sample=input('Enter your sample size?');
a=sample/(Max-Min);
% t=Min:a:Max;
% array_part=0;
L1 = length(0:(1/a):Max);
L2 = length(1:1:length(n));
matrix = zeros(L1,L2);
yfinal=zeros(L1,L2);
zfinal=zeros(L1,L2);
l=zeros(L1,L2);
it=0;i=1; 
yy = 0;
qw=100;

% asd=1:1:100;
as=zeros(L1,L2);
%% ploting both function
for v=0:(1/a):(Max+500)
   def=(Delay - (v)-Min);
   Def=delay-Min;
   for x=1:1:length(n)
   if x>def*a && x<=def*a+(stop1)*a
       y(x)=Amp1*exp(Decay1*(x+v*a)*a/sample);
   else
    y(x)=0;
   end
  
   if x>Def*a && x<=Def*a+(stop2)*a
       z(x)=Amp2*exp(Decay2*x*a/sample);
   else
    z(x)=0;
   end
   end
   %% Convolution:
       it = it +1;
     cv=fliplr(y);
     Yfinal(it,:)=cv; 
    Zfinal(it,:)=z;
     v=Yfinal(it,:).*Zfinal(it,:);
     as(it,:)=sum(v)+as(it,:);
    %% plotting ther graph
    subplot(2,2,1)
plot(as)
subplot(2,2,2)
plot(n,cv,'linewidth',3)% plot the graph
% pause(0.1)
% hold on
% plot(-n,y,'linewidth',3)
subplot(2,2,3)
plot(n,z,'r')
% pause(0.1);

subplot(2,2,1)
  plot(as,'linewidth',3)
  xlabel('time');
  ylabel('Convolution function');
  title('Convolution');
  subplot(2,2,2)
  plot(n,cv,'linewidth',3)% plot the graph
  xlabel('time');
  ylabel('function h(T-t)');
  title('h(T-t)');
  pause(0.1)
% hold on
% plot(-n,y,'linewidth',3)
  subplot(2,2,3)
  plot(n,z,'r','linewidth',3)
  xlabel('time');
  ylabel('function f(t)');
  title('f(t)');
  pause(0.1);
  subplot(2,2,4)
  plot(n,cv,'linewidth',3)% plot the graph
  hold on
  plot(n,z,'r')
  hold off
  xlabel('time');
  ylabel('Convolution function');
  title('f(t)*h(T-t)');
  pause(0.1)
% asdf=1:1:100;
% subplot(3,2,4)
% hold on
% plot(as)
% pause(0.1);
% hold off
% xlabel('Time (in sec)');
% ylabel('Amplitude')
% hold off
% pause(0.1)
% o=conv(y,z)
% close
% figure()
% plot(t,l)
% pause(1);
% clear(sprintf('%v',v))
% grid on;
% xlabel('time');
% ylabel('Amplitude');
  p=repmat(y,L1,L2);
end
    case 12
        %% exponential function h(t)
clc
clear all
close all
disp ('This is exponential function h(t)')
f=0;
% tstart=1;
% tstop=10;
minrange='What is the start point of the graph? Enter your Start point?';
Min=input(minrange);
% Min=-5;
maxrange='what is the end point of the graph?';
Max=input(maxrange);
% Max=5;
Sample=('what is your sample size?');
sample=input(Sample);
% sample=1000;
amp='What is your amplitude?';
Amp1=input(amp);
% Amp1=1;
% amp=1;
decay='What is your decay rate?';
Decay1=input(decay);
% Decay1=1;
Delay=input('What is your delay?');
% Delay=0;
stop1=input('From where your function should stop?');
% stop1=5;
def=Delay-Min;
%% Ramp Signal
disp 'This is ramp Function f(t)'
Min=input('What is your Minimum value of your plot?')
% Min=-10;
Max=input('What is your Maximum value of your plot?')
% Max=10;
% def=Max-Min;
sample=input('How mnay samples do you need?');
% sample=100;
a=sample/(Max-Min);
n=Min:1/a:Max;
start=0;
End=input('From where your ramp should end?')
% End=5;
delay=input('What is your Dealy?');
% delay=0;
decay=input('What is ypur decay rate? Should be 1 or -1');
% decay=1;
Def=0-Min;
% Delay=-1;
%% Intilization of variables for convolution
array_part=0;
L1 = length(0:(1/a):Max);
L2 = length(1:1:length(n));
matrix = zeros(L1,L2);
yfinal=zeros(L1,L2);
zfinal=zeros(L1,L2);
l=zeros(L1,L2);
it=0;i=1;
yy = 0;
qw=100;
% asd=1:1:100;
as=zeros(L1,L2);
%% Calculation of F(t) & h(T-t)
for v=0:(1/a):(Max+500)
   def=Delay-(v)-Min;
   Def=0-Min;
   for x=1:1:length(n)
   if x>def*a && x<=def*a+(stop1)*a
       y(x)=Amp1*exp(Decay1*(x+v*a)*a/sample);
   else
    y(x)=0;
   end
   if x>(Def+delay)*a && x<=(Def+End+delay)*a
        z(x)=decay*n(x-(delay*a));
    else
       z(x)=0;
    end
   end
   it = it +1;
 
   
    cv=fliplr(y);
    Yfinal(it,:)=cv; 
    Zfinal(it,:)=z;
    v=Yfinal(it,:).*Zfinal(it,:);
     as(it,:)=sum(v)+as(it,:);
  % end
%   ed=Min:Sample:Max;
   subplot(2,2,1)
  plot(as,'linewidth',3)
  xlabel('time');
  ylabel('Convolution function');
  title('Convolution');
  subplot(2,2,2)
  plot(n,cv,'linewidth',3)% plot the graph
  xlabel('time');
  ylabel('function h(T-t)');
  title('h(T-t)');
  pause(0.1)
% hold on
% plot(-n,y,'linewidth',3)
  subplot(2,2,3)
  plot(n,z,'r')
  xlabel('time');
  ylabel('function f(t)');
  title('f(t)');
  pause(0.1);
  subplot(2,2,4)
  plot(n,cv,'linewidth',3)% plot the graph
  hold on
  plot(n,z,'r')
  hold off
  xlabel('time');
  ylabel('Convolution function');
  title('f(t)*h(T-t)');
  pause(0.1)
% asdf=1:1:100;
% subplot(3,2,4)
% hold on
% plot(as)
% pause(0.1);
% hold off
% xlabel('Time (in sec)');
% ylabel('Amplitude')
% hold off
% pause(0.1)
% o=conv(y,z)
% close
% figure()
% plot(t,l)
% pause(1);
% clear(sprintf('%v',v))
% grid on;
% xlabel('time');
% ylabel('Amplitude');
%   p=repmat(y,L1,L2);
end
    case 13
       clc
close all
clear all
%% Ramp function:
disp 'This is Ramp Function h(t)'
Min=input('What is your Minimum value of your plot?')
% Min=-10;
Max=input('What is your Maximum value of your plot?')
% Max=10;
% def=Max-Min;
Sample=input('How mnay samples do you need?');
% Sample=1000;
a=Sample/(Max-Min);
n=Min:1/a:Max;
start=0;
End=input('From where your ramp should end?')
% End=5;
Delay=input('What is your Dealy?');
% Delay=0;
decay=input('What is ypur decay rate? Should be 1 or -1');
% decay=-1;
def=0-Min;
% Delay=-1;
%% step function:
disp('This is Step function:f(t)')
min='what is your min range?';
Min=input(min);
% Min=-10;
max='What is your max range?';
Max=input(max);
% Max=10;
Sample='how many sample you want?';
sample=input(Sample);
% sample=1000;
a=sample/(Max-Min);
Delay='How much time shift do you need?';
delay=input(Delay);
% delay=0;
Def=delay-Min;
% def=def-(def-Def);
Stepwidth2=input('How much is the step width?');
% Stepwidth2=5;
n=Min:(Max-Min)/sample:Max;% Scalling of the x axis according to sample
Amp2=10;
Amp2=input('What is your Amplitude?');
%% Intilization of variables for convolution
array_part=0;
L1 = length(0:(1/a):Max);
L2 = length(1:1:length(n));
matrix = zeros(L1,L2);
yfinal=zeros(L1,L2);
zfinal=zeros(L1,L2);
l=zeros(L1,L2);
it=0;i=1;
yy = 0;
qw=100;
% asd=1:1:100;
as=zeros(L1,L2);
%% Calculation:
for v=0:(1/a):(Max+500)
   def=0-(v)-Min;
   Def=delay-Min;
   for x=1:1:length(n)
    if x>(def+Delay)*a && x<=(def+End+Delay)*a
        y(x)=decay*n(x-(Delay*a)+(v*a));
    else
       y(x)=0;
    end
    

    if x>=Def*a && x<Def*a+(Stepwidth2)*a% for selecting the Rabge where Step is heppening according to the step width
        
        z(x)=Amp2;% take the value pf the amplitude till the step width
        
    else 
        
        z(x)=0;
        
    end
  
   end
    it = it +1;
    cv=fliplr(y);
    Yfinal(it,:)=cv; 
    Zfinal(it,:)=z;
    v=Yfinal(it,:).*Zfinal(it,:);
   as(it,:)=sum(v)+as(it,:);
  % end
  ed=Min:Sample:Max;
   subplot(2,2,1)
  plot(as,'linewidth',3)
  xlabel('time');
  ylabel('Convolution function');
  title('Convolution');
  subplot(2,2,2)
  plot(n,cv,'linewidth',3)% plot the graph
  xlabel('time');
  ylabel('function h(T-t)');
  title('h(T-t)');
  pause(0.1)
% hold on
% plot(-n,y,'linewidth',3)
  subplot(2,2,3)
  plot(n,z,'r')
  xlabel('time');
  ylabel('function f(t)');
  title('f(t)');
  pause(0.1);
  subplot(2,2,4)
  plot(n,cv,'linewidth',3)% plot the graph
  hold on
  plot(n,z,'r')
  hold off
  xlabel('time');
  ylabel('Convolution function');
  title('f(t)*h(T-t)');
  pause(0.1)
% asdf=1:1:100;
% subplot(3,2,4)
% hold on
% plot(as)
% pause(0.1);
% hold off
% xlabel('Time (in sec)');
% ylabel('Amplitude')
% hold off
% pause(0.1)
% o=conv(y,z)
% close
% figure()
% plot(t,l)
% pause(1);
% clear(sprintf('%v',v))
% grid on;
% xlabel('time');
% ylabel('Amplitude');
%   p=repmat(y,L1,L2);
     end
  
    case 14
        %% ramp function:
clc
close all
clear all
disp 'This is Ramp Function h(t)'
Min=input('What is your Minimum value of your plot?')
% Min=-10;
Max=input('What is your Maximum value of your plot?')
% Max=10;
def=Max-Min;
Sample=input('How mnay samples do you need?');
% Sample=1000;
a=Sample/(Max-Min);
n=Min:1/a:Max;
start=0;
End=input('From where your ramp should end?')
% End=5;
Delay=input('What is your Dealy?');
% Delay=5;
decay=input('What is ypur decay rate? Should be 1 or -1');
% decay=-1;
def=0-Min;
% Delay=-1;
%% Sinosudial function:
disp('This is sin function f(t)');
Minrange=input('What is your minimum range of the plot?');
% Minrange=-10;
Maxrange=input(maxrange);
% Maxrange=10;
sample=('How much samples do you need? ');
Sample=input(sample);
% Sample=1000;
s=Sample/(Maxrange-Minrange);
t=Minrange:1/s:Maxrange;
frequency=('What is your frequency?');
Frequency=input(frequency);
% Frequency=1;
% frequency=2;
amp=('what is your Amplitude?');
Amp=input(amp);
% Amp=5;
% Stepwidth=2;
Stepwidth=input('From where your function should stop plotting your function?');
% amp=4;
phase=('What is your phase shift? (in rad)');
Phase=input(phase);
% Phase=0;
delay=input('What is your delay? ');
% delay=0;
def=Delay-Minrange;
Def=delay-Min;
% close
Delay=Delay-(def-Def);
array_part=0;
 %% Intilization of variables for convolution
L1 = length(0:(1/a):Max);
L2 = length(1:1:length(n));
matrix = zeros(L1,L2);
yfinal=zeros(L1,L2);
zfinal=zeros(L1,L2);
l=zeros(L1,L2);
it=0;i=1;
yy = 0;
qw=100;
% asd=1:1:100;
as=zeros(L1,L2);
input=input('Choose from cos and sin function ? 1 for sin & 0 for cos.');
er=input('Press 1 for sin function & 0 for cosine function');
for v=0:(1/a):(Max+500)
   def=0-(v)-Minrange;
   Def=delay-Min;
   for x=1:1:length(n)
    if x>(def+Delay)*a && x<=(def+End+Delay)*a
        y(x)=decay*n(x-(Delay*a)+(v*a));
    else
       y(x)=0;
    end
   end
   for x=1:1:length(n)
    if er==1
    
         if x>=Def*a && x<Def*a+(Stepwidth)*a
              z(x)=Amp*sin(2*pi*Frequency*x*(1/s)+Phase);
         else
              z(x)=0;
         end
         if er==0
             if x>=Def*a && x<=Def*a+(Stepwidth)*a
                 z(x)=Amp*cos(2*pi*Frequency*x*(1/s)+Phase);
             else
                 z(x)=0;
             end
         end
    end
   end
  
    %% Calculation for convolution:
% figure
% eval(sprintf('%v = y(x)', v))
% save('file.mat', sprintf('%v',v)) 
     it = it +1;
    cv=fliplr(y);
    Yfinal(it,:)=cv; 
     Zfinal(it,:)=z;
% g=Yfinal(it,it).*Zfinal(it,it)
% l(1:it,:)=Yfinal.*Zfinal; 
% v=Yfinal(i,round(def*a)+1).*Zfinal(i,round(Def*a)+1);
% for i=1:1:qw
% v=Yfinal(i,round(def*a)+1).*Zfinal(i,round(Def*a)+1)
     v=Yfinal(it,:).*Zfinal(it,:);
    as(it,:)=sum(v)+as(it,:);
% end
%% Plotting the convolution:
    subplot(2,2,1)
     plot(as,'linewidth',3)
    xlabel('time');
    ylabel('Convolution function');
    title('Convolution');
  subplot(2,2,2)
  plot(n,cv,'linewidth',3)% plot the graph
  xlabel('time');
  ylabel('function h(T-t)');
  title('h(T-t)');
  pause(0.1)
% hold on
% plot(-n,y,'linewidth',3)
  subplot(2,2,3)
  plot(n,z,'r')
  xlabel('time');
  ylabel('function f(t)');
  title('f(t)');
  pause(0.1);
  subplot(2,2,4)
  plot(n,cv,'linewidth',3)% plot the graph
  hold on
  plot(n,z,'r')
  hold off
  xlabel('time');
  ylabel('Convolution function');
  title('f(t)*h(T-t)');
  pause(0.1)
% asdf=1:1:100;
% subplot(3,2,4)
% hold on
% plot(as)
% pause(0.1);
% hold off
% xlabel('Time (in sec)');
% ylabel('Amplitude')
% hold off
% pause(0.1)
% o=conv(y,z)
% close
% figure()
% plot(t,l)
% pause(1);
% clear(sprintf('%v',v))
% grid on;
% xlabel('time');
% ylabel('Amplitude');
  p=repmat(y,L1,L2);
end
case 15
        
%% Ramp Function h(t)
clc
close all
clear all
disp 'This is ramp Function'
Min=input('What is your Minimum value of your plot?')
% Min=-10;
Max=input('What is your Maximum value of your plot?')
% Max=10;
def=Max-Min;
Sample=input('How mnay samples do you need?');
% Sample=1000;
a=Sample/(Max-Min);
n=Min:1/a:Max;
start=0;
End=input('From where your ramp should end?')
% End1=5;
Delay=input('What is your Dealy?');
% Delay=0;
decay1=input('What is ypur decay rate? Should be 1 or -1');
% decay=-1;
% decay1=1;
def=0-Min;
%% Ramp Function hft)
disp 'This is Ramp Function f(t)'
clc
close all
clear all
Min=input('What is your Minimum value of your plot?')
% Min=-10;
Max=input('What is your Maximum value of your plot?')
% Max=10;
def=Max-Min;
Sample=input('How mnay samples do you need?');
% Sample=1000;
a=Sample/(Max-Min);
n=Min:1/a:Max;
start=0;
End2=input('From where your ramp should end?')
% End2=3;
delay=input('What is your Dealy?');
% delay=0;
decay2=input('What is ypur decay rate? Should be 1 or -1');
% decay2=-1;
Def=0-Min;
% def=Delay-(v)-Min;
%    Def=delay-Min;
%% Variables
L1 = length(0:(1/a):Max);
L2 = length(1:1:length(n));
matrix = zeros(L1,L2);
yfinal=zeros(L1,L2);
zfinal=zeros(L1,L2);
l=zeros(L1,L2);
it=0;i=1;
yy = 0;
qw=100;
% asd=1:1:100;
as=zeros(L1,L2);

%% Calculation of both functions
for v=0:(1/a):(Max+500)
   Delay=5;
   End1=5;
   decay1=-1;
   def=0-(v)-Min;
   Def=0-Min;
   for x=1:1:length(n)
    if x>(def+Delay)*a && x<=(def+End1+Delay)*a
        y(x)=decay1*n(x-(Delay*a)+(v*a));
    else
       y(x)=0;
    end
    if x>(Def+delay)*a && x<=(Def+End2+delay)*a
        z(x)=decay2*n(x-(delay*a));
    else
       z(x)=0;
    end
   end
   %% Calculation for convolution:
% figure
% eval(sprintf('%v = y(x)', v))
% save('file.mat', sprintf('%v',v)) 
     it = it +1;
    cv=fliplr(y);
    Yfinal(it,:)=cv; 
     Zfinal(it,:)=z;
% g=Yfinal(it,it).*Zfinal(it,it)
% l(1:it,:)=Yfinal.*Zfinal; 
% v=Yfinal(i,round(def*a)+1).*Zfinal(i,round(Def*a)+1);
% for i=1:1:qw
% v=Yfinal(i,round(def*a)+1).*Zfinal(i,round(Def*a)+1)
     v=Yfinal(it,:).*Zfinal(it,:);
    as(it,:)=sum(v)+as(it,:);
% end
%% Plotting the convolution:
    subplot(2,2,1)
     plot(as,'linewidth',3)
    xlabel('time');
    ylabel('Convolution function');
    title('Convolution');
  subplot(2,2,2)
  plot(cv,'linewidth',3)% plot the graph
  xlabel('time');
  ylabel('function h(T-t)');
  title('h(T-t)');
  pause(0.1)
% hold on
% plot(-n,y,'linewidth',3)
  subplot(2,2,3)
  plot(n,z,'r')
  xlabel('time');
  ylabel('function f(t)');
  title('f(t)');
  pause(0.1);
  subplot(2,2,4)
  plot(n,cv,'linewidth',3)% plot the graph
  hold on
  plot(n,z,'r')
  hold off
  xlabel('time');
  ylabel('Convolution function');
  title('f(t)*h(T-t)');
  pause(0.1)
% asdf=1:1:100;
% subplot(3,2,4)
% hold on
% plot(as)
% pause(0.1);
% hold off
% xlabel('Time (in sec)');
% ylabel('Amplitude')
% hold off
% pause(0.1)
% o=conv(y,z)
% close
% figure()
% plot(t,l)
% pause(1);
% clear(sprintf('%v',v))
% grid on;
% xlabel('time');
% ylabel('Amplitude');
  p=repmat(y,L1,L2);
end

    case 16
        clc
close all
clear all
%% Ramp function:
disp 'This is Ramp Function h(f)'
Min=input('What is your Minimum value of your plot?')
% Min=-10;
Max=input('What is your Maximum value of your plot?')
% Max=10;
% def=Max-Min;
Sample=input('How mnay samples do you need?');
% Sample=1000;
a=Sample/(Max-Min);
n=Min:1/a:Max;
start=0;
End=input('From where your ramp should end?')
% End=5;
Delay=input('What is your Dealy?');
% Delay=0;
decay=input('What is ypur decay rate? Should be 1 or -1');
% decay=1;
def=0-Min;
% Delay=-1;
%% exponential function:
%% Exponential Function f(t)
disp 'This is exponential function f(t)'
f=0;
% tstart=1;
% tstop=10;
min='What is the start point of the graph? Enter your Start point?';
Min=input(minrange);
% Min=-5;
max='what is the end point of the graph?';
Max=input(maxrange);
% Max=5;
sample=('what is your sample size?');
Sample=input(sample);
% sample=1000;
amp='What is your amplitude?';
Amp2=input(amp);
% Amp2=1;
% amp=1;
decay='What is your decay rate?';
Decay=input(decay);
% Decay=1;
delay=Input('What is your delay?');
% delay=0;
stop=input('From where your graph will stop?');
% stop=1;
Def=delay-Min;
%% initialization of variables for convolution:
% Delayx=Delay-(def-Def);
% delay=input('What is your delay?');
% decay=1;
% include time shiftsamp
% sample=input('Enter your sample size?');
a=sample/(Max-Min);
% t=Min:a:Max;
array_part=0;
L1 = length(0:(1/a):Max);
L2 = length(1:1:length(n));
matrix = zeros(L1,L2);
yfinal=zeros(L1,L2);
zfinal=zeros(L1,L2);
l=zeros(L1,L2);
it=0;i=1;
yy = 0;
qw=100;
% asd=1:1:100;
as=zeros(L1,L2);
% Calculation:
for v=0:(1/a):(Max+500)
   def=-(v)-Min;
   Def=delay-Min;
   for x=1:1:length(n)
    if x>(def+Delay)*a && x<=(def+End+Delay)*a
        y(x)=decay*n(x-(Delay*a)+(v*a));
    else
       y(x)=0;
    end
    if x>Def*a && x<=Def*a+(stop)*a
       z(x)=Amp2*exp(Decay*x/a);
   else
    z(x)=0;
    end
   end
   it = it +1;
    cv=fliplr(y);
    Yfinal(it,:)=cv; 
    Zfinal(it,:)=z;
    v=Yfinal(it,:).*Zfinal(it,:);
   as(it,:)=sum(v)+as(it,:);
  % end
  ed=Min:Sample:Max;
   subplot(2,2,1)
  plot(as,'linewidth',3)
  xlabel('time');
  ylabel('Convolution function');
  title('Convolution');
  subplot(2,2,2)
  plot(n,cv,'linewidth',3)% plot the graph
  xlabel('time');
  ylabel('function h(T-t)');
  title('h(T-t)');
  pause(0.1)
% hold on
% plot(-n,y,'linewidth',3)
  subplot(2,2,3)
  plot(n,z,'r')
  xlabel('time');
  ylabel('function f(t)');
  title('f(t)');
  pause(0.1);
  subplot(2,2,4)
  plot(n,cv,'linewidth',3)% plot the graph
  hold on
  plot(n,z,'r')
  hold off
  xlabel('time');
  ylabel('Convolution function');
  title('f(t)*h(T-t)');
  pause(0.1)
% asdf=1:1:100;
% subplot(3,2,4)
% hold on
% plot(as)
% pause(0.1);
% hold off
% xlabel('Time (in sec)');
% ylabel('Amplitude')
% hold off
% pause(0.1)
% o=conv(y,z)
% close
% figure()
% plot(t,l)
% pause(1);
% clear(sprintf('%v',v))
% grid on;
% xlabel('time');
% ylabel('Amplitude');
%   p=repmat(y,L1,L2);
     end
 
    otherwise
        disp('Renter the value again')
end

