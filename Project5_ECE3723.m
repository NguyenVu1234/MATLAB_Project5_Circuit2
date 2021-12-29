clc;
close all;
clear all;

%Part 1
% Define the signal parameters
Vm = 1;
T = 1;
w0 = 2*pi/T;

% Define the symbolic variables
syms n t;

% Define the signal
v1 = Vm*sin(4*pi*t/T);
v2 = 3*Vm*sin(4*pi*t/T);

% Evaluate the fourier series integral
ak = 2/T*int(v1*cos(n*w0*t),0,T/2) + 2/T*int(v2*cos(n*w0*t),T/2,T);
bk = 2/T*int(v1*sin(n*w0*t),0,T/2) + 2/T*int(v2*sin(n*w0*t),T/2,T);
av = 1/T*int(v1,0,T/2) + 1/T*int(v2,T/2,T);

% Declare the range of terms
nMax = 100;
n = 1:nMax;
a = subs(ak);
b = subs(bk);

% define the time vector
ts = 1e-2; % ts is sampling the
t = -1:ts:3*T-ts;

% directly plot the signal x(t)
t1 = -1:ts:0-ts;
v1 = Vm*sin(4*pi*t1/T).*(t1<=-T/2);
v2 = 3*Vm*sin(4*pi*t1/T).*(t1>-T/2).*(t1<0);
v = v1+v2;
x = repmat(v,1,4);


plot(t,x,'r','linewidth',2);
grid on;
title('Original signal v(t)')
xlabel('time(s)');
ylabel('Amplitude');
ylim([-3 3]);

% %Part2
% %First 4 terms Harmonic n = 1,2,3,4
% for i = 1:length(t)
% for k = 1
% x(k,i) = a(k)*cos(k*w0*t(i)) + b(k)*sin(k*w0*t(i));
% end
% y(i) = av+sum(x(:,i)); % Add DC term
% end
% subplot(2,2,1)
% plot(t,y)
% xlabel('time(s)');ylabel('Amplitude');
% ylim([-3 3]);
% title('n = 1')
% 
% for i = 1:length(t)
% for k = 1:2
% x(k,i) = a(k)*cos(k*w0*t(i)) + b(k)*sin(k*w0*t(i));
% end
% y(i) = av+sum(x(:,i)); % Add DC term
% end
% subplot(2,2,2)
% plot(t,y)
% xlabel('time(s)');ylabel('Amplitude');
% ylim([-3 3]);
% title('n = 2')
% for i = 1:length(t)
% for k = 1:3
% x(k,i) = a(k)*cos(k*w0*t(i)) + b(k)*sin(k*w0*t(i));
% end
% y(i) = av+sum(x(:,i)); % Add DC term
% end
% subplot(2,2,3)
% plot(t,y)
% xlabel('time(s)');ylabel('Amplitude');
% ylim([-3 3]);
% title('n = 3')
% 
% for i = 1:length(t)
% for k = 1:4
% x(k,i) = a(k)*cos(k*w0*t(i)) + b(k)*sin(k*w0*t(i));
% end
% y(i) = av+sum(x(:,i)); % Add DC term
% end
% subplot(2,2,4)
% plot(t,y)
% xlabel('time(s)');ylabel('Amplitude');
% ylim([-3 3]);
% title('n = 4')

% %Part 3
% %Harmonic n = 3,5,10,50
% for i = 1:length(t)
% for k = 1:3
% x(k,i) = a(k)*cos(k*w0*t(i)) + b(k)*sin(k*w0*t(i));
% end
% y(i) = av+sum(x(:,i)); % Add DC term
% end
% subplot(2,2,1)
% plot(t,y)
% xlabel('time(s)');ylabel('Amplitude');
% ylim([-3 3]);
% title('n = 3')
% 
% for i = 1:length(t)
% for k = 1:5
% x(k,i) = a(k)*cos(k*w0*t(i)) + b(k)*sin(k*w0*t(i));
% end
% y(i) = av+sum(x(:,i)); % Add DC term
% end
% subplot(2,2,2)
% plot(t,y)
% xlabel('time(s)');ylabel('Amplitude');
% ylim([-3 3]);
% title('n = 5')
% for i = 1:length(t)
% for k = 1:10
% x(k,i) = a(k)*cos(k*w0*t(i)) + b(k)*sin(k*w0*t(i));
% end
% y(i) = av+sum(x(:,i)); % Add DC term
% end
% subplot(2,2,3)
% plot(t,y)
% xlabel('time(s)');ylabel('Amplitude');
% ylim([-3 3]);
% title('n = 10')
% 
% for i = 1:length(t)
% for k = 1:50 
% x(k,i) = a(k)*cos(k*w0*t(i)) + b(k)*sin(k*w0*t(i));
% end
% y(i) = av+sum(x(:,i)); % Add DC term
% end
% subplot(2,2,4)
% plot(t,y)
% xlabel('time(s)');ylabel('Amplitude');
% ylim([-3 3]);
% title('n = 50')

% %Part 4
% % Error vs Time
% for i = 1:length(t)
% for numTerm = 1:3 
% %p preresents for Vf
% p(numTerm,i) = a(numTerm)*cos(numTerm*w0*t(i)) + b(numTerm)*sin(numTerm*w0*t(i));
% end
% 
% y(i) = av+sum(p(:,i)); 
% end
% % fnc preresent for error
% fnc = ((abs(y-x))/3)*100;
% subplot(221)
% plot(t,fnc)
% title('n = 3')
% xlabel('time(s)')
% ylabel('error(%)')
% 
% for i = 1:length(t)
% for numTerm = 1:5 
% p(numTerm,i) = a(numTerm)*cos(numTerm*w0*t(i)) + b(numTerm)*sin(numTerm*w0*t(i));
% end
% y(i) = av+sum(p(:,i)); 
% end
% fnc = ((abs(y-x))/3)*100;
% subplot(222)
% plot(t,fnc)
% title('n = 5')
% xlabel('time(s)')
% ylabel('error(%)')
% 
% for i = 1:length(t)
% for numTerm = 1:10 
% p(numTerm,i) = a(numTerm)*cos(numTerm*w0*t(i)) + b(numTerm)*sin(numTerm*w0*t(i));
% end
% y(i) = av+sum(p(:,i)); 
% end
% fnc = ((abs(y-x))/3)*100;
% subplot(223)
% plot(t,fnc)
% title ('n = 10')
% xlabel('time(s)')
% ylabel('error(%)')
% 
% for i = 1:length(t)
% for numTerm = 1:50 
% p(numTerm,i) = a(numTerm)*cos(numTerm*w0*t(i)) + b(numTerm)*sin(numTerm*w0*t(i));
% end
% y(i) = av+sum(p(:,i)); 
% end
% fnc = ((abs(y-x))/3)*100;
% subplot(224)
% plot(t,fnc)
% title ('n = 50')
% xlabel('time(s)')
% ylabel('error(%)')

% %Part 5
% % Create the vector that contain maximum error 11.3177, 7.27565, 4.287,
% % 0.849166
% maxError = [11.3177 7.27565 4.287 0.849166];
% numTerm = linspace(3,50,4.1);
% plot(numTerm, maxError);
% grid on;
% xlabel('number of terms');
% ylabel('%Error');
% title('%Max Error vs number of terms');






