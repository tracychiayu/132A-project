function h = rrc(t,alpha,T)
% function def
%         Generates the impulse response, h(t), of a RRC filter

% INPUT: 
%       - t: Time vector at which to evaluate the RRC filter impulse response, h(t).
%       - alpha: Roll-off factor of RRC filter (0<=alpha<=1)
%       - T: Symbol period (T=1/Rs)
% OUTPUT: 
%       - P: Impulse response of RCC filter, denoted as h(t) in time domain

% t(t==0) = eps; % T/1e5; % gets rid of t==0, which is not defined
h = zeros(size(t));
h(t~=0) = 4*alpha/(pi*sqrt(T)) * ( cos((1+alpha)*pi*t(t~=0)/T) +...
    sin((1-alpha)*pi*t(t~=0)/T)./(4*alpha*t(t~=0)/T) ) ./...
    (1 - (4*alpha*t(t~=0)/T).^2);
h(t==0) = (pi+alpha*(4-pi))/(pi*sqrt(T));
h((t==T/4/alpha)|(t==-T/4/alpha)) = (alpha/sqrt(2*T))*...
    ((1+2/pi)*sin(pi/4/alpha) + (1-2/pi)*cos(pi/4/alpha));