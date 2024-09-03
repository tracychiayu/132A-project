function h = rrc(t,alpha,T)
% function def
%         Generates the impulse response, h(t), of a RRC filter

% INPUT: 
%       - t: Time vector at which to evaluate the RRC filter impulse response, h(t).
%       - alpha: Roll-off factor of RRC filter (0<=alpha<=1)
%       - T: Symbol period (T=1/Rs)
% OUTPUT: 
%       - h: Impulse response of RCC filter, denoted as h(t) in time domain

% To see the detailed formula for impulse response h(t) of RRC filter,
% please refer to the corresponding formula in the text section of the
% MATLAB live script.
h = zeros(size(t));
for i=1:length(t)
    if t(i)==0
        h(i)=(1+alpha*(4/pi-1))/T;
    elseif t(i)==T/(4*alpha) || t(i)==-T/(4*alpha)
        h(i)=(alpha/T/sqrt(2))*((1+2/pi)*sin(pi/4/alpha)+...
              (1-2/pi)*cos(pi/4/alpha));
    else
        h(i)=(sin((pi*t(i)/T)*(1-alpha))+(4*alpha*t(i)/T)*cos((pi*t(i)/T)*(1+alpha)))/...
             ((pi*t(i))*(1-(4*alpha*t(i)/T)^2));
    end
end