function [value, isterminal, directions] = event_delays( t, xi)
%EVENTS Summary of this function goes here
%   Detailed explanation goes here
global N ieta il itau delay

value       = [-min(xi(ieta));sum(xi(il)==1&kron(ones(N,1),xi(itau))>=delay)==sum(xi(il)==1)];
isterminal  = [1;1];
directions   = [1; 0];


end

