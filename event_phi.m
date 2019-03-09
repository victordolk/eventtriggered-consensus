function [value, isterminal, directions] = event_phi( t, phi,gamma,lambda,phi10)
%EVENTS Summary of this function goes here
%   Detailed explanation goes here

value       = [1/lambda*phi(2)-phi(1);
               phi(1)-lambda*phi10];
isterminal  = [1;1];
directions  = [-1;-1];

end

