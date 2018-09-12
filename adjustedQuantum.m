% Ion number starts from 1.
% For the neutral argon atoms put in 1.
%
function [ y ] = adjustedQuantum( Z,iP )

y=(Z)*(2*iP)^(-1/2);
end