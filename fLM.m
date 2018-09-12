function [ y ] = fLM( l,m )

y=(2*l+1)*factorial(l+abs(m))/((2^(abs(m))*factorial(abs(m))*factorial(l-abs(m))));
end