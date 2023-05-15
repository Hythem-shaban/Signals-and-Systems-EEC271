function x = delta(time_initial,time_final,to,t)
time = abs(time_final-time_initial);
x = [zeros(1,round(length(t)*abs(to-time_initial)/time)) 1 zeros(1,round(length(t)*abs(time_final-to)/(time))-1)] ;
end