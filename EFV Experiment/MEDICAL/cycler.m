
%% TESTING CYCLE
%Filename save
clear all
close all
clc

s = connectDAQ(1,1000);
%input('\nSTART PUMP & Press Enter to Start Test: \n');
%Create Listener
figure(1)
hold on
lh = addlistener(s,'DataAvailable',@stashDatab);

s.startBackground()
