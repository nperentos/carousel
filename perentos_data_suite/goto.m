function [out] = goto(fileBase)
% [out] = goto(fileBase) changes the current path to the fileBase's
% processed path and returns the path to the 'out' variable
cd(getfullpath(fileBase));
pwd
out = pwd;