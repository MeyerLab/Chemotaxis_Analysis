function subdir=getSubdirectories(path,regExp)
% returns a list of the subdirectories in a specified directory

d0=dir(path);
ind=[d0.isdir];
subdir={d0(ind).name};
subdir=vect(setdiff(subdir,{'.','..'}));

if nargin>1
   subdir=subdir(boolRegExp(subdir,regExp));
end
