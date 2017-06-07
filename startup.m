curpath = pwd;
addpath(curpath);
dirs = {'MAT', 'MAT_HADA', 'TT', 'TT_HADA'};
for s = dirs
  addpath(sprintf('%s/%s',curpath,s{:}))
end
