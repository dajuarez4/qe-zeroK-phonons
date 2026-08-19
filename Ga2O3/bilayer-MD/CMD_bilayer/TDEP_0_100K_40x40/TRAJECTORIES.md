# Trajectory storage

The complete 16,000-atom trajectories (about 32 MB each) are retained in `/home/dajuarez4/bilayer_ga2o3/tdep_<T>K_40x40/` and are intentionally not duplicated in Git. Each finite-temperature trajectory contains 26 frames sampled every 1 ps over the 25 ps production window. The five frames used for each TDEP fit are represented by `infile.positions`, `infile.forces`, `infile.stat`, and `infile.meta` in the local calculation; the large position/force matrices are also omitted here to keep the report repository compact.
