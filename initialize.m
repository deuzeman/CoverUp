global almanac
global opts

comp = readdata('compact.dat'); %load dataset
b190 = comp(comp.beta==1.90,:);
b195 = comp(comp.beta==1.95,:);
b210 = comp(comp.beta==2.10,:);

nf211 = [b190; b195; b210];
nf211_un = nf211;
nf211_un([1, 2, 8, 9, 16],:) = [];

almanac.mpi     = 139.6;
almanac.fpi     = 130.7;
almanac.mK      = 493.7;
almanac.fK      = 159.8;
almanac.mN      = 939;
almanac.mDelta  = 1232;
almanac.gA      = 1.2695;
almanac.r0      = 0.44;
almanac.w0      = 0.1755;
almanac.mev_fm  = 197.326;
almanac.fm_mev  = almanac.mev_fm;

almanac.l1.ave = -0.4;
almanac.l1.std =  0.6; 
almanac.l2.ave =  4.3;
almanac.l2.std =  0.1;
almanac.l3.ave =  2.9;
almanac.l3.std =  2.4;
almanac.l4.ave =  4.4;
almanac.l4.std =  0.2;

opts.asq  = 'ON'; % Other option: 'OFF'
opts.iso = 'OFF'; % Oter option: 'ON'
opts.fvol = 'CDH'; % Other options:  'CWW', 'BMW', 'GL', 'OFF'
opts.plot = 'ON'; % Other option: 'OFF'
opts.priors = 'MIN'; % Other options: 'ON', 'OFF', 'NOZP'
opts.nnlo = 'OFF'; % Other options: 'BAER', 'ON'
opts.wipe = 'ON'; % Other options: 'OFF'
opts.nboot = 100; % Number of samples
