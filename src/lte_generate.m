function [wave, grid, info] = lte_generate(nrb, rho, ncellid, nsubframes)
%LTE_GENERATE  Synthesise an LTE downlink baseband waveform (module M1).
%
%   [WAVE, GRID, INFO] = LTE_GENERATE(NRB, RHO, NCELLID, NSUBFRAMES)
%
%   Generates a 3GPP-style LTE FDD downlink signal using base MATLAB only.
%   No LTE Toolbox function is called; the resource grid is populated by
%   hand and the OFDM modulation is a plain IFFT plus cyclic prefix.
%
%   Inputs
%     NRB        number of downlink resource blocks: 6, 25, 50 or 100
%     RHO        resource occupancy (load factor) of the PDSCH, 0.0 .. 1.0
%     NCELLID    physical cell identity, 0 .. 503
%     NSUBFRAMES number of 1 ms subframes to generate (>= 1)
%
%   Outputs
%     WAVE   complex column vector, time domain baseband samples,
%            normalised to unit mean power
%     GRID   complex matrix, 12*NRB by 14*NSUBFRAMES resource grid,
%            row 1 is the lowest occupied subcarrier, the DC subcarrier is
%            not a row of the grid because it is never transmitted
%     INFO   struct with fields
%              nfft, fs, cplen, occupied_bw          (required by the PRD)
%              nrb, rho, ncellid, nsubframes, nsc, scs,
%              cp_first, cp_other, samples_per_subframe, nsamples,
%              papr_db, rms_scale, n_pdsch_avail, n_pdsch_used, rho_actual
%
%   Bandwidth table implemented (sample rate is always nfft * 15000):
%
%     nrb  BW        nfft   fs (MHz)  occupied SC  occupied BW (MHz)
%       6  1.4 MHz    128     1.92         72            1.08
%      25  5   MHz    512     7.68        300            4.50
%      50  10  MHz   1024    15.36        600            9.00
%     100  20  MHz   2048    30.72       1200           18.00
%
%   Physical signals inserted, in this order:
%     1. CRS  - OFDM symbols 0 and 4 of every slot, every 6th subcarrier,
%               frequency shift mod(ncellid,6), QPSK values from the LTE
%               length-31 Gold (pseudo random) sequence
%     2. PSS  - last OFDM symbol of slots 0 and 10, central 62 subcarriers,
%               Zadoff-Chu with root 25 / 29 / 34 chosen by mod(ncellid,3)
%     3. SSS  - second to last OFDM symbol of slots 0 and 10, central 62
%               subcarriers, the 36.211 m-sequence construction
%     4. PDSCH- a random fraction RHO of every remaining occupied resource
%               element, filled with unit power QPSK, the rest left at zero
%
%   Note on the CRS lattice: 36.211 staggers the second CRS symbol of each
%   slot by three subcarriers.  The PRD for this project specifies a single
%   frequency shift mod(ncellid,6) for every CRS symbol, so that is what is
%   implemented here.  The consequence is a strictly periodic 6-subcarrier
%   comb, which is what the delay ambiguity prediction of M2 assumes.
%
%   A note on the name of the second output.  It is called grid because the
%   PRD's contract calls it that, and it shadows the built-in grid function
%   in whatever workspace receives it.  Scripts that plot should therefore
%   ask for [wave, ~, info] unless they actually want the resource grid.
%
%   Example
%     [w, g, info] = lte_generate(50, 1.0, 17, 10);
%
%   See also CAF, SYNTH_RX.

% ---------------------------------------------------------------- checks
if nargin < 4
    error('lte_generate:nargin', 'Four input arguments are required.');
end
if ~any(nrb == [6 25 50 100])
    error('lte_generate:nrb', 'nrb must be one of 6, 25, 50, 100 (got %g).', nrb);
end
if ~isscalar(rho) || rho < 0 || rho > 1
    error('lte_generate:rho', 'rho must be a scalar in [0, 1] (got %g).', rho);
end
if ~isscalar(ncellid) || ncellid < 0 || ncellid > 503 || ncellid ~= fix(ncellid)
    error('lte_generate:ncellid', 'ncellid must be an integer in [0, 503].');
end
if ~isscalar(nsubframes) || nsubframes < 1 || nsubframes ~= fix(nsubframes)
    error('lte_generate:nsubframes', 'nsubframes must be a positive integer.');
end

% ------------------------------------------------------- fixed constants
SCS        = 15000;            % subcarrier spacing, Hz
SYM_PER_SF = 14;               % OFDM symbols per 1 ms subframe, normal CP
SYM_PER_SL = 7;                % OFDM symbols per 0.5 ms slot, normal CP
NMAXRB     = 110;              % largest LTE bandwidth, used by the CRS spec

nfft   = nfft_for_nrb(nrb);
fs     = nfft * SCS;
scale  = nfft / 2048;
cp_first = 160 * scale;        % first symbol of every slot
cp_other = 144 * scale;        % remaining six symbols of every slot
cplen  = [cp_first, repmat(cp_other, 1, SYM_PER_SL - 1)];
cplen  = [cplen, cplen];       % one subframe = two slots

nsc    = 12 * nrb;             % occupied subcarriers, DC excluded
nsym   = SYM_PER_SF * nsubframes;
spsf   = SYM_PER_SF * nfft + sum(cplen);   % samples per subframe

% --------------------------------------------------------- resource grid
grid = zeros(nsc, nsym);
used = false(nsc, nsym);       % resource elements taken by CRS / PSS / SSS

half   = nsc / 2;
centre = (half - 30) : (half + 31);        % central 62 subcarriers, no DC

x1 = gold_x1(1600 + 4 * NMAXRB + 32);      % same for every symbol, hoist it
[crs_rows, crs_m] = crs_positions(nsc, nrb, ncellid, NMAXRB);

for sf = 0:(nsubframes - 1)
    sfn = mod(sf, 10);                     % subframe number in the frame
    col0 = sf * SYM_PER_SF;                % 0-based column offset

    % --- 1. cell specific reference signals -----------------------------
    for slot = 0:1
        ns = 2 * sfn + slot;               % slot number in the radio frame
        for l = [0 4]
            sym = col0 + slot * SYM_PER_SL + l + 1;    % 1-based column
            r = crs_values(ns, l, ncellid, crs_m, x1, NMAXRB);
            grid(crs_rows, sym) = r;
            used(crs_rows, sym) = true;
        end
    end

    % --- 2. primary synchronisation signal ------------------------------
    % --- 3. secondary synchronisation signal ----------------------------
    if sfn == 0 || sfn == 5
        pss_sym = col0 + SYM_PER_SL;       % 1-based: last symbol of slot 0
        sss_sym = pss_sym - 1;             % second to last symbol of slot 0

        grid(centre, pss_sym) = pss_sequence(ncellid);
        used(centre, pss_sym) = true;

        grid(centre, sss_sym) = sss_sequence(ncellid, sfn);
        used(centre, sss_sym) = true;
    end
end

% --- 4. physical downlink shared channel --------------------------------
avail   = find(~used);                     % linear indices of free REs
n_avail = numel(avail);
n_fill  = round(rho * n_avail);
if n_fill > 0
    pick  = avail(randperm(n_avail, n_fill));
    bits  = randi([0 1], 2, n_fill);
    grid(pick) = ((1 - 2 * bits(1, :)) + 1j * (1 - 2 * bits(2, :))).' / sqrt(2);
end

% ------------------------------------------------------ OFDM modulation
wave = zeros(nsubframes * spsf, 1);
lo_rows = 1:half;                          % negative frequency half
hi_rows = (half + 1):nsc;                  % positive frequency half
lo_bins = (nfft - half + 1):nfft;          % bins for offsets -half .. -1
hi_bins = 2:(half + 1);                    % bins for offsets +1 .. +half

ptr = 1;
for s = 1:nsym
    X = zeros(nfft, 1);
    X(lo_bins) = grid(lo_rows, s);
    X(hi_bins) = grid(hi_rows, s);         % X(1) is DC and stays zero

    x  = ifft(X, nfft) * sqrt(nfft);       % unit power per active subcarrier
    cp = cplen(mod(s - 1, SYM_PER_SF) + 1);
    blk = [x(end - cp + 1:end); x];        % cyclic prefix

    wave(ptr:ptr + cp + nfft - 1) = blk;
    ptr = ptr + cp + nfft;
end

% ------------------------------------------------------- normalise, info
pmean = mean(abs(wave).^2);
if pmean > 0
    rms_scale = 1 / sqrt(pmean);
else
    rms_scale = 1;                         % degenerate: empty grid
end
wave = wave * rms_scale;

info = struct();
info.nfft        = nfft;
info.fs          = fs;
info.cplen       = cplen;
info.occupied_bw = nsc * SCS;
info.nrb         = nrb;
info.rho         = rho;
info.ncellid     = ncellid;
info.nsubframes  = nsubframes;
info.nsc         = nsc;
info.scs         = SCS;
info.cp_first    = cp_first;
info.cp_other    = cp_other;
info.samples_per_subframe = spsf;
info.nsamples    = numel(wave);
info.papr_db     = 10 * log10(max(abs(wave).^2) / mean(abs(wave).^2));
info.rms_scale   = rms_scale;
info.n_pdsch_avail = n_avail;
info.n_pdsch_used  = n_fill;
if n_avail > 0
    info.rho_actual = n_fill / n_avail;
else
    info.rho_actual = 0;
end
end

% ======================================================================
% local functions
% ======================================================================

function n = nfft_for_nrb(nrb)
switch nrb
    case 6,   n = 128;
    case 25,  n = 512;
    case 50,  n = 1024;
    case 100, n = 2048;
end
end

function [rows, mvals] = crs_positions(nsc, nrb, ncellid, nmaxrb)
%CRS_POSITIONS  Grid rows carrying CRS and the matching sequence indices.
%   The comb has period 6 with offset mod(ncellid,6).  m indexes the CRS
%   within the transmitted band, m' indexes it within the 110 RB reference
%   band that the 36.211 sequence is defined over.
v     = mod(ncellid, 6);
m     = (0:(2 * nrb - 1)).';               % 2 CRS per resource block
rows  = 6 * m + v + 1;                     % 1-based grid rows
rows  = rows(rows <= nsc);
m     = m(1:numel(rows));
mvals = m + (nmaxrb - nrb);                % m' of 36.211 section 6.10.1.2
end

function r = crs_values(ns, l, ncellid, mvals, x1, nmaxrb)
%CRS_VALUES  QPSK reference symbols from the LTE length-31 Gold sequence.
%   36.211 section 6.10.1.1, normal cyclic prefix (N_CP = 1).
cinit = 2^10 * (7 * (ns + 1) + l + 1) * (2 * ncellid + 1) + 2 * ncellid + 1;
c     = gold_sequence(cinit, 4 * nmaxrb, x1);
idx   = 2 * mvals;                         % c(2m') and c(2m'+1), 0-based
r     = ((1 - 2 * c(idx + 1)) + 1j * (1 - 2 * c(idx + 2))) / sqrt(2);
r     = r(:);
end

function x1 = gold_x1(n)
%GOLD_X1  First constituent m-sequence of the LTE Gold sequence.
x1 = zeros(n, 1);
x1(1) = 1;                                 % x1(0)=1, x1(1..30)=0
for i = 1:(n - 31)
    x1(i + 31) = mod(x1(i + 3) + x1(i), 2);
end
end

function c = gold_sequence(cinit, len, x1)
%GOLD_SEQUENCE  Length-31 Gold sequence of 36.211 section 7.2.
NC = 1600;
n  = NC + len + 31;
x2 = zeros(n, 1);
b  = cinit;
for i = 1:31                               % x2(i) = bit i of cinit
    x2(i) = mod(b, 2);
    b = floor(b / 2);
end
for i = 1:(n - 31)
    x2(i + 31) = mod(x2(i + 3) + x2(i + 2) + x2(i + 1) + x2(i), 2);
end
c = mod(x1(NC + (1:len)) + x2(NC + (1:len)), 2);
end

function d = pss_sequence(ncellid)
%PSS_SEQUENCE  Length-62 Zadoff-Chu primary synchronisation signal.
roots = [25 29 34];
u = roots(mod(ncellid, 3) + 1);
n = (0:30).';
d = [exp(-1j * pi * u * n .* (n + 1) / 63); ...
     exp(-1j * pi * u * (n + 32) .* (n + 33) / 63)];
end

function d = sss_sequence(ncellid, sfn)
%SSS_SEQUENCE  Length-62 secondary synchronisation signal, 36.211 6.11.2.1.
nid1 = floor(ncellid / 3);
nid2 = mod(ncellid, 3);

qp = floor(nid1 / 30);
q  = floor((nid1 + qp * (qp + 1) / 2) / 30);
mp = nid1 + q * (q + 1) / 2;
m0 = mod(mp, 31);
m1 = mod(m0 + floor(mp / 31) + 1, 31);

st = mseq31([0 0 0 0 1], [2 0]);           % x(i+5) = x(i+2) + x(i)
ct = mseq31([0 0 0 0 1], [3 0]);           % x(i+5) = x(i+3) + x(i)
zt = mseq31([0 0 0 0 1], [4 2 1 0]);       % x(i+5) = x(i+4)+x(i+2)+x(i+1)+x(i)

n  = (0:30).';
s0 = st(mod(n + m0, 31) + 1);
s1 = st(mod(n + m1, 31) + 1);
c0 = ct(mod(n + nid2, 31) + 1);
c1 = ct(mod(n + nid2 + 3, 31) + 1);
z0 = zt(mod(n + mod(m0, 8), 31) + 1);
z1 = zt(mod(n + mod(m1, 8), 31) + 1);

d = zeros(62, 1);
if sfn == 0                                % subframe 0
    d(1:2:end) = s0 .* c0;
    d(2:2:end) = s1 .* c1 .* z0;
else                                       % subframe 5
    d(1:2:end) = s1 .* c0;
    d(2:2:end) = s0 .* c1 .* z1;
end
end

function s = mseq31(init, taps)
%MSEQ31  Length-31 binary m-sequence mapped to +-1.
x = zeros(31, 1);
x(1:5) = init(:);
for i = 1:26
    x(i + 5) = mod(sum(x(i + taps)), 2);
end
s = 1 - 2 * x;
end
