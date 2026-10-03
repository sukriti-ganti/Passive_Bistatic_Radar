function [wave, grid, info] = wave_cache(nrb, rho, ncellid, nsubframes)
%WAVE_CACHE  Generate an LTE waveform, or reload it from data/ if it is there.
%
%   [WAVE, GRID, INFO] = WAVE_CACHE(NRB, RHO, NCELLID, NSUBFRAMES)
%
%   Follows the data contract of the PRD: waveforms live in
%   data/wave_nrb{N}_rho{RRR}.mat, where RRR is 100*rho written with three
%   digits.  If the stored file was made with the same cell identity and
%   length it is reloaded, otherwise it is regenerated and overwritten.
%
%   data/ is in .gitignore.  Waveform files are large and can always be
%   regenerated from this function, so they are never committed.  Reloading
%   rather than regenerating is also what keeps the whole team inside the
%   20 hour monthly budget of the free MATLAB Online tier.
%
%   See also LTE_GENERATE.

root = fileparts(fileparts(mfilename('fullpath')));
d    = fullfile(root, 'data');
if ~exist(d, 'dir')
    mkdir(d);
end
f = fullfile(d, sprintf('wave_nrb%d_rho%03d.mat', nrb, round(rho * 100)));

if exist(f, 'file')
    s = load(f);
    if isfield(s, 'info') && s.info.ncellid == ncellid ...
            && s.info.nsubframes == nsubframes && s.info.nrb == nrb
        wave = s.wave;
        grid = s.grid;
        info = s.info;
        return
    end
end

[wave, grid, info] = lte_generate(nrb, rho, ncellid, nsubframes);
save(f, 'wave', 'grid', 'info', '-v7');
end
