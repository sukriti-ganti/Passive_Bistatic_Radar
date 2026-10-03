function out = save_fig(h, name)
%SAVE_FIG  Export a figure to figures/<name>.pdf.
%
%   OUT = SAVE_FIG(H, NAME) writes figure H to figures/NAME.pdf and returns
%   the path.  PDF, not PNG: the report is typeset in LaTeX and a vector
%   figure stays sharp at any zoom, which matters for a range-Doppler map
%   with a 60 dB dynamic range.
%
%   The '-bestfit' option of PRINT only exists in MATLAB, so it is applied
%   inside a try and the call falls back to a plain '-dpdf' elsewhere.

root = fileparts(fileparts(mfilename('fullpath')));
d    = fullfile(root, 'figures');
if ~exist(d, 'dir')
    mkdir(d);
end
out = fullfile(d, [name '.pdf']);

set(h, 'PaperOrientation', 'landscape');
try
    print(h, out, '-dpdf', '-bestfit');
catch
    print(h, out, '-dpdf');
end
fprintf('  figure -> %s\n', out);
end
