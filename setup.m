function setup(force)
%SETUP Download and configure external dependencies for the PHIMV experiments.
%
%   setup
%   setup(true)
%
% External packages are downloaded into the local "external" directory.
% Only the directory containing each package's main routine is added to the
% MATLAB path.  In particular, this avoids recursively adding bundled copies
% of other solvers (for example, the copy of KIOPS contained in BAMPHI's
% numerical_experiments directory).
%
% The optional input FORCE=true removes and redownloads the external
% packages.
%
% Installed packages:
%   - phi_funm  (Al-Mohy and Liu)
%   - BAMPHI    (Caliari, Cassini, and Zivcovich)
%   - KIOPS     (Gaudreault, Rainwater, and Tokman)
%   - PHIPM     (Niesen and Wright, ACM Algorithm 919)
%   - Anymatrix (Higham and Mikaitis)
%   - matrices-expm collection (Xiaobo Liu), installed as the Anymatrix
%     group expm
%
% MATLAB R2022b or later is recommended.
%
% Some high-precision reference computations additionally require the
% Advanpix Multiprecision Computing Toolbox.  This is commercial software
% and is therefore not installed by this script.

if nargin < 1
    force = false;
end

if ~(islogical(force) || (isnumeric(force) && isscalar(force)))
    error('The optional input FORCE must be a logical scalar.');
end
force = logical(force);

rootdir     = fileparts(mfilename('fullpath'));
externaldir = fullfile(rootdir,'external');

if ~exist(externaldir,'dir')
    mkdir(externaldir);
end

fprintf('\nPHIMV dependency setup\n');
fprintf('Repository root : %s\n',rootdir);
fprintf('External folder : %s\n\n',externaldir);

% -------------------------------------------------------------------------
% Package definitions
% -------------------------------------------------------------------------

packages = struct( ...
    'name', { ...
        'phi_funm', ...
        'bamphi', ...
        'kiops', ...
        'phipm', ...
        'anymatrix'}, ...
    'url', { ...
        'https://github.com/xiaobo-liu/phi_funm/archive/refs/heads/main.zip', ...
        'https://github.com/francozivcovich/bamphi/archive/refs/heads/main.zip', ...
        'https://gitlab.com/stephane.gaudreault/kiops/-/archive/master/kiops-master.zip', ...
        'https://calgo.acm.org/919.zip', ...
        'https://github.com/north-numerical-computing/anymatrix/archive/refs/heads/master.zip'}, ...
    'marker', { ...
        'phi_funm.m', ...
        'bamphi.m', ...
        'kiops.m', ...
        'phipm.m', ...
        'anymatrix.m'} );

% Remove any previously added paths lying inside external/.  This is
% important if an older setup used ADDPATH(GENPATH(externaldir)), since
% that may have put BAMPHI's bundled KIOPS ahead of the standalone KIOPS.
remove_external_paths(externaldir);

% Download/install each package.
for k = 1:numel(packages)
    install_package(packages(k),externaldir,force);
end

% -------------------------------------------------------------------------
% Install Xiaobo Liu's matrices-expm collection inside Anymatrix.
% -------------------------------------------------------------------------
%
% The matrices-expm README prescribes the following Anymatrix integration:
%   1. download the matrices-expm repository;
%   2. rename the downloaded collection directory to "expm";
%   3. remove LICENSE and README files from that collection directory;
%   4. place the resulting "expm" directory at the root of Anymatrix;
%   5. run "anymatrix scan".
%
% We perform steps 1--4 here. The scan is performed below after Anymatrix
% itself has been added to the MATLAB path.

anymatrix_pkgdir = fullfile(externaldir,'anymatrix');
anymatrix_marker = dir(fullfile(anymatrix_pkgdir,'**','anymatrix.m'));

if isempty(anymatrix_marker)
    error('Could not locate anymatrix.m after installing Anymatrix.');
end

if numel(anymatrix_marker) > 1
    paths = arrayfun(@(d) fullfile(d.folder,d.name), ...
                     anymatrix_marker,'UniformOutput',false);
    depths = cellfun(@path_depth,paths);
    [~,idx] = min(depths);
    anymatrix_marker = anymatrix_marker(idx);
end

anymatrix_root = anymatrix_marker.folder;

install_matrices_expm(anymatrix_root,force);

% -------------------------------------------------------------------------
% Add ONLY the intended package directories.
% -------------------------------------------------------------------------

package_dirs = cell(size(packages));

fprintf('\nConfiguring MATLAB path:\n');

for k = 1:numel(packages)

    pkgdir = fullfile(externaldir,packages(k).name);

    marker = dir(fullfile(pkgdir,'**',packages(k).marker));

    if isempty(marker)
        error('Could not locate %s after installing %s.', ...
              packages(k).marker,packages(k).name);
    end

    if numel(marker) > 1
        % Prefer the shallowest occurrence of the main routine.
        paths = arrayfun(@(d) fullfile(d.folder,d.name), ...
                         marker,'UniformOutput',false);
        depths = cellfun(@path_depth,paths);
        [~,idx] = min(depths);
        marker = marker(idx);
    end

    package_dirs{k} = marker.folder;

    % Add only the folder containing the package's main routine.
    addpath(package_dirs{k},'-begin');

    fprintf('  %-12s %s\n',packages(k).name,package_dirs{k});
end

rehash;

% Register the newly installed expm collection with Anymatrix.
fprintf('\nScanning Anymatrix collections...\n');

try
    anymatrix('scan');
catch ME
    error(['Anymatrix was installed, but the scan failed:\n%s'],ME.message);
end

expm_collection = fullfile(anymatrix_root,'expm');

if ~exist(fullfile(expm_collection,'anymatrix_expm.m'),'file')
    error(['The matrices-expm collection was not installed correctly at:\n' ...
           '  %s'],expm_collection);
end

fprintf('  expm collection: %s\n',expm_collection);

% -------------------------------------------------------------------------
% Verify that MATLAB resolves every solver to the intended package.
% -------------------------------------------------------------------------

fprintf('\nChecking installed functions:\n');

all_ok = true;

for k = 1:numel(packages)

    fname    = erase(packages(k).marker,'.m');
    resolved = which(fname);

    expected_root = canonical_path(fullfile(externaldir,packages(k).name));
    actual_file   = canonical_path(resolved);

    if isempty(resolved)
        fprintf('  %-12s NOT FOUND\n',fname);
        all_ok = false;

    elseif ~starts_with_path(actual_file,expected_root)
        fprintf('  %-12s WRONG COPY:\n      %s\n',fname,resolved);
        all_ok = false;

    else
        fprintf('  %-12s %s\n',fname,resolved);
    end
end

% Explicitly check that KIOPS is the standalone downloaded copy rather
% than BAMPHI's bundled numerical-experiment copy.
kiops_path = which('kiops');

if ~isempty(kiops_path) && contains(lower(kiops_path), ...
        [filesep 'bamphi' filesep])
    all_ok = false;
    fprintf(['\nERROR: MATLAB is resolving KIOPS through the BAMPHI tree:\n' ...
             '  %s\n'],kiops_path);
end

% -------------------------------------------------------------------------
% Optional high-precision dependency
% -------------------------------------------------------------------------

has_mp = (exist('mp','class') == 8) || ...
         (exist('mp','file') == 2) || ...
         (exist('mp','file') == 6);

fprintf('\nOptional high-precision dependency:\n');

if has_mp
    fprintf('  Advanpix Multiprecision Computing Toolbox: found.\n');
else
    fprintf(['  Advanpix Multiprecision Computing Toolbox: not found.\n' ...
             '  This commercial toolbox is not installed automatically.\n' ...
             '  It is required only for scripts that regenerate the\n' ...
             '  multiprecision reference solutions.\n']);
end

if all_ok
    fprintf('\nSetup completed successfully.\n');
    fprintf(['Only the intended package directories were added to the\n' ...
             'current MATLAB session; bundled solver copies were excluded.\n\n']);
else
    error(['Setup finished, but one or more external routines resolve to ' ...
           'an unintended copy. See the messages above.']);
end

end


% =========================================================================
% Local functions
% =========================================================================

function install_package(pkg,externaldir,force)
%INSTALL_PACKAGE Download and extract one external package.

dest = fullfile(externaldir,pkg.name);

marker = [];

if exist(dest,'dir')
    marker = dir(fullfile(dest,'**',pkg.marker));
end

if ~force && ~isempty(marker)
    fprintf('%-12s already installed; skipping.\n',pkg.name);
    return
end

if exist(dest,'dir')
    fprintf('%-12s removing existing local copy...\n',pkg.name);
    rmdir(dest,'s');
end

mkdir(dest);

fprintf('%-12s downloading...\n',pkg.name);

zipfile = [tempname,'.zip'];

try
    websave(zipfile,pkg.url);
catch ME
    clean_failed_install(zipfile,dest);
    error('Failed to download %s:\n%s',pkg.name,ME.message);
end

fprintf('%-12s extracting...\n',pkg.name);

try
    unzip(zipfile,dest);
catch ME
    clean_failed_install(zipfile,dest);
    error('Failed to extract %s:\n%s',pkg.name,ME.message);
end

if exist(zipfile,'file')
    delete(zipfile);
end

marker = dir(fullfile(dest,'**',pkg.marker));

if isempty(marker)
    error(['Package %s was downloaded, but %s was not found in the ' ...
           'archive.'],pkg.name,pkg.marker);
end

fprintf('%-12s installed.\n',pkg.name);

end


% -------------------------------------------------------------------------

function install_matrices_expm(anymatrix_root,force)
%INSTALL_MATRICES_EXPM Install Xiaobo Liu's matrices-expm as Anymatrix/expm.
%
% The upstream repository instructs users to rename the repository
% directory to "expm", remove LICENSE and README from that directory, and
% place it at the root of Anymatrix.

    url = ...
        'https://github.com/xiaobo-liu/matrices-expm/archive/refs/heads/main.zip';

    dest = fullfile(anymatrix_root,'expm');
    marker = fullfile(dest,'anymatrix_expm.m');

    if ~force && exist(marker,'file')
        fprintf('%-12s already installed in Anymatrix; skipping.\n', ...
                'matrices-expm');
        return
    end

    if exist(dest,'dir')
        fprintf('%-12s removing existing Anymatrix collection...\n', ...
                'matrices-expm');
        rmdir(dest,'s');
    end

    zipfile = [tempname,'.zip'];
    stagedir = tempname;
    mkdir(stagedir);

    fprintf('%-12s downloading...\n','matrices-expm');

    try
        websave(zipfile,url);
    catch ME
        cleanup_staging(zipfile,stagedir);
        error('Failed to download matrices-expm:\n%s',ME.message);
    end

    fprintf('%-12s extracting...\n','matrices-expm');

    try
        unzip(zipfile,stagedir);
    catch ME
        cleanup_staging(zipfile,stagedir);
        error('Failed to extract matrices-expm:\n%s',ME.message);
    end

    if exist(zipfile,'file')
        delete(zipfile);
    end

    % Locate the repository root using the Anymatrix bridge file.
    d = dir(fullfile(stagedir,'**','anymatrix_expm.m'));

    if isempty(d)
        cleanup_staging('',stagedir);
        error(['matrices-expm was downloaded, but anymatrix_expm.m ' ...
               'was not found.']);
    end

    if numel(d) > 1
        paths = arrayfun(@(x) fullfile(x.folder,x.name), ...
                         d,'UniformOutput',false);
        depths = cellfun(@path_depth,paths);
        [~,idx] = min(depths);
        d = d(idx);
    end

    source_root = d.folder;

    % Required by the matrices-expm README for Anymatrix integration.
    remove_if_exists(fullfile(source_root,'LICENSE'));
    remove_if_exists(fullfile(source_root,'LICENSE.txt'));
    remove_if_exists(fullfile(source_root,'README'));
    remove_if_exists(fullfile(source_root,'README.md'));
    remove_if_exists(fullfile(source_root,'readme'));
    remove_if_exists(fullfile(source_root,'readme.md'));

    % Moving the repository root to this destination simultaneously renames
    % the collection directory to "expm" and places it in Anymatrix.
    [ok,msg] = movefile(source_root,dest);

    if ~ok
        cleanup_staging('',stagedir);
        error('Failed to install matrices-expm as Anymatrix/expm:\n%s',msg);
    end

    if exist(stagedir,'dir')
        rmdir(stagedir,'s');
    end

    if ~exist(marker,'file')
        error(['matrices-expm installation completed, but the marker file ' ...
               'is missing:\n  %s'],marker);
    end

    fprintf('%-12s installed as Anymatrix group "expm".\n', ...
            'matrices-expm');

end


% -------------------------------------------------------------------------

function remove_if_exists(filename)
%REMOVE_IF_EXISTS Delete a file if it exists.

    if exist(filename,'file')
        delete(filename);
    end

end


% -------------------------------------------------------------------------

function cleanup_staging(zipfile,stagedir)
%CLEANUP_STAGING Remove temporary matrices-expm download files.

    if ~isempty(zipfile) && exist(zipfile,'file')
        delete(zipfile);
    end

    if ~isempty(stagedir) && exist(stagedir,'dir')
        rmdir(stagedir,'s');
    end

end


% -------------------------------------------------------------------------

function remove_external_paths(externaldir)
%REMOVE_EXTERNAL_PATHS Remove all current MATLAB path entries under externaldir.

entries = strsplit(path,pathsep);

external_canon = canonical_path(externaldir);

to_remove = {};

for k = 1:numel(entries)

    entry = entries{k};

    if isempty(entry)
        continue
    end

    entry_canon = canonical_path(entry);

    if starts_with_path(entry_canon,external_canon)
        to_remove{end+1} = entry; %#ok<AGROW>
    end
end

if ~isempty(to_remove)
    rmpath(to_remove{:});
end

end


% -------------------------------------------------------------------------

function clean_failed_install(zipfile,dest)

if exist(zipfile,'file')
    delete(zipfile);
end

if exist(dest,'dir')
    rmdir(dest,'s');
end

end


% -------------------------------------------------------------------------

function d = path_depth(filename)
%PATH_DEPTH Number of path components; used to prefer the main package file.

parts = strsplit(char(filename),filesep);
d = numel(parts);

end


% -------------------------------------------------------------------------

function p = canonical_path(p)
%CANONICAL_PATH Normalize a path for reliable prefix comparisons.

if isempty(p)
    return
end

p = char(p);

% Normalize slash direction.
p = strrep(p,'/',filesep);
p = strrep(p,'\',filesep);

% Remove a trailing file separator.
while numel(p) > 1 && p(end) == filesep
    p(end) = [];
end

% Windows paths are case-insensitive.
if ispc
    p = lower(p);
end

end


% -------------------------------------------------------------------------

function tf = starts_with_path(candidate,parent)
%STARTS_WITH_PATH True if CANDIDATE is PARENT or lies below PARENT.

if isempty(candidate) || isempty(parent)
    tf = false;
    return
end

tf = strcmp(candidate,parent) || ...
     startsWith(candidate,[parent filesep]);

end
