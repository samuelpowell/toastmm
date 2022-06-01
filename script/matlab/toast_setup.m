% This script configures the toast installation by adding the appropriate
% directories to the MATLAB path.

function toast_setup()

    
    fprintf(1,'\nConfiguring search paths for the toast toolbox.\n\n')
    
    % Get the installation root
    [toastdir,name,ext] = fileparts(which('toast_setup.m'));
    
    % Remove current toast references from path
    fprintf('Searching for existing toast path entries ...\n');
    p = path;
    p(find(p==' ')) = '^';
    p(find(p==pathsep)) = ' ';
    p = textscan(p,'%s');
    p = p{1};
    for i = 1:size(p,1)        % restore spaces
        p{i}(find(p{i}=='^')) = ' ';
    end
    
    k = strfind(p,'toast');
    nfound = 0;
    for i=1:length(k)
        if length(k{i}) > 0
            fprintf('Removing search path %s\n', p{i});
            rmpath(p{i});
            nfound = nfound+1;
        end
    end
    if nfound > 0
        fprintf('Removed %d existing toast path entries\n', nfound);
    else
        fprintf('No existing toast paths found.\n');
    end
    
    fprintf ('\nToast root directory: %s\n', toastdir);
    
    % Sanity checks: assert valid file structure
    assertfile(fullfile(toastdir, 'toast_setup.m'));
    assertdir (fullfile(toastdir, 'toast'));
    assertdir (fullfile(toastdir, 'utilities'));
    assertdir (fullfile(toastdir, 'demos'));
    assertdir (fullfile(toastdir, 'html'));
    assertfile(fullfile(toastdir, 'toast', ['toastmex.' mexext]));
    
    % Add all directories under the script/matlab node   
    toastdir_core = fullfile(toastdir, 'toast');
    toastdir_util = fullfile(toastdir, 'utilities');
    toastdir_demo = fullfile(toastdir, 'demos');
    toast_paths = {toastdir_core, toastdir_util, toastdir_demo};
    
    for i = 1:length(toast_paths)   
        disp(['Adding search path ' toast_paths{i}]);
        addpath(toast_paths{i});
    end
    
    % Report success
    if usejava('desktop')
        fprintf(1,'\nToast toolbox configured, save the path to make available permanently\n');
        pathtool
    else
        fprintf(1,'\nToast toolbox configured, use savepath to store configuration\n');
    end
    
      
    function assertfile(pathstr)
    if exist(pathstr,'file') ~= 2
        error('\nToast toolbox file structure mismatch. Expected file not found:\n%s\n\nDetected toast root folder:   %s\n\nPlease check your toast installation.', pathstr, toastdir);
    end
    end
    
    function assertdir(pathstr)
    if exist(pathstr,'dir') ~= 7
        error('\nToast toolbox file structure mismatch. Expected directory not found:\n%s\n\nDetected toast root folder:   %s\n\nPlease check your toast installation.', pathstr, toastdir);
    end
    end
    
    end