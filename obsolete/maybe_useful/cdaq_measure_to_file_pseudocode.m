% This is a "code snippet" for saving from the daq to a file. The idea is that the results
% are saved to the working DorB_operation directory, read out into workspace and then into
% the meas struct, and then the file is deleted so it won't get copied with the
% "save_everything" command. This will need to be adapted a bit to handle a read that's
% broken into multiple files. 
% Some definitions in the cdaq objects are fairly inconsistent (durations in seconds or
% minutes), or just plain stupid (listenered not being defined as property of the object).

cdaq_filename = 'cdaq_meas.mat';

cDAQ.setDuration(duration_in_sec + 2);
lh = cDAQ.measureToFile(cdaq_filename,duration_in_min + 2/60,false);

% Have the experiment run
...
Something.turn.on
pause(3)
Something.turn.off
pause(duration_in_sec - 10)

% After the experiment is finished, read file. A bit of a brute force approach to finding
% relevant delay until file is read
tries = 0; success = 0;
while tries < 10 && ~success
    try
        cdaq_data = load(cdaq_filename);
        success = 1;
    end
    if ~success
        tries = tries + 1;
        pause(1);
    end
end
if ~success; error('Didn''t get cdaq file'); end

% Clear listener and reset daq (obviously the listener needs to be a property of the
% object and this needs to be written in a smarter way)
cDAQ.clearListener(lh);
get_daq;

% Get rid of the file (we don't want it always there)
delete(cdaq_filename);
