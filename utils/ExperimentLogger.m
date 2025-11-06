classdef ExperimentLogger < handle
    % ExperimentLogger is a LOGGER class that prints to a text file
    %
    %   class  ExperimentLogger  is inhereted from class handle and
    %   implements a logger that prints all information to
    %   a (temporary) file in text mode.
    %
    % Methods:
    %   oTL = ExperimentLogger(filename, overwrite)  instantiate a 
    %         ExperimentLogger object
    %
    %         filename is the name of the logging file, and defaults to a
    %         temporary name if it is not given
    %
    %         overwrite denotes whether the potentially existing content 
    %         of the logging file should be overwritten, defaults to false.
    %
    %   oTL.log(str)    print datetime + log information to file.
    %   oTL.DEBUG(str)    print datetime + 'DEBUG' + log information to file.
    %   oTL.INFO(str)    print datetime + 'INFO' + log information to file.
    %   oTL.WARNING(str)    print datetime + 'WARNING' + log information to file.
    %   oTL.ERROR(str)    print datetime + 'ERROR' + log information to file.
    %
    %   close_log(oTL)  delete ExperimentLogger object and closes file
    %
    % Properties:
    %   oTL.filename  the name of the logging file
    %   oTL.level     the current line in file
    % 
    
    % +++++++++++++++++++++++++++++++++++++++++++++++++++++++
    properties (SetAccess = public, GetAccess = public)
    end
    
    % +++++++++++++++++++++++++++++++++++++++++++++++++++++++
    properties (SetAccess = public, GetAccess = public)
        filename = '';  % name of logging file
    end
    
    % +++++++++++++++++++++++++++++++++++++++++++++++++++++++
    properties (SetAccess = public, GetAccess = public)
        level = 0;
    end
    
    % +++++++++++++++++++++++++++++++++++++++++++++++++++++++
    properties (SetAccess = protected, GetAccess = public)
        filehandle = [];
    end
        % +++++++++++++++++++++++++++++++++++++++++++++++++++++++
    properties (SetAccess = protected, GetAccess = public)
        also_disp = 1;
    end

    % +++++++++++++++++++++++++++++++++++++++++++++++++++++++
    methods        
        function obj = ExperimentLogger(filename, overwrite, also_disp)
            if nargin < 1
                filename = tempname();
            end
            if nargin < 2 || isempty(overwrite)
                overwrite = 1;
            end
            if nargin >= 3 && ~isempty(also_disp)
                obj.also_disp = also_disp;
            end
            
            obj.filename = [filename, '.txt'];
            if ~overwrite 
                obj.filehandle = fopen(obj.filename, 'at');
            else
                obj.filehandle = fopen(obj.filename, 'wt');
            end
            obj.heading()
        end
        function heading(obj)
            fprintf(obj.filehandle, '%s\n', '####################################################');
            fprintf(obj.filehandle, '%s\n', '############      Experiment Logger     ############');
            fprintf(obj.filehandle, '%s\n', '####################################################');
            fprintf(obj.filehandle, '%s\n', ' ');
            fprintf(obj.filehandle, '%s\n', ' ');
            obj.level = obj.level + 5;
        end
        function log(obj, str, also_disp)
            if nargin < 3
                also_disp = obj.also_disp;
            end

            dated_str = [datestr(datetime) '  ' str];
            fprintf(obj.filehandle, '[+] %s\n', dated_str);
            if also_disp
                disp(['[+] ' dated_str]);
            end
            obj.level = obj.level + 1;
        end
        function DEBUG(obj, str)
            new_str = ['DEBUG' str];
            obj.log(new_str)    
        end
        function INFO(obj, str)
            new_str = ['INFO' str];
            obj.log(new_str)    
        end
        function WARNING(obj, str)
            new_str = ['WARNING' str];
            obj.log(new_str)    
        end
        function ERROR(obj, str)
            new_str = ['ERROR' str];
            obj.log(new_str)    
        end
        function close_log(obj)
            fclose(obj.filehandle);
        end
        
    end
    
end
