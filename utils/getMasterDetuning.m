function detuningGHz = getMasterDetuning(isotope, F, Fp)
    % getMasterDetuning Returns detuning (GHz) of an Rb D1 transition relative to Rb-85 F=3→F'=3.
    % 
    % detuning = getMasterDetuning(isotope, F, Fp) computes the detuning for 
    % the given isotope (85 or 87), ground state F, and excited state F'. 
    % Fp can be numeric (e.g. 2) or a string like '2/3' for crossovers.
    % Crossovers are the average of the two involved transitions.
    %
    % Example: getMasterDetuning(87, 1, '1/2') → 5.385
    
    
    % Validate isotope
    if ~ismember(isotope, [85, 87])
        error('Isotope must be 85 or 87.');
    end
    
    % Define transitions [isotope, F, F', detuning (GHz)]
    % Detunings relative to Rb-85 F=3 → F'=3 = 0.000 GHz
    transitions = [
        85, 3, 3,  0;
        85, 3, 2, -0.36;
        85, 2, 3, 3.04;
        85, 2, 2, 2.68;
        87, 2, 2, -1.04;
        87, 2, 1, -1.85;
        87, 1, 2, 5.80;
        87, 1, 1, 4.97
    ];
    
    % Parse Fp input
    if ischar(Fp) || isstring(Fp)
        crossoverStr = strrep(char(Fp), ' ', '');
        parts = split(crossoverStr, '/');
        if length(parts) ~= 2
            error('Invalid crossover format. Use e.g. "2/3" or "2 / 3".');
        end
        Fp1 = str2double(parts{1});
        Fp2 = str2double(parts{2});
    
        if isnan(Fp1) || isnan(Fp2)
            error('Crossover string must contain valid numbers.');
        end
    
        % Find the two transitions to average
        rows = (transitions(:,1) == isotope) & (transitions(:,2) == F) & ...
               ((transitions(:,3) == Fp1) | (transitions(:,3) == Fp2));
    
        if sum(rows) ~= 2
            error('Crossover transition not found for isotope %d, F=%d between F''=%d and F''=%d.', isotope, F, Fp1, Fp2);
        end
    
        detuningGHz = mean(transitions(rows, 4));
    else
        % Standard transition
        if ~isscalar(Fp) || ~ismember(Fp, [1, 2, 3])
            error('F'' must be 1, 2, or 3 for standard transitions.');
        end
    
        row = (transitions(:,1) == isotope) & (transitions(:,2) == F) & (transitions(:,3) == Fp);
    
        if ~any(row)
            error('Transition not found for isotope %d, F=%d → F''=%d.', isotope, F, Fp);
        end
    
        detuningGHz = transitions(row, 4);
    end
end