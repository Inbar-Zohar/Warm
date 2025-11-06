function srPSD = NMRG_phase_freq_srPSD(f, srSphi, T2, varargin)

srPSD = abs(NMRG_phase_freq_transfer_func(2*pi*1i*f, T2, varargin{:})) * 3600 * srSphi;
