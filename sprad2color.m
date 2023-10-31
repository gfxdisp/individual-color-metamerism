classdef sprad2color
    % Convert spectral radiance to trichromatic values, such as CIE XYZ
    % 
    % All data files are form http://www.cvrl.org/
    
    properties
        lambda = [];
        CMF = [];
        weights = [];
    end
    
    methods
        function sc = sprad2color( col_match_func )
            
            switch col_match_func
                case { 'xyz1931', 'xyz31', 'XYZ1931', 'XYZ31' } % 2 deg CIE 1931 colour matching functions
                    fname = 'ciexyz31.csv';
                    sc.weights = [1 1 1] * 683.002;
                case { 'xyz1964', 'xyz64', 'XYZ1964', 'XYZ64' } % 10 Deg CIE 1964 colour matching functions
                    fname = 'ciexyz64.csv';
                    sc.weights = [1 1 1] * 683.002;
                case { 'xyz2006', 'xyz2006-2deg' } % 2 deg CIE 2006 (proposed) colour matching functions
                    fname = 'lin2012xyz2e_1_7sf.csv';
                    sc.weights = [1 1 1] * 683.002;
                case { 'cie2006lms', 'lms' }
                    fname = 'linss2_10e_1.csv';                    
                    sc.weights = [0.689903 0.348322 0.0371597] * 683.002;
                case { 'smith&pokorny', 'sp' } % Cone fundamentals from http://www.cvrl.org/
                    fname = 'cone_smith_pokorny_1975.csv';                    
                    sc.weights = [1 1 1] * 683.002;
                case { 'lum2006', 'vl2006' } % CIE 2006 Luminous efficiency function
                    fname = 'vl_linCIE2008v2e_1.csv';                    
                    sc.weights = 683.002;
                case { 'lum1924', 'vl1924' } % CIE 1924 Luminous efficiency function
                    fname = 'vl1924e_1.csv';                    
                    sc.weights = 683.002; 
                case { 'sclum1951', 'sclum' } % CIE 1951 Scotopic Luminous efficiency function
                    fname = 'scvle.csv';                    
                    sc.weights = 1700; 
                otherwise
                    error( 'Unknown color matching functions' );
            end
            
            D = dlmread( fname );
            sc.lambda = D(:,1);
            sc.CMF = D(:,2:end);            
        end
        
        function col = get_color( sc, lambda, L )
            % lambda is a vector with wavelengths in nm
            % L is a vector with spectral radiance
            
            assert( all(size(lambda)==size(L)) ); % both vectors must be the same length
            
            % Resample CMF to match the measurement sampling 
            CMF_resampl = zeros(numel(L),size(sc.CMF,2));
            for cc=1:size(sc.CMF,2)
                CMF_resampl(:,cc) = interp1( sc.lambda, sc.CMF(:,cc), lambda, 'linear', 0 ); % set values outside the range to 0
            end

            col = zeros(1,size(sc.CMF,2));
            for cc=1:size(sc.CMF,2)
                col(cc) = sc.weights(cc) * trapz( lambda, CMF_resampl(:,cc).*L(:) );
            end
        end
        
        function plot( sc, line_style )
            
            if ~exist( 'line_style', 'var' )
                line_style = '-';
            end
            
            COLORs = { 'r', 'g', 'b' };
            for cc=1:size(sc.CMF,2)
                plot( sc.lambda, sc.CMF(:,cc)*sc.weights(cc), strcat( line_style, COLORs{cc} ) );
                hold on
            end
            
        end
        
        function plot_normalized( sc, line_style, COLORs )
            
            if ~exist( 'line_style', 'var' )
                line_style = '-';
            end
            
            if ~exist( 'COLORs', 'var')
                COLORs = { 'r', 'g', 'b' };
            end
            
            for cc=1:size(sc.CMF,2)
                plot( sc.lambda, sc.CMF(:,cc)./max(sc.CMF(:,cc)),... 
                    strcat( line_style, COLORs{cc} ) );
                hold on
            end
            
        end
        
    end    
    
end