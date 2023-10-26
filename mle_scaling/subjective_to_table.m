% clean up subjective data

data_path = '../ObserverMetamerism/ObserverMetamerism/Data/ObsResponses/';

files = dir( fullfile( data_path, 'Corrected_*.csv' ) );

D = struct( 'observer_id', [], 'display_pair', [], 'color_id', [], 'similarity', [] );

N = 6*22;

display_pairs = [12, 13, 14, 23, 24, 34]';
color_ids = [1:11 1:11];

display_pairs_mat = repmat( display_pairs, 1, 22 );
color_ids_mat = repmat( color_ids, 6, 1 );

for kk=1:length(files)
    
    if strcmp( files(kk).name, 'Corrected_junwoo_test.csv' )
        continue;
    end

    fprintf( 1, "processing %s\n", files(kk).name );

    SS = dlmread( fullfile( data_path, files(kk).name ) );

    obs_id = files(kk).name((end-6):(end-4));

    D.observer_id = cat( 1, D.observer_id, repmat( { obs_id }, N, 1 ) );
    D.display_pair = cat( 1, D.display_pair, display_pairs_mat(:) );
    D.color_id = cat( 1, D.color_id, color_ids_mat(:) );
    D.similarity = cat( 1, D.similarity, SS(:) );

end

T = struct2table( D );
writetable( T, 'subjective_color_similarity.csv' );