% Scale colour difference similarity scores

if ~exist( 'mos_scale', 'file' )
    addpath( fullfile( pwd, '../../../pwcmp/') );
end

T = readtable( 'subjective_color_similarity.csv' );

for kk=1:height(T)
    T.condition_id{kk} = strcat( num2str(T.display_pair(kk)), '-', num2str(T.color_id(kk)) );
end

[phi_rec, delta_rec, v_rec, w_rec, OBSs, CONDs] = scale_inter_intra( T.observer_id, T.condition_id, T.similarity, "no_delta", true );

%%
Ts = grpstats( T, "condition_id", { @mean, @std }, 'DataVars', "similarity" );

st = struct();
st.mle_mean_similarity = phi_rec';
st.mle_std_similarity = sqrt(w_rec)';
st.condition_id = CONDs;

ST = struct2table(st);

writetable( ST, "color_similarities_scaled.csv" );

J = join(ST, Ts, "Keys","condition_id");

%%
clf;
subplot( 1, 2, 1 );
plot( [0 6], [0 6], '--k' );
hold on
scatter( J.mean_similarity, J.mle_mean_similarity );
hold off
xlabel( 'Mean similarity' )
ylabel( 'MLE similarity' )

subplot( 1, 2, 2 );
plot( [0 1.5], [0 1.5], '--k' );
hold on
scatter( J.std_similarity, J.mle_std_similarity );
xlabel( 'Std similarity' )
ylabel( 'MLE std similarity' )
