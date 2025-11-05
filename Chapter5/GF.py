Genomic Data Analysis Script
Author: Shaul Sapielkin
Date: Thu Jan 30 2025

This script aggregates various analysis functions for genomic and recombination data:
- Merge genomic and recombination datasets based on marker intervals
- Correlation analysis between genomic features and recombination parameters
- Interaction analysis using machine learning models
- Differential effect analysis between mutant and wild types
- Genome context profiling via clustering
- Comprehensive report combining multiple analyses

Usage:
- Load your data into variables `genomic_df`, `mutant_df`, and `wild_df`.
- Run this script in Spyder or other IDE.
- The code includes all functions; call `main()` to execute the full analysis.

Dependencies:
- pandas, numpy, scipy, sklearn, re
"""

import pandas as pd
import numpy as np
import re
from scipy import stats
from sklearn.ensemble import RandomForestRegressor
from sklearn.linear_model import LinearRegression
from sklearn.metrics import r2_score
from sklearn.cluster import KMeans
from sklearn.preprocessing import StandardScaler

def parse_marker(marker):
    """
    Parses a marker string in the format: chromosome number + genome letter + position.
    Example: '3A45' returns (3, 'A', 45)
    """
    match = re.match(r'(\d+)([A-Za-z]+)(\d+)', marker)
    if match:
        chrom, genome, pos = match.groups()
        return int(chrom), genome, int(pos)
    return None

def is_interval_match(interval_name, markers_trio):
    """
    Checks if a trio of markers corresponds to a genomic interval.
    """
    if not isinstance(markers_trio, str):
        return False
    trio = markers_trio.split(' ')
    if len(trio) != 3:
        return False
    interval_start, interval_end = interval_name.split('_')
    start_parse = parse_marker(interval_start)
    end_parse = parse_marker(interval_end)
    trio_parses = [parse_marker(m) for m in trio]
    if not all([start_parse, end_parse] + trio_parses):
        return False
    # Check if chromosomes and genomes match
    if not all(p[0] == start_parse[0] and p[1] == start_parse[1] for p in trio_parses):
        return False
    start_pos, end_pos = start_parse[2], end_parse[2]
    trio_positions = [p[2] for p in trio_parses]
    # Check if trio markers are within interval bounds
    return min(start_pos, end_pos) <= min(trio_positions) and max(trio_positions) <= max(start_pos, end_pos)

def merge_genomic_and_recomb_data(genomic_df, mutant_df, wild_df):
    """
    Merge genomic data with mutant and wild type data based on interval matches.
    """
    merged_rows = []

    for _, genomic_row in genomic_df.iterrows():
        int_name = genomic_row['int_name']
        mutant_matches = mutant_df[mutant_df['markers'].apply(lambda x: is_interval_match(int_name, x))]
        wild_matches = wild_df[wild_df['markers'].apply(lambda x: is_interval_match(int_name, x))]

        if not mutant_matches.empty and not wild_matches.empty:
            for _, mut_row in mutant_matches.iterrows():
                for _, wild_row in wild_matches.iterrows():
                    merged_row = {
                        'int_name': int_name,
                        'markers': mut_row['markers'],
                        # Genomic features
                        **{f"{k}_gen": v for k, v in genomic_row.items() if k != 'int_name'},
                        # Mutant features
                        **{f"{k}_mut": v for k, v in mut_row.items() if k != 'markers'},
                        # Wild features
                        **{f"{k}_wild": v for k, v in wild_row.items() if k != 'markers'}
                    }
                    merged_rows.append(merged_row)

    final_df = pd.DataFrame(merged_rows)
    print(f"Found {len(merged_rows)} matched intervals.")
    return final_df

def analyze_correlations(merged_data):
    """
    Analyze Pearson correlations among genomic features and recombination parameters.
    Returns a dictionary of correlation DataFrames.
    """
    # Define features
    def get_numeric_features(df, suffix=''):
        return [col for col in df.columns if col.endswith(suffix) and pd.api.types.is_numeric_dtype(df[col])]

    # Detect features
    genomic_features = get_numeric_features(merged_data, '_gen')
    mut_features = [f + '_mut' for f in ['interference', 'coincidence', 'r1', 'r2']]
    wild_features = [f + '_wild' for f in ['interference', 'coincidence', 'r1', 'r2']]

    correlations = {
        'genomic_mut': pd.DataFrame(index=genomic_features, columns=mut_features),
        'genomic_wild': pd.DataFrame(index=genomic_features, columns=wild_features),
        'mut_vs_wild': pd.DataFrame(index=mut_features, columns=wild_features)
    }

    # Compute correlations
    for df_key, feature_rows in correlations.items():
        for feat1 in feature_rows.index:
            for feat2 in feature_rows.columns:
                try:
                    data1 = pd.to_numeric(merged_data[feat1], errors='coerce')
                    data2 = pd.to_numeric(merged_data[feat2], errors='coerce')
                    mask = ~(data1.isna() | data2.isna())
                    if mask.sum() > 1:
                        corr, p_value = stats.pearsonr(data1[mask], data2[mask])
                        correlations[df_key].loc[feat1, feat2] = f"{corr:.3f} (p={p_value:.3e})"
                except:
                    correlations[df_key].loc[feat1, feat2] = None

    return correlations

def save_correlations_to_excel(correlations, filename='correlations.xlsx'):
    """
    Save the correlation DataFrames in an Excel file with separate sheets.
    """
    with pd.ExcelWriter(filename) as writer:
        for sheet_name, df in correlations.items():
            df.to_excel(writer, sheet_name=sheet_name[:31])  # Excel sheet name limit

def analyze_feature_interactions(merged_data):
    """
    Analyze feature interactions using Random Forest and linear regression.
    """
    # Define features
    genomic_features = [col for col in merged_data.columns if col.endswith('_gen')]
    # Example target variables (modify as appropriate)
    target_vars = ['interference_mut', 'coincidence_mut', 'r1_mut', 'r2_mut']
    results = {}

    for target in target_vars:
        if target not in merged_data.columns:
            continue
        # Prepare data
        X = merged_data[genomic_features].apply(pd.to_numeric, errors='coerce')
        y = pd.to_numeric(merged_data[target], errors='coerce')
        mask = ~(X.isna().any(axis=1) | y.isna())
        X_clean = X[mask]
        y_clean = y[mask]

        if len(X_clean) < 10:
            print(f"Not enough data for {target}")
            continue

        # Random Forest for feature importance
        rf = RandomForestRegressor(n_estimators=100, random_state=42)
        rf.fit(X_clean, y_clean)
        importance_df = pd.DataFrame({'feature': genomic_features, 'importance': rf.feature_importances_})
        importance_df = importance_df.sort_values('importance', ascending=False)
        print(f"\nTop important features for {target}")
        print(importance_df.head(10))

        # Interaction analysis with top features
        top_features = importance_df.head(5)['feature'].tolist()
        interaction_data = X_clean[top_features].copy()

        # Generate pairwise interaction terms
        for i, feat1 in enumerate(top_features):
            for feat2 in top_features[i+1:]:
                interaction_name = f"{feat1}_x_{feat2}"
                interaction_data[interaction_name] = X_clean[feat1] * X_clean[feat2]

        # Regression with interaction terms
        scaler = StandardScaler()
        X_scaled = scaler.fit_transform(interaction_data)
        reg = LinearRegression()
        reg.fit(X_scaled, y_clean)

        r2 = r2_score(y_clean, reg.predict(X_scaled))
        print(f"R^2 with interactions for {target}: {r2:.3f}")

        coef_df = pd.DataFrame({
            'feature': list(interaction_data.columns),
            'coefficient': reg.coef_,
        }).sort_values('abs_coefficient', ascending=False)
        results[target] = {
            'importance': importance_df,
            'interaction_coefficients': coef_df,
            'r2': r2
        }
    return results

def analyze_differential_response(merged_data):
    """
    Analyze how mutations affect the response of features compared to wild type.
    """
    diff_results = {}
    comparison_pairs = [
        ('r1_hald_mut', 'r1_hald_wild'),
        ('r2_hald_mut', 'r2_hald_wild'),
        ('c_null_mut', 'c_null_wild'),
        ('c_free_mut', 'c_free_wild')
    ]

    for mut_var, wild_var in comparison_pairs:
        if mut_var not in merged_data.columns or wild_var not in merged_data.columns:
            continue
        mut_data = pd.to_numeric(merged_data[mut_var], errors='coerce')
        wild_data = pd.to_numeric(merged_data[wild_var], errors='coerce')
        delta = mut_data - wild_data
        merged_data[f'delta_{mut_var}'] = delta
        # Try correlations with genomic features
        correlations_list = []
        genomic_features = [col for col in merged_data.columns if col.endswith('_gen')]
        for feature in genomic_features:
            feature_data = pd.to_numeric(merged_data[feature], errors='coerce')
            mask = ~(feature_data.isna() | delta.isna())
            if mask.sum() > 5:
                corr, p_value = stats.pearsonr(feature_data[mask], delta[mask])
                correlations_list.append({'feature': feature, 'corr': corr, 'p': p_value})
        # Save significant predictors
        significant = [item for item in correlations_list if abs(item['corr']) > 0.3 and item['p'] < 0.05]
        print(f"\nDifferential response analysis for {mut_var} vs {wild_var}")
        print(f"Top predictors: {significant[:10]}")
        diff_results[f'{mut_var}_vs_{wild_var}'] = {
            'correlations': correlations_list,
            'significant_predictors': significant
        }
    return diff_results

def analyze_genomic_context_profiles(merged_data):
    """
    Cluster genomic intervals based on key features to identify contexts.
    """
    key_features = [
        'DNaseI_av_gen', 'MNase_score_gen',
        'L_TE_cs_gen', 'L_TE_rel_cs_gen',
        'N_HC_sv_gen', 'L_HC_sv_gen',
        'H3K4me1_av_gen', 'H3K36me3_av_gen'
    ]
    available_features = [f for f in key_features if f in merged_data.columns]
    if len(available_features) < 3:
        print("Not enough features for clustering.")
        return None
    X = merged_data[available_features].apply(pd.to_numeric, errors='coerce')
    mask = ~X.isna().any(axis=1)
    X_clean = X[mask]
    if len(X_clean) < 20:
        print("Not enough data points for clustering.")
        return None
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X_clean)
    n_clusters = 5
    kmeans = KMeans(n_clusters=n_clusters, random_state=42)
    clusters = kmeans.fit_predict(X_scaled)
    merged_clean = merged_data[mask].copy()
    merged_clean['genomic_cluster'] = clusters
    print(f"Generated {n_clusters} clusters")
    # Profile clusters
    cluster_profiles = []
    for cluster_id in np.unique(clusters):
        sel = merged_clean[merged_clean['genomic_cluster'] == cluster_id]
        profile = {
            'cluster_id': cluster_id,
            'size': len(sel)
        }
        for feature in available_features:
            profile[f'mean_{feature}'] = sel[feature].mean()
        # Optional: profile other features like interference etc.
        cluster_profiles.append(profile)
    profile_df = pd.DataFrame(cluster_profiles)
    print("\nGenomic context profiles:")
    print(profile_df)
    # Find clusters with max/min interference
    if 'mean_interference_mut' in profile_df.columns:
        max_c = profile_df.loc[profile_df['mean_interference_mut'].idxmax()]
        min_c = profile_df.loc[profile_df['mean_interference_mut'].idxmin()]
        print(f"Max interference in cluster {max_c['cluster_id']}")
        print(f"Min interference in cluster {min_c['cluster_id']}")
    return {
        'cluster_profiles': profile_df,
        'clustered_data': merged_clean,
        'features': available_features
    }

def save_comprehensive_results(results, filename='comprehensive_analysis.xlsx'):
    """
    Save all analysis results into an Excel file.
    """
    with pd.ExcelWriter(filename) as writer:
        # Save feature importance and interaction coefficients
        if 'interactions' in results:
            for target, res in results['interactions'].items():
                if 'importance' in res:
                    res['importance'].to_excel(writer, sheet_name=f'{target}_importance')
                if 'interaction_coefficients' in res:
                    res['interaction_coefficients'].to_excel(writer, sheet_name=f'{target}_interactions')
        # Save differential analysis
        if 'differential' in results:
            for comp, res in results['differential'].items():
                if 'significant_predictors' in res:
                    df = pd.DataFrame(res['significant_predictors'])
                    df.to_excel(writer, sheet_name=comp[:31])
        # Save genomic context profiles
        if 'contexts' in results and results['contexts']:
            results['contexts']['cluster_profiles'].to_excel(writer, sheet_name='genomic_contexts')

def main():
    """
    Main execution function.
    Load data and run various analyses.
    """
    # Load your data here:
    # Example if data stored in Data.spydata, load variables accordingly.
    # Alternatively, load from CSV/Excel files if available.
    # For example:
    # import scipy.io
    # data = scipy.io.loadmat('Data.spydata')
    # Replace with your actual data loading procedure.

    # Assuming variables: genomic_df, mutant_df, wild_df are in your workspace.
    # For this example, if in Spyder environment, ensure these are loaded previously.

    # Example:
    # from your_data_module import genomic_df, mutant_df, wild_df
    # Or load from files:
    # genomic_df = pd.read_csv('genomic_data.csv')
    # mutant_df = pd.read_csv('mutant_data.csv')
    # wild_df = pd.read_csv('wild_data.csv')

    # For demonstration, assuming variables already loaded in workspace

    # 1. Merge datasets
    merged_data = merge_genomic_and_recomb_data(genomic_df, mutant_df, wild_df)

    # 2. Correlation analysis
    correlations = analyze_correlations(merged_data)
    save_correlations_to_excel(correlations)

    # 3. Feature interactions
    interaction_results = analyze_feature_interactions(merged_data)
    print("\nFeature Interaction Analysis Complete.")

    # 4. Differential response analysis
    differential_results = analyze_differential_response(merged_data)
    print("\nDifferential Response Analysis Complete.")

    # 5. Genome context profiling
    context_results = analyze_genomic_context_profiles(merged_data)
    print("\nGenome Context Profiling Complete.")

    # 6. Comprehensive report (if needed)
    results_grid = {
        'interactions': interaction_results,
        'differential': differential_results,
        'contexts': context_results
    }
    save_comprehensive_results(results_grid)

    print("Analysis complete!")

# Execute main if run as script
if __name__ == "__main__":
    main()
