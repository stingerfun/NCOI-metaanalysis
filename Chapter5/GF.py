# -*- coding: utf-8 -*-
"""
Created on Thu Jan 30 10:26:00 2025

@author: user
"""

def merge_genomic_and_recomb_data(genomic_df, mutant_df, wild_df):
    def parse_marker(marker):
        """Парсит маркер в формате: номер хромосомы + буква генома + порядковый номер"""
        import re
        match = re.match(r'(\d+)([A-Za-z]+)(\d+)', marker)
        if match:
            chrom, genome, pos = match.groups()
            return int(chrom), genome, int(pos)
        return None

    def is_interval_match(interval_name, markers_trio):
        """Проверяет, соответствует ли тройка маркеров интервалу"""
        if not isinstance(markers_trio, str):
            return False
            
        # Разбиваем тройку маркеров и интервал
        trio = markers_trio.split(' ')
        if len(trio) != 3:
            return False
            
        interval_start, interval_end = interval_name.split('_')
        
        # Парсим маркеры
        start_parse = parse_marker(interval_start)
        end_parse = parse_marker(interval_end)
        trio_parses = [parse_marker(m) for m in trio]
        
        if not all([start_parse, end_parse] + trio_parses):
            return False
            
        # Проверяем, что хромосома и геном совпадают
        if not all(p[0] == start_parse[0] and p[1] == start_parse[1] for p in trio_parses):
            return False
            
        # Получаем позиции
        start_pos = start_parse[2]
        end_pos = end_parse[2]
        trio_pos = [p[2] for p in trio_parses]
        
        # Проверяем, что тройка маркеров находится внутри интервала
        return min(start_pos, end_pos) <= min(trio_pos) and max(trio_pos) <= max(start_pos, end_pos)

    # Создаем списки для хранения результатов
    merged_rows = []
    
    # Для каждого интервала в геномных данных
    for _, genomic_row in genomic_df.iterrows():
        int_name = genomic_row['int_name']
        
        # Ищем соответствующие записи в данных мутантного и дикого типа
        mutant_matches = mutant_df[mutant_df['markers'].apply(lambda x: is_interval_match(int_name, x))]
        wild_matches = wild_df[wild_df['markers'].apply(lambda x: is_interval_match(int_name, x))]
        
        # Если нашли соответствия в обоих наборах данных
        if not mutant_matches.empty and not wild_matches.empty:
            for _, mut_row in mutant_matches.iterrows():
                for _, wild_row in wild_matches.iterrows():
                    merged_row = {
                        'int_name': int_name,
                        'markers': mut_row['markers'],
                        # Геномные характеристики
                        **{f"{k}_gen": v for k, v in genomic_row.items() if k != 'int_name'},
                        # Характеристики мутантного типа
                        **{f"{k}_mut": v for k, v in mut_row.items() if k not in ['markers']},
                        # Характеристики дикого типа
                        **{f"{k}_wild": v for k, v in wild_row.items() if k not in ['markers']}
                    }
                    merged_rows.append(merged_row)
    
    import pandas as pd
    final_df = pd.DataFrame(merged_rows)
    
    # Добавляем проверку результатов
    print(f"Найдено {len(merged_rows)} соответствий")
    
    return final_df

# Функция для проверки работы парсинга
def test_marker_parsing(genomic_df, mutant_df, wild_df):
    def parse_marker(marker):
        import re
        match = re.match(r'(\d+)([A-Za-z]+)(\d+)', marker)
        if match:
            chrom, genome, pos = match.groups()
            return int(chrom), genome, int(pos)
        return None
    
    print("Тестирование парсинга маркеров:")
    
    # Проверяем несколько интервалов
    print("\nПримеры парсинга интервалов:")
    for int_name in genomic_df['int_name'].head():
        start, end = int_name.split('_')
        print(f"{int_name}: {parse_marker(start)} - {parse_marker(end)}")
    
    # Проверяем несколько троек маркеров
    print("\nПримеры парсинга троек маркеров:")
    for markers in mutant_df['markers'].head():
        if isinstance(markers, str):
            trio = markers.split(' ')
            print(f"{markers}: {[parse_marker(m) for m in trio]}")

# Пример использования:
# merged_data = merge_genomic_and_recomb_data(genomic_df, mutant_df, wild_df)
# test_marker_parsing(genomic_df, mutant_df, wild_df)
# split_markers = data_001['markers'].str.replace(' - ', ' ').str.split()
def replace_with_order(marker_list):
    return ' '.join(data['order'][data['Marker'] == marker].values[0] for marker in marker_list)
# data_001['markers'] = split_markers.apply(replace_with_order)
test_marker_parsing(genomic_df, mutant_df, wild_df)

# Если парсинг работает корректно, объединяем данные
merged_data = merge_genomic_and_recomb_data(genomic_df, mutant_df, wild_df)

# Сохраняем результат
merged_data.to_csv('merged_genomic_recomb_data.csv', index=False)

def analyze_correlations(merged_data):
   import pandas as pd
   import numpy as np
   from scipy import stats
   
   # Получаем геномные характеристики (колонки с суффиксом _gen)
   genomic_features = [col for col in merged_data.columns if col.endswith('_gen')]
   
   # Характеристики рекомбинации и интерференции для обоих типов
   recomb_features = ['interference', 'coincidence', 'r1', 'r2']
   mut_features = [f + '_mut' for f in recomb_features]
   wild_features = [f + '_wild' for f in recomb_features]
   
   # Создаем словари для хранения корреляций
   correlations = {
       'genomic_mut': pd.DataFrame(index=genomic_features, columns=mut_features),
       'genomic_wild': pd.DataFrame(index=genomic_features, columns=wild_features),
       'mut_vs_wild': pd.DataFrame(index=mut_features, columns=wild_features)
   }
   
   # Вычисляем корреляции Пирсона и p-value для геномных характеристик с мутантным типом
   for gen in genomic_features:
       for mut in mut_features:
           corr, p_value = stats.pearsonr(merged_data[gen].astype(float), 
                                        merged_data[mut].astype(float))
           correlations['genomic_mut'].loc[gen, mut] = f"{corr:.3f} (p={p_value:.3e})"
           
   # Вычисляем корреляции для геномных характеристик с диким типом
   for gen in genomic_features:
       for wild in wild_features:
           corr, p_value = stats.pearsonr(merged_data[gen].astype(float), 
                                        merged_data[wild].astype(float))
           correlations['genomic_wild'].loc[gen, wild] = f"{corr:.3f} (p={p_value:.3e})"
           
   # Вычисляем корреляции между мутантным и диким типом
   for mut in mut_features:
       for wild in wild_features:
           corr, p_value = stats.pearsonr(merged_data[mut].astype(float), 
                                        merged_data[wild].astype(float))
           correlations['mut_vs_wild'].loc[mut, wild] = f"{corr:.3f} (p={p_value:.3e})"
   
   # Функция для выделения значимых корреляций
   def get_significant_correlations(corr_df, threshold=0.5):
       significant = []
       for idx in corr_df.index:
           for col in corr_df.columns:
               value = corr_df.loc[idx, col]
               if value:
                   corr = float(value.split()[0])
                   p_value = float(value.split('p=')[1].rstrip(')'))
                   if abs(corr) >= threshold and p_value < 0.05:
                       significant.append((idx, col, corr, p_value))
       return significant
   
   # Выводим результаты
   print("Значимые корреляции (|r| >= 0.5, p < 0.05):")
   
   print("\n1. Геномные характеристики vs Мутантный тип:")
   significant = get_significant_correlations(correlations['genomic_mut'])
   for feat1, feat2, corr, p_value in sorted(significant, key=lambda x: abs(x[2]), reverse=True):
       print(f"{feat1} vs {feat2}: r={corr:.3f} (p={p_value:.3e})")
       
   print("\n2. Геномные характеристики vs Дикий тип:")
   significant = get_significant_correlations(correlations['genomic_wild'])
   for feat1, feat2, corr, p_value in sorted(significant, key=lambda x: abs(x[2]), reverse=True):
       print(f"{feat1} vs {feat2}: r={corr:.3f} (p={p_value:.3e})")
       
   print("\n3. Мутантный vs Дикий тип:")
   significant = get_significant_correlations(correlations['mut_vs_wild'])
   for feat1, feat2, corr, p_value in sorted(significant, key=lambda x: abs(x[2]), reverse=True):
       print(f"{feat1} vs {feat2}: r={corr:.3f} (p={p_value:.3e})")
   
   return correlations

# Пример использования:
# correlations = analyze_correlations(merged_data)

# Сохранение результатов в Excel
def save_correlations_to_excel(correlations, filename='correlations.xlsx'):
   with pd.ExcelWriter(filename) as writer:
       correlations['genomic_mut'].to_excel(writer, sheet_name='Genomic_vs_Mutant')
       correlations['genomic_wild'].to_excel(writer, sheet_name='Genomic_vs_Wild')
       correlations['mut_vs_wild'].to_excel(writer, sheet_name='Mutant_vs_Wild')

# Использование:
# save_correlations_to_excel(correlations)

# Проводим анализ корреляций
correlations = analyze_correlations(merged_data)

# Сохраняем результаты
save_correlations_to_excel(correlations)


def analyze_correlations(merged_data):
    import pandas as pd
    import numpy as np
    from scipy import stats
    
    # Функция для проверки, является ли колонка числовой
    def is_numeric_column(column):
        try:
            pd.to_numeric(merged_data[column])
            return True
        except:
            return False
    
    # Получаем только числовые геномные характеристики
    genomic_features = [col for col in merged_data.columns if col.endswith('_gen') and is_numeric_column(col)]
    
    # Характеристики рекомбинации и интерференции для обоих типов
    recomb_features = ['interference', 'coincidence', 'r1', 'r2']
    mut_features = [f + '_mut' for f in recomb_features]
    wild_features = [f + '_wild' for f in recomb_features]
    
    print("Анализируемые геномные характеристики:", len(genomic_features))
    print("Первые 5 характеристик:", genomic_features[:5])
    
    # Создаем словари для хранения корреляций
    correlations = {
        'genomic_mut': pd.DataFrame(index=genomic_features, columns=mut_features),
        'genomic_wild': pd.DataFrame(index=genomic_features, columns=wild_features),
        'mut_vs_wild': pd.DataFrame(index=mut_features, columns=wild_features)
    }
    
    # Вычисляем корреляции
    for gen in genomic_features:
        gen_data = pd.to_numeric(merged_data[gen])
        for mut in mut_features:
            mut_data = pd.to_numeric(merged_data[mut])
            mask = ~(gen_data.isna() | mut_data.isna())
            if mask.sum() > 1:  # Проверяем, что есть хотя бы две пары значений
                corr, p_value = stats.pearsonr(gen_data[mask], mut_data[mask])
                correlations['genomic_mut'].loc[gen, mut] = f"{corr:.3f} (p={p_value:.3e})"
            
    for gen in genomic_features:
        gen_data = pd.to_numeric(merged_data[gen])
        for wild in wild_features:
            wild_data = pd.to_numeric(merged_data[wild])
            mask = ~(gen_data.isna() | wild_data.isna())
            if mask.sum() > 1:
                corr, p_value = stats.pearsonr(gen_data[mask], wild_data[mask])
                correlations['genomic_wild'].loc[gen, wild] = f"{corr:.3f} (p={p_value:.3e})"
            
    for mut in mut_features:
        mut_data = pd.to_numeric(merged_data[mut])
        for wild in wild_features:
            wild_data = pd.to_numeric(merged_data[wild])
            mask = ~(mut_data.isna() | wild_data.isna())
            if mask.sum() > 1:
                corr, p_value = stats.pearsonr(mut_data[mask], wild_data[mask])
                correlations['mut_vs_wild'].loc[mut, wild] = f"{corr:.3f} (p={p_value:.3e})"
    
    # Функция для выделения значимых корреляций
    def get_significant_correlations(corr_df, threshold=0.5):
        significant = []
        for idx in corr_df.index:
            for col in corr_df.columns:
                value = corr_df.loc[idx, col]
                if pd.notna(value):
                    try:
                        corr = float(value.split()[0])
                        p_value = float(value.split('p=')[1].rstrip(')'))
                        if abs(corr) >= threshold and p_value < 0.05:
                            significant.append((idx, col, corr, p_value))
                    except:
                        continue
        return significant
    
    # Выводим результаты
    print("\nЗначимые корреляции (|r| >= 0.5, p < 0.05):")
    
    print("\n1. Геномные характеристики vs Мутантный тип:")
    significant = get_significant_correlations(correlations['genomic_mut'])
    for feat1, feat2, corr, p_value in sorted(significant, key=lambda x: abs(x[2]), reverse=True):
        print(f"{feat1} vs {feat2}: r={corr:.3f} (p={p_value:.3e})")
        
    print("\n2. Геномные характеристики vs Дикий тип:")
    significant = get_significant_correlations(correlations['genomic_wild'])
    for feat1, feat2, corr, p_value in sorted(significant, key=lambda x: abs(x[2]), reverse=True):
        print(f"{feat1} vs {feat2}: r={corr:.3f} (p={p_value:.3e})")
        
    print("\n3. Мутантный vs Дикий тип:")
    significant = get_significant_correlations(correlations['mut_vs_wild'])
    for feat1, feat2, corr, p_value in sorted(significant, key=lambda x: abs(x[2]), reverse=True):
        print(f"{feat1} vs {feat2}: r={corr:.3f} (p={p_value:.3e})")
    
    return correlations

def save_correlations_to_excel(correlations, filename='correlations.xlsx'):
    with pd.ExcelWriter(filename) as writer:
        for name, df in correlations.items():
            df.to_excel(writer, sheet_name=name[:31])  # Excel ограничивает длину имени листа

# Использование:
correlations = analyze_correlations(merged_data)
save_correlations_to_excel(correlations)

import pandas as pd

def merge_consecutive_intervals(df):
    """
    Merge pairs of consecutive intervals if they share end/start points
    """
    df_result = df.copy()
    
    # Convert 'no' to 0 in specific columns
    cols_to_convert = ['W1', 'M1', 'M2', 'W2', 'av_W', 'av_M']
    for col in cols_to_convert:
        df_result[col] = df_result[col].replace('no', 0)
        df_result[col] = pd.to_numeric(df_result[col], errors='coerce')
    
    preserve_cols = ['int_name', 'genome', 'chr', 'arm', 'rec_reg_cs']
    rows_to_drop = []
    
    # Step through pairs of rows
    for i in range(0, len(df_result)-1, 2):
        current_interval = df_result.iloc[i]['int_name']
        next_interval = df_result.iloc[i+1]['int_name']
        
        current_start, current_end = current_interval.split('_')
        next_start, next_end = next_interval.split('_')
        
        # Check if intervals are consecutive
        if current_end == next_start:
            # Create new merged interval name
            new_interval = f"{current_start}_{next_end}"
            df_result.iloc[i, df_result.columns.get_loc('int_name')] = new_interval
            
            # Sum all numerical columns
            for column in df_result.columns:
                if column not in preserve_cols:
                    try:
                        value1 = df_result.iloc[i][column]
                        value2 = df_result.iloc[i+1][column]
                        
                        if pd.notnull(value1) and pd.notnull(value2):
                            df_result.iloc[i, df_result.columns.get_loc(column)] = value1 + value2
                    except:
                        continue
            
            rows_to_drop.append(i+1)
    
    # Drop all merged rows at once
    df_result = df_result.drop(df_result.index[rows_to_drop])
    
    return df_result

def merge_duplicate_intervals(df):
    """
    Merge rows that have the same interval name in the 'int_name' column.
    For numerical columns: calculate the average
    For string columns: concatenate unique values with space delimiter
    
    Args:
        df: DataFrame with 'int_name' column and various data types
    Returns:
        DataFrame with merged rows
    """
    # Function to determine aggregation method for each column
    def agg_method(col_name):
        # Get the column data
        col_data = df[col_name]
        
        # Check if all values are numeric (allowing NaN)
        numeric_mask = pd.to_numeric(col_data, errors='coerce').notna() | col_data.isna()
        all_numeric = numeric_mask.all()
        
        if all_numeric:
            return 'mean'  # Average for numeric columns
        else:
            # For string columns, concatenate unique values
            return lambda x: ' '.join(sorted(set([str(i) for i in x if pd.notna(i) and str(i) != 'nan'])))
    
    # Create a dictionary of aggregation methods for each column
    agg_dict = {col: agg_method(col) for col in df.columns if col != 'int_name'}
    
    # Group by interval name and apply the aggregation
    result = df.groupby('int_name', as_index=False).agg(agg_dict)
    
    
    return result

def analyze_correlations(merged_data):
    import pandas as pd
    import numpy as np
    from scipy import stats
    
    # Function to check if column is numeric
    def is_numeric_column(column):
        try:
            pd.to_numeric(merged_data[column])
            return True
        except:
            return False
    
    # Get only numeric genomic features
    genomic_features = [col for col in merged_data.columns if col.endswith('_gen') and is_numeric_column(col)]
    
    # Define recombination features that exist in your data
    mut_features = [
        'r1_hald_mut', 'r2_hald_mut', 
        'c_null_mut', 'c_free_mut', 'c_equal_mut', 'c1_diff_mut', 'c2_diff_mut',
        'r1_group1_mut', 'r2_group1_mut', 'r1_group2_mut', 'r2_group2_mut'
    ]
    
    wild_features = [
        'r1_hald_wild', 'r2_hald_wild', 
        'c_null_wild', 'c_free_wild', 'c_equal_wild', 'c1_diff_wild', 'c2_diff_wild',
        'r1_group1_wild', 'r2_group1_wild', 'r1_group2_wild', 'r2_group2_wild'
    ]
    
    # Filter to ensure all features exist in the dataframe
    mut_features = [f for f in mut_features if f in merged_data.columns]
    wild_features = [f for f in wild_features if f in merged_data.columns]
    
    print("Analyzing genomic features:", len(genomic_features))
    print("First 5 genomic features:", genomic_features[:5])
    print("Mutant features:", mut_features)
    print("Wild features:", wild_features)
    
    # Create dictionaries for storing correlations
    correlations = {
        'genomic_mut': pd.DataFrame(index=genomic_features, columns=mut_features),
        'genomic_wild': pd.DataFrame(index=genomic_features, columns=wild_features),
        'mut_vs_wild': pd.DataFrame(index=mut_features, columns=wild_features)
    }
    
    # Calculate correlations
    for gen in genomic_features:
        gen_data = pd.to_numeric(merged_data[gen], errors='coerce')
        for mut in mut_features:
            mut_data = pd.to_numeric(merged_data[mut], errors='coerce')
            mask = ~(gen_data.isna() | mut_data.isna())
            if mask.sum() > 1:  # Ensure we have at least two pairs of values
                corr, p_value = stats.pearsonr(gen_data[mask], mut_data[mask])
                correlations['genomic_mut'].loc[gen, mut] = f"{corr:.3f} (p={p_value:.3e})"
            
    for gen in genomic_features:
        gen_data = pd.to_numeric(merged_data[gen], errors='coerce')
        for wild in wild_features:
            wild_data = pd.to_numeric(merged_data[wild], errors='coerce')
            mask = ~(gen_data.isna() | wild_data.isna())
            if mask.sum() > 1:
                corr, p_value = stats.pearsonr(gen_data[mask], wild_data[mask])
                correlations['genomic_wild'].loc[gen, wild] = f"{corr:.3f} (p={p_value:.3e})"
            
    for mut in mut_features:
        mut_data = pd.to_numeric(merged_data[mut], errors='coerce')
        for wild in wild_features:
            wild_data = pd.to_numeric(merged_data[wild], errors='coerce')
            mask = ~(mut_data.isna() | wild_data.isna())
            if mask.sum() > 1:
                corr, p_value = stats.pearsonr(mut_data[mask], wild_data[mask])
                correlations['mut_vs_wild'].loc[mut, wild] = f"{corr:.3f} (p={p_value:.3e})"
    
    # Function to extract significant correlations
    def get_significant_correlations(corr_df, threshold=0.5):
        significant = []
        for idx in corr_df.index:
            for col in corr_df.columns:
                value = corr_df.loc[idx, col]
                if pd.notna(value):
                    try:
                        corr = float(value.split()[0])
                        p_value = float(value.split('p=')[1].rstrip(')'))
                        if abs(corr) >= threshold and p_value < 0.05:
                            significant.append((idx, col, corr, p_value))
                    except:
                        continue
        return significant
    
    # Print results
    print("\nSignificant correlations (|r| >= 0.5, p < 0.05):")
    
    print("\n1. Genomic features vs Mutant type:")
    significant = get_significant_correlations(correlations['genomic_mut'])
    for feat1, feat2, corr, p_value in sorted(significant, key=lambda x: abs(x[2]), reverse=True):
        print(f"{feat1} vs {feat2}: r={corr:.3f} (p={p_value:.3e})")
        
    print("\n2. Genomic features vs Wild type:")
    significant = get_significant_correlations(correlations['genomic_wild'])
    for feat1, feat2, corr, p_value in sorted(significant, key=lambda x: abs(x[2]), reverse=True):
        print(f"{feat1} vs {feat2}: r={corr:.3f} (p={p_value:.3e})")
        
    print("\n3. Mutant vs Wild type:")
    significant = get_significant_correlations(correlations['mut_vs_wild'])
    for feat1, feat2, corr, p_value in sorted(significant, key=lambda x: abs(x[2]), reverse=True):
        print(f"{feat1} vs {feat2}: r={corr:.3f} (p={p_value:.3e})")
    
    return correlations

def analyze_feature_interactions(merged_data):
    """
    Анализ взаимодействий между геномными features и их комбинированных эффектов
    Это то, что требует профессор: не просто одномерные корреляции, а взаимодействия
    """
    import pandas as pd
    import numpy as np
    from sklearn.ensemble import RandomForestRegressor
    from sklearn.feature_selection import SelectKBest, f_regression
    from sklearn.preprocessing import StandardScaler
    from sklearn.linear_model import LinearRegression
    from sklearn.metrics import r2_score
    import itertools
    
    # Получаем геномные features
    genomic_features = [col for col in merged_data.columns if col.endswith('_gen')]
    
    # Основные target variables (то что профессор хочет предсказывать)
    target_vars = ['interference_mut', 'coincidence_mut', 'r1_mut', 'r2_mut']
    
    # Убираем non-numeric features
    numeric_features = []
    for feature in genomic_features:
        try:
            pd.to_numeric(merged_data[feature], errors='raise')
            numeric_features.append(feature)
        except:
            continue
    
    print(f"Анализируем {len(numeric_features)} числовых геномных характеристик")
    
    results = {}
    
    for target in target_vars:
        if target not in merged_data.columns:
            continue
            
        print(f"\n=== Анализ для {target} ===")
        
        # Подготавливаем данные
        X = merged_data[numeric_features].apply(pd.to_numeric, errors='coerce')
        y = pd.to_numeric(merged_data[target], errors='coerce')
        
        # Убираем NaN
        mask = ~(X.isna().any(axis=1) | y.isna())
        X_clean = X[mask]
        y_clean = y[mask]
        
        if len(X_clean) < 10:  # Минимум для анализа
            print(f"Недостаточно данных для {target}")
            continue
        
        # 1. RANDOM FOREST для feature importance
        rf = RandomForestRegressor(n_estimators=100, random_state=42)
        rf.fit(X_clean, y_clean)
        
        # Feature importance
        importance_df = pd.DataFrame({
            'feature': numeric_features,
            'importance': rf.feature_importances_
        }).sort_values('importance', ascending=False)
        
        print("Top 10 важных features (Random Forest):")
        print(importance_df.head(10))
        
        # 2. INTERACTION TERMS анализ
        top_features = importance_df.head(5)['feature'].tolist()
        
        print(f"\nАнализ взаимодействий между топ-5 features...")
        
        # Создаем interaction terms для топ features
        interaction_data = X_clean[top_features].copy()
        interaction_names = top_features.copy()
        
        # Добавляем попарные взаимодействия
        for i, feat1 in enumerate(top_features):
            for feat2 in top_features[i+1:]:
                interaction_name = f"{feat1}_x_{feat2}"
                interaction_data[interaction_name] = X_clean[feat1] * X_clean[feat2]
                interaction_names.append(interaction_name)
        
        # Multiple regression с interaction terms
        scaler = StandardScaler()
        X_scaled = scaler.fit_transform(interaction_data)
        
        reg = LinearRegression()
        reg.fit(X_scaled, y_clean)
        
        r2 = r2_score(y_clean, reg.predict(X_scaled))
        
        print(f"R² с interaction terms: {r2:.3f}")
        
        # Коэффициенты модели
        coef_df = pd.DataFrame({
            'feature': interaction_names,
            'coefficient': reg.coef_,
            'abs_coefficient': np.abs(reg.coef_)
        }).sort_values('abs_coefficient', ascending=False)
        
        print("Топ interaction effects:")
        print(coef_df.head(10))
        
        results[target] = {
            'feature_importance': importance_df,
            'interaction_coefficients': coef_df,
            'r2_with_interactions': r2
        }
    
    return results

def analyze_differential_response(merged_data):
    """
    Анализ дифференциального ответа: как мутация меняет эффекты features
    Это ключевое требование профессора: Δ(effect) = effect_mutant - effect_wildtype
    """
    import pandas as pd
    import numpy as np
    from scipy import stats
    
    # Получаем геномные features
    genomic_features = [col for col in merged_data.columns if col.endswith('_gen')]
    
    # Парные сравнения мутант vs wild
    comparison_pairs = [
        ('r1_hald_mut', 'r1_hald_wild'),
        ('r2_hald_mut', 'r2_hald_wild'),
        ('c_null_mut', 'c_null_wild'),
        ('c_free_mut', 'c_free_wild')
    ]
    
    results = {}
    
    for mut_var, wild_var in comparison_pairs:
        if mut_var not in merged_data.columns or wild_var not in merged_data.columns:
            continue
            
        print(f"\n=== Дифференциальный анализ: {mut_var} vs {wild_var} ===")
        
        # Вычисляем дифференциальный ответ
        mut_data = pd.to_numeric(merged_data[mut_var], errors='coerce')
        wild_data = pd.to_numeric(merged_data[wild_var], errors='coerce')
        
        differential_response = mut_data - wild_data
        merged_data[f'delta_{mut_var}'] = differential_response
        
        # Анализируем корреляции с геномными features
        diff_correlations = []
        
        for feature in genomic_features:
            feature_data = pd.to_numeric(merged_data[feature], errors='coerce')
            mask = ~(feature_data.isna() | differential_response.isna())
            
            if mask.sum() > 5:
                corr, p_value = stats.pearsonr(feature_data[mask], differential_response[mask])
                diff_correlations.append({
                    'feature': feature,
                    'correlation_with_delta': corr,
                    'p_value': p_value,
                    'significant': p_value < 0.05 and abs(corr) > 0.3
                })
        
        diff_df = pd.DataFrame(diff_correlations).sort_values('correlation_with_delta', 
                                                             key=abs, ascending=False)
        
        print("Топ correlations с дифференциальным ответом:")
        significant_df = diff_df[diff_df['significant']]
        print(significant_df.head(10))
        
        results[f'{mut_var}_vs_{wild_var}'] = {
            'differential_correlations': diff_df,
            'significant_predictors': significant_df
        }
    
    return results

def analyze_genomic_context_profiles(merged_data):
    """
    Анализ комбинированных геномных профилей
    Профессор хочет видеть: High accessibility + Low TE + High gene density = ?
    """
    import pandas as pd
    import numpy as np
    from sklearn.cluster import KMeans
    from sklearn.preprocessing import StandardScaler
    
    print("=== Анализ комбинированных геномных контекстов ===")
    
    # Ключевые features для профилирования
    key_features = [
        'DNaseI_av_gen', 'MNase_score_gen',  # Chromatin accessibility
        'L_TE_cs_gen', 'L_TE_rel_cs_gen',    # Transposable elements
        'N_HC_sv_gen', 'L_HC_sv_gen',        # Gene content
        'H3K4me1_av_gen', 'H3K36me3_av_gen'  # Active histone marks
    ]
    
    # Отбираем features, которые есть в данных
    available_features = [f for f in key_features if f in merged_data.columns]
    
    if len(available_features) < 3:
        print("Недостаточно ключевых features для профилирования")
        return None
    
    # Подготавливаем данные
    X = merged_data[available_features].apply(pd.to_numeric, errors='coerce')
    mask = ~X.isna().any(axis=1)
    X_clean = X[mask]
    
    if len(X_clean) < 20:
        print("Недостаточно данных для кластеризации")
        return None
    
    # Нормализуем данные
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X_clean)
    
    # Кластеризация для выделения геномных контекстов
    kmeans = KMeans(n_clusters=5, random_state=42)
    clusters = kmeans.fit_predict(X_scaled)
    
    # Добавляем кластеры к данным
    merged_clean = merged_data[mask].copy()
    merged_clean['genomic_cluster'] = clusters
    
    print(f"Выделено {len(np.unique(clusters))} геномных контекстов")
    
    # Анализируем характеристики каждого кластера
    cluster_profiles = []
    
    for cluster_id in np.unique(clusters):
        cluster_mask = merged_clean['genomic_cluster'] == cluster_id
        cluster_data = merged_clean[cluster_mask]
        
        profile = {
            'cluster_id': cluster_id,
            'n_intervals': len(cluster_data)
        }
        
        # Средние значения features
        for feature in available_features:
            profile[f'mean_{feature}'] = cluster_data[feature].mean()
        
        # Средние значения recombination response
        target_vars = ['interference_mut', 'r1_mut', 'r2_mut']
        for target in target_vars:
            if target in merged_clean.columns:
                profile[f'mean_{target}'] = cluster_data[target].mean()
        
        cluster_profiles.append(profile)
    
    profile_df = pd.DataFrame(cluster_profiles)
    
    print("\nПрофили геномных контекстов:")
    print(profile_df)
    
    # Находим наиболее контрастные кластеры
    if 'mean_interference_mut' in profile_df.columns:
        max_interference_cluster = profile_df.loc[profile_df['mean_interference_mut'].idxmax()]
        min_interference_cluster = profile_df.loc[profile_df['mean_interference_mut'].idxmin()]
        
        print(f"\nКластер с максимальной интерференцией (кластер {max_interference_cluster['cluster_id']}):")
        print(f"Средняя интерференция: {max_interference_cluster['mean_interference_mut']:.3f}")
        
        print(f"\nКластер с минимальной интерференцией (кластер {min_interference_cluster['cluster_id']}):")
        print(f"Средняя интерференция: {min_interference_cluster['mean_interference_mut']:.3f}")
    
    return {
        'cluster_profiles': profile_df,
        'clustered_data': merged_clean,
        'feature_names': available_features
    }

def comprehensive_analysis_report(merged_data):
    """
    Комплексный анализ, который отвечает на требования профессора
    """
    print("="*60)
    print("КОМПЛЕКСНЫЙ АНАЛИЗ ГЕНОМНЫХ FEATURES И РЕКОМБИНАЦИИ")
    print("="*60)
    
    # 1. Анализ взаимодействий features
    print("\n1. АНАЛИЗ ВЗАИМОДЕЙСТВИЙ МЕЖДУ ГЕНОМНЫМИ ХАРАКТЕРИСТИКАМИ")
    interaction_results = analyze_feature_interactions(merged_data)
    
    # 2. Дифференциальный анализ мутант vs wild
    print("\n2. ДИФФЕРЕНЦИАЛЬНЫЙ АНАЛИЗ (МУТАНТ vs ДИКИЙ ТИП)")
    differential_results = analyze_differential_response(merged_data)
    
    # 3. Геномные контексты
    print("\n3. АНАЛИЗ КОМБИНИРОВАННЫХ ГЕНОМНЫХ КОНТЕКСТОВ")
    context_results = analyze_genomic_context_profiles(merged_data)
    
    return {
        'interactions': interaction_results,
        'differential': differential_results,
        'contexts': context_results
    }

# ИСПОЛЬЗОВАНИЕ:
# Добавьте эти функции к вашему существующему коду и запустите:

# comprehensive_results = comprehensive_analysis_report(merged_data)

# Сохранение расширенных результатов
def save_comprehensive_results(comprehensive_results, filename='comprehensive_analysis.xlsx'):
    """Сохраняет все результаты комплексного анализа"""
    with pd.ExcelWriter(filename) as writer:
        
        # Interaction analysis results
        if 'interactions' in comprehensive_results:
            for target, results in comprehensive_results['interactions'].items():
                if 'feature_importance' in results:
                    results['feature_importance'].to_excel(writer, 
                                                         sheet_name=f'{target}_importance')
                if 'interaction_coefficients' in results:
                    results['interaction_coefficients'].to_excel(writer, 
                                                               sheet_name=f'{target}_interactions')
        
        # Differential analysis results
        if 'differential' in comprehensive_results:
            for comparison, results in comprehensive_results['differential'].items():
                if 'significant_predictors' in results:
                    results['significant_predictors'].to_excel(writer, 
                                                             sheet_name=f'{comparison[:20]}_diff')
        
        # Genomic context results
        if 'contexts' in comprehensive_results and comprehensive_results['contexts']:
            comprehensive_results['contexts']['cluster_profiles'].to_excel(writer, 
                                                                         sheet_name='genomic_contexts')

# comprehensive_results = comprehensive_analysis_report(merged_data)
# save_comprehensive_results(comprehensive_results)
