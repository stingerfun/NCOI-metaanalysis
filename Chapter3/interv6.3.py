import numpy as np
import pandas as pd
from scipy.optimize import differential_evolution, minimize
from scipy.stats import chi2
from typing import List, Dict, Tuple, Optional
from math import log
import os
from datetime import datetime

class GeneticDataAnalyzer:
    def __init__(self, min_distance: float = 7.0, max_distance: float = 25.0):
        self.min_distance = min_distance
        self.max_distance = max_distance
        
    def read_genetic_data(self, file_path: str) -> pd.DataFrame:
        """
        Reads and preprocesses genetic data from file
        """
        try:
            with open(file_path, 'r') as f:
                lines = f.readlines()
            
            processed_lines = []
            for line in lines:
                # Remove 'T' or space at start
                if line.startswith('T') or line.startswith(' '):
                    line = line.lstrip(' T')
                
                fields = [x.strip() for x in line.split() if x.strip()]
                
                # Remove single 'S' in second position
                if len(fields) > 1 and fields[1] == 'S':
                    fields.pop(1)
                
                processed_lines.append(fields)
            
            # Find maximum field count
            max_fields = max(len(line) for line in processed_lines)
            
            # Pad shorter lines with '-'
            data = []
            for line in processed_lines:
                if len(line) < max_fields:
                    line.extend(['-'] * (max_fields - len(line)))
                data.append(line)
            
            df = pd.DataFrame(data)
            df[0] = pd.to_numeric(df[0], errors='coerce')
            
            return df
            
        except Exception as e:
            raise ValueError(f"Error reading genetic data: {str(e)}")

    def find_marker_triplets(self, positions: np.ndarray) -> List[Tuple[int, int, int]]:
        """
        Identifies valid marker triplets for interference analysis
        """
        triplets = []
        n = len(positions)
        
        for i in range(n-2):
            for j in range(i+1, n-1):
                dist1 = positions[j] - positions[i]
                if dist1 < self.min_distance or dist1 > self.max_distance:
                    continue
                    
                for k in range(j+1, n):
                    dist2 = abs(positions[k] - positions[j])
                    if dist2 < self.min_distance or dist2 > self.max_distance:
                        continue
                        
                    total_dist = abs(positions[k] - positions[i])
                    if total_dist > 2 * self.max_distance:
                        continue
                        
                    triplets.append((i, j, k))
        
        return triplets

    def calculate_genotype_frequencies(self, triple: Tuple[List[str], List[str], List[str]], 
                                     population_size: int) -> Tuple[Dict[str, int], int]:
        """
        Calculates genotype frequencies from genetic data
        """
        ia1, ia2, ia3 = triple
        counts = {f"n{i}{j}{k}": 0 for i in range(1,4) for j in range(1,4) for k in range(1,4)}
        
        genotype_map = {
            ('H','H','H'): 'n333', ('A','A','A'): 'n111', ('B','B','B'): 'n222',
            ('A','B','B'): 'n122', ('A','A','B'): 'n112', ('B','B','A'): 'n221',
            ('B','A','A'): 'n211', ('A','B','A'): 'n121', ('B','A','B'): 'n212',
            ('H','A','A'): 'n311', ('A','H','A'): 'n131', ('A','A','H'): 'n113',
            ('B','B','H'): 'n223', ('B','H','B'): 'n232', ('H','B','B'): 'n322',
            ('H','H','A'): 'n331', ('H','A','H'): 'n313', ('A','H','H'): 'n133',
            ('H','H','B'): 'n332', ('H','B','H'): 'n323', ('B','H','H'): 'n233',
            ('A','B','H'): 'n123', ('B','A','H'): 'n213', ('A','H','B'): 'n132',
            ('B','H','A'): 'n231', ('H','A','B'): 'n312', ('H','B','A'): 'n321'
        }
        
        valid_count = 0
        
        for j in range(min(len(ia1), len(ia2), len(ia3), population_size)):
            if ia1[j] == '-' or ia2[j] == '-' or ia3[j] == '-':
                continue
                
            valid_count += 1
            genotype = (ia1[j], ia2[j], ia3[j])
            
            if genotype in genotype_map:
                counts[genotype_map[genotype]] += 1
        
        return counts, valid_count

    def calculate_recombination_rates(self, 
                                    observed_freqs: Dict[str, int]) -> Tuple[Tuple[float, float], Tuple[float, float]]:
        """
        Calculates recombination rates for both sex-specific groups
        """
        total = sum(observed_freqs.values())
        if total == 0:
            return (0.0, 0.0), (0.0, 0.0)
        
        # Pattern definitions for recombination rates
        patterns = {
            'group1': {
                'r1': {
                    'n121': 1.0, 'n122': 1.0, 'n123': 1.0,
                    'n311': 0.5, 'n312': 0.5, 'n313': 0.5,
                    'n331': 0.5, 'n332': 0.5, 'n333': 0.25
                },
                'r2': {
                    'n111': 1.0, 'n112': 1.0, 'n113': 1.0,
                    'n311': 0.5, 'n312': 0.5, 'n313': 0.5,
                    'n331': 0.5, 'n332': 0.5, 'n333': 0.25
                }
            },
            'group2': {
                'r1': {
                    'n211': 1.0, 'n212': 1.0, 'n213': 1.0,
                    'n321': 0.5, 'n322': 0.5, 'n323': 0.5,
                    'n331': 0.5, 'n332': 0.5, 'n333': 0.25
                },
                'r2': {
                    'n221': 1.0, 'n222': 1.0, 'n223': 1.0,
                    'n321': 0.5, 'n322': 0.5, 'n323': 0.5,
                    'n331': 0.5, 'n332': 0.5, 'n333': 0.25
                }
            }
        }
        
        # Calculate recombination rates
        rates = {}
        for group in ['group1', 'group2']:
            for rate_type in ['r1', 'r2']:
                pattern = patterns[group][rate_type]
                rate_sum = sum(observed_freqs.get(k, 0) * v for k, v in pattern.items())
                rates[f"{group}_{rate_type}"] = min(0.5, rate_sum / (2 * total))
        
        return ((rates['group1_r1'], rates['group1_r2']), 
                (rates['group2_r1'], rates['group2_r2']))

    def optimize_model(self, observed_freqs: Dict[str, int], 
                  initial_rates: Tuple[Tuple[float, float], Tuple[float, float]],
                  model_type: str = 'free') -> Dict[str, float]:
        """
        Optimizes the genetic interference model using Nelder-Mead method.
        
        Parameters:
        -----------
        observed_freqs : Dict[str, int]
            Observed genotype frequencies
        initial_rates : Tuple[Tuple[float, float], Tuple[float, float]]
            Initial recombination rates ((r1_1, r2_1), (r1_2, r2_2))
        model_type : str
            Type of model to fit:
            - 'null': No sex difference, c = 1
            - 'free': No sex difference, c free
            - 'equal_c': Sex difference with c1 = c2
            - 'diff_c': Sex difference with c1 ? c2
        """
        total = sum(observed_freqs.values())
        (r1_1, r2_1), (r1_2, r2_2) = initial_rates
        max_r = max(r1_1, r2_1, r1_2, r2_2)
        if max_r < 0.001:
            max_r = 0.5
    
        if model_type == 'null':
            # Model 1: No sex difference, c = 1
            return {
                'r1_group1': r1_1,
                'r2_group1': r2_1,
                'r1_group2': r1_1,
                'r2_group2': r2_1,
                'c1': 1.0,
                'c2': 1.0,
                'log_likelihood': self._calculate_likelihood([1.0, 1.0], 
                                                           observed_freqs, r1_1, r2_1, r1_1, r2_1, total)
            }
        
        elif model_type == 'free':
            # Model 2: No sex difference, c free (current implementation)
            def objective_with_bounds(params):
                c = params[0]
                if c < 0 or c > 1/max_r:
                    return float('inf')
                return self._calculate_likelihood([c, c], observed_freqs, r1_1, r2_1, r1_1, r2_1, total)
                
            initial_points = [[0.01], [1.0], [1.5]]
            
        elif model_type == 'equal_c':
            # Model 3: Sex difference with c1 = c2 (current implementation)
            def objective_with_bounds(params):
                c = params[0]
                if c < 0 or c > 1/max_r:
                    return float('inf')
                return self._calculate_likelihood([c, c], observed_freqs, r1_1, r2_1, r1_2, r2_2, total)
                
            initial_points = [[0.01], [1.0], [1.5]]
            
        else:  # model_type == 'diff_c'
            # Model 4: Sex difference with c1 ? c2
            def objective_with_bounds(params):
                c1, c2 = params
                if c1 < 0 or c2 < 0 or c1 > 1/max_r or c2 > 1/max_r:
                    return float('inf')
                return self._calculate_likelihood([c1, c2], observed_freqs, r1_1, r2_1, r1_2, r2_2, total)
                
            initial_points = [[0.01, 0.01], [1.0, 1.0], [1.5, 1.5]]
    
        # Optimization for models 2-4
        best_result = None
        best_likelihood = float('inf')
        
        for x0 in initial_points:
            try:
                result = minimize(
                    objective_with_bounds,
                    x0=x0,
                    method='Nelder-Mead',
                    options={'maxiter': 200, 'xatol': 1e-8, 'fatol': 1e-8}
                )
                
                if result.success and result.fun < best_likelihood:
                    if all(0 <= x <= 1/max_r for x in result.x):
                        best_likelihood = result.fun
                        best_result = result
                        
            except Exception as e:
                continue
        
        if best_result is None:
            raise ValueError("Optimization failed for all initial points")
            
        if model_type in ['free', 'equal_c']:
            c_opt = min(max(0, best_result.x[0]), 1/max_r)
            return {
                'r1_group1': r1_1,
                'r2_group1': r2_1,
                'r1_group2': r1_2 if model_type == 'equal_c' else r1_1,
                'r2_group2': r2_2 if model_type == 'equal_c' else r2_1,
                'c1': c_opt,
                'c2': c_opt,
                'log_likelihood': -best_likelihood
            }
        else:  # diff_c model
            c1_opt = min(max(0, best_result.x[0]), 1/max_r)
            c2_opt = min(max(0, best_result.x[1]), 1/max_r)
            return {
                'r1_group1': r1_1,
                'r2_group1': r2_1,
                'r1_group2': r1_2,
                'r2_group2': r2_2,
                'c1': c1_opt,
                'c2': c2_opt,
                'log_likelihood': -best_likelihood
            }  
    def _calculate_theoretical_frequencies(self, params: np.ndarray, r1_f: float, r2_f: float,
                                    r1_m: float, r2_m: float, total: int) -> Dict[str, float]:
        """
        Calculates theoretical genotype frequencies for F2 intercross model.
        """
        c_f, c_m = params
        
        # Calculate basic probabilities for female parent
        i111 = (1 - r1_f - r2_f + c_f*r1_f*r2_f)/2    # No recombination
        i112 = (r2_f - c_f*r1_f*r2_f)/2                # Second interval recombination
        i121 = c_f*r1_f*r2_f/2                         # Double recombination
        i122 = (r1_f - c_f*r1_f*r2_f)/2                # First interval recombination
        i211 = i122                                     # Symmetric classes
        i212 = i121
        i221 = i112
        i222 = i111
        
        # Calculate basic probabilities for male parent
        a111 = (1 - r1_m - r2_m + c_m*r1_m*r2_m)/2    # No recombination
        a112 = (r2_m - c_m*r1_m*r2_m)/2                # Second interval recombination
        a121 = c_m*r1_m*r2_m/2                         # Double recombination
        a122 = (r1_m - c_m*r1_m*r2_m)/2                # First interval recombination
        a211 = a122                                     # Symmetric classes
        a212 = a121
        a221 = a112
        a222 = a111
        
        # Calculate F2 genotype frequencies
        freq = {
            # Triple homozygotes
            'n111': (i111*a111) * total,                # AABBCC
            'n222': (i222*a222) * total,                # aabbcc
            
            # Double homozygotes with one heterozygote
            'n112': (i112*a112) * total,                # AABBcc
            'n121': (i121*a121) * total,                # AAbbCC
            'n122': (i122*a122) * total,                # AAbbcc
            'n211': (i211*a211) * total,                # aaBBCC
            'n212': (i212*a212) * total,                # aaBBcc
            'n221': (i221*a221) * total,                # aabbCC
            
            # Triple heterozygote
            'n333': (i111*a222 + i222*a111 + i112*a221 + i221*a112 + 
                     i122*a211 + i211*a122 + i121*a212 + i212*a121) * total,  # AaBbCc
            
            # Double heterozygotes
            'n311': (i111*a211 + i211*a111) * total,    # AABbCC
            'n113': (i111*a112 + i112*a111) * total,    # AABBCc
            'n223': (i221*a222 + i222*a221) * total,    # aabbCc
            'n232': (i212*a222 + i222*a212) * total,    # aaBbcc
            'n322': (i122*a222 + i222*a122) * total,    # aaBbcc
            
            # Complex heterozygous combinations
            'n331': (i111*a221 + i221*a111 + i121*a211 + i211*a121) * total,  # AaBBCC
            'n313': (i111*a212 + i212*a111 + i211*a112 + i112*a211) * total,  # AABbCc
            'n133': (i111*a122 + i112*a121 + i121*a112 + i122*a111) * total,  # AABCc
            'n332': (i112*a222 + i222*a112 + i122*a212 + i212*a122) * total,  # AaBBcc
            'n323': (i121*a222 + i222*a121 + i221*a122 + i122*a221) * total,  # aaBbCc
            'n233': (i211*a222 + i212*a221 + i221*a212 + i222*a211) * total,  # aaBCc
            
            # Additional combinations
            'n123': (i121*a122 + i122*a121) * total,    # AAbbCc
            'n213': (i211*a212 + i212*a211) * total,    # aaBBCc
            'n132': (i112*a122 + i122*a112) * total,    # AABbcc
            'n231': (i211*a221 + i221*a211) * total,    # aaBbCC
            'n312': (i112*a212 + i212*a112) * total,    # AABbcc
            'n321': (i121*a221 + i221*a121) * total,    # aaBbCC
            
            # Single heterozygotes
            'n131': (i111*a121 + i121*a111) * total     # AABbCC
        }
        
        # Normalize frequencies
        freq_sum = sum(freq.values())
        if freq_sum > 0:
            for key in freq:
                freq[key] = (freq[key] / freq_sum) * total
                
        return freq

        
    def _calculate_likelihood(self, params: np.ndarray, observed_freqs: Dict[str, int], 
                            r1_1: float, r2_1: float, r1_2: float, r2_2: float, 
                            total: int) -> float:
        """
        Calculates the log-likelihood for the genetic model
        """
        theoretical = self._calculate_theoretical_frequencies(
            params, r1_1, r2_1, r1_2, r2_2, total
        )
        
        ll = 0
        for genotype, freq in theoretical.items():
            if genotype in observed_freqs:
                if freq <= 0:
                    return float('inf')
                ll += observed_freqs[genotype] * np.log(freq + 1e-10)
        
        return -ll

    def analyze_interference(self, data: pd.DataFrame) -> pd.DataFrame:
        """
        Performs complete interference analysis on genetic data with four different models:
        1. No sex difference, c = 1
        2. No sex difference, c free
        3. Sex difference with c1 = c2
        4. Sex difference with c1 ? c2
        """
        results = []
        positions = data.iloc[:, 0].values
        markers = data.iloc[:, 1].values
        genotype_data = data.iloc[:, 2:].values
        
        triplets = self.find_marker_triplets(positions)
        
        for i, j, k in triplets:
            try:
                # Calculate intervals
                interval1 = abs(positions[j] - positions[i])
                interval2 = abs(positions[k] - positions[j])
                
                if interval1 > self.max_distance or interval2 > self.max_distance:
                    continue
                
                # Process genotype data
                triple = (
                    [str(x) for x in genotype_data[i]],
                    [str(x) for x in genotype_data[j]],
                    [str(x) for x in genotype_data[k]]
                )
                
                observed_freqs, valid_count = self.calculate_genotype_frequencies(
                    triple, len(triple[0])
                )
                
                if valid_count < 30:
                    continue
                
                # Calculate recombination rates
                recombination_rates = self.calculate_recombination_rates(observed_freqs)
                
                # Fit all four models
                model1 = self.optimize_model(
                    observed_freqs,
                    (recombination_rates[0], recombination_rates[0]),
                    model_type='null'
                )
                
                model2 = self.optimize_model(
                    observed_freqs,
                    (recombination_rates[0], recombination_rates[0]),
                    model_type='free'
                )
                
                model3 = self.optimize_model(
                    observed_freqs,
                    recombination_rates,
                    model_type='equal_c'
                )
                
                model4 = self.optimize_model(
                    observed_freqs,
                    recombination_rates,
                    model_type='diff_c'
                )
                
                # Calculate likelihood ratio tests
                lr_1vs2 = 2 * (model2['log_likelihood'] - (-model1['log_likelihood']))
                lr_2vs3 = 2 * (model3['log_likelihood'] - model2['log_likelihood'])
                lr_3vs4 = 2 * (model4['log_likelihood'] - model3['log_likelihood'])
                
                results.append({
                    'start_position': positions[i],
                    'middle_position': positions[j],
                    'end_position': positions[k],
                    'markers': f"{markers[i]} - {markers[j]} - {markers[k]}",
                    'interval1_length': interval1,
                    'interval2_length': interval2,
                    'r1_hald': model1['r1_group1'],
                    'r2_hald': model1['r2_group1'],
                    # Model likelihoods
                    'model1_ll': model1['log_likelihood']*-1,
                    'model2_ll': model2['log_likelihood'],
                    'model3_ll': model3['log_likelihood'],
                    'model4_ll': model4['log_likelihood'],
                    
                    # Likelihood ratio tests
                    'lr_1vs2': lr_1vs2,
                    'lr_2vs3': lr_2vs3,
                    'lr_3vs4': lr_3vs4,
                    
                    # P-values
                    'p_1vs2': 1 - chi2.cdf(lr_1vs2, df=1),
                    'p_2vs3': 1 - chi2.cdf(lr_2vs3, df=3),
                    'p_3vs4': 1 - chi2.cdf(lr_3vs4, df=1),
                    
                    # Model parameters
                    'c_null': model1['c1'],
                    'c_free': model2['c1'],
                    'c_equal': model3['c1'],
                    'c1_diff': model4['c1'],
                    'c2_diff': model4['c2'],
                    
                    # Basic recombination rates
                    'r1_group1': model4['r1_group1'],
                    'r2_group1': model4['r2_group1'],
                    'r1_group2': model4['r1_group2'],
                    'r2_group2': model4['r2_group2'],
                    
                    # Additional information
                    'valid_samples': valid_count,
                    'observed_freqs': observed_freqs
                })
                
            except Exception as e:
                print(f"Error processing triplet {i},{j},{k}: {str(e)}")
                continue
        
        return pd.DataFrame(results)
                              

# Initialize analyzer
analyzer = GeneticDataAnalyzer(min_distance=5.0, max_distance=25.0)

# Directory containing genetic data files
directory = r"/lustre1/home/fahima/ssapielki1/F2INTERF/data"

# Get list of text files
files = [f for f in os.listdir(directory) if f.endswith('.txt')]

print(f"Found {len(files)} files to process:")
for f in files:
    print(f"  - {f}")

all_results = []

for filename in files:
    file_path = os.path.join(directory, filename)
    
    try:
        # Read data
        data = analyzer.read_genetic_data(file_path)
        print(f"\nProcessing {filename}")
        print(f"Read {data.shape[0]} markers and {data.shape[1]-2} individuals")
        
        # Analyze interference
        results = analyzer.analyze_interference(data)
        
        if not results.empty:
            results['source_file'] = filename
            all_results.append(results)
            print(f"Successfully processed {filename}")
            print(f"Found {len(results)} marker triplets")
        else:
            print(f"No valid results found for {filename}")
            
    except Exception as e:
        print(f"Error processing file {filename}: {str(e)}")
        continue

# Combine and save results if any were found
if all_results:
    combined_results = pd.concat(all_results, ignore_index=True)
    #combined_results = combined_results.sort_values('log_likelihood_sex_specific', ascending=False)
    
    print("\nFINAL RESULTS:")
    print(f"Total files processed: {len(all_results)}")
    print(f"Total marker triplets found: {len(combined_results)}")
    
    # Create timestamped filename
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    output_filename = f"combined_results_{timestamp}.xlsx"
    output_path = os.path.join(directory, output_filename)
    
    # Save results
    combined_results.to_excel(output_path, index=False)
    print(f"\nResults saved to: {output_path}")
    
    # Print detailed results
    for _, row in combined_results.iterrows():
        print("\n" + "="*80)
        print(f"\nInterval: {row['start_position']:.2f} - {row['middle_position']:.2f} - {row['end_position']:.2f} cM")
        print(f"Interval lengths: {row['interval1_length']:.2f}, {row['interval2_length']:.2f} cM")
        print(f"Markers: {row['markers']}")
        
        print("\nNull Model (No Sex Differences):")
        print(f"Recombination rate r1: {row['r1']:.3f}")
        print(f"Recombination rate r2: {row['r2']:.3f}")
        print(f"Coincidence coefficient: {row['coincidence']:.3f}")
        print(f"Interference: {row['interference']:.3f}")
        print(f"Log-likelihood: {row['log_likelihood_null']:.3f}")
        
        print("\nSex-Specific Model:")
        print(f"Group 1: r1 = {row['r1_group1']:.3f}, r2 = {row['r2_group1']:.3f}, c1 = {row['c1']:.3f}")
        print(f"Group 2: r1 = {row['r1_group2']:.3f}, r2 = {row['r2_group2']:.3f}, c2 = {row['c2']:.3f}")
        print(f"Log-likelihood: {row['log_likelihood_sex_specific']:.3f}")
        print(f"LR statistic: {row['lr_statistic']:.3f}")
        print(f"P-value: {row['p_value']:.3e}")
        print(f"Valid samples: {row['valid_samples']}")
        
        if row['observed_freqs']:
            print("\nGenotype frequencies:")
            total = sum(row['observed_freqs'].values())
            for genotype, count in row['observed_freqs'].items():
                if count > 0:
                    print(f"{genotype}: {count} ({count/total:.3f})")
else:
    print("No valid results found for analysis")