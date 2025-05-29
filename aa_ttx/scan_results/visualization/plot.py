import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.ticker import LogLocator, NullFormatter, MaxNLocator
import sys
import random

# -----------------------------------------------------------------------------
# CONFIGURATION
# -----------------------------------------------------------------------------
SORT_MODE = 'global' # 'global' or 'local'
GLOBAL_SORT_ANCHOR_SQRT_S = 1600.0 # e.g., 1600.0 or None (for mean sort)
SQRT_S_MIN_DISPLAY = None
SQRT_S_MAX_DISPLAY = 9001.0
LOWER_INSET_WIDTH_SCALE = 1.0 # Adjust to align x-axes data areas
TOP_PERCENTILES_TO_PLOT = [0.1, 0.5]

# Heatmap Y-axis labels for Global Sort mode:
#    - Provide a list of graph names, e.g., ["GL258", "GL34", "GL292"]
#    - Set to None or [] to show all graph names with HEATMAP_ALL_GRAPHS_FONTSIZE
#    - Set to a non-empty list to show only those specific graph names.
HEATMAP_GLOBAL_SORT_Y_LABELS = None # Example: ["GL258", "GL34", "GL292"], or None for all
HEATMAP_ALL_GRAPHS_FONTSIZE = 3  # Font size for y-tick labels when all graphs are shown (if HEATMAP_GLOBAL_SORT_Y_LABELS is None/empty)
# -----------------------------------------------------------------------------

def generate_mock_raw_data_file(filename="raw_data.py", N_supergraphs=124, N_energies=32):
    """Generates a mock raw_data.py file for testing, trying to mimic image features."""
    sqrt_s_values = np.logspace(np.log10(200), np.log10(10000), N_energies)
    graph_names = []
    predefined_names = ["GL258", "GL386", "GL212", "GL34", "GL80", "GL246", "GL158", "GL190", "GL376"]
    for i in range(N_supergraphs):
        if i < len(predefined_names) and predefined_names[i] not in graph_names:
            graph_names.append(predefined_names[i])
        else: 
            new_name = f"GL{i+N_supergraphs:03d}" 
            while new_name in graph_names: 
                new_name = f"GL{random.randint(N_supergraphs*2, N_supergraphs*10):03d}"
            graph_names.append(new_name)
    
    if "GL292" not in graph_names and N_supergraphs > len(predefined_names):
        # Try to place GL292 if it's not there and we have space (e.g. as last element for variety)
        if N_supergraphs > 0 : graph_names[min(N_supergraphs-1, len(predefined_names))] = "GL292" 
    elif "GL292" in graph_names and graph_names[min(N_supergraphs-1, len(predefined_names))] != "GL292": 
        # If GL292 is already present but not at this specific position, swap it
        try:
            idx_gl292 = graph_names.index("GL292")
            target_idx = min(N_supergraphs-1, len(predefined_names))
            if idx_gl292 != target_idx:
                 graph_names[idx_gl292], graph_names[target_idx] = graph_names[target_idx], graph_names[idx_gl292]
        except ValueError: pass    
    
    # Ensure uniqueness if any accidental duplicates (less likely with new logic but good safeguard)
    if len(graph_names) != len(set(graph_names)):
        print("Warning: Duplicate graph names generated in mock data. Re-generating generic names.")
        graph_names = [f"MockGL{i+1:03d}" for i in range(N_supergraphs)]

    output_data = []
    # More refined mock data generation (conceptual)
    for s_idx, s_val in enumerate(sqrt_s_values):
        energy_point_data = []
        # Simplified total cross-section shape for mock data
        base_val = 0.002 * (s_val/200.0)**0.5 
        if s_val < 350: # Below ttbar threshold (approx)
            total_xs_at_s = base_val * np.random.uniform(0.001, 0.01) * (s_val/350.0)**3 # Suppressed
        elif 350 <= s_val < 700: # Rise
            total_xs_at_s = base_val * ( (s_val-350)/350 * 1.0 + 0.01) * np.random.uniform(0.5, 1.0)
        elif 700 <= s_val < 2000: # Peak region
            total_xs_at_s = base_val * (1.0 - ((s_val-1350)/650)**2 *0.7 ) * np.random.uniform(0.7, 1.3) # Parabolic peak around 1350
            total_xs_at_s = max(total_xs_at_s, base_val * 0.1) # ensure it's not too small
        elif 2000 <= s_val < 5000: # Fall-off start
            total_xs_at_s = base_val * (0.5 + (5000-s_val)/3000 * 0.5) * np.random.uniform(0.8,1.2)
        else: # High energy fall-off
            total_xs_at_s = base_val * (1.0 * (2000.0/s_val)**0.8) * np.random.uniform(0.8,1.2)
        total_xs_at_s = max(total_xs_at_s, 1e-7) # Floor value

        # Define cancellation characteristics
        cancellation_factor = np.random.uniform(10.0, 250.0); 
        small_xs_mag = 1e-4 * base_val 
        if abs(total_xs_at_s) < small_xs_mag : cancellation_factor = np.random.uniform(1.5, 15.0)
        
        if total_xs_at_s > 0:
            target_sum_pos = total_xs_at_s * cancellation_factor
            target_sum_neg = total_xs_at_s * (1.0 - cancellation_factor)
        elif total_xs_at_s < 0: # Should not happen with current total_xs_at_s logic
            target_sum_neg = total_xs_at_s * cancellation_factor
            target_sum_pos = total_xs_at_s * (1.0 - cancellation_factor)
        else: # total_xs_at_s is zero or very near
            target_sum_pos = np.random.uniform(1e-5, 1e-4) * (1000.0/s_val)*base_val
            target_sum_neg = -target_sum_pos * np.random.uniform(0.8, 1.2)

        # Generate contributions
        contributions = np.random.randn(N_supergraphs) # Start with random numbers
        # Make some graphs dominant
        dominant_indices = np.random.choice(N_supergraphs, N_supergraphs // 8, replace=False)
        contributions[dominant_indices] *= np.random.uniform(5, 15, size=len(dominant_indices))
        
        # Specific behavior for some graphs based on name (example)
        if "GL258" in graph_names: contributions[graph_names.index("GL258")] *= (s_val/1000)**0.5 * (1 if s_val < 3000 else -1) # Rises, then flips
        if "GL34" in graph_names: contributions[graph_names.index("GL34")] *= (1 - s_val/10000) # Decreases with energy
        if "GL292" in graph_names: contributions[graph_names.index("GL292")] *= -1.5 * (2000/s_val)**0.3 # Negative and falls slowly

        # At high energy, some non-dominant graphs get suppressed further
        if s_val > 3000:
            zero_out_fraction = (s_val - 3000) / 7000 
            num_to_zero = int(N_supergraphs * zero_out_fraction * 0.7)
            non_dominant_indices = [i for i in range(N_supergraphs) if i not in dominant_indices]
            if num_to_zero > 0 and len(non_dominant_indices) >= num_to_zero:
                indices_to_zero = np.random.choice(non_dominant_indices, num_to_zero, replace=False)
                contributions[indices_to_zero] *= 1e-5 # Effectively zero them out

        # Normalize to achieve target_sum_pos and target_sum_neg
        pos_mask = contributions > 0
        neg_mask = contributions <= 0
        current_sum_pos = np.sum(contributions[pos_mask])
        current_sum_neg = np.sum(contributions[neg_mask])

        if np.sum(pos_mask) > 0:
            if abs(current_sum_pos) > 1e-12: contributions[pos_mask] *= (target_sum_pos / current_sum_pos)
            elif target_sum_pos != 0: contributions[pos_mask] = target_sum_pos / np.sum(pos_mask) # Distribute if all were zero
        elif target_sum_pos > 0 and N_supergraphs > 0: contributions[np.random.randint(0,N_supergraphs)] += target_sum_pos # Add if no positive graphs

        if np.sum(neg_mask) > 0:
            if abs(current_sum_neg) > 1e-12: contributions[neg_mask] *= (target_sum_neg / current_sum_neg)
            elif target_sum_neg != 0: contributions[neg_mask] = target_sum_neg / np.sum(neg_mask) # Distribute if all were zero
        elif target_sum_neg < 0 and N_supergraphs > 0: contributions[np.random.randint(0,N_supergraphs)] += target_sum_neg # Add if no negative graphs
        
        for i in range(N_supergraphs):
            energy_point_data.append([graph_names[i], contributions[i]])
        
        # Shuffle the order of graphs within this energy point to simulate unordered real data
        random.shuffle(energy_point_data)
        output_data.append([s_val, energy_point_data])

    # Write to file
    with open(filename, "w") as f:
        f.write("[\n")
        for i, entry in enumerate(output_data):
            s_val_str = f"{entry[0]:.8e}"
            # IMPORTANT: Sort graphs by name before writing for consistent parsing later
            sorted_graph_entries = sorted(entry[1], key=lambda e: e[0])
            graphs_str_list = [f'  ["{g[0]}", {g[1]:.8e}]' for g in sorted_graph_entries]
            graphs_str = ",\n  ".join(graphs_str_list)
            f.write(f"  [{s_val_str}, [\n  {graphs_str}\n  ]]")
            f.write(",\n" if i < len(output_data) - 1 else "\n")
        f.write("]\n")
    print(f"Generated mock data in {filename}")

# --- Call once if needed ---
## import os
# if not os.path.exists("raw_data.py"): # Generate only if it doesn't exist
#    generate_mock_raw_data_file(N_supergraphs=124, N_energies=32) 
# Forcing regeneration for testing this script:
## generate_mock_raw_data_file(N_supergraphs=124, N_energies=32)
# ---

# --- 1. Load and Parse Data ---
try:
    with open("./raw_data.py", 'r') as f:
        raw_data_content = f.read()
    # Using eval is generally unsafe with untrusted input. 
    # For this specific case where we control the file generation, it's acceptable.
    # A safer method for arbitrary files would be ast.literal_eval or a JSON parser.
    loaded_data_list_full = eval(raw_data_content)
    # Ensure graphs within each energy point are sorted by name for consistent processing
    loaded_data_list_full = [[l_entry[0], sorted(l_entry[1], key=lambda e: e[0])] for l_entry in loaded_data_list_full]

except Exception as e:
    print(f"Error loading or parsing data from raw_data.py: {e}")
    sys.exit(1)

# Filter by sqrt(s) range
sqrt_s_all = np.array([entry[0] for entry in loaded_data_list_full])
s_min_actual = np.min(sqrt_s_all) if len(sqrt_s_all) > 0 else 0.0
s_max_actual = np.max(sqrt_s_all) if len(sqrt_s_all) > 0 else 0.0

s_filter_min = SQRT_S_MIN_DISPLAY if SQRT_S_MIN_DISPLAY is not None else s_min_actual
s_filter_max = SQRT_S_MAX_DISPLAY if SQRT_S_MAX_DISPLAY is not None else s_max_actual

loaded_data_list = [entry for entry in loaded_data_list_full if s_filter_min <= entry[0] <= s_filter_max]

if not loaded_data_list:
    print(f"No data points in the range [{s_filter_min:.1f}, {s_filter_max:.1f}]. Original range: [{s_min_actual:.1f}, {s_max_actual:.1f}]. Exiting.")
    sys.exit(0)

sqrt_s_values_loaded = np.array([entry[0] for entry in loaded_data_list])
# Graph names are taken from the first energy point, assuming consistency after sorting
graph_names_original_order = [graph_entry[0] for graph_entry in loaded_data_list[0][1]]
N_supergraphs_loaded = len(graph_names_original_order)
N_energies_loaded = len(sqrt_s_values_loaded)

# Populate the data matrix
data_matrix_loaded = np.zeros((N_supergraphs_loaded, N_energies_loaded))
for j_idx, entry in enumerate(loaded_data_list):
    if len(entry[1]) != N_supergraphs_loaded:
        raise ValueError(f"Inconsistent number of supergraphs at sqrt(s)={entry[0]}. Expected {N_supergraphs_loaded}, got {len(entry[1])}")
    for i_idx, graph_entry in enumerate(entry[1]):
        # This check relies on graphs being sorted by name within each energy point
        if graph_entry[0] != graph_names_original_order[i_idx]:
            raise ValueError(f"Inconsistent graph name or order at sqrt(s)={entry[0]}. Expected {graph_names_original_order[i_idx]} at index {i_idx}, got {graph_entry[0]}")
        data_matrix_loaded[i_idx, j_idx] = graph_entry[1]
print(f"Parsed filtered data: {N_supergraphs_loaded} SGs from {sqrt_s_values_loaded[0]:.1f} to {sqrt_s_values_loaded[-1]:.1f} GeV ({N_energies_loaded} energy points).")


# --- 2. Process Data for Plotting ---
heatmap_display_data = np.zeros((N_supergraphs_loaded, N_energies_loaded)) # This will store normalized values for heatmap display
epsilon_norm = 1e-9 # To avoid division by zero or tiny numbers in normalization
actual_anchor_s_val_for_title = None # For global sort title
sorted_indices_global = np.arange(N_supergraphs_loaded) # Default: original order
graph_names_for_heatmap_y_axis = list(graph_names_original_order) # Default graph names for y-axis

# Determine the data to be shown in the heatmap based on SORT_MODE
if SORT_MODE == 'global':
    if GLOBAL_SORT_ANCHOR_SQRT_S is not None and N_energies_loaded > 0:
        # Sort by absolute contribution at a specific anchor sqrt(s)
        anchor_idx = np.argmin(np.abs(sqrt_s_values_loaded - GLOBAL_SORT_ANCHOR_SQRT_S))
        anchor_contributions = np.abs(data_matrix_loaded[:, anchor_idx])
        sorted_indices_global = np.argsort(anchor_contributions)[::-1]
        actual_anchor_s_val_for_title = sqrt_s_values_loaded[anchor_idx]
    elif N_supergraphs_loaded > 0 and N_energies_loaded > 0: # Default global sort: by mean absolute contribution
        mean_abs_contributions_global = np.mean(np.abs(data_matrix_loaded), axis=1)
        sorted_indices_global = np.argsort(mean_abs_contributions_global)[::-1]
    
    # Apply the global sort to the data used for normalization and to graph names
    data_for_heatmap_normalization = data_matrix_loaded[sorted_indices_global, :]
    graph_names_for_heatmap_y_axis = [graph_names_original_order[i] for i in sorted_indices_global]
else: # 'local' sort mode
    # For local sort, normalization happens on locally sorted columns,
    # but heatmap_display_data rows still correspond to original/global order if displayed so.
    # Here, data_for_heatmap_normalization will be locally sorted inside the loop.
    data_for_heatmap_normalization = np.copy(data_matrix_loaded) # Will be modified column by column if local sort


# Normalize data for the heatmap display
for j in range(N_energies_loaded): # For each energy column
    if SORT_MODE == 'local':
        # For local sort, get the original column, sort it, normalize, then place back into heatmap_display_data
        # The rows of heatmap_display_data will represent ranks for that column
        current_column_original_data = data_matrix_loaded[:, j]
        sort_indices_local_j = np.argsort(np.abs(current_column_original_data))[::-1]
        column_data_to_normalize_this_iter = current_column_original_data[sort_indices_local_j]
    else: # 'global' sort
        # Data is already globally sorted by rows
        column_data_to_normalize_this_iter = data_for_heatmap_normalization[:, j]

    # Perform normalization on column_data_to_normalize_this_iter
    significant_pos_in_col = column_data_to_normalize_this_iter[column_data_to_normalize_this_iter > epsilon_norm]
    significant_neg_in_col = column_data_to_normalize_this_iter[column_data_to_normalize_this_iter < -epsilon_norm]
    
    normalizer_pos = np.max(significant_pos_in_col) if len(significant_pos_in_col) > 0 else 0.0
    normalizer_neg = np.min(significant_neg_in_col) if len(significant_neg_in_col) > 0 else 0.0 # This is negative

    for i in range(N_supergraphs_loaded): # For each graph/row in the (potentially locally sorted) column
        val = column_data_to_normalize_this_iter[i]
        scaled_val = 0.0
        if val > epsilon_norm:
            if normalizer_pos > epsilon_norm: # Check normalizer_pos itself
                scaled_val = val / normalizer_pos
        elif val < -epsilon_norm:
            if abs(normalizer_neg) > epsilon_norm: # Check abs(normalizer_neg)
                scaled_val = val / abs(normalizer_neg) # val is negative, abs(normalizer_neg) is positive -> scaled_val is negative
        
        # Store in heatmap_display_data. For 'local' sort, i is the rank.
        # For 'global' sort, i is the index in the globally sorted list.
        heatmap_display_data[i, j] = scaled_val


# Calculate aggregate sums for the second panel (always on original, unsorted data)
sum_positive_abs = np.sum(np.maximum(0, data_matrix_loaded), axis=0)
sum_negative_abs = np.sum(np.minimum(0, data_matrix_loaded), axis=0) # These are negative or zero
total_cross_section = np.sum(data_matrix_loaded, axis=0)

# Calculate relative contributions for the second panel
small_total_xs_threshold_for_ratio = 1e-15 # Avoid division by zero for tiny total_xs
rel_sum_positive = np.divide(sum_positive_abs, total_cross_section, 
                              out=np.zeros_like(sum_positive_abs), 
                              where=np.abs(total_cross_section) > small_total_xs_threshold_for_ratio)
rel_sum_negative = np.divide(sum_negative_abs, total_cross_section, 
                              out=np.zeros_like(sum_negative_abs), 
                              where=np.abs(total_cross_section) > small_total_xs_threshold_for_ratio)


# Calculate contributions from top percentiles for the third panel
sum_contrib_top_percentiles = {p: np.zeros(N_energies_loaded) for p in TOP_PERCENTILES_TO_PLOT}
rel_contrib_top_percentiles = {p: np.zeros(N_energies_loaded) for p in TOP_PERCENTILES_TO_PLOT}

if N_supergraphs_loaded > 0: # Ensure there are graphs to process
    for j in range(N_energies_loaded): # For each energy point
        current_s_contributions = data_matrix_loaded[:, j] # Use original data column
        abs_mags = np.abs(current_s_contributions)
        sorted_indices_this_s = np.argsort(abs_mags)[::-1] # Sort graphs by |contrib| at this s

        for p_val in TOP_PERCENTILES_TO_PLOT:
            N_top_p = max(1, int(p_val * N_supergraphs_loaded))
            # Sum actual contributions of these top graphs
            sum_contrib_top_percentiles[p_val][j] = np.sum(current_s_contributions[sorted_indices_this_s[:N_top_p]])
            
            # Calculate relative to total_cross_section at this s
            if np.abs(total_cross_section[j]) > small_total_xs_threshold_for_ratio:
                rel_contrib_top_percentiles[p_val][j] = sum_contrib_top_percentiles[p_val][j] / total_cross_section[j]
            else: # Handle cases where total_xs is zero or tiny
                rel_contrib_top_percentiles[p_val][j] = 0.0 # Or np.nan, or based on sum_contrib itself


# --- 3. Create the Visualization (4 panels) ---
fig, axes = plt.subplots(4, 1, figsize=(15, 22), sharex=True, 
                         gridspec_kw={'height_ratios': [3, 1.5, 1.5, 1.5]})
ax_heatmap, ax_lineplot, ax_percent_contrib, ax_total_xs_loglog = axes.flatten()

fig_title_sort_part = f"Sort: {SORT_MODE}"
if SORT_MODE == 'global' and actual_anchor_s_val_for_title is not None:
    fig_title_sort_part += f" (Anchor ~{actual_anchor_s_val_for_title:.0f} GeV)"
elif SORT_MODE == 'global':
    fig_title_sort_part += " (Mean |Contrib.|)"
# Main figure title is set later after all subplots are defined


# Panel 1: Heatmap
cmap_heatmap = plt.get_cmap('coolwarm_r') # Red positive, Blue negative, White zero

# Define edges for pcolormesh to center cells on data points
y_edges = np.arange(N_supergraphs_loaded + 1) - 0.5 # Centered on integer indices

if N_energies_loaded > 1:
    log_s_values = np.log(sqrt_s_values_loaded)
    # Midpoints between log(s) values
    log_s_edges_mid = (log_s_values[:-1] + log_s_values[1:]) / 2.0
    # Extrapolate first and last edges
    first_log_edge = log_s_values[0] - (log_s_edges_mid[0] - log_s_values[0]) if N_energies_loaded > 2 else log_s_values[0] - (log_s_values[1]-log_s_values[0])/2 if N_energies_loaded == 2 else log_s_values[0] - 0.1
    last_log_edge = log_s_values[-1] + (log_s_values[-1] - log_s_edges_mid[-1]) if N_energies_loaded > 2 else log_s_values[-1] + (log_s_values[-1]-log_s_values[-2])/2 if N_energies_loaded == 2 else log_s_values[-1] + 0.1
    x_log_edges = np.concatenate(([first_log_edge], log_s_edges_mid, [last_log_edge]))
elif N_energies_loaded == 1: # Single energy point
    x_log_edges = np.array([np.log(sqrt_s_values_loaded[0]) - 0.1, np.log(sqrt_s_values_loaded[0]) + 0.1]) # Arbitrary small width
else: # No energy points
    x_log_edges = np.array([-0.5, 0.5]) # Default for empty case
x_edges = np.exp(x_log_edges)


im = ax_heatmap.pcolormesh(x_edges, y_edges, heatmap_display_data, 
                           cmap=cmap_heatmap, vmin=-1.0, vmax=1.0, shading='flat')
ax_heatmap.set_xscale('log') # Set x-axis to log scale for the heatmap

heatmap_title_str = r'Individual Supergraph Contributions (Per-$\sqrt{s}$ Slice Normalization' # Base title
heatmap_ylabel_str = ''
ytick_positions = []
ytick_labels_display = []
current_ytick_label_fontsize = plt.rcParams['ytick.labelsize'] # Get default matplotlib fontsize

if SORT_MODE == 'local':
    heatmap_ylabel_str = f'Supergraph Rank (within each $\sqrt{{s}}$ slice)'
    heatmap_title_str += ' & Sort)'
    if N_supergraphs_loaded > 0:
        num_ytick_labels = min(10, N_supergraphs_loaded)
        ytick_positions = np.linspace(0, N_supergraphs_loaded - 1, num_ytick_labels, dtype=int)
        ytick_labels_display = [f"Rank {pos + 1}" for pos in ytick_positions]
else: # 'global' sort mode
    if actual_anchor_s_val_for_title is not None:
        sort_desc = fr'by |Contrib.| at $\sqrt{{s}} \approx {actual_anchor_s_val_for_title:.0f}$ GeV'
    else: # Default global: sort by mean abs contribution
        sort_desc = 'by Avg. |Contrib.| Across Displayed $\sqrt{s}$ Range'
    heatmap_ylabel_str = f'Supergraph (Order {sort_desc})'
    heatmap_title_str += f', Globally Sorted {sort_desc})'

    if (HEATMAP_GLOBAL_SORT_Y_LABELS is None or not HEATMAP_GLOBAL_SORT_Y_LABELS) and N_supergraphs_loaded > 0:
        # User wants ALL graph names labeled
        ytick_positions = np.arange(N_supergraphs_loaded) # A tick for every graph row
        ytick_labels_display = graph_names_for_heatmap_y_axis # Use the globally sorted names
        current_ytick_label_fontsize = HEATMAP_ALL_GRAPHS_FONTSIZE # Apply custom font size
    elif HEATMAP_GLOBAL_SORT_Y_LABELS and N_supergraphs_loaded > 0: # Non-empty list provided
        # User provided a specific list of names to label
        temp_ytick_pos = []
        temp_ytick_labels = []
        for name_to_label in HEATMAP_GLOBAL_SORT_Y_LABELS:
            try:
                # graph_names_for_heatmap_y_axis is already sorted according to sorted_indices_global
                idx = graph_names_for_heatmap_y_axis.index(name_to_label)
                temp_ytick_pos.append(idx) # idx is the row number in the heatmap
                temp_ytick_labels.append(name_to_label)
            except ValueError:
                print(f"Warning: Graph name '{name_to_label}' for Y-label not found in sorted list.")
        
        if temp_ytick_pos: # If any specified labels were found
            sorted_label_info = sorted(zip(temp_ytick_pos, temp_ytick_labels)) # Sort by row index
            ytick_positions = [info[0] for info in sorted_label_info]
            ytick_labels_display = [info[1] for info in sorted_label_info]
        else: # No specified labels found, fallback to Rank #1 and #N
            ytick_positions = [0, N_supergraphs_loaded - 1] if N_supergraphs_loaded > 1 else [0] if N_supergraphs_loaded == 1 else []
            ytick_labels_display = ["Rank #1", f"Rank #{N_supergraphs_loaded}"] if N_supergraphs_loaded > 1 else ["Rank #1"] if N_supergraphs_loaded == 1 else []
            if ytick_positions: # Add graph names if possible for Rank #1 and #N
                ytick_labels_display[0] += f" ({graph_names_for_heatmap_y_axis[ytick_positions[0]]})"
                if len(ytick_positions) > 1:
                    ytick_labels_display[1] += f" ({graph_names_for_heatmap_y_axis[ytick_positions[1]]})"

    elif N_supergraphs_loaded > 0: # Default for global if HEATMAP_GLOBAL_SORT_Y_LABELS is an empty list (but not None)
        ytick_positions = [0, N_supergraphs_loaded - 1] if N_supergraphs_loaded > 1 else [0] if N_supergraphs_loaded == 1 else []
        ytick_labels_display = ["Rank #1", f"Rank #{N_supergraphs_loaded}"] if N_supergraphs_loaded > 1 else ["Rank #1"] if N_supergraphs_loaded == 1 else []
        if ytick_positions: # Add graph names if possible for Rank #1 and #N
            ytick_labels_display[0] += f" ({graph_names_for_heatmap_y_axis[ytick_positions[0]]})"
            if len(ytick_positions) > 1:
                ytick_labels_display[1] += f" ({graph_names_for_heatmap_y_axis[ytick_positions[1]]})"


ax_heatmap.set_ylabel(heatmap_ylabel_str, fontsize=12)
if len(ytick_positions) > 0:
    ax_heatmap.set_yticks(ytick_positions) # Positions are 0 to N_supergraphs-1
    ax_heatmap.set_yticklabels(ytick_labels_display, fontsize=current_ytick_label_fontsize)

ax_heatmap.invert_yaxis() # Show rank 1 / top contributor at the top
ax_heatmap.set_title(heatmap_title_str, fontsize=14)
# ax_heatmap.tick_params(axis='y', labelsize=9) # Default, overridden by set_yticklabels fontsize

cbar = fig.colorbar(im, ax=ax_heatmap, orientation='vertical', fraction=0.046, pad=0.02) # Use default pad for multi-panel
cbar.set_label(r'Normalized Contribution (scaled to max +ve/-ve in each $\sqrt{s}$ slice)', fontsize=10)


# Panel 2: Relative Aggregate Contributions
ax_lineplot.plot(sqrt_s_values_loaded, rel_sum_positive, label='Sum of Positive Contrib. / Total', color='crimson', linestyle='-', lw=2)
ax_lineplot.plot(sqrt_s_values_loaded, rel_sum_negative, label='Sum of Negative Contrib. / Total', color='royalblue', linestyle='-', lw=2)
ax_lineplot.fill_between(sqrt_s_values_loaded, rel_sum_positive, 0, color='lightcoral', alpha=0.5, interpolate=True)
ax_lineplot.fill_between(sqrt_s_values_loaded, rel_sum_negative, 0, color='lightskyblue', alpha=0.5, interpolate=True)
ax_lineplot.set_ylabel('Relative Contribution to Total', fontsize=10)
ax_lineplot.grid(True, which='both', linestyle=':', linewidth=0.7)
ax_lineplot.legend(loc='best', fontsize=8)
ax_lineplot.set_title('Relative Aggregate Contributions & Cancellations', fontsize=12)
ax_lineplot.tick_params(axis='both', which='major', labelsize=9)
ax_lineplot.axhline(0, color='black', linewidth=0.7, linestyle='-')


# Panel 3: Relative Contribution from Top % Supergraphs
percentile_colors = ['darkgreen', 'darkorchid', 'sienna', 'teal', 'gold'] # Colors for different percentiles
if N_supergraphs_loaded > 0: # Only plot if there are supergraphs
    for i, p_val in enumerate(TOP_PERCENTILES_TO_PLOT):
        label = f'Sum Top {int(p_val*100)}% / Total' # e.g., "Sum Top 10% / Total"
        color = percentile_colors[i % len(percentile_colors)]
        ax_percent_contrib.plot(sqrt_s_values_loaded, rel_contrib_top_percentiles[p_val], label=label, color=color, linestyle='-')

ax_percent_contrib.axhline(1.0, label='Net Total / Total (=1)', color='black', linestyle='-', linewidth=1.5) # Line at y=1 for reference
ax_percent_contrib.set_ylabel('Sum / Total Net XS', fontsize=10)
ax_percent_contrib.grid(True, which='both', linestyle=':', linewidth=0.7)
ax_percent_contrib.legend(loc='best', fontsize=8)
ax_percent_contrib.set_title(r'Relative Contribution from Top % Supergraphs (Ranked per $\sqrt{s}$ by |Magnitude|)', fontsize=12)
ax_percent_contrib.tick_params(axis='both', which='major', labelsize=9)
ax_percent_contrib.axhline(0, color='black', linewidth=0.7, linestyle='-') # Zero line


# Panel 4: Total Cross-Section (Log-Log)
log_plot_threshold = 1e-9 # Threshold for masking zero or near-zero values for log plot

# Mask data for log plotting to avoid errors with non-positive values
has_plottable_positive_total_xs = np.any(total_cross_section > log_plot_threshold)
has_plottable_negative_total_xs = np.any(total_cross_section < -log_plot_threshold) 

# Prepare masked arrays for positive and absolute negative parts
total_xs_positive_part_masked = np.ma.masked_where(total_cross_section <= log_plot_threshold, total_cross_section)
total_xs_negative_part_abs_masked = np.ma.masked_where(total_cross_section >= -log_plot_threshold, np.abs(total_cross_section))

if has_plottable_positive_total_xs:
    ax_total_xs_loglog.plot(sqrt_s_values_loaded, total_xs_positive_part_masked, label='Total Net XS (Positive)', color='darkgreen', linestyle='-', lw=2)
if has_plottable_negative_total_xs: # If there are significant negative contributions
    ax_total_xs_loglog.plot(sqrt_s_values_loaded, total_xs_negative_part_abs_masked, label='|Total Net XS (Negative)|', color='darkred', linestyle='--', lw=2)

ax_total_xs_loglog.set_xscale('log')
ax_total_xs_loglog.set_yscale('log')
ax_total_xs_loglog.set_xlabel(r'$\sqrt{s}$ (GeV)', fontsize=12)
ax_total_xs_loglog.set_ylabel('Total Cross-Section [arb. units]', fontsize=10) # Adjust unit if known
if has_plottable_positive_total_xs or has_plottable_negative_total_xs: # Show legend only if something is plotted
    ax_total_xs_loglog.legend(loc='best', fontsize=8)
ax_total_xs_loglog.grid(True, which='both', linestyle=':', linewidth=0.7)
ax_total_xs_loglog.set_title('Total Net Cross-Section (Log-Log Scale)', fontsize=12)
ax_total_xs_loglog.tick_params(axis='both', which='major', labelsize=9)


# X-axis limits and ticks for the bottom-most plot (which controls shared x-axis)
if N_energies_loaded > 0 :
    # Set xlim for heatmap based on calculated x_edges to match pcolormesh extent
    if x_edges is not None and len(x_edges) >= 2:
        ax_heatmap.set_xlim(x_edges[0], x_edges[-1])
        # Since x-axes are shared, setting xlim on one should affect others if they are truly shared.
        # However, explicit setting for the bottom plot ensures it's correct.
        for ax_shared in [ax_lineplot, ax_percent_contrib, ax_total_xs_loglog]:
            ax_shared.set_xlim(x_edges[0], x_edges[-1])


    min_s_disp, max_s_disp = np.nanmin(sqrt_s_values_loaded), np.nanmax(sqrt_s_values_loaded)
    if N_energies_loaded > 1 and max_s_disp / min_s_disp > 10 : # Check if log scale is meaningful
        # Use LogLocator for major ticks on the bottom-most plot
        ax_total_xs_loglog.xaxis.set_major_locator(LogLocator(base=10.0, numticks=15)) # Auto numticks
        # Define minor ticks based on major tick range
        subs_major_range = np.log10(max_s_disp) - np.log10(min_s_disp)
        subs_for_minor = np.arange(2,10,1) * 0.1 if subs_major_range > 1.0 else None # More minor ticks for wider range
        if subs_for_minor is not None:
            ax_total_xs_loglog.xaxis.set_minor_locator(LogLocator(base=10.0, subs=subs_for_minor))
        else: # Fallback if range is small
            ax_total_xs_loglog.xaxis.set_minor_locator(LogLocator(base=10.0, subs=(0.2,0.4,0.6,0.8)))
        ax_total_xs_loglog.xaxis.set_minor_formatter(NullFormatter()) # Hide minor tick labels
    else: # Fallback for linear scale or very few points
        ax_total_xs_loglog.xaxis.set_major_locator(MaxNLocator(nbins=max(2, min(5, N_energies_loaded)), prune='both'))


# Overall Figure Title and Layout Adjustments
fig.suptitle(fr'$\gamma\gamma \to t\bar{{t}}$ Supergraph Contributions ({fig_title_sort_part})', fontsize=18)
plt.tight_layout(rect=[0, 0.01, 1, 0.96]) # rect to make space for suptitle

# Fine-tune subplot widths to align data areas (experimental)
# This needs to be done after initial tight_layout and drawing to get correct positions
fig.canvas.draw() # Ensure positions are calculated

heatmap_bbox = ax_heatmap.get_position()
h_left = heatmap_bbox.x0
h_width = heatmap_bbox.width

lower_plot_target_width = h_width * LOWER_INSET_WIDTH_SCALE
# Calculate the left position to center the scaled lower plots relative to the heatmap's data area
lower_plot_left_centering_offset = (h_width * (1.0 - LOWER_INSET_WIDTH_SCALE)) / 2.0
lower_plot_actual_target_left = h_left + lower_plot_left_centering_offset

for ax_lower in [ax_lineplot, ax_percent_contrib, ax_total_xs_loglog]:
    current_lower_bbox = ax_lower.get_position()
    ax_lower.set_position([
        lower_plot_actual_target_left, 
        current_lower_bbox.y0, 
        lower_plot_target_width, 
        current_lower_bbox.height
    ])


# Save the figure
save_name_suffix = f"sort_{SORT_MODE}"
if SORT_MODE == 'global' and GLOBAL_SORT_ANCHOR_SQRT_S is not None:
    save_name_suffix += f"_anchor_{GLOBAL_SORT_ANCHOR_SQRT_S:.0f}"
elif SORT_MODE == 'global': # Default global sort by mean
    save_name_suffix += "_anchor_mean"
save_name_suffix += f"_insetScale_{LOWER_INSET_WIDTH_SCALE:.2f}"

plt.savefig(f"supergraph_visualization_final_{save_name_suffix}.png", dpi=300)
print(f"Saved plot to supergraph_visualization_final_{save_name_suffix}.png")
plt.show()
