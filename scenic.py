#!/usr/bin/env python3

import argparse
import os
import subprocess
import sys
import pandas as pd
import scanpy as sc
import anndata as ad
import pickle
import scipy.sparse

try:
    from pyscenic.export import export2loom
except ImportError:
    print("Error: pyscenic.export module not found. Please ensure pySCENIC is installed correctly.")
    sys.exit(1)


def run_command(command, log_file_path):
    """Executes a shell command and logs its output to a file."""
    print(f"Executing: {' '.join(command)}")
    os.makedirs(os.path.dirname(log_file_path), exist_ok=True)
    try:
        with open(log_file_path, 'w') as log_file:
            process = subprocess.Popen(command, stdout=log_file, stderr=subprocess.STDOUT, text=True)
            process.wait()
        
        if process.returncode != 0:
            print(f"Error: Command failed with return code {process.returncode}. Check log: {log_file_path}")
            sys.exit(process.returncode)
    except Exception as e:
        print(f"An exception occurred while running command: {e}. Check log: {log_file_path}")
        sys.exit(1)


def preprocess_anndata(input_h5ad_path, output_dir, min_genes, min_cells, norm_target_sum_str, log_dir):
    """Preprocesses AnnData and extracts expression matrix."""
    print("\n--- Starting Preprocessing ---")
    output_h5ad_path = os.path.join(output_dir, "preprocessed_data.h5ad")
    output_tsv_path = os.path.join(output_dir, "expr_matrix_for_scenic.tsv")
    py_log_file_path = os.path.join(log_dir, "preprocess_anndata_script.log")

    os.makedirs(output_dir, exist_ok=True)
    os.makedirs(log_dir, exist_ok=True)

    with open(py_log_file_path, 'w') as log_f:
        try:
            adata = ad.read_h5ad(input_h5ad_path)
            log_f.write(f"Read AnnData from {input_h5ad_path}\n")

            sc.pp.filter_cells(adata, min_genes=min_genes)
            sc.pp.filter_genes(adata, min_cells=min_cells)
            log_f.write(f"Filtered cells (min_genes={min_genes}) and genes (min_cells={min_cells})\n")

            adata_for_tsv = adata.copy()
            
            actual_target_sum = None
            if norm_target_sum_str and norm_target_sum_str.lower() != 'none':
                try:
                    actual_target_sum = float(norm_target_sum_str)
                except ValueError:
                    warning_msg = f"Warning: could not parse norm_target_sum '{norm_target_sum_str}'. Skipping normalization for SCENIC TSV."
                    log_f.write(warning_msg + "\n")
                    print(warning_msg)
            
            if actual_target_sum is not None:
                sc.pp.normalize_total(adata_for_tsv, target_sum=actual_target_sum)
                log_f.write(f"Normalized total for TSV to target_sum={actual_target_sum}\n")

            if scipy.sparse.issparse(adata_for_tsv.X):
                expr_matrix_for_tsv = pd.DataFrame(adata_for_tsv.X.toarray(), index=adata_for_tsv.obs_names, columns=adata_for_tsv.var_names).T
            else:
                expr_matrix_for_tsv = pd.DataFrame(adata_for_tsv.X, index=adata_for_tsv.obs_names, columns=adata_for_tsv.var_names).T
            
            expr_matrix_for_tsv.to_csv(output_tsv_path, sep='\t')
            log_f.write(f"Expression matrix for SCENIC written to {output_tsv_path}\n")

            if actual_target_sum is not None and 'normalize_total' not in adata.uns_keys():
                 sc.pp.normalize_total(adata, target_sum=actual_target_sum)
            sc.pp.log1p(adata)
            adata.write_h5ad(output_h5ad_path)
            log_f.write(f"Preprocessed AnnData written to {output_h5ad_path}\n")

        except Exception as e:
            error_msg = f"Error during preprocessing: {e}. Detailed log: {py_log_file_path}"
            log_f.write(f"Error: {e}\n")
            print(error_msg)
            sys.exit(1)
            
    print(f"--- Preprocessing Complete ---")
    return output_h5ad_path, output_tsv_path


def run_pyscenic_grn(exp_matrix_tsv_path, tf_file_path, output_dir, grn_method, num_workers, log_dir):
    print("\n--- Starting pySCENIC GRN ---")
    adj_file_path = os.path.join(output_dir, "adj.tsv")
    cli_log_file_path = os.path.join(log_dir, "pyscenic_grn_cli.log")

    command = [
        "pyscenic", "grn",
        exp_matrix_tsv_path,
        tf_file_path,
        "--output", adj_file_path,
        "--method", grn_method,
        "--num_workers", str(num_workers)
    ]
    run_command(command, cli_log_file_path)
    print(f"--- pySCENIC GRN Complete ---")
    return adj_file_path


def run_pyscenic_ctx(adj_file_path, exp_matrix_tsv_path, motif_annotations_file, motif_database_paths, output_dir, num_workers, log_dir):
    print("\n--- Starting pySCENIC CTX (Regulon Prediction) ---")
    regulons_pickle_path = os.path.join(output_dir, "regulons.p")
    cli_log_file_path = os.path.join(log_dir, "pyscenic_ctx_cli.log")
    
    command = [
        "pyscenic", "ctx",
        adj_file_path
    ]
    command.extend(motif_database_paths)
    command.extend([
        "--annotations_fname", motif_annotations_file,
        "--expression_mtx_fname", exp_matrix_tsv_path,
        "--output", regulons_pickle_path,
        "--num_workers", str(num_workers),
        "--mode", "custom_multiprocessing",
        "--mask_dropouts"
    ])
    run_command(command, cli_log_file_path)
    print(f"--- pySCENIC CTX Complete ---")
    return regulons_pickle_path


def run_pyscenic_aucell(exp_matrix_tsv_path, regulons_pickle_path, output_dir, num_workers, seed, log_dir):
    print("\n--- Starting pySCENIC AUCell ---")
    auc_mtx_path = os.path.join(output_dir, "auc_mtx.csv")
    cli_log_file_path = os.path.join(log_dir, "pyscenic_aucell_cli.log")

    command = [
        "pyscenic", "aucell",
        exp_matrix_tsv_path,
        regulons_pickle_path,
        "--output", auc_mtx_path,
        "--num_workers", str(num_workers),
        "--seed", str(seed)
    ]
    run_command(command, cli_log_file_path)
    print(f"--- pySCENIC AUCell Complete ---")
    return auc_mtx_path


def create_final_loom(preprocessed_h5ad_path, auc_mtx_csv_path, regulons_pickle_path, output_dir, log_dir):
    """Creates the final loom file."""
    print("\n--- Starting Loom File Creation ---")
    final_loom_path = os.path.join(output_dir, "pyscenic_output.loom")
    py_log_file_path = os.path.join(log_dir, "create_loom_script.log")

    with open(py_log_file_path, 'w') as log_f:
        try:
            adata = ad.read_h5ad(preprocessed_h5ad_path)
            log_f.write(f"Read preprocessed AnnData from {preprocessed_h5ad_path}\n")

            auc_mtx_df = pd.read_csv(auc_mtx_csv_path, index_col=0)
            log_f.write(f"Read AUCell matrix from {auc_mtx_csv_path}\n")

            with open(regulons_pickle_path, 'rb') as f:
                regulons_obj = pickle.load(f)
            log_f.write(f"Read regulons from {regulons_pickle_path}\n")
            
            export2loom(
                ex_mtx=adata,
                auc_mtx=auc_mtx_df,
                regulons=regulons_obj,
                out_fname=final_loom_path
            )
            log_f.write(f"Loom file created successfully at {final_loom_path}\n")

        except Exception as e:
            error_msg = f"Error during loom creation: {e}. Detailed log: {py_log_file_path}"
            log_f.write(f"Error: {e}\n")
            print(error_msg)
            sys.exit(1)
            
    print(f"--- Loom File Creation Complete ---")
    return final_loom_path


def main():
    parser = argparse.ArgumentParser(description="Run a streamlined SCENIC pipeline.")

    # Input files
    parser.add_argument("--input_h5ad", type=str, required=True, help="Path to the input AnnData h5ad file.")
    parser.add_argument("--tf_file", type=str, required=True, help="Path to the list of transcription factors.")
    parser.add_argument("--motif_annotations", type=str, required=True, help="Path to motif annotations file.")
    parser.add_argument("--motif_databases", type=str, nargs='+', required=True, help="Path(s) to feather motif database(s).")

    # Output directory
    parser.add_argument("--output_dir", type=str, required=True, help="Directory to store all outputs and logs.")

    # Preprocessing parameters
    parser.add_argument("--min_genes", type=int, default=200, help="Min genes per cell for filtering.")
    parser.add_argument("--min_cells", type=int, default=3, help="Min cells per gene for filtering.")
    parser.add_argument("--norm_target_sum", type=str, default="10000", help="Target sum for normalization.")

    # SCENIC tool parameters
    parser.add_argument("--grn_method", type=str, default="grnboost2", choices=["grnboost2", "genie3"], help="Method for GRN inference.")
    parser.add_argument("--grn_workers", type=int, default=4, help="Number of workers for pySCENIC GRN.")
    parser.add_argument("--ctx_workers", type=int, default=4, help="Number of workers for pySCENIC CTX.")
    parser.add_argument("--aucell_workers", type=int, default=4, help="Number of workers for pySCENIC AUCell.")
    parser.add_argument("--aucell_seed", type=int, default=123, help="Seed for AUCell reproducibility.")

    args = parser.parse_args()

    log_dir = os.path.join(args.output_dir, "pipeline_logs")
    os.makedirs(log_dir, exist_ok=True)

    print(f"Starting SCENIC pipeline. Outputs will be in: {args.output_dir}")
    
    preprocessed_h5ad_path, scenic_exp_matrix_tsv = preprocess_anndata(
        args.input_h5ad, args.output_dir, args.min_genes, args.min_cells, args.norm_target_sum, log_dir
    )

    adj_file_path = run_pyscenic_grn(
        scenic_exp_matrix_tsv, args.tf_file, args.output_dir, args.grn_method, args.grn_workers, log_dir
    )

    regulons_pickle_path = run_pyscenic_ctx(
        adj_file_path, scenic_exp_matrix_tsv, args.motif_annotations, args.motif_databases, 
        args.output_dir, args.ctx_workers, log_dir
    )

    auc_mtx_path = run_pyscenic_aucell(
        scenic_exp_matrix_tsv, regulons_pickle_path, args.output_dir, args.aucell_workers, args.aucell_seed, log_dir
    )

    final_loom_path = create_final_loom(
        preprocessed_h5ad_path, auc_mtx_path, regulons_pickle_path, args.output_dir, log_dir
    )

    print(f"\nSCENIC Pipeline finished successfully!")
    print(f"Main outputs in: {args.output_dir}")
    print(f"  Preprocessed AnnData: {preprocessed_h5ad_path}")
    print(f"  Final Loom file: {final_loom_path}")

if __name__ == "__main__":
    main()