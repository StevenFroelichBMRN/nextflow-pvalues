#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.input_dir = "s3://r6333-pep-nppc-oi-bmn333-dev/gwas_testing_nwk/Gene_Editing/arc-virtual-cell-atlas/STATE"
params.output_dir = "s3://r6333-pep-nppc-oi-bmn333-dev/gwas_testing_nwk/Gene_Editing/arc-virtual-cell-atlas/STATE/pvalues"
params.log2fc_threshold = 0.25
params.min_cells_treated = 5
params.min_cells_control = 10

process extract_pvalues {
    tag "${h5ad_file.simpleName}"
    memory { 64.GB * task.attempt }
    cpus 4
    time '6h'
    maxRetries 2
    errorStrategy { task.exitStatus in [137, 143] ? 'retry' : 'finish' }

    publishDir params.output_dir, mode: 'copy'

    input:
    path h5ad_file

    output:
    path "${h5ad_file.simpleName}_pvalues.parquet", emit: results
    path "${h5ad_file.simpleName}_summary.json", emit: summary

    script:
    """
    extract_pvalues.py \
        --input ${h5ad_file} \
        --output ${h5ad_file.simpleName}_pvalues.parquet \
        --summary ${h5ad_file.simpleName}_summary.json \
        --log2fc-threshold ${params.log2fc_threshold} \
        --min-cells-treated ${params.min_cells_treated} \
        --min-cells-control ${params.min_cells_control} \
        --threads ${task.cpus}
    """
}

process merge_summaries {
    memory '4 GB'
    cpus 1
    time '30m'

    publishDir params.output_dir, mode: 'copy'

    input:
    path summaries

    output:
    path "run_summary.json"

    script:
    """
    #!/usr/bin/env python3
    import json, glob

    all_summaries = []
    for f in sorted(glob.glob("*_summary.json")):
        with open(f) as fh:
            all_summaries.append(json.load(fh))

    total_rows = sum(s.get("total_results", 0) for s in all_summaries)
    total_errors = sum(1 for s in all_summaries if s.get("status") == "error")

    combined = {
        "total_files": len(all_summaries),
        "total_result_rows": total_rows,
        "total_errors": total_errors,
        "per_file": all_summaries
    }

    with open("run_summary.json", "w") as out:
        json.dump(combined, out, indent=2)
    """
}

workflow {
    h5ad_files = Channel.fromPath("${params.input_dir}/c*.h5ad")

    extract_pvalues(h5ad_files)

    merge_summaries(extract_pvalues.out.summary.collect())
}
