# Neoantigen-Prioritization-Pipeline

[![Python 3.10+](https://img.shields.io/badge/python-3.10%2B-blue.svg)](https://www.python.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](https://opensource.org/licenses/MIT)

NeoScorer is a research prototype for ranking candidate HLA Class I peptides from pre-annotated somatic mutations and gene-level expression values. It includes a command-line pipeline and a local dashboard for exploring reports from multiple patients.

The workflow filters missense variants, resolves protein sequences, generates mutant and wild-type 9-mers, runs selectable predictors, and applies a heuristic ranking score. It does not process raw WES/RNA-seq reads or produce validated vaccine formulations. Predictions and tiers are exploratory; no experimental or clinical validation of this pipeline is documented in the repository.

## Contents

- [Implemented workflow](#implemented-workflow)
- [Scoring and interpretation](#scoring-and-interpretation)
- [Bundled TCGA example](#bundled-tcga-example)
- [Local dashboard](#local-dashboard)
- [Repository structure](#repository-structure)
- [Quickstart and CLI](#quickstart-and-cli)
- [Limitations](#limitations)
- [Verification](#verification)
- [Data provenance and references](#data-provenance-and-references)

## Implemented workflow

1. **Read and validate JSON inputs.** Mutation records use a cBioPortal-compatible schema with gene identifiers, mutation type, protein change, and tumor alternate/reference read counts. Expression maps Entrez gene IDs to finite, non-negative values. Patient identity is inferred from mutations or supplied explicitly.
2. **Filter missense variants by expression.** The default cutoff is `1.0`. Missing gene-expression entries are treated as zero and filtered at this cutoff. The code calls this quantity TPM, but users must verify input units; see [limitations](#limitations).
3. **Resolve protein context.** A local cache or UniProt reviewed human sequence search supplies a gene-level reference sequence. Resolution checks the one-based protein position and reference residue, then extracts up to ten residues on each side. Unresolved sequences, unsupported protein-change strings, out-of-range positions, and residue mismatches share an exclusion counter. Uncached genes are prefetched with up to 16 workers.
4. **Generate and score peptides.** Sliding windows retain 9-mers spanning the altered residue, paired with wild-type sequences. One HLA allele is evaluated per run, defaulting to `HLA-A*02:01`.
5. **Select and rank.** The lowest predicted mutant IC50 selects one peptide per screened variant. The composite score then ranks these representatives. Reports contain all selected variants, including deprioritized rows; `--top` only limits console output.
6. **Export and explore.** Current runs write a TSV and JSON metadata sidecar. The dashboard discovers reports, plots metrics, filters/sorts candidates, exports CSV/TSV, and runs one background pipeline job at a time.

### Predictor modes

| Mode | Current implementation |
| --- | --- |
| `ensemble` (default) | MHCflurry scores all unique peptides. NetMHCpan is requested for each variant's best MHCflurry mutant peptide and its wild-type partner, plus pairs with mutant MHCflurry IC50 ≤2,500 nM. Available paired predictions are combined. Other peptides retain MHCflurry results with source `MHCflurry_Filter_NonBinder`. |
| `mhcflurry` | Local `Class1PresentationPredictor`, returning affinity, integrated presentation score, and presentation percentile. Requires separately installed model weights. |
| `netmhcpan` | IEDB MHC-I API requests using `method=netmhcpan_ba`, with a persistent cache. The wrapper labels results NetMHCpan-4.1 but does not pin the remote service version. |
| `pwm` | Hand-specified position weight matrix and heuristic affinity/rank mappings for HLA-A*02:01 9-mers. An illustrative baseline, without documented empirical calibration. |

MHCflurry's presentation score combines affinity and processing predictions ([official documentation](https://openvax.github.io/mhcflurry/api.html)). This pipeline supplies peptides without flanking context to the model and does not export separate cleavage or TAP transport probabilities. Combining predictors here is a project-defined heuristic, not independent validation of presentation or immunogenicity.

## Scoring and interpretation

The CLI uses functions in `src/prioritizer.py` to calculate:

$$S_{neo} = \operatorname{round}\left(\min\left(100, C_{AF} E_{TPM} P_{MHC} A_{bonus} \times 20\right),1\right)$$

These are heuristic weights, not fitted probabilities of T-cell response. Scores can saturate at 100 and are not calibrated for comparisons across patients.

### Allele-frequency term

$$AF = \frac{\text{tumor alternate reads}}{\max(1,\text{alternate reads}+\text{reference reads})}$$
$$C_{AF} = \min(1,\max(0,2AF))$$

This is a VAF-based proxy. The pipeline does not correct for tumor purity or copy number and does not estimate cancer-cell fraction. A zero-read record yields AF zero.

### Expression term

$$E_{TPM} = \log_2(1+\max(0,x))$$

Here $x$ is the supplied gene-expression value. The historical variable/column name is TPM; the bundled fetcher does not verify or convert expression units. Gene-level expression does not establish expression of the mutant allele.

### Affinity-derived term

$$P_{MHC} = \max\left(0,\min\left(1,1-\frac{\log_{10}(\min(50000,\max(1,IC_{50})))}{\log_{10}(50000)}\right)\right)$$

This is an affinity transformation, not a calibrated presentation probability. It is separate from the MHCflurry `Presentation_Score` exported in the TSV. Predicted IC50 values are in nM; lower values indicate stronger predicted binding. Affinity-only classes use ≤50 nM for strong, >50 to ≤500 nM for weak, and >500 nM for non-binders.

### Mutant/wild-type affinity ratio

$$AI = \frac{IC_{50}^{WT}}{\max(0.01,IC_{50}^{MT})}$$

For $AI \ge 1$:

$$A_{bonus} = \min(1.5,1+0.15\log_2(AI+1))$$

For $AI < 1$:

$$A_{bonus} = \max(0.75,0.90+0.10AI)$$

The ratio measures predicted differential binding. It does not measure TCR recognition, self-tolerance, or autoimmune risk.

### Combined predictions and tiers

When paired model results are available, the ensemble uses:

$$IC_{50}^{combined}=\sqrt{\max(0.01,IC_{50}^{MHCflurry})\max(0.01,IC_{50}^{NetMHCpan})}$$
$$Rank^{combined}=\frac{Rank^{MHCflurry}+Rank^{NetMHCpan}}{2}$$

MHCflurry presentation percentiles and NetMHCpan binding-affinity percentiles describe different outputs. Their arithmetic mean has no documented calibration. The MHCflurry integrated presentation score is retained, rather than averaged.

The active CLI in `run_pipeline.py` assigns Tier 1 when score ≥50 **and** mutant IC50 ≤50 nM, Tier 2 when score ≥25 **and** mutant IC50 ≤500 nM (after Tier 1), and Tier 3 otherwise. These are implementation thresholds, not clinical actionability categories. The separate `score_and_rank_candidates` helper uses a different Tier 2 rule (`or`); it is not the CLI's tier-assignment path.

Existing filenames and output fields retain terms such as `clinical_vaccine_prescription`, `ClinicalTier`, `Vaccine Payload`, and `INCLUDE`/`REJECT`. These are legacy software labels, not prescriptions or treatment recommendations. A high score alone does not imply a high priority tier: expression can raise the score of a weakly binding peptide.

## Bundled TCGA example

The repository includes mutation and expression JSONs for sample `TCGA-D3-A1Q1-06` from the cBioPortal TCGA SKCM PanCancer Atlas study, and a saved report for `TCGA-D3-A1Q1`. This is a worked example using public cancer data, not a validation cohort or evidence of vaccine efficacy. The default HLA-A*02:01 setting is an analysis assumption; the inputs do not establish this patient's HLA genotype.

These counts can be read from the bundled inputs and TSV without rerunning prediction:

| Quantity | Count | Interpretation |
| --- | ---: | --- |
| Input mutation records | 486 | Supplied annotated variants, not a measured mutation burden per megabase |
| Missense records | 410 | Variants eligible for this workflow |
| Missense records below expression cutoff 1.0 | 66 | Includes missing entries treated as zero |
| Saved report rows | 326 | One representative peptide per exported variant |
| Affinity-only strong binders | 28 | Mutant IC50 ≤50 nM |
| Affinity-only weak binders | 51 | Mutant IC50 >50 and ≤500 nM |
| Dashboard strong-binder count | 29 | IC50 ≤50 nM **or** rank in 0–0.5%; a different definition |
| Saved Tier 1 / Tier 2 / Tier 3 rows | 6 / 35 / 285 | Recorded heuristic assignments |
| Rows labeled combined / MHCflurry-only filtered | 211 / 115 | Not every row contains two model predictions |

The saved TSV has no JSON sidecar or explicit `HLA_Allele` column. Current runs add both. Historical claims of 18 QC exclusions and 2,912 candidate windows are not recorded in this report's metadata, so they are not presented as verified funnel statistics.

These six rows carry Tier 1 labels in the saved report. Rank is their position in the full score-sorted report, not a new ranking of this subset:

| Saved rank | Gene | Change | Mutant 9-mer | Score |
| ---: | --- | --- | --- | ---: |
| 1 | AP4B1 | G525V | LLLVVIDEV | 100.0 |
| 4 | TNFRSF10B | R145L | FLEEDSPEM | 82.9 |
| 7 | TMEM245 | G754W | FLWTYWAAV | 64.3 |
| 8 | PCNX3 | L1607I | SLEPFIYGL | 61.2 |
| 9 | GPR158 | L535F | RMLAVILFV | 59.8 |
| 12 | PML | S388Y | RLQDLSYCI | 52.3 |

For example, the saved TNFRSF10B row has a predicted wild-type/mutant affinity ratio of 340.78. This is a model output, not measured binding or demonstrated T-cell recognition. The repository does not establish structural mechanisms or prove that excluded peptides are biologically unpresentable.

## Local dashboard

The FastAPI application uses vanilla JavaScript and bundled Chart.js. Existing reports can be explored without a frontend build or internet connection. Features include patient discovery/switching, expression/affinity and mutant/wild-type charts, search, tier filters, multi-column sorting, model evidence cards, CSV/TSV export, and background runs from selected or uploaded JSON inputs.

It binds to `127.0.0.1:8520` and supports one active job per server. It is a local single-operator application without the authentication, audit infrastructure, or validation of a deployed clinical service. See [DASHBOARD.md](DASHBOARD.md) for input schemas, report interpretation, job behavior, and verification commands.

## Repository structure

```text
data/                    # Bundled example inputs and protein/prediction caches
reports/                 # Saved TSV; new runs also write JSON metadata
src/
  data_models.py         # Variant/context dataclasses
  protein_resolver.py    # Gene-level protein lookup and residue checks
  hla_scorer.py          # 9-mer generation and HLA-A*02:01 PWM baseline
  dnn_predictor.py        # MHCflurry and IEDB predictor wrappers
  scale_hla_screener.py  # Expression filtering, peptide screening and selection
  prioritizer.py         # Heuristic scoring functions
scripts/                 # cBioPortal fetcher and input inspection utility
dashboard_backend/       # Report parsing, discovery and background jobs
dashboard_static/        # Frontend and bundled Chart.js
tests/test_dashboard.py  # API, identity, exports and job regression tests
patient_io.py            # Shared input validation and report naming
run_pipeline.py          # CLI scoring, tier assignment and output generation
dashboard.py             # Local web entrypoint
requirements.txt         # Runtime dependencies
requirements-dev.txt     # Development/test dependencies
DASHBOARD.md             # Dashboard operation and technical details
```

## Quickstart and CLI

Use Python 3.10+ and an environment compatible with the predictor dependencies. In this `C:\Tools` workspace, use the shared environment at `C:\Tools\.venv` by default. From the project directory:

```powershell
& "C:\Tools\.venv\Scripts\python.exe" -m pip install -r requirements.txt
& "C:\Tools\.venv\Scripts\python.exe" dashboard.py --open-browser
```

On Windows the existing workspace launcher is `C:\Tools\LAUNCHERS\Launch_NeoScorer_Dashboard.bat`. Outside this workspace, substitute your Python interpreter. Installing packages alone does not install MHCflurry model weights; follow the [MHCflurry installation instructions](https://openvax.github.io/mhcflurry/intro.html) for the presentation bundle required by `Class1PresentationPredictor`.

Run the example through the CLI:

```powershell
& "C:\Tools\.venv\Scripts\python.exe" run_pipeline.py `
  --mutations data/real_patient_TCGA-D3-A1Q1-06_mutations.json `
  --expression data/real_patient_TCGA-D3-A1Q1-06_expression.json `
  --hla "HLA-A*02:01" --predictor ensemble `
  --tpm-threshold 1.0 --top 12 --output-dir reports/example_rerun
```

The separate output directory preserves the bundled report. The CLI replaces existing outputs with the same patient-derived name in its chosen directory; the dashboard refuses to overwrite existing patient reports.

| Option | Meaning and default |
| --- | --- |
| `--predictor` | `ensemble`, `mhcflurry`, `netmhcpan`, or `pwm`; default `ensemble` |
| `--hla` | One target allele; default `HLA-A*02:01`. PWM supports only this allele. |
| `--tpm-threshold` | Cutoff on the supplied expression scale; default `1.0` |
| `--top` | Console display count, not an export limit; default `12` |
| `--mutations`, `--expression` | Input JSON paths; default bundled example files |
| `--patient-id` | Optional ID, inferred from mutations if omitted |
| `--output-dir` | Report directory; default `reports/` |

## Limitations

- **Scientific scope:** There is no documented held-out benchmark, experimental peptide-binding/presentation assay, T-cell response validation, or clinical evaluation of this ranking. Predictor publications do not validate the pipeline's scoring weights or tiers. It does not model vaccine manufacture, delivery, dose, efficacy, or safety.
- **Expression units and missingness:** `scripts/fetch_patient_data.py` copies values from `skcm_tcga_pan_can_atlas_2018_rna_seq_v2_mrna` without a TPM conversion. cBioPortal describes TCGA RNASeqV2 values as RSEM normalized expression ([FAQ](https://docs.cbioportal.org/user-guide/faq/)). The bundled values should not be described as verified TPM. Their scale affects the gate, scores, and chart labels. Missing entries are conflated with zero; gene-level abundance does not establish mutant-transcript abundance. Expression JSON has no patient identity, so sample matching remains the user's responsibility.
- **Variants and reference mapping:** Only parseable single-amino-acid missense changes and 9-mers are screened. Indels, frameshifts, fusions, splice-derived peptides, other lengths, and HLA Class II are outside this implementation. UniProt lookup selects one reviewed human search result per gene, without transcript/isoform-aware mapping or pinned sequence releases. Residue matching is a local check, not full annotation validation; excluded variants may still be biologically relevant.
- **HLA and tumor context:** One user-specified allele is used per run. There is no HLA typing, multi-allele patient aggregation, purity/copy-number correction, HLA-loss assessment, or model of antigen-processing defects. VAF and expression terms are proxies, not evidence of tumor-wide target coverage.
- **Predictor coverage and failures:** The ensemble prefilter leaves some peptides with MHCflurry-only results. IEDB exceptions can produce fallback IC50 50,000 nM/rank 50%; missing or malformed responses can also lead to defaults downstream. Such values can influence scores and are not observed negative binding results. Combined source labels do not reliably preserve all fallback provenance. Cache keys lack model-version metadata, and dependency/service changes can change reruns.
- **Selection and traceability:** The report saves only the best-affinity peptide per variant, not every screened window. Exclusions are aggregate counts rather than a per-variant reason log. New sidecars record filenames, parameters, timestamp and funnel statistics, but not input hashes, sequence releases, or full model/dependency versions. The bundled legacy report has less provenance.
- **Local operation:** Viewing is local; prediction runs can send peptide sequences/HLA to IEDB and gene queries to UniProt. Inputs and logs persist under `.dashboard/jobs/`. The current protein/data fetchers disable TLS certificate verification, which limits remote-data integrity assurance. The dashboard is not designed for public exposure or a shared clinical workflow.

## Verification

Install development dependencies if needed and run:

```powershell
& "C:\Tools\.venv\Scripts\python.exe" -m pip install -r requirements-dev.txt
& "C:\Tools\.venv\Scripts\python.exe" -m pytest tests/ -q
```

Existing tests cover dashboard/API behavior, patient identity, validation, filtering/export, overwrite protection, failed jobs, and a subprocess CLI run with a synthetic variant excluded by expression. They do not establish full ensemble inference correctness, scientific ranking performance, or clinical validity. See [DASHBOARD.md](DASHBOARD.md) for additional static checks and manual browser scenarios.

## Data provenance and references

The fetch utility specifies cBioPortal study `skcm_tcga_pan_can_atlas_2018`, sample `TCGA-D3-A1Q1-06`, mutation profile `skcm_tcga_pan_can_atlas_2018_mutations`, and expression profile `skcm_tcga_pan_can_atlas_2018_rna_seq_v2_mrna`. Cached sequences and predictions support reuse but are not a complete reproducibility manifest.

- Hoadley et al. *Cell-of-Origin Patterns Dominate the Molecular Classification of 10,000 Tumors from 33 Types of Cancer*. Cell (2018). [DOI](https://doi.org/10.1016/j.cell.2018.03.022).
- O'Donnell et al. *MHCflurry 2.0: Improved Pan-Allele Prediction of MHC Class I-Presented Peptides by Incorporating Antigen Processing*. Cell Systems (2020). [DOI](https://doi.org/10.1016/j.cels.2020.06.010).
- Reynisson et al. *NetMHCpan-4.1 and NetMHCIIpan-4.0: Improved predictions of MHC antigen presentation by concurrent motif deconvolution and integration of MS eluted ligand data*. Nucleic Acids Research (2020). [DOI](https://doi.org/10.1093/nar/gkaa379).
- Cerami et al. *The cBio Cancer Genomics Portal: An Open Platform for Exploring Multidimensional Cancer Genomics Data*. Cancer Discovery (2012). [DOI](https://doi.org/10.1158/2159-8290.CD-12-0095).
- The UniProt Consortium. *UniProt: the Universal Protein Knowledgebase in 2023*. Nucleic Acids Research (2023). [DOI](https://doi.org/10.1093/nar/gkac1052).

These references describe data resources and external predictors. They are not evidence that NeoScorer produces clinically actionable vaccine candidates.
