# Tutorial output provenance

These figures and tables were generated from the public, non-human mothur
MiSeq SOP mouse-gut and mock-community teaching dataset:

- Description: https://mothur.org/wiki/miseq_sop/
- Public archive: https://mothur.s3.us-east-2.amazonaws.com/wiki/miseqsopdata.zip
- SHA-256: `728f7544824f8ca0223618b18bb60fce6a5dc20ae3d476267dbd5cf5594a35d4`

On 15 November 2026 the existing input archive's SHA-256 was checked against a
fresh download from that public URL; the hashes matched exactly. Sample names
and input counts are also published in the mothur SOP. No private study data
were used. Raw FASTQ files are not included in this repository.

Generation: `docs/tutorial/run.jl`, Julia 1.11.1, DenoiseDNA 0.2.0 committed
release source on Athena (batch job 4052432). The full processing run completed
successfully. Taxonomy classification is a separately documented optional step
and was not used to generate these figures.

Please credit the dataset authors: Kozich JJ, Westcott SL, Baxter NT,
Highlander SK, Schloss PD (2013), Development of a dual-index sequencing
strategy and curation pipeline for analyzing amplicon sequence data on the
MiSeq Illumina sequencing platform, Applied and Environmental Microbiology
79(17):5112–5120.
