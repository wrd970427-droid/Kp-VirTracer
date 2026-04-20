This folder stores reference assets only.

Runtime input for Kp-VirTracer is NOT a FASTA file.
You must provide a BLAST DB prefix created by makeblastdb -out.

Required preparation:
  makeblastdb -in /path/to/virulence_genes.fasta -dbtype nucl -out /path/to/db/virulence_db

Then run Kp-VirTracer with:
  --virulence-db /path/to/db/virulence_db

Expected virulence genes usually include:
  iroB, iucA, peg344, rmpA, rmpA2
