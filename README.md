# README

- A collection of useful scripts

## Join and split table / sequences / other entities

- assign multiple files into multiple chunks: `scripts/chunkify-fasta.py`  
- combined multple fasta file in a single one: `scripts/combine-fasta.py`  
- chunkify a fasta file: `scripts/split-fasta.py`  
- chunkify a text file: `scripts/split-text.py`
- split stockholm format alignment: `scripts/split-stockholm.py`
- split covarience models: `scripts/split-cm-models.py`
- combined table: `scripts/concatenate-table.py`
- group table by customized group id extractor
  - `scripts/group-text.py`
  - Line correspond to the same group should be consecutive
  - `scripts/group-text.py -i pairwise.0309.txt -e 'lambda x:x.split(":GCF")[0]' -r 'lambda x,n:x.split(":")[0] +"-" + str(n).zfill(5)' -wh -od pairwise.0309`

## Tabular data reformatting
- reformat fragGeneScan results: `scripts/fgs2bed.py` 
- nhmmer tabuler to gff format: `scripts/nhmmer-tbl-to-gff.py`
- infernal tabuler to gff format: `scripts/infernal-tbl-to-gff.py`
- :TODO: protein hmmsearch
- convert gff format to bed format: `scripts/gff2bed.py`


## fasta processing
- subsetting fasta by sequence id: `scripts/subsetting-fasta.py`
  - similar to `seqtk subseq seq-ids.txt`
- remove sequence very similar to that present in another set:`scripts/extract-non-representative-segments.py`
- remove nearly identical sequence by minimap2 search: `scripts/reduce-redundancy.py`
- drop sequence with the same id: `scripts/drop-sequence-with-same-id.py`
- rename sequence in input fasta file `scripts/rename-fasta.py`

## Homolog search
- Homolog search with mmseqs
  - RNA search: `scripts/RNA-homolog-search.py`
  - protein search: `scripts/protein-homolog-search.py`

- filter mmseqs hits in blast+ format
  - `scripts/filter-hits.py`

## Sequence cluster processing

- clustering pairwise search results with MCL algorithm
  - `scripts/mcl-clustering.py`
 
- clustering pairwise search results with leiden algorithm
  - `scripts/leiden-partitioning.py`
 
- reformat cd-hit/cd-hit-est clustering results
  - `scripts/cd-hit-to-clustering-table.py`
 
- reformat MCL clustering results
  - `scripts/mcl-to-clustering-table.py`

- group sequence by clustering table
  - `scripts/group-sequences.py`


## Interval manipulation

- Annotate genomic context of specified interval in a prokaryote genome given gene bed file
  - `scripts/annotate-intervals.py`

- Pick local max of overlapping intervals
  - `scripts/pick-local-max.py`

- Select mutual closest features in 2 bed files, now used for transcription unit determination
  - `scripts/extract-transcription-unit-from-genome.py` 
 

## k-mer related analysis

- count kmer frequency: `scripts/kmer-frequency-fitter.py`
- simulate sequence from kmer frequency: `scripts/kmer-emitter.py`

- classify sequences with kmer composition: 
  - `kmer-profile-classification.py`
  - `scripts/kmer-profile-inference.py`
 
## Sequence shuffling
- shuffle sequence while preserve kmer frequency: `scripts/kmer-preserved-shuffling.py`
- shuffled homolog sequence while preserve phylogenetic signal: `scripts/phylogeny-preserved-shuffling.py`


