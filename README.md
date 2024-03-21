# README

- A collection of useful scripts

## Join and split table / sequences / other entities

- assign multiple files into multiple chunks: scripts/chunkify-fasta.py  
- combined multple fasta file in a single one: scripts/combine-fasta.py  
- chunkify a fasta file: scripts/split-fasta.py  
- chunkify a text file: scripts/split-text.py
- split stockholm format alignment: scripts/split-stockholm.py
- split covarience models: scripts/split-cm-models.py 
- combined table: scripts/concatenate-table.py


## Tabular data reformatting
- reformat fraggenescan results: scripts/fgs2bed.py 
- nhmmer tabuler to gff format: scripts/nhmmer-tbl-to-gff.py
- infernal tabuler to gff format: scripts/infernal-tbl-to-gff.py
- :TODO: protein hmmsearch
- convert gff format to bed format: scripts/gff2bed.py


## fasta processing
- subsetting fasta by sequence id: scripts/subsetting-fasta.py
  - similar to `seqtk subseq seq-ids.txt`
- remove nearly identical sequence by minimap2 search: scripts/reduce-redundancy.py
- remove sequence very similar to that present in another set:scripts/extract-non-representative-segments.py
- drop sequence with the same id: scripts/drop-sequence-with-same-id.py


- group-sequences.py  
- cmfinder-search.py  

