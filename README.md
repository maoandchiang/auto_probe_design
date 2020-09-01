_Updated: September. 01, 2020_
## Automatic Probe Design

In capture-based short read sequencing, the sequencing probes are designed from target sequences.
Without loss of generality, we take T-cell receptor (TCR) V genes probe design as an example.
From IMGT database, we downloaded all the TCRV alleles, which mostly range from 192 to 352 bp (with an outlier pseudo gene TRAV8-5\*01 length 1355 bp).
Hence, we set three 60 pb probes for each TCRV allele.
Ideally, three probes are situated in the start, middle, and end of the target allele sequence.
However, we want to avoid the part with high GC content, which can decrease the capture ability of our probes.

```
python3 probe_design.py -fa TCRV_alleles.fasta -pn 3 -pl 60 -fo TCRV_probes.csv -mr 20 -gct 60
``` 

The program ```probe_design.py``` read target sequences' fasta file ```TCRV_alleles.fasta``` and generate designed probes file ```TCRV_probes.csv```.
- pn is the probe number
- pl is the probe length
- mr is the search range for lowest GC content probe
- gct is the GC content threshold, the program will report probes with GC content above the threshold
