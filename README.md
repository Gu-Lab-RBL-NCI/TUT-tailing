# TUT-tailing

## Classification of templated and non-templated isomiRs
### Example 1: Templated isomiR (miR-148b-3p)
```
     TCAGTGCATCACAGAACTTTGTc      read obtained
     |||||||||||||||||||||||
GAAAGTCAGTGCATCACAGAACTTTGTCTCGA  genomic reference
 ```

### Example 2: Untemplated  isomiR (miR-148b-3p)
 ```
     TCAGTGCATCACAGAACTTTGTt      read obtained
     ||||||||||||||||||||||0
GAAAGTCAGTGCATCACAGAACTTTGTCTCGA  genomic reference
 ```

### Example 3: Templated  isomiR with multiple paralogs (miR-1-3p) 
```
    TGGAATGTAAAGAAGTATGTATt      read obtained
    |||||||||||||||||||||||
GCTATGGAATGTAAAGAAGTATGTATCTCA  genomic reference paralog 1
GCTATGGAATGTAAAGAAGTATGTATTTTT  genomic reference paralog 2
```
In order to avoid confounding effects, if a read could be templated at least for one of the paralogs generating that miRNA, it was considered as a templated read.
