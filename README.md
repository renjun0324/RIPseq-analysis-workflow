# RIPseq-analysis-workflow
Summarized existing analysis methods for RIP-seq, including alignment, peak calling, and annotation.

<a id='section1'></a>
## analysis-workflow
1. fastqc
2. create index & sequence
3. duplicated & make index
4. peak calling
5. annotation

<a id='section1'></a>
## peak calling methods

### 1. CLAM
https://github.com/Xinglab/CLAM

#### installation
```
conda activate python2.7
pip install --index-url https://test.pypi.org/simple/ --no-deps CLAM
```

### 2. Piranha
https://github.com/smithlabcode/piranha

#### installation
```
conda activate python2.7
conda install -c bioconda piranha
```
### 3. RIPSeeker
https://bioconductor.riken.jp/packages/3.4/bioc/html/RIPSeeker.html
https://github.com/yueli-compbio/RIPSeeker
