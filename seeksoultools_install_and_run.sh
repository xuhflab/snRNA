# seeksoultools download and install
## seeksoultools software are avaliable on https://github.com/seekgenebio/seeksoultools
git clone git@github.com:seekgenebio/seeksoultools.git
cd seeksoultools
conda env create -n seeksoultools -f conda_dependencies.yml
conda activate seeksoultools
pip install .


# run seeksoultools
## tutorials of seeksoultools are avaliable on http://seeksoul.seekgene.com/en/v1.2.0/2.tutorial/1.rna/2.run.html
seeksoultools rna run \
--fq1 /path/to/demo_dd/demo_dd_S39_L001_R1_001.fastq.gz \
--fq2 /path/to/demo_dd/demo_dd_S39_L001_R2_001.fastq.gz \
--samplename demo_dd \
--genomeDir /path/to/silkbase/star \
--gtf /path/to/silkbase/genes/genes.gtf \
--chemistry DDV2 \
--core 4 \
--include-introns


