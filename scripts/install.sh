#!/usr/bin/env bash

# Install Python environment
if [ -d "venv" ]; then
  echo "Python virtual venv already exists"
else
  echo "Create new Python venv"
  python3 -m venv venv
  source venv/bin/activate
  echo "Using python from $(which python3), $(python3 -V)"

  pip install -r requirements.txt --ignore-installed
fi

# Install R environment
if [ -f "renv.lock" ]; then
  echo "File renv.lock exits"
else
  Rscript -e '
    if (!"renv" %in% rownames(installed.packages())) {
      install.packages("renv", lib = Sys.getenv("R_LIBS_USER"), repos = "https://cloud.r-project.org")
    } else {
      packageVersion("renv")
    }
    renv::init(bioconductor = T)
  '

  Rscript -e '
    renv::install(c(
      "Matrix@1.6-5",
      "Seurat@4.3.0",
      "SeuratObject@4.1.3",
      "getopt@1.20.4",
      "limma@3.58.1",
      "xgboost@1.7.7.1",
      "cmapR@1.14.0",
      "immunogenomics/presto",
      "jinworks/CellChat"
    ))
  '
fi

# Download drug reference data
if [ -d data/reference ]; then
  echo "Reference data downloaded"
else
  mkdir -p data/reference
  mkdir -p data/reference/drug_validation
  mkdir -p data/reference/l1000
  mkdir -p data/reference/tahoe
  wget -O data/reference/gdsc.csv 'https://www.dropbox.com/scl/fi/quaaz4afke3f053nz6wse/gdsc.csv?rlkey=dzeavs1d4ejglh72wjuzp40gz'
  wget -O data/reference/sider.csv 'https://www.dropbox.com/scl/fi/lx76j27uxah10e9kr1pag/sider.csv?rlkey=aoertif7c48zq9vj1vm8r1ud4'
  wget -O data/reference/drug_validation/liver.csv 'https://www.dropbox.com/scl/fi/ow0gk1aydpjuuk3662li1/liver.csv?rlkey=ywgd0mukuy6u91sro025q9ct6'
  wget -O data/reference/drug_validation/prostate.csv 'https://www.dropbox.com/scl/fi/ljhth4hh70xru1bnnip2t/prostate.csv?rlkey=kq7x9t975g4t0xa9uek0qb6id'
  wget -O data/reference/l1000/GSE70138.tar.gz 'https://www.dropbox.com/scl/fi/vn6lhq00yku71lv1725ib/GSE70138.tar.gz?rlkey=cdhbpziwcgc3xwafj7mf3vlc7'
  wget -O data/reference/l1000/GSE92742.tar.gz 'https://www.dropbox.com/scl/fi/qv51o6ldwwokyxn5tacqp/GSE92742.tar.gz?rlkey=8q5wom3h8rv6ilca0t3qby1wt'
  wget -O data/reference/tahoe/drug_ref.rds 'https://www.dropbox.com/scl/fi/e7pv5jyxs1nqzgqwndeex/drug_ref.rds?rlkey=5g95qul6lu83b78gbph6w8erj'

  # De-compress tar gz
  tar -xzvf data/reference/l1000/GSE70138.tar.gz -C data/reference/l1000
  tar -xzvf data/reference/l1000/GSE92742.tar.gz -C data/reference/l1000
  rm data/reference/l1000/GSE70138.tar.gz
  rm data/reference/l1000/GSE92742.tar.gz
fi
