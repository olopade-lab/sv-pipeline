# Installing Conda: https://developers.google.com/earth-engine/python_install-conda#mac

curl https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-x86_64.sh -o ~/miniconda.sh
bash ~/miniconda.sh -b -p
rm ~/miniconda.sh
source $HOME/miniconda3/bin/activate
printf '\n# add path to conda\nexport PATH="$HOME/miniconda3/bin:$PATH"\n' >> ~/.bashrc
export PATH="$HOME/miniconda3/bin:$PATH"
conda activate

# Installing Delly: https://anaconda.org/bioconda/delly

conda install -c bioconda delly
brew install dos2unix samtools bcftools
pip install parsl

# Installing Survivor

wget https://github.com/fritzsedlazeck/SURVIVOR/archive/master.tar.gz -O SURVIVOR.tar.gz
tar xzvf SURVIVOR.tar.gz
cd SURVIVOR-master/Debug/
make
./SURVIVOR
