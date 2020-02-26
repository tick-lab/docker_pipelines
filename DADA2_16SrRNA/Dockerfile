FROM bioconductor/release_core2
MAINTAINER Allison E. Mann (allison.e.mann@gmail.com)

RUN mkdir -p /01_data
RUN mkdir -p /02_code
RUN mkdir -p /03_output

COPY /02_code/install_packages.R /02_code/install_packages.R
COPY /02_code/script.R /02_code/script.R

# Install background packages
RUN sudo apt-get update \
	&& sudo apt-get install \
	&& sudo apt-get -y install python3-pip

# Install R packages
RUN Rscript /02_code/install_packages.R

# Install system packages
# Install and update conda
RUN wget -q https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -O /tmp/miniconda3.sh \
	&& /bin/bash /tmp/miniconda3.sh -bp /opt/miniconda3 \
	&& rm /tmp/miniconda3.sh \
	&& /opt/miniconda3/bin/conda update -yq conda \
	&& /opt/miniconda3/bin/conda install -yq conda-build

# Install python packages
RUN pip3 install --user --upgrade cutadapt \
	&& export PATH=/root/.local/bin:$PATH

# Install git packages
RUN git clone https://github.com/lh3/seqtk.git \
	&& cd seqtk \
	&& make

# Set conda environment
RUN export PATH=/opt/miniconda3/bin:$PATH
ENV PATH /opt/miniconda3/bin:$PATH

# Install conda packages
RUN /opt/miniconda3/bin/conda install -c bioconda/label/cf201901 vsearch
RUN wget http://data.qiime2.org/distro/core/qiime2-2019.10-py36-linux-conda.yml \
	&& /opt/miniconda3/bin/conda env create -n qiime2-2019.10 --file qiime2-2019.10-py36-linux-conda.yml \
	&& rm qiime2-2019.10-py36-linux-conda.yml \
	&& conda init bash

CMD Rscript /02_code/script.R