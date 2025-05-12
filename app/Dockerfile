FROM --platform=linux/amd64 r-base:4.3.0

ENV TZ=Asia/Tokyo
USER root

# apt install
RUN apt-get update &&\
    apt-get -y install \
	bc \
	build-essential \
    make \
    cmake \
    wget \
    git \
	default-jre \
    unzip \
	sed \
    gawk \
    rename \
    rsync &&\
	apt-get clean &&\
    rm -rf /var/lib/apt/lists/*

# install miniforge3
ENV CONDA_DIR=/opt/conda
RUN wget -O Miniforge3-Linux-x86_64.sh "https://github.com/conda-forge/miniforge/releases/download/24.9.0-0/Miniforge3-Linux-x86_64.sh" && \
    sh Miniforge3-Linux-x86_64.sh -b -p ${CONDA_DIR} && \
    rm -f Miniforge3-Linux-x86_64.sh
# set path
ENV PATH=${CONDA_DIR}/bin:${PATH}


# # ----------    Install dependencies of SNPcaster    ----------
# Install BactSNP ver. 1.1.0
WORKDIR /usr/local/src

# install as website instructions
# write the location of BactSNP in line 58 of bactsnp_snippy_gubbins.sh
# e.g. DIR1=${HOME}/programs/bactsnp-1.1.0.linux64
RUN wget https://github.com/IEkAdN/BactSNP/releases/download/v1.1.0/bactsnp-1.1.0.linux64.tgz && \
    tar xf bactsnp-1.1.0.linux64.tgz && \
    rm bactsnp-1.1.0.linux64.tgz
RUN cd bactsnp-1.1.0.linux64 && \
    make && \
    cp bactsnp /usr/local/bin

# # ----------    Install BLAST    ----------
RUN wget https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.15.0/ncbi-blast-2.15.0+-x64-linux.tar.gz &&\
    tar -xzf ncbi-blast-2.15.0+-x64-linux.tar.gz &&\
    ln -s /usr/local/src/ncbi-blast-2.15.0+/bin/* /usr/local/bin/

# # ----------    Install libraries from conda    ----------
COPY setup/conda/ /usr/local/src/snpcaster/setup/conda/
WORKDIR /usr/local/src/snpcaster/setup/conda
RUN mamba env update -f=conda_env_base.yml
RUN chmod +x /usr/local/src/snpcaster/setup/conda/conda_env_create.sh
RUN ./conda_env_create.sh conda_env
    
# # ----------    Set an entrypoint    ----------
RUN apt-get update && apt-get -y install gosu
COPY setup/entrypoint.sh /usr/local/bin/entrypoint.sh
RUN chmod +x /usr/local/bin/entrypoint.sh
ENTRYPOINT ["/usr/local/bin/entrypoint.sh"]
                
# # ----------    Install grape & snpcaster   ----------
ENV BAC1=/usr/local/src/snpcaster
COPY src/ /usr/local/src/snpcaster/
RUN ln -s /usr/local/src/snpcaster/grape_qc_assembly/find_strain_pairs.py /usr/local/bin/find_strain_pairs.py && \
    ln -s /usr/local/src/snpcaster/grape_qc_assembly/grape_qc_assembly.sh /usr/local/bin/grape_qc_assembly.sh && \
    ln -s /usr/local/src/snpcaster/grape_qc_assembly/grape_qc_assembly_reduced_tree.sh /usr/local/bin/grape_qc_assembly_reduced_tree.sh && \
    ln -s /usr/local/src/snpcaster/snpcaster/snpcaster.sh /usr/local/bin/snpcaster.sh && \
    ln -s /usr/local/src/snpcaster/snpcaster/snpclipper.sh /usr/local/bin/snpclipper.sh && \
    ln -s /usr/local/src/snpcaster/snpcaster/snpcaster_gubbins.sh /usr/local/bin/snpcaster_gubbins.sh && \
    ln -s /usr/local/src/snpcaster/snpcaster/raxml-ng.sh /usr/local/bin/raxml-ng.sh && \
    ln -s /usr/local/src/snpcaster/snpcaster/art_fq_generator.sh /usr/local/bin/art_fq_generator.sh && \
    ln -s /usr/local/src/snpcaster/snpcaster/n_combine.sh /usr/local/bin/n_combine.sh && \
    ln -s /usr/local/src/snpcaster/snpcaster/repeat_detect.sh /usr/local/bin/repeat_detect.sh && \
    ln -s /usr/local/src/snpcaster/snpcaster/iqtree.sh /usr/local/bin/iqtree.sh

# # ----------    Create notebook folder   ----------
COPY notebook/ /home/jovyan/snpcaster/notebook/

# # ----------    Grant execution permission   ----------
RUN find /usr/local/src -type f -exec chmod +x {} \;
RUN find /usr/local/bin -type f -exec chmod +x {} \;

# # ----------    Create mounted folder   ----------
RUN mkdir /home/jovyan/snpcaster/notebook/project
WORKDIR /home/jovyan/snpcaster/notebook/project
RUN useradd -m -u 1001 snpcaster
RUN chown -R snpcaster:snpcaster /home/jovyan/snpcaster/notebook/project
