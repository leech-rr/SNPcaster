FROM --platform=linux/amd64 r-base:4.3.0

ENV TZ=Asia/Tokyo
USER root

# ARG & ENV
# Default values for ARGs are set in docker-compose.yml using variable substitution.
# Default values are also set here for building without Docker Compose.
ARG HTTP_PROXY=""
ARG HTTPS_PROXY=""
ARG APT_PIPELINING=""
ARG APT_NO_CACHE="1"
ARG APT_HTTP_SETTING_FILE="/etc/apt/apt.conf.d/http_setting"
ENV HTTP_PROXY=$HTTP_PROXY
ENV http_proxy=$HTTP_PROXY
ENV HTTPS_PROXY=$HTTPS_PROXY
ENV https_proxy=$HTTPS_PROXY

# apt setting
RUN if [ -n "${HTTP_PROXY}" ]; then \
        echo "Acquire::http::Proxy \""${HTTP_PROXY}"\";" >> "${APT_HTTP_SETTING_FILE}"; \
    fi; \
    if [ -n ""${HTTPS_PROXY}"" ]; then \
        echo "Acquire::https::Proxy \""${HTTPS_PROXY}"\";" >> "${APT_HTTP_SETTING_FILE}"; \
    fi;
RUN if [ -n "${APT_PIPELINING}" ]; then \
        if [ "${APT_PIPELINING}" -ge 0 ]; then \
            echo "Acquire::http::Pipeline-Depth ${APT_PIPELINING};" >> "${APT_HTTP_SETTING_FILE}"; \
        else \
            echo "[Error] APT_PIPELINING [${APT_PIPELINING}] must be a non-negative integer. Exit..." >&2; \
            exit 1; \
        fi; \
    else \
        if [ -n "${HTTP_PROXY}" ] || [ -n "${HTTPS_PROXY}" ]; then \
            echo "Acquire::http::Pipeline-Depth 0;" >> "${APT_HTTP_SETTING_FILE}"; \
        fi; \
    fi;

RUN if [ "${APT_NO_CACHE}" -eq 0 ]; then \
            echo "Acquire::http::No-Cache false;" >> "${APT_HTTP_SETTING_FILE}"; \
    elif [ "${APT_NO_CACHE}" -eq 1 ]; then \
        echo "Acquire::http::No-Cache true;" >> "${APT_HTTP_SETTING_FILE}"; \
    else \
        echo "[Error] APT_NO_CACHE [${APT_NO_CACHE}] must be 0 or 1. Exit..." >&2; \
        exit 1; \
    fi;

RUN cat "${APT_HTTP_SETTING_FILE}"
RUN exit 1

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

# # ----------    Grant execution permission   ----------
RUN find /usr/local/src -type f -exec chmod +x {} \;
RUN find /usr/local/bin -type f -exec chmod +x {} \;
    
# # ----------    Create SNPcaster user   ----------
RUN useradd -m -u 1001 snpcaster

# # ----------    Create notebook folder   ----------
COPY notebook/ /home/snpcaster/notebook/

# # ----------    Create mounted folder   ----------
RUN mkdir /home/snpcaster/notebook/project
RUN chown -R snpcaster:snpcaster /home/snpcaster/notebook/project
WORKDIR /home/snpcaster/notebook/project
