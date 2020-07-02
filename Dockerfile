FROM ssidk/bifrost-base:2.0.5

ARG version="v1.0.3"
ARG last_updated="02/07/2020"
ARG name="qcquickie"
ARG full_name="bifrost-${name}"

LABEL \
    name=${name} \
    description="Docker environment for ${full_name}" \
    version=${version} \
    resource_version=${last_updated} \
    maintainer="stca@ssi.dk;"

#- Tools to install:start---------------------------------------------------------------------------
RUN \
    conda install -yq -c conda-forge -c bioconda -c defaults bbmap==38.58; \
    conda install -yq -c conda-forge -c bioconda -c defaults samtools==1.9; \
    conda install -yq -c conda-forge -c bioconda -c defaults cyvcf2==0.11.4; \
    conda install -yq -c conda-forge -c bioconda -c defaults fastqc; \
    # Don't use conda for Quast they cap the python version which causes issues with install
    # Note prokka has a 1 year deadline due to tbl2asn. 1.14.6 was made available Feb 20th
#- Tools to install:end ----------------------------------------------------------------------------

#- Additional resources (files/DBs): start ---------------------------------------------------------
RUN cd /bifrost_resources && \
    wget -q https://raw.githubusercontent.com/ssi-dk/bifrost/master/setup/adapters.fasta && \
    chmod +r adapters.fasta;
#- Additional resources (files/DBs): end -----------------------------------------------------------

#- Source code:start -------------------------------------------------------------------------------
RUN cd /bifrost && \
    git clone --branch ${version} https://github.com/stefanocardinale/${full_name}.git ${name};
#- Source code:end ---------------------------------------------------------------------------------

#- Set up entry point:start ------------------------------------------------------------------------
ENV PATH /bifrost/${name}/:$PATH
ENTRYPOINT ["launcher.py"]
CMD ["launcher.py", "--help"]
#- Set up entry point:end --------------------------------------------------------------------------