FROM mambaorg/micromamba:1.5.10

COPY --chown=$MAMBA_USER:$MAMBA_USER environment.yml /tmp/environment.yml

RUN micromamba install -y -n base -f /tmp/environment.yml && \
    micromamba clean --all --yes

WORKDIR /work

ENTRYPOINT ["/usr/local/bin/_entrypoint.sh"]
CMD ["snakemake", "--help"]
