# CLASSEQ: A clade-informed sequence identification tool

Literature:

CLASSEQ: A clade-informed sequence identification tool. Leveraging phylogenetic insights for biological sequences classification.

## Install

## Quick start

See quick [start](docs/mds/README.md).

## Serve Classeq as API

```bash

export CONTEXT_MODELS="docs/manuscript/assays/bacillus-subtilis-group/supplementary-material/05-config.yaml"

uvicorn classeq.ports.api.main:app \
    --reload \
    --host 0.0.0.0 \
    --port 8080

```
