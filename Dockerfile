FROM python:3.11-slim

RUN apt-get update && apt-get install -y --no-install-recommends procps && rm -rf /var/lib/apt/lists/*

RUN pip install --no-cache-dir \
    anndata==0.10.8 \
    h5py \
    scipy \
    statsmodels \
    numpy \
    pandas==2.1.4 \
    pyarrow

COPY bin/extract_pvalues.py /usr/local/bin/extract_pvalues.py
RUN chmod +x /usr/local/bin/extract_pvalues.py
