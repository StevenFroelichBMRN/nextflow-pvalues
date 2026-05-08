FROM python:3.11-slim

RUN pip install --no-cache-dir \
    anndata \
    h5py \
    scipy \
    statsmodels \
    numpy \
    pandas \
    pyarrow

COPY bin/extract_pvalues.py /usr/local/bin/extract_pvalues.py
RUN chmod +x /usr/local/bin/extract_pvalues.py
