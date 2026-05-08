FROM python:3.11-slim

RUN pip install --no-cache-dir \
    anndata==0.10.8 \
    h5py \
    scipy \
    statsmodels \
    numpy \
    pandas==2.1.4 \
    pyarrow
