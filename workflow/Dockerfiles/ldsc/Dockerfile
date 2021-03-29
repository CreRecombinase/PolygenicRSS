FROM python:2.7.18-stretch
RUN apt-get update \
   && apt-get install -y --no-install-recommends build-essential gfortran libgfortran3 git curl ca-certificates bcftools libopenblas-base libopenblas-dev
RUN pip install bitarray==0.8 nose==1.3.7 pybedtools==0.7 scipy==0.18 pandas==0.20 numpy==1.16
RUN pip install git+https://github.com/xinhe-lab/ldsc.git && which ldsc.py
