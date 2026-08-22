FROM python:3.11-slim

LABEL maintainer="Fabiano Menegidio <labiomics@bioinformatica.com.br>"
LABEL description="Container image for ampliQC: Context-aware Quality Control engine for Amplicon sequencing data (Short & Long Reads)"

WORKDIR /app

RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    zlib1g-dev \
    git \
    && rm -rf /var/lib/apt/lists/*

COPY . /app

RUN pip install --no-cache-dir .

ENTRYPOINT ["ampliqc"]
CMD ["--help"]
