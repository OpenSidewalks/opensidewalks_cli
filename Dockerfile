FROM python:3.8.5-buster
MAINTAINER Nick Bolten <nbolten@gmail.com>

RUN apt-get update && \
    apt-get install -y \
      libspatialindex-dev \
      gdal-bin \
      osmium-tool

RUN pip3 install --upgrade pip

RUN mkdir -p /work
WORKDIR /work

COPY opensidewalks_cli /work/opensidewalks_cli
COPY ./pyproject.toml /work

RUN pip3 install /work
