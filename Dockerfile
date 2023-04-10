FROM docker.artifactory01.ghdna.io/ghbibaseimage:CentOS8.5.2111
# temporary build, replace with final before release

ARG ARTIFACTORY_USERNAME
ARG ARTIFACTORY_PASSWORD

RUN dnf update -y && dnf install -y python38

# -- Install bcl2fastq -- #
RUN curl -SL -u "${ARTIFACTORY_USERNAME}:${ARTIFACTORY_PASSWORD}" \
        https://docker.artifactory01.ghdna.io:443/artifactory/biInternalSoftware/bcl2fastq-2.20.0.422.tar.gz \
        | tar -xzC /usr/bin/

WORKDIR /app
COPY requirements.txt .

RUN pip3.8 install -r requirements.txt --index-url https://${ARTIFACTORY_USERNAME}:${ARTIFACTORY_PASSWORD}@artifactory01.ghdna.io/artifactory/api/pypi/pypi/simple

# -- Install dependencies -- #
COPY g360_eioqc.py .
COPY configs /app/configs
COPY scripts /app/scripts
COPY test /app/test
COPY VERSION.txt ./
