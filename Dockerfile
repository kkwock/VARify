RUN dnf update -y && dnf install -y python38

# -- Install sam2pairwise -- #
RUN wget https://github.com/mlafave/sam2pairwise/archive/refs/tags/v1.0.0.tar.gz -O - | tar -xz

WORKDIR /app
COPY requirements.txt .

RUN pip3.8 install -r requirements.txt

# -- Install dependencies -- #
COPY extract_codons.py .
COPY configs /app/configs
COPY scripts /app/scripts
COPY test /app/test

