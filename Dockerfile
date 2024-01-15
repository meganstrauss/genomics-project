FROM python:3.12-slim

WORKDIR /app

# install the shared library first so its layer caches independently
COPY packages/biocore /app/packages/biocore
RUN pip install --no-cache-dir -e /app/packages/biocore

COPY requirements.txt /app/requirements.txt
RUN pip install --no-cache-dir -r /app/requirements.txt

COPY . /app

ENV PYTHONUNBUFFERED=1
EXPOSE 8000 8001

# default: serve the phylogeny app. Override the command to serve AMR.
CMD ["python", "cli.py", "serve", "phylo", "--host", "0.0.0.0", "--port", "8000"]
