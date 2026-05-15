#!/bin/bash
# ---------------------------------------------------------
# Test Dataset Fetcher for EasyRNA-Seq Benchmark
# Target: GSE52778 (Airway smooth muscle cells) - Standard DESeq2 test dataset
# ---------------------------------------------------------

mkdir -p test_dataset/fastq
cd test_dataset/fastq

echo "🔥 Starting Download of Test FASTQ files (GSE52778)..."
echo "We will fetch the first 1,000,000 reads to ensure test runs finish quickly."

# SRA Tools Docker Image
IMG="quay.io/biocontainers/sra-tools:3.0.3--h87f3376_0"

# Sample 1: SRR1039508 (Dexamethasone treated, rep 1)
echo "📥 Downloading SRR1039508..."
docker run --rm -v $(pwd):/data -w /data $IMG \
    fasterq-dump -X 1000000 --split-files -p SRR1039508

# Sample 2: SRR1039509 (Untreated, rep 1)
echo "📥 Downloading SRR1039509..."
docker run --rm -v $(pwd):/data -w /data $IMG \
    fasterq-dump -X 1000000 --split-files -p SRR1039509

echo "📦 Compressing FASTQ files..."
gzip -f *.fastq

echo "✅ Test dataset download complete! Files are saved in test_dataset/fastq/"
