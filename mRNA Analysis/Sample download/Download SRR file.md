
# Download SRR file 
*The following step was performed using a laptop. The files could preferably be downloaded directly to an external disk or USB drive, since SRA files can be large.* 

## Requirements
SRA Toolkit must be installed and available in the terminal:
```
prefetch --version
```

## 1. Go to the dataset directory
Create and/or move to a directory dedicated to the dataset:
```
cd /your_directory/sample_download/your_dataset
```
## 2. Create a sample sheet
Create a tab-separated samplesheet:
```
nano samples.tsv
```
Add the SRAIDs and the sample names. The file must be tab-separated:
```
SRA_ID SAMPLE_NAME
SRRxxxxxxxxx	xxx_Control
SRRxxxxxxxxx	xxx_OE
SRRxxxxxxxxx	xxx_OE
SRRxxxxxxxxx	xxx_OE
SRRxxxxxxxxx	xxx_OE
............ ......
```
Save the file with ```Ctrl + O```, press ```Enter```, and exit with ```Ctrl+X```.

## 3. Create a download script
Create a script for downloading all samples listed in the sample sheet:
```
nano prefetch_only.sh
```

Paste the following code:
```
#!/bin/bash
set -euo pipefail

PROJECT="/your_directory/sample_download/your_dataset"
SAMPLESHEET="$PROJECT/samples.tsv"

mkdir -p "$PROJECT/sra"
cd "$PROJECT/sra"

# Skip header, read only SRA_ID
tail -n +2 "$SAMPLESHEET" | while IFS=$'\t' read -r SRA_ID SAMPLE_NAME
do
    echo "================================"
    echo "Downloading $SRA_ID ($SAMPLE_NAME)"
    echo "================================"

    prefetch "$SRA_ID"

    echo "Finished: $SRA_ID"
done

echo "All SRR files are downloaded."
```
Save the file with ```Ctrl + O```, press ```Enter```, and exit with ```Ctrl+X```.

## 4. Make the script runnable
```
chmod +x prefetch_only.sh
```
## 5. Run the download
```
./prefetch_only.sh | tee prefetch.log
```
Downloaded ```.sra``` files will be saved in :
/your_directory/sample_download/your_dataset/sra
