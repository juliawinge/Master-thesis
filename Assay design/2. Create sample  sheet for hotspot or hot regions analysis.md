# Create a sample sheet from the downloaded TSV file
*This step describes how to create a simplified sample sheet from the downloaded COSMIC `.tsv` file. The output file will be used as input for the hotspot and hot region analysis scripts.*

 ## 1. Open the TSV file in Excel
Open a new Excel document/ new page. 
Go to **Data** in the top menu, and click **"From text/CSV"**.
<img width="936" height="305" alt="image" src="https://github.com/user-attachments/assets/01cd2a02-9de5-48cd-ac4e-addb5a8b9d7a" />

## 2. Select the downloaded TSV file
In the file selection window, open the file type menu and select **All files** to make the `.tsv` files visible
Select the TSV file of interest, and click **Import**.
<img width="955" height="597" alt="image" src="https://github.com/user-attachments/assets/6c7ce6ec-73e5-4253-847b-0e982a857e85" />

## 3. Import the file
 Use the import settings shown below, and click **Read file**.
 <img width="917" height="664" alt="image" src="https://github.com/user-attachments/assets/05a07bc6-f587-400e-946c-040e900710c9" />

## 4. Create a simplified sample sheet
Create a new Excel file/ page.
Copy the following columns from the imported COSMIC file:
```text
Sample ID
Sample Name
CDS Mutations
```

Rename the columns to:
```text
sample_id
sample_name
cds_mutations
```
These column names are required as input for the downstream scripts. 
<img width="536" height="168" alt="image" src="https://github.com/user-attachments/assets/9b7b8d58-5a7b-423d-9a1f-622faf7b08ad" />

## 5. Save the sample sheet
Save the file in `.tsv` or `.txt` format in the following folder
```
/your_directory/hot_spot_region_analysis
```
Name the file using the gene name followed by tissue or histology type.
<img width="465" height="91" alt="image" src="https://github.com/user-attachments/assets/a086b365-e734-411f-982a-614bfeb7b865" />
