# Design target assays
**

## 1. Retrieve genomic coordinates from COSMIC
Go to COSMIC and open the selected gene and tissue as mentioned in section 1. 
Copy the **Genomic Coordinates GRCh38** for the target site or region. 
<img width="495" height="236" alt="image" src="https://github.com/user-attachments/assets/c0ba2e5f-2d43-4b97-8808-55cc4f9193f8" />

## 2. Open the region in UCSC Genome Browser
Go to the UCSC Genome Browser using the **Human GRCh38/hg38** assembly. 
Enter the genomic coordinate in the search bar and click **Search**.
<img width="1211" height="131" alt="image" src="https://github.com/user-attachments/assets/d1fc838f-ed75-42b0-a515-26ac3d86f961" />

## 3. Get the DNA sequence
Hover over **View** found in the menu, and select **DNA sequence**.
<img width="1203" height="192" alt="image" src="https://github.com/user-attachments/assets/78391cf5-9483-4472-bf3c-01427ae83500" />

Add **105 bases** upstream and downstream of the target site. 
Click **Get DNA**.
<img width="1375" height="467" alt="image" src="https://github.com/user-attachments/assets/11945372-1275-4036-91b3-ab5b41d3dbdf" />

Copy the DNA sequence. 
<img width="795" height="159" alt="image" src="https://github.com/user-attachments/assets/415b41eb-e237-4f2c-8153-1dbe8b32ed15" />

## 4. 
Go to:https://www.ncbi.nlm.nih.gov/tools/primer-blast/
Paste the copied DNA sequence.

Set the range according to the 
Needs to take avalible regions into consideration, target site is located at position 106-107. With a safety margin of at least 5 bases on each side of the target, the target interval including safety margins was 100-113. Forward <100 Reverse >113

Then the identified highly similar regions needs to be taken into consideration, as in previus example in part 7.
```
1-86
94-99
127-129
152-168
```
Apply the appropriate ranges: 
<img width="1119" height="253" alt="image" src="https://github.com/user-attachments/assets/04768f38-60c6-4c62-83a9-88fc9509021b" />


Apply the following settings:
PCR product size: Min= 60, max to 80
<img width="710" height="268" alt="image" src="https://github.com/user-attachments/assets/6c6c767f-727f-41a6-841e-3e9614fe8ad4" />

Database: Genomes for selected eukaryotic organisms
<img width="1164" height="400" alt="image" src="https://github.com/user-attachments/assets/8fd2e86c-f662-4c75-8f4e-11d5bd6b9168" />


Scroll down and open the advanced parameters:
<img width="313" height="43" alt="image" src="https://github.com/user-attachments/assets/fd7ee89f-b547-4019-9353-4781018e53c1" />

Set Max Self complimentary: Any = 5, 3´=3
Max pair complamantary:  Any = 5, 3´=3
<img width="772" height="108" alt="image" src="https://github.com/user-attachments/assets/1b0fbaac-e67a-465c-be56-5922297a71c3" />


Click "Get parameters"
<img width="605" height="61" alt="image" src="https://github.com/user-attachments/assets/73070795-6832-49ac-b193-a4e99ca6eb5e" />


Click Check
<img width="605" height="61" alt="image" src="https://github.com/user-attachments/assets/e1ce37e5-a7a9-4027-9711-a76baa2b531d" />


check the box and click submit

click submit
<img width="1490" height="305" alt="image" src="https://github.com/user-attachments/assets/5ec2b887-eee8-42c2-a902-196c5aed9e28" />


Hower over the diffrent results and denote the range where the forvard and reverse primers are within. 
<img width="610" height="63" alt="image" src="https://github.com/user-attachments/assets/f6d574d4-d3e8-4d0e-9673-30ac26a3cac1" />

<img width="1465" height="222" alt="image" src="https://github.com/user-attachments/assets/bf2a4c8b-4a2b-4682-9258-df865cd4d095" />

Forward är 79-96
94, 95, 96 Last 3 is uniqe (3´anchering) also the start 79-86 unuiq 

126-145
127-129 är unik, 126 är inte de så inte 3´anchoring

<img width="1415" height="139" alt="image" src="https://github.com/user-attachments/assets/258e4f7e-7909-4b77-8198-a13858a11831" />

Condition is 3´anchering where at lest the 3 last bases are uniqe.

The primers were designed for multiplex PCR, where several targets are amplified at the same time. 

This requiers that the primers need similar conditions, especially GC content, since it affects the melting temperature. The amplicons were kept short to make them suitable for fragmented DNA. 

Prioritize primers that have low self 3´complimentrary , then also self complimentary. 

Select 3 primers per site that shows the most promesing conditions. 
