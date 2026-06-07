# GC content analysis
*This section describes how to retrieve the target DNA sequence using UCSC In Silico PCR and calculate the GC content using the VectorBuilder GC content calculator.*

## 1. Run In Silico PCR
Go to UCSC In-Silico PCR:
https://genome.ucsc.edu/cgi-bin/hgPcr

Paste the forward and reverse sequences into respective fields. 
Click submit. 
<img width="1671" height="218" alt="image" src="https://github.com/user-attachments/assets/946778d8-fcca-4ae5-937a-ece5c48c95cc" />


## 2. Open the PCR product
Click the result link to open the predicted PCR product. 
If a **fix** option appears, do not click it. 
<img width="954" height="167" alt="image" src="https://github.com/user-attachments/assets/f938ab5f-8914-49bf-979f-874ac2d3bce1" />

## 3. Get the DNA sequence
Hover over **View** found in the menu, and select **DNA sequence**.
<img width="1551" height="313" alt="image" src="https://github.com/user-attachments/assets/2f7aea88-f314-4828-a6b1-d6934e501bb2" />


Set both upstream and downstream sequences to **0** bp. 
Click **Get DNA**.
<img width="1575" height="571" alt="image" src="https://github.com/user-attachments/assets/b9d96367-d44a-4097-84a4-6b8377191f8c" />

Copy the DNA sequence. 
<img width="644" height="97" alt="image" src="https://github.com/user-attachments/assets/11a47c41-9dc0-4048-8ab4-6a61d55414ad" />


## 4. Calculate CG content
Go to the VectorBuilder GC content calculator:
https://en.vectorbuilder.com/tool/gc-content-calculator.html

Paste the copied DNA sequence into the sequence box. 
Click **Submit**.
<img width="1443" height="665" alt="image" src="https://github.com/user-attachments/assets/d54bc155-84cc-4883-b2eb-675daa884fca" />


## 5. Output
The calculator will return the GC content of the target sequence. This could help understand some variations in qPCR results.  
<img width="1236" height="446" alt="image" src="https://github.com/user-attachments/assets/eb7f34a8-3c08-4b0a-83cf-56bf7f6ad842" />


