# Manual validation of hotspots and hot regions
*This step describes how to manually validate hotspots and hot regions identified in the previous analysis.*
*Manual validation was required because sample annotation in COSMIC was not always consistent. In some cases, samples from the same individual had different sample names. This could make one patient appear as several unique samples in the downloaded data.*

## 1. Go back to COSMIC
Open COSMIC and go to the gene and tissue or histology category where a hotspot or hot region was identified. 

## 2. Search for the identified site
Use the identified cDNA position or region from the hotspot or hot region output file. 
Enter it in the COSMIC search bar, as seen below. 
This will filter the mutation table to samples with mutations at the specific site. 
<img width="1056" height="528" alt="image" src="https://github.com/user-attachments/assets/a6efdc63-e12b-425f-9a8a-fd6e5c6ae5a0" />

## 3. Check linked samples from the same individual
Click on a sample name, as seen above.
Scroll down to **Additional data**, check the section **Other samples linked to the same individual**. 
This section shows whether other samples are linked to the same patient. 
<img width="701" height="347" alt="image" src="https://github.com/user-attachments/assets/35771b88-458d-4f43-9bda-25523302e4a2" />

## 4. Interpret the validation
If several samples contributing to the hotspot or hotregion are linked to the same individual, the finding should not be considered independent across patients.
In the example above, four samples were identified at the same site. However, the samples were linked to the same individual despite having different sample names. 
This means that the identified hotspot was a false positive caused by inconsistent sample annotation. 
If some samples are linked to the same individual, but the required number of unique individuals remains above the limit after manual validation, the finding can still be considered a hotspot or hot region. However, the correct number of unique individuals should be noted. 

## 5. Use of the previous analysis
The hotspot and hot region analysis are still useful for identifying candidate target sites efficiently. However, the candidate site must be manually validated in COSMIC before being used for target selection or downstream analysis.
