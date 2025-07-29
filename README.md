# GeneLab Methylation Sequencing Data Processing Workflow

> GeneLab, part of [NASA's Open Science Data Repository (OSDR)](https://www.nasa.gov/osdr), has wrapped each step of the Methylation Sequencing (Methyl-Seq) consensus pipeline ([MethylSeq](https://github.com/nasa/GeneLab_Data_Processing/tree/master/Methyl-Seq)) into a Nextflow workflow with validation and verification of output files built in after each step. This repository contains the Nextflow workflow code (NF_MethylSeq) along with instructions for installation and usage. Exact workflow run info and MethylCP version used to process specific datasets that have been released are available in the \*nextflow\_processing\_info.txt file on the [Open Science Data Repository (OSDR)](https://osdr.nasa.gov/bio/repo/), which can be found under 'Files' -> 'GeneLab Processed Methyl-Seq Files' -> 'Supplemental Materials'.

# NF\_MethylSeq Workflow Information and Usage Instructions

## General workflow info
The current GeneLab Methyl-Seq sequencing data processing pipeline (MethylSeq), [GL-DPPD-7113.md](https://github.com/nasa/GeneLab_Data_Processing/tree/master/Methyl-Seq/Pipeline_GL-DPPD-7113_Versions/GL-DPPD-7113.md), is implemented as a [Nextflow](https://nextflow.io/) DSL2 workflow and can utilize images managed through either [Docker](https://docs.docker.com/get-started/) or [Singularity](https://docs.sylabs.io/guides/3.10/user-guide/introduction.html), or [conda])https://docs.conda.io/en/latest/) for environments (GeneLab uses Singularity when processing GLDS datasets). This workflow (NF_MethylSeq) is run using the command line interface (CLI) of any unix-based system. While knowledge of creating workflows in Nextflow is not required to run the workflow as is, [the Nextflow documentation](https://nextflow.io/docs/latest/index.html) is a useful resource for users who want to modify and/or extend this workflow.   

## Utilizing the workflow

coming soon.....

<br>

---

## Licenses

The software for the Methyl-Seq pipeline and workflow is released under the [NASA Open Source Agreement (NOSA) Version 1.3](License/Methylation_Sequencing_NOSA_License.pdf).


### 3rd Party Software Licenses

Licenses for the 3rd party open source software utilized in the Methyl-Seq pipeline and workflow can be found in the [3rd_Party_Licenses sub-directory](License/3rd_Party_Licenses/). 

<br>

---

## Notices

Copyright © 2024 United States Government as represented by the Administrator of the National Aeronautics and Space Administration.  All Rights Reserved. 

### Disclaimers

No Warranty: THE SUBJECT SOFTWARE IS PROVIDED "AS IS" WITHOUT ANY WARRANTY OF ANY KIND, EITHER EXPRESSED, IMPLIED, OR STATUTORY, INCLUDING, BUT NOT LIMITED TO, ANY WARRANTY THAT THE SUBJECT SOFTWARE WILL CONFORM TO SPECIFICATIONS, ANY IMPLIED WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, OR FREEDOM FROM INFRINGEMENT, ANY WARRANTY THAT THE SUBJECT SOFTWARE WILL BE ERROR FREE, OR ANY WARRANTY THAT DOCUMENTATION, IF PROVIDED, WILL CONFORM TO THE SUBJECT SOFTWARE. THIS AGREEMENT DOES NOT, IN ANY MANNER, CONSTITUTE AN ENDORSEMENT BY GOVERNMENT AGENCY OR ANY PRIOR RECIPIENT OF ANY RESULTS, RESULTING DESIGNS, HARDWARE, SOFTWARE PRODUCTS OR ANY OTHER APPLICATIONS RESULTING FROM USE OF THE SUBJECT SOFTWARE.  FURTHER, GOVERNMENT AGENCY DISCLAIMS ALL WARRANTIES AND LIABILITIES REGARDING THIRD-PARTY SOFTWARE, IF PRESENT IN THE ORIGINAL SOFTWARE, AND DISTRIBUTES IT "AS IS."

Waiver and Indemnity: RECIPIENT AGREES TO WAIVE ANY AND ALL CLAIMS AGAINST THE UNITED STATES GOVERNMENT, ITS CONTRACTORS AND SUBCONTRACTORS, AS WELL AS ANY PRIOR RECIPIENT.  IF RECIPIENT'S USE OF THE SUBJECT SOFTWARE RESULTS IN ANY LIABILITIES, DEMANDS, DAMAGES, EXPENSES OR LOSSES ARISING FROM SUCH USE, INCLUDING ANY DAMAGES FROM PRODUCTS BASED ON, OR RESULTING FROM, RECIPIENT'S USE OF THE SUBJECT SOFTWARE, RECIPIENT SHALL INDEMNIFY AND HOLD HARMLESS THE UNITED STATES GOVERNMENT, ITS CONTRACTORS AND SUBCONTRACTORS, AS WELL AS ANY PRIOR RECIPIENT, TO THE EXTENT PERMITTED BY LAW.  RECIPIENT'S SOLE REMEDY FOR ANY SUCH MATTER SHALL BE THE IMMEDIATE, UNILATERAL TERMINATION OF THIS AGREEMENT. 

The “GeneLab Methylation Sequencing Processing Pipeline and Workflow” software also makes use of the following 3rd party Open Source software, released under the licenses indicated above.  A complete listing of 3rd Party software notices and licenses made use of in “GeneLab Methylation Sequencing Processing Pipeline and Workflow” can be found in the [Methyl-Seq_3rd_Party_Software.md](License/3rd_Party_Licenses/Methyl-Seq_3rd_Party_Software.md) file as indicated above. 

<br>

---
**Developed by:

Michael Lee
Crystal Han

**Maintained by:

