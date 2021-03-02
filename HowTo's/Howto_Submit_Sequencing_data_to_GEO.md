# What is GEO? :earth_americas:

The **Gene Expression Omnibus (GEO)** is a public repository that archives and freely distributes comprehensive sets of microarray, NGS, and other forms of high-throughput functional genomic data submitted by the scientific community. In addition to data storage, a collection of web-based interfaces and applications are available to help users query and download the studies and gene expression patterns stored in GEO. For more information about various aspects of GEO, please see our [FAQs page](https://www.ncbi.nlm.nih.gov/geo/info/faq.html#what).

Most journals require deposit of sequence data to a MIAME- or MINSEQE-compliant public repository like GEO for publication; on the contrary <ins>journal publication is not a requirement for data submission to GEO</ins>. It is recommended to upload data previous to submission, as many journals require accession numbers for the data before acceptance of a paper for publication. 

Records may remain private until your manuscript is publicly available. During the submission process, you are prompted to specify a release date for your records. Although the maximum limit is 3 years, this date may be brought forward or pushed back at any time. This feature allows a submitter to deposit data and receive a GEO accession number to quote in a manuscript before the data become public. If you need editors or reviewers to access your data, you can generate a reviewer access token to access your private GEO records during the review process. 

Once a GEO accession number is quoted in a manuscript (including publicly posted unpublished preprints through servers like bioRxiv), the records must be released so that the data are accessible to the scientific community. <ins>Omit GEO accession numbers from preprint manuscripts if you don’t want your data released before the official publication</ins>. If GEO accession numbers are found to be quoted in any publication before the scheduled release date, GEO staff are obligated to release those records, even if a second manuscript describing the same data is pending.  
</br>

## Submitting high-throughput sequence data to GEO
Previous to submission, you need to create a _My NCBI_ account (you can create one [here](https://www.ncbi.nlm.nih.gov/account/register/?back_url=/geo/submitter/)) and complete a _My GEO Profile_ form with the contact information to be displayed on the GEO records.

Every submission needs to contain raw data, processed data, and descriptive information about the samples, protocols and overall study in a supported deposit format. You can find detailed information on how to submit High Throughput Data here (https://www.ncbi.nlm.nih.gov/geo/info/seq.html).

### Metadata spreadsheet 
Download and fill in this [metadata spreadsheet](https://www.ncbi.nlm.nih.gov/geo/info/examples/seq_template.xlsx). It contains the descriptive information about the overall study, individual samples, all protocols, and references to processed and raw data file names.



GEO accepts next generation sequence data that examine quantitative gene expression, gene regulation, epigenomics or other aspects of functional genomics using methods such as RNA-seq, miRNA-seq, ChIP-seq, RIP-seq, HiC-seq, methyl-seq, etc. We process all components of your study, including the samples, project description, processed data files, and we submit the raw data files to the Sequence Read Archive (SRA) on your behalf.  
</br>

### Processed data files
Processed data are a required part of GEO submissions. The final processed data are defined as the data on which the conclusions in the related manuscript are based. We do not expect standard alignment files (e.g., BAM, SAM, BED) as processed data since conclusions are expected to be based on further-processed data. When standard alignments are the only processed data available, please write to us to inquire about whether your data are suitable for submission to GEO. Requirements for processed data files are not fully standardized and will depend on the nature of the experiment:

Expression profiling analysis usually generates quantitative data for features of interest. Features of interest may be genes, transcripts, exons, miRNA, or some other genetic entity. Two levels of data are often generated:
raw counts of sequencing reads for the features of interest, and/or
normalized abundance measurements, e.g., output from Cufflinks, Cuffdiff, DESeq, edgeR, etc.
Either or both of these data types may be supplied as processed data. They may be formatted either as a matrix table or individual files for each sample. Provide complete data with values for all features (e.g., genes) and all samples, not only lists of differentially-expressed genes.
ChIP-Seq data might include peak files with quantitative data, tag density files, etc. Common formats include WIG, bigWig, bedGraph.
Features (e.g., genes, transcripts) in processed data files should be traceable using public accession numbers or chromosome coordinates. The reference assembly used (e.g., hg19, mm9, GCF_000001405.13) should be provided in the metadata spreadsheet.

A description of the format and content of processed data files should be provided in the metadata spreadsheet data processing fields.

If you provide WIG, bedGraph, GFF, or GTF files, please refer to the UCSC file format FAQ for format requirements.  
</br>

### Raw data files
Raw data are a required part of GEO submissions. We will submit raw data files to SRA for you. The raw data files should be the original files containing reads and quality scores, as generated by the sequencing instrument (unless the raw files are barcoded/multiplexed, see below for further instructions).

Raw Data File Formats: Acceptable file formats include FASTQ, as well as other formats described in the SRA File Format Guide. Files that do not conform to supported format requirements will be deleted from our systems.

Barcode/Multiplexed Data: Whenever possible, we do require that files be demultiplexed so that each barcoded sample ends up with a dedicated run file. However, for single-cell sequencing studies (e.g. 10x Genomics, Drop-Seq, InDrops), we can support the submission of multiplexed files in cases where these files are required for reanalysis in your pipeline, or when demultiplexing would create an unmanageable number of files.

Paired-end Experiments: We usually expect 2 files per run (4 files per run when sequences and qualities are included in separate files). If submitting FASTQ files, please submit the original unedited files from the Illumina pipeline. Edited files may not be processed correctly by SRA.

MD5 Checksums: We recommend that submitters provide MD5 checksums for their raw data files. The checksums are used to verify file integrity. Checksums can be calculated using the following methods:

Unix: md5sum <file>
OS X: md5 <file>
Windows: Application required. Many are available for free download.
Data File Compression: Individual files can be compressed to speed transfer, but this is not required. Acceptable compression formats are gzip and bzip2 (i.e. files ending with a .gz or .bz2 extension). Never compress binary files (e.g., BAM, bigWig, bigBed), and DO NOT upload ZIP archives (files with a .zip extension).


Uploading your submission Back to top
There are two steps for submission:

1. Transfer all your files to the GEO FTP server. 
If your files are >1 terabyte, please contact us and do not transfer your files until you hear back from us. Transfer Files	2. After the FTP transfer is complete, notify GEO using the Submit to GEO web formNotify GEO

