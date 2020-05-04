# How to download files from Google drive using `wget` :arrow_down:

Files can be downloaded from google drive using `wget`. Before that you need to know if your files are small or large sized in Google Drive (Files less than 100MB are regarded as small files, whereas files greater than 100MB are regarded as large files).

Before the file can be downloaded it has to be shared publicly or to everyone with the sharing link.

Steps:

1. Select a file that has to be downloaded and do right click.
2. Click **Share** and a dialog box will appear.
3. Click Advance in the right bottom corner.
4. Click on the Change.. under who has access.
5. Make it Public or accessible with the link.
6. Click Save button.
7. Copy the link for sharing, e.g. https://drive.google.com/file/d/1-uiB7SntziHO7ItAxZ2PuQGcsadG8ZJ1/view?usp=sharing
8. Extract FILEID, in this case: **1-uiB7SntziHO7ItAxZ2PuQGcsadG8ZJ1**

For small files, run the following command on your terminal:
````
wget --no-check-certificate 'https://docs.google.com/uc?export=download&id=FILEID' -O FILENAME
````

In the above command change the **FILEID** by above id extracted and rename **FILENAME** as the name you want to save the file with.  
<br/>

For large files, run the following command with accordinly changes in FILEID and FILENAME:
````
wget --load-cookies /tmp/cookies.txt "https://docs.google.com/uc?export=download&confirm=$(wget --quiet --save-cookies /tmp/cookies.txt --keep-session-cookies --no-check-certificate 'https://docs.google.com/uc?export=download&id=FILEID' -O- | sed -rn 's/.*confirm=([0-9A-Za-z_]+).*/\1\n/p')&id=FILEID" -O FILENAME && rm -rf /tmp/cookies.txt
````
