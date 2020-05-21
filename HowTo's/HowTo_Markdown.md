# How to use Markdown

Markdown is a simple formatting syntax for authoring HTML, PDF, and MS Word documents. Its design allows it to be converted to many output formats, but the original tool by the same name only supports HTML. Markdown is often used to format readme files, for writing messages in online discussion forums, and to create rich text using a plain text editor.<br/>

EXPLAIN EXTENSION (.md)

EXPLAIN WHERE YOU CAN CREATE A MD. For example: 
- You can create mds in GitHub to make nice-readable files, with format (vs. an unformatted txt file)
- You can create mds in an app??
- You can create mds in R to make nice reports including text, code and the figures all in one file.


## SYNTAX: How to format a Markdown file

### Headings
To create a heading, add a hash symbol (#) in front of a word or phrase. The number hash symbols you use should correspond to the heading level. For each extra (#) that you use your heading will be smaller.<br/>

### Paragraphs
To create paragraphs, use a blank line to separate one or more lines of text.<br/>

### Line breaks
To create a line break, end a line with two or more spaces, and then type return.<br/>
````
<br/>
````

### Emphasis

#### Bold
To bold text, add two asterisks or underscores before and after a word or phrase. To bold the middle of a word for emphasis, add two asterisks without spaces around the letters. 
````
**example**
__example__
````

### Italic
To italicize text, add one asterisk or underscore before and after a word or phrase. To italicize the middle of a word for emphasis, add one asterisk without spaces around the letters.
````
*example*
_example_
````


### Bold and Italic
To emphasize text with bold and italics at the same time, add three asterisks or underscores before and after a word or phrase. To bold and italicize the middle of a word for emphasis, add three asterisks without spaces around the letters.
````
***example***
___example___
````

## Ordering
### Bullet Points
To create bullet points add a hypen and a space before the word or sentence to write.
````
* Example 1
* Example 2
````

### List
To create an ordered list, add line items with numbers followed by periods. The numbers don’t have to be in numerical order, but the list should start with the number one.

````
1. Example 1
2. Example 2
````

For more information you can visit this web pages: 
- https://www.markdownguide.org/basic-syntax/
- https://rstudio.com/resources/cheatsheets/
- https://rstudio.com/resources/webinars/getting-started-with-r-markdown/

# How to use R Markdown
R Markdown provides an authoring framework for data science. You can use a single R Markdown file to both:

- save and execute code.
- generate high quality reports that can be shared with an audience.

R Markdown documents are fully reproducible and support dozens of static and dynamic output formats. (more info: https://rmarkdown.rstudio.com/)

### Install and open Markdown
You can install it directly from RStudio:
Tools > Install Packages
Choose "Repository CRAN" and write "Markdown" in "Packages (separate multiple with space or comma:)"

To open a new file, click File > New File > R Markdown in the RStudio menu bar. A window will pop up where you can name the project an the author ("Title" and "Author") and select the specific type of output that you wish to build (HTML, PDF or Word). Rememeber that you can swith to the other output formats anytime.

A template will be opened to generete your Markdown report. Here you have a summary of the template in order of appearance:

**1. Important information of your document such as the title, author, date (will apear in the final document) and output format (you can change it at any moment):**

````
title: "Test"
author: "Mireia"
date: "11/05/2020"
output: pdf_document
````

**2. An R code in a grey square that you should leave by default:**

````
{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
````

**3. Text about Markdown. Note that this is for you to edit with info about your specific project. After the two hash symbols you can write the title of the project that you want to show. In the paragraph write the description of the project, what you are going to show or the results of your analysis (or everything).
You can use this paragraph with different information in multiple ocasions during the report.**

````
## R Markdown

This is an R Markdown document. Markdown is a simple formatting syntax for authoring HTML, PDF, and MS Word documents. For more details on using R Markdown see <http://rmarkdown.rstudio.com>.

When you click the **Knit** button a document will be generated that includes both content as well as the output of any embedded R code chunks within the document. You can embed an R code chunk like this:

````
**4. R code that you want to show. Note that you can run that code directly from that window and you can see the progress in the "Console" window of R Studio.**

There are two way to include code to the document:

**1. Embeded code:** Insert a chunck of code as the example below. When you compile, R markdown will run the code and include its results. R markdown will also remove the ```{r} and ```. Yo can hide the code by clicking the triangle botton that there are in the line.
````
```{r}
# some code
```
````

Useful tips options in the brackets after r:

* echo = FALSE --> hides the code (useful if you want to show only a plot).
* eval = FALSE --> prevent the code from being run. As a result, no results will be displayed with the code.
* message = FALSE --> supresses messages from appearing in the output, for example warnings.
* engine = 'python' --> to embed non R code, set the engine option to the language you want to embed. 
````
```{r echo = FALSE}
# some code
```
````
````
```{r eval = FALSE}
# some code
```
````
````
```{r message = FALSE}
# some code
```
````
````
```{r engine = 'python'}
# some code
```
````
More chunck options in: https://rstudio.com/resources/cheatsheets/

**2. Inline code:** Place code in a sentence with `r #code`. R markdown will replace the code with its results.
````
Today is `r Sys.date()` --> Today is 2020-05-20
````

**5. A place to embed your plots with echo = FALSO to avoid showing the code:**
```{r name of the plot, echo=FALSE}
plot(name of the plot)
```

Finally, to get your repot you can click on "Knit" botton and choose your favovourite format: HTML, PDF or Word.

I stongly recommend you this webinar to start using R Markdown:
https://rstudio.com/resources/webinars/getting-started-with-r-markdown/
