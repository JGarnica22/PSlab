# How to use Markdown :page_facing_up: :pencil2:

Markdown is a lightweight markup language with plain-text-formatting syntax. Its design allows it to be converted to many output formats: HTML, PDF, Word documents and slides. Markdown is often used to format readme files, and to create rich text using a plain text editor creating an .md file (MD).<br/>

Implementations of Markdown are available for over a dozen programming languages; in addition, many platforms and frameworks support Markdown. For example, Markdown plugins exist for every major blogging platform.<br/>

Examples:
- **GitHub** Flavored Markdown (GFM) ignores underscores in words, and adds syntax highlighting, task lists, and tables.
- **RStudio**, an IDE for R, provides a C++ wrapper function for a markdown variant called sundown.
- The sourcecode documentation generator Doxygen supports Markdown with extra features.
- Discount – a C implementation.
- MarkAPL is a converter written in Dyalog APL. It supports fenced blocks, smart typography, link references, and special attributes, and can generate a table of contents.
- PHP Markdown - a library package that includes the PHP Markdown parser and its sibling PHP Markdown Extra with additional features.
- Markdig – A .NET library that follows the CommonMark specifications, and includes a collection of extensions and the ability for the user to create their own.

Markdown is often used <ins>to format readme files</ins>, and to create rich text using a plain text editor creating an **.md file (MD)**. An MD file is a text file created using one of several possible dialects of the Markdown language. It is saved in plain text format but includes inline text symbols that define how to format the text (e.g., bold, indentations, headers, table formatting).<br/>

MD files are not only useful for HTML documentation systems, but also for source code version control. This is because the files can be compared against historical revisions in human-readable text (whereas a binary format cannot be easily compared).  

Projects created with GitHub, often use a file named README.md, which contains the readme for the project.  

To sum up, it's a very simple language used to create beautiful and presentable readme files that summarize your project where you can include text but also code and figures in the same file.  
<br/>


## SYNTAX: How to format a Markdown file

### Headings
To create a heading, add a hash symbol (#) in front of a word or phrase. The number of hash symbols you use should correspond to the heading level. For each extra (#) that you use your heading will be smaller.<br/>

### Paragraphs
To create paragraphs, use a blank line to separate one or more lines of text.<br/>

### Line breaks
To create a line break, end a line with two or more spaces, and then type \<br/> (return).<br/>

### Emphasis

#### Bold
To bold text, add two asterisks or underscores before and after a word or phrase. To bold the middle of a word for emphasis, add two asterisks without spaces around the letters. 

   **example**
````
**example**
__example__
````

#### Italic
To italicize text, add one single asterisk or underscore before and after a word or phrase. To italicize the middle of a word for emphasis, add one asterisk without spaces around the letters.

   *example*
````
*example*
_example_
````

#### Bold and Italic
To emphasize text with bold and italics at the same time, add three asterisks or underscores before and after a word or phrase. To bold and italicize the middle of a word for emphasis, add three asterisks without spaces around the letters.

   ***example***
````
***example***
___example___
````

### Ordering
#### Bullet Points
To create bullet points add a hyphen or an asterisk and a space before the word or sentence to write.

* Example 1
- Example 2
````
* Example 1
- Example 2
````

#### List
To create an ordered list, add line items with numbers followed by periods. The numbers don’t have to be in numerical order, but the list should start with the number one.

1. Example 1
2. Example 2
````
1. Example 1
2. Example 2
````

### Special characters without effect
To produce a literal character (an asterisk or underscore, for example) at a position where it would otherwise be used as a formatting character, you can backslash escape it.  

\*this text is surrounded by literal asterisks\*
````
\*this text is surrounded by literal asterisks\*
````


### Naming tools with backticks
The backtick (also known as the grave accent or backquote) is used to mention a tool or package. Include a backtick after and before the world that you want to highlight.<br/>

`Example`
````
`Example`
````


### Introduce chunks of code
To insert chunks of code into your file type 4 grave accents \`\`\`\` before and after the code. In a conventional .md file (to read in GitHub or export to html or PDF), the code chunks will appear within a blue background.<br/>

\`\`\`\`<br/>
Code Example <br/>
\`\`\`\`<br/>

````
Code Example
````


### Add links
To embed links into a topic you can either add the link (will be seen as the whole link text) or use the markdown syntax as below where the word "link" will be converted to a clickable link (use [ ] around the linkable word, followed without space for the link in brackets.
````
Go to this [link](https://write_here_your_link.com)
````

### Embed images
To create an inline image link, enter an exclamation point ( ! ), wrap the alt text in brackets ( [ ] ), and then wrap the link in parenthesis ( ( ) ). (\*_Alt text_ is a phrase or sentence that describes the image.)
````
![Alt Text](url)
````

### Embed emojis :blush:
You can also add emojis to your text or report copying its code (\:emojiX\:). You can find the code for different emojis in multiple web pages, for example this one: https://gist.github.com/rxaviers/7360908  
<br/>

For more information on Markdown you can visit these web pages:<br/>
- https://www.markdownguide.org/basic-syntax/  
- https://rstudio.com/resources/cheatsheets/  
- https://rstudio.com/resources/webinars/getting-started-with-r-markdown/  
<br/>

## How to use R Markdown
R Markdown provides an authoring framework for data science. You can use a single R Markdown file to both:

- Save and execute code.
- Generate high quality reports that can be shared with an audience.

R Markdown documents are fully reproducible and support dozens of static and dynamic output formats ([more info](https://rmarkdown.rstudio.com/)).

### Install and open Markdown
You can install it directly from RStudio:  

Tools > Install Packages <br/>
Choose "Repository CRAN" and write "Markdown" in "Packages (separate multiple with space or comma:)"

To open a new file, click File > New File > R Markdown in the RStudio menu bar. 
A window will pop up where you can name the project an the author ("Title" and "Author") and select the specific type of output that you wish to build (HTML, PDF or Word). Remember that you can switch to the other output formats anytime.  


A template will be opened to generate your Markdown report. The Markdown template contains:

**1. A description of your document:**
It includes info on the title, author, date and output format. This description will apear in the final document; you can change it at any moment.

````
title: "Test"
author: "Mireia"
date: "11/05/2020"
output: pdf_document
````  

**2. A first R code chunk that you should leave by default:**

````
{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
````  

**3. Text.** 
You can add text along your file to describe your project or explain each code chunk or results. Use Markdown syntax as explained above.

````
## R Markdown

This is an R Markdown document. Markdown is a simple formatting syntax for authoring HTML, PDF, and MS Word documents. For more details on using R Markdown see <http://rmarkdown.rstudio.com>.

When you click the **Knit** button a document will be generated that includes both content as well as the output of any embedded R code chunks within the document. You can embed an R code chunk like this:

````  

**4. R code chunks.** 
When you render your .Rmd file, R Markdown will run each code chunk and embed the results beneath the code chunk in your final report (unless you use special chunk options described below). Note that you can run that code directly from that window and you can see the progress in the "Console" window of R Studio.  

There are two ways to include code to the document:

**A) Embeded code:** 
Similarly to what was explained above for conventional markdown, you can add chunks of code in R Markdown by using 3 backticks before the code \`\`\` followed by {r _info to contextualize this chunk_ and _code chunk options_}. Close the chunk by using again 3 backticks after code \`\`\`. When you compile, R markdown will run the code and by default include the code and its results in the output. R markdown will remove the first \`\`\`{r} and last \`\`\` lines.  

````
```{r}
Some code
```
````

Options in the brackets after {r}:

* echo = FALSE --> hides the code (useful if you want to show only the results/plot).
* eval = FALSE --> prevents the code from being run. As a result, no results will be displayed with the code.
* message = FALSE --> supresses messages from appearing in the output, for example warnings.
* engine = 'python' --> to embed non R code, set the engine option to the language you want to embed.  

See examples below:

````
```{r echo = FALSE}
Some code
```
````

![](https://github.com/patriciasolesanchez/PSlab/blob/master/HowTo's/Screenshots/echo_example.png) 

I WOULD SHOW ONE IMAGE WITH echo=FALSE AND ONE WITH echo=TRUE SO THAT THE DIFFERENCE IS CLEAR. SAME FOR THE OTHER OPTIONS!!

ALTERNATIVELY, SHOW A FIRST IMAGE (AS A REFERENCE) WITH NO OPTIONS, SO THAT THE OPTIONS BELOW ARE COMPARED TO THIS ONE. BUT THEN ALL EXAMPLES OF CODE WOULD NEED TO BE THE SAME... I THINK THE FIRST OPTION MIGHT BE EASIER :)


````
```{r eval = FALSE}
Some code
```
````

![](https://github.com/patriciasolesanchez/PSlab/blob/master/HowTo's/Screenshots/eval_example.png) 

````
```{r message = FALSE}
Some code
```
````

![](https://github.com/patriciasolesanchez/PSlab/blob/master/HowTo's/Screenshots/message_example.png) 

````
```{r engine = 'python'}
Some code
```
````

![](https://github.com/patriciasolesanchez/PSlab/blob/master/HowTo's/Screenshots/engine_example.png) 


More chunk options in: https://rstudio.com/resources/cheatsheets/  
<br/>

**B) Inline code:** 
Place code in a sentence with \`r _followed by your code_\`. R markdown will replace the code with its results.
````
Today is `r Sys.date()`
````
_Today is 2020-05-20_
<br/>

I DONT' UNDERSTAND WELL HOW THIS INLINE CODE WORKS. I TRIED IN A MARKDOWN AND RESULTS APPEAR WEIRD. LET'S DISCUSS.  
<br/>

Finally, to get your report you can click on "Knit" botton and choose your favourite format: HTML, PDF or Word.

I strongly recommend you this webinar to start using R Markdown:
https://rstudio.com/resources/webinars/getting-started-with-r-markdown/

