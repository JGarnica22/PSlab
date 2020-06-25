# How to replace commas with newlines in Excel: :bar_chart:

This tutorial will show you how to replace commas within a list of elements in a single Excel cell by newlines, as represented in the image:

![](https://cdn.extendoffice.com/images/stories/doc-excel/dc-replace-commas-to-newline/doc-replace-commas-to-newline-1.png)

This can be achieved either by using **replace function** or **VBA coding**. \<ins>*In Mac, only the VBA option works!</ins>

## Replace function

Select the cells you will replace all commas with newlines. Then click **Home > Find & Select > Replace**.

![](https://cdn.extendoffice.com/images/stories/doc-excel/dc-replace-commas-to-newline/doc-replace-commas-to-newline-2.png)

In the opening Find and Replace dialog box and under the Replace tab, you need to:

1. Type a comma into the **Find what** box
2. Click on the **Replace** with box, then press the **Ctrl + Shift + J** keys simultaneously;
3. Click the **Replace All** button.

![](https://cdn.extendoffice.com/images/stories/doc-excel/dc-replace-commas-to-newline/doc-replace-commas-to-newline-3.png)

Then click OK when prompt box pops up and finally close dialog box.

## VBA Code

Select the cells containing the commas you need to replace with newlines, then press the **Alt + F11** keys simultaneously to open the **Microsoft Visual Basic** for Applications window.

**Click Insert > Module.** Then copy and paste VBA code into the Code window.

````
Sub ReplaceComma()
Dim rngCell As Range
  For Each rngCell In Selection
    rngCell.Value = Replace(rngCell, ",", vbLf)
  Next
End Sub
````

![](https://cdn.extendoffice.com/images/stories/doc-excel/dc-replace-commas-to-newline/doc-replace-commas-to-newline-5.png)

Press the **F5** key or click the **Run** button to run the code.


This tutorial was extracted from this [website.](https://www.extendoffice.com/documents/excel/4859-excel-replace-comma-with-newline-alt-enter.html)
