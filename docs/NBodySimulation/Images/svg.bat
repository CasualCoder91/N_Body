latex cone.tex
dvips -E cone.dvi
ps2pdf -dAutoRotatePages#/None cone.ps
"E:\Programms\pdf2svg\dist-64bits\pdf2svg.exe" cone.pdf cone.svg
pause
