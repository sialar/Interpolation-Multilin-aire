
pdflatex main.tex

mkdir temp
mv main.tex main.pdf temp
rm -rf main.* ressources/body.aux
mv temp/* .
rm -rf temp

evince main.pdf
