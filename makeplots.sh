#!/bin/sh

for i in "p" "phi" "chi" "pip" "pic" "lnF(p)" "lnF(phi)" "lnF(chi)" "lnF(pip)" "lnF(pic)"; do
    python plot.py "$i" "$i.pdf"
done
