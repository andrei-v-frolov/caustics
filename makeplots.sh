#!/bin/sh

for i in "phi" "chi" "pip" "pic" "lna" "p" "lnF(p)" "lnF(phi)" "lnF(chi)" "lnF(pip)" "lnF(pic)"; do
    python plot.py "$i" "$i.pdf"
done
