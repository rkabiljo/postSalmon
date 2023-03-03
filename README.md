# postSalmon
analyse Salmon results

Sulev's annotation had gene ID appended to transcript ID, like this <br>

```
Name	Length	EffectiveLength	TPM	NumReads
ENST00000456328.2|ENSG00000223972.5|OTTHUMG00000000961.2|OTTHUMT00000362751.1|DDX11L1-202|DDX11L1|1657|processed_transcript|	1657	1307.000	0.000000	0.000
```
and I needed it to be like this:
```
Name	Length	EffectiveLength	TPM	NumReads
ENST00000456328.2 1657 1356.000 0.000000 0.000
```
Run this to change it, just change the path:
```
 find /Users/renatakabiljo/Documents/MNDquant/pawsey0360/MNDsalmon/MNDquant/ -type f -name '*.sf' -exec sh -c 'awk -F "\t" "{gsub(/\|.*/, \"\", \$1)}1" "$1" > tmp && mv tmp "$1"' _ {} \; 
```
