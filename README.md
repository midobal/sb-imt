# Segment-Based IMT
Segment-based Interactive Machine Translation implementation based on Moses XML scheme.

## Requirements ##
Python 2.

## Usage ##
```
python simulation.py -s source_file -r reference_file -m moses_ini -a alignments [options]
```

### Options ###
```
-v             verbose mode on.
-XML           show XML Markup.
-p theshold    probability threshold. (Default: 0.)
```

## Citation ##
On using this software, please cite the following paper:

  Domingo, M., Peris, Á., and Casacuberta, F. (2016). Interactive-predictive translation based on multiple word-segments. In Proceedings of the Annual Meeting of the European Association for Machine Translation, pages 282–291.
