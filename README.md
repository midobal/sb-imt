# Segment-Based IMT
Segment-based Interactive Machine Translation implementation based on Moses XML scheme.

## Requirements
[Moses](https://github.com/moses-smt/mosesdecoder) and [mgiza](https://github.com/moses-smt/mgiza) are needed in order to run this software. Alternatively, you can [run it through Docker](#run-through-docker).

## Usage
### Alignments generation
This software requires of an alignment model. Therefore, prior to its use, alignments must be trained by doing:
```
tools/alignments.sh -s src_file -t tgt_file -o output_file -m mgiza_bin {options}

options:
     -b: use IBM Model 1. (Default: Use HMM.)
```

where `src_file` and `tgt_file` are the source and target of the training dataset; `output_file` is the file in which to store the alignments (which will be required for using the segment-based software); `mgiza_bin` is the path to mgiza's bin folder (e.g., */opt/moses/mgiza/mgizapp/bin*); and the `-b` flags switches from using *Hidden Markov Models* to using *IBM Model 1*.

### User simulation
You can simulate a user working on a segment-based IMT framework by doing:
```
simulation.py [-h] -s source_file -r reference_file -c moses_ini -a
                     alignments_file [-v] [-x] [-p threshold] [-m moses_bin]
                     -l log_file

optional arguments:
  -h, --help            show this help message and exit
  -s source_file, --sources source_file
                        file containing the source segments.
  -r reference_file, --references reference_file
                        file containing the reference segments.
  -c moses_ini, --config moses_ini
                        file containing moses configuration.
  -a alignments_file, --alignments alignments_file
                        file containing the alignments generated by alignments.sh
  -v, --verbose         Activate verbose mode.
  -x, --xml             Show XML markup.
  -p threshold, --probability threshold
                        probability threshold. (Default 0.)
  -m moses_bin, --moses moses_bin
                        Path to moses bin. (Default: /opt/moses/bin/moses.)
  -l log_file, --log log_file
                        File to store the log.
```

## Run through Docker
If you want to run this software through Docker you can have a look at this [repo](https://github.com/midobal/dockerfiles/tree/master/sb-imt).

## Citation
On using this software, please cite the following paper:

>Miguel Domingo and Álvaro Peris and Francisco Casacuberta. Segment-Based Interactive-Predictive Machine Translation. Machine Translation
Journal, 31:163–185, 2017
